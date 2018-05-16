import numpy as np
from config.my_config import parse_key_as_numpy_array
from config.simulation_type import get_mesh_symmetry_factor
from config.elements import load_elements
from class_and_config_conversion import set_attributes_with_config_section

def photoemission_parameters(config,args,top,w3d,f3d):
  """
  Sets up the parameters in the warp objects, top, w3d, and f3d, according to
  what's in the config file.  Returns these objects and some of the parameters that
  are used elsewhere.
  
  Args:
    config:  A config parser object cointaining the values specified in the configuration file.
    args: An argparser object containing values specified at the command line interface.
    top:  One of the 3 warp objects.
    w3d:  One of the 3 warp objects.
    f3d:  One of the 3 warp objects.
  Returns:
    top:  One of the 3 warp objects (edited).
    w3d:  One of the 3 warp objects (edited).
    f3d:  One of the 3 warp objects (edited).
    adv_dt:  A numpy array containing the different dt's used during simulation.
    adv_steps:  A numpy array containing the number of steps associated with each dt.
  """

  adv_dt = np.array(config.get("Simulation parameters", "adv_dt"))# Numpy array of dts
  adv_steps = np.array(config.get("Simulation parameters", "adv_steps")) #Number of steps for each dt 
  try:
    adv_dt.size
  except:
    adv_dt = np.array([adv_dt])
    adv_steps = np.array([adv_steps])

  #Mesh size
  dx = config.get("Simulation parameters","dx")
  dz = config.get("Simulation parameters","dz")
  xmax = config.get("Simulation parameters","xmax")
  zmin = config.get("Simulation parameters","zmin")
  zmax = config.get("Simulation parameters","zmax")


  #Load the parameters for the w3d and top objects from the config.
  set_attributes_with_config_section(top, config, "top parameters", {",":parse_key_as_numpy_array})
  set_attributes_with_config_section(w3d, config, "w3d parameters")
  set_attributes_with_config_section(f3d, config, "f3d parameters")

  sym_factor = get_mesh_symmetry_factor("wrz", top,w3d)

  #Derived top and w3d attribute modifications
  top.dt = adv_dt[0]     # Specify initial dt advance for generate etc 
  top.prwall = xmax  # cylinder radius to absorb particles if they exceed this distance
  w3d.nx = sym_factor*int(xmax/dx)
  w3d.ny = sym_factor*int(xmax/dx)
  w3d.nz = int(zmax/dz)
  w3d.xmmin = -xmax
  w3d.xmmax =  xmax
  w3d.ymmin = -xmax
  w3d.ymmax =  xmax 
  w3d.zmmin =  zmin
  w3d.zmmax =  zmax 

  return top, w3d, f3d, adv_dt, adv_steps

def handle_cathode_and_anode(config,Ea,V=None):
  """
  Loads the cathode and anode (actually all conductor elements) from the configuration file
  and applies the correct potential to them to create the
  desired extraction field.
  
  Args:
    config:  A config parser object cointaining the values specified in the configuration file.
    Ea: The extraction field strength in V.
    V:  (Optional)  The potential of the anode.  Overwrites the potential calculation from Ea if specifed.
  Returns:
    conductor_elements: A list of conductor elements containing the cathode and anode.
  """
  conductor_elements = load_elements(config,"Conductor elements")
  try:
    conductor_elements.cathode
    conductor_elements.anode
  except:
    raise Exception("Either the cathode or anode are not present.")
  if V is None:
    anode_z = conductor_elements.anode.zmin #The placement of the inner surface of the anode.
    cathode_z = conductor_elements.cathode.zcent + 0.5*conductor_elements.cathode.length#The placement of the inner surface of the cathode.
    d = anode_z - cathode_z #The separation distance between the diodes.
    V = Ea*d*1e6 #Assumed in units of MV/m
  conductor_elements.cathode.voltage = 0.  # Cathode bias [V]
  conductor_elements.anode.voltage = V # Anode bias [V]
  print("Cathode voltage = " + str(conductor_elements.cathode.voltage) + " V")
  print("Anode voltage = " + str(conductor_elements.anode.voltage) + " V")
  return conductor_elements
