import numpy as np
from config.my_config import parse_key_as_numpy_array
from config.simulation_type import get_mesh_symmetry_factor
from config.elements import load_elements
from class_and_config_conversion import set_attributes_with_config_section
from injectors.io import get_distribution_average

def moving_bunch_parameters(config,args,top,w3d,f3d,particle_coordinate_file):
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
    particle_coordinate_file:  The filepath to the input file containing particle coordinates.  Each row is a particle and columns correspond to x,y,z,px,py,pz.
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
    adv_dt[0]
  except:
    adv_dt = np.array([adv_dt])
    adv_steps = np.array([adv_steps])

  #Get the average particle z-coordinate to establish where the center of the grid should be.
  zavg = get_distribution_average(particle_coordinate_file)

  #Mesh size
  dx = config.get("Simulation parameters","dx")
  dz = config.get("Simulation parameters","dz")
  xmax = config.get("Simulation parameters","xmax")
  z_extent = config.get("Simulation parameters","z_extent")

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
  w3d.nz = int(2*z_extent/dz)
  w3d.nz = int(z_extent/dz)
  w3d.xmmin = -xmax
  w3d.xmmax =  xmax
  w3d.ymmin = -xmax
  w3d.ymmax =  xmax 
  w3d.zmmin =  zavg - z_extent
  w3d.zmmax =  zavg + z_extent

  return top, w3d, f3d, adv_dt, adv_steps
