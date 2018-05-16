# Load Warp 
description="""
Python input script for a Warp simulation of 
the MSU UEM experiment using input photoemission data. 

For help and/or additional information contact:

     Steve Lund     lund@frib.msu.edu    (517) 908-7291

For documentation see the Warp web-site:

     https://warp.lbl.gov

All code inputs are mks except where noted.  

To run this Warp script:
  Interactive (return to interpreter after executing):

    % python -i uem.py [options] input_file 

  Non-interactive (exit interpreter after executing):

    % python uem.py [options] input_file

To see a list and desciption of options:
   
    % python uem.py -h
    % python uem.py --help
    % python uem.py
"""
from warpoptions import *
#Handle command line arguments and default values with argparse.
parser.description = description
parser.formatter_class=argparse.RawDescriptionHelpFormatter
parser.add_argument('input_file', type=str, 
                    help='The path to the file containing the intial ' +
                    'conditions pickled dictionary - specified with the ' +
                    'keys time, x, y, z, px, py, and pz.')
parser.add_argument('config_file', type=str, 
                    help='The config file.  Contains a bunch of parameters ' +
                    'and references to other config files.')
parser.add_argument('-e','--extraction_field',  
                    dest="extraction_field", type=float,
                    help='Specifies the extraction electric field gradient in MV.', 
                    default=1.)
parser.add_argument('--potential',
                    dest='potential', type=float,  
                    help="specified the potential of the anode assuming the cathode is at 0 V.  This is an alternative aproach to specifying the extraction field, and this overwrites that approach.  Units in V.",
                    default=None)
parser.add_argument('-m','--electrons_per_macroparticle', 
                    dest="electrons_per_macroparticle", type=float, 
                    help='The number of electrons per macroparticle for the ' +
                    'simulation.  The defaults to 100', default=100.)
parser.add_argument('--turn_off_field_solver', dest='field_solver_off',
                    action='store_true', help='Turn off field solver (for runs ' +
                    'with applied field but no self fields).  Default has the fieldsolver ' + 
                    'on.',default=False)
parser.add_argument('--iterative_phase_space_dump', dest="iterative_dump", type=int,
                    help='Tells the program to print the x, y, z, px, py, pz coordinates ' +
                    'of all the macroparticles at every iterative_dump steps.  Default is ' + 
                    'to skip this step.', default=None)
parser.add_argument('--phase_space_dump_step_list', dest="dump_list", type=str,
                    help='Tells the program to print the x, y, z, px, py, pz coordinates ' +
                    'of all the macroparticles at the listed steps separated by a comma.  ' + 
                    'For example, 35,46,72 will tell the program to dump the  coordinats after ' + 
                    'the 35th, 46th, and 72nd steps.  Default is ' + 
                    'to skip this step.', default=None)

args = parser.parse_args()

print "Argument dictionary: " 
print "\t" + "\n\t".join([k + " = " + str(v) for k, v in vars(args).iteritems()])

from config.my_config import MyConfigParser, handle_rz_1D_mesh_refinements
from config.simulation_type import get_solver
from diagnostics.diagnostic_classes import DiagnosticsByTimes, DumpBySteps
from diagnostics.phase_volume import dump_phase_volume
from diagnostics.steves_uem_diagnostics import steves_plots, electric_potential_plots
from discrete_fourspace.mesh import get_supremum_index
from injectors.injector_classes import ElectronInjector
from injectors.steves_uem_injection import steves_injectelectrons
from injectors.io import photoemission_loader
from photoemission_support import photoemission_parameters, handle_cathode_and_anode
from warp import *
#from histplot import *

# Invoke setup routine: needed to created a cgm file for plots
setup()

#Load the config file and parameters from config file.
config = MyConfigParser()
config.read(args.config_file)

#Load the parameters into the warp objects and register solver
top, w3d, f3d, adv_dt, adv_steps = photoemission_parameters(config,args,top,w3d,f3d)
solver = get_solver("wrz", top, w3d)
registersolver(solver)

#Install the conductor elements
conductor_elements = handle_cathode_and_anode(config,args.extraction_field,args.potential)
for conductor_element in conductor_elements:
   installconductor(conductor_element)  # install conductors into the warp framework
   print "Installed "+str(conductor_element)

#Add grid refinement near diodes
handle_rz_1D_mesh_refinements(config,solver,conductor_elements,config.get("Simulation parameters","xmax"))

# Create the electron beam species.  This function reads in input coordinates
# and saves them in particle objects.  The time coordinate is used to specify 
# at what step the program will ``create'' the relavent particles. 
electron_injector = ElectronInjector(steves_injectelectrons,top, args.input_file,
                    top.echarge/top.emass, args.electrons_per_macroparticle,loader=photoemission_loader)
installuserinjection(electron_injector.callFunction)  # install injection function in timestep 

#Diagnostics to the cgm file.
diagnostic_time_interval = config.get("Simulation parameters","diagnostic_time_interval")
t_final   = sum(adv_steps*adv_dt)  # Total time advance [s] 
diagnostic_times = arange(diagnostic_time_interval,t_final,diagnostic_time_interval) #List of times to run diagnostics.
diagnostics = DiagnosticsByTimes(steves_plots,top,top,diagnostic_times)
installafterstep(diagnostics.callFunction) # install function steves_plots() to be called after each timestep.  This goes to the cgm file.

#The following two if statements handle writing (dumping) the phase coordinates of all of the macroparticles in the simulation.
#These are completely optional.
#The iterative dump writes a file every specified step.
#The dump_list requires the specification of certain iteration steps.
#They have the same interface as they feed the required steps for dumping to the function;
#so the if statement first defines the steps for dumping.  The function will check the current
#step against this list, and if it is the correct step, it will dump the data.
if args.iterative_dump is not None: 
  steps_tot = sum(adv_steps)         # Total steps to take
  dump_steps = range(args.iterative_dump,int(steps_tot),args.iterative_dump) #Every iterative_dump step.
  phase_volume_dump = DumpBySteps(dump_phase_volume,electron_injector.getElectronContainer(),
        args.electrons_per_macroparticle*top.emass,top,dump_steps)
  installafterstep(phase_volume_dump.callFunction) # install function steves_plots() to be called after each timestep.  This goes to the cgm file. 
if args.dump_list is not None: 
  dump_steps = [int(s) for s in args.dump_list.split(",")]
  phase_volume_dump = DumpBySteps(dump_phase_volume,electron_injector.getElectronContainer(),
        args.electrons_per_macroparticle*top.emass,top,dump_steps)
  installafterstep(phase_volume_dump.callFunction)

#Set up the simulation.
package("w3d") 
generate() 

# Generate initial plots
solver.drawboxzr()#Draws mesh refinement boxes.
ix_cen = get_supremum_index(w3d.xmesh,0)
iy_cen = get_supremum_index(w3d.ymesh,0)
electric_potential_plots(ix_cen,iy_cen)

if args.field_solver_off:
  print "Turning field solver off."
  solver.ldosolve = False

# Advance simulation through each step interval set  
for ii in range(len(adv_steps)):
  top.dt = adv_dt[ii]
  step(adv_steps[ii])  

# Print out timing statistics of run 
printtimers() 

# Make sure that last plot is flushed from buffer
fma() 
