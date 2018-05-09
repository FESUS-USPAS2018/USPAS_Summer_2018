In this directory are a number of python scripts with basic functionality.  All scripts have been written with argparse support; this allows a common command line interface that will be described shortly.

To run any of the scripts, you need to first [install the correct modules](https://github.com/FESUS-USPAS2018/USPAS_Summer_2018/tree/master/python_sandbox).  If you are using a sandbox, you also need to [turn it on](https://github.com/FESUS-USPAS2018/USPAS_Summer_2018/tree/master/python_sandbox#turning-the-sandbox-on-and-off).  Then navigate to the directory where the scipt is stored and run

```
python script_name.py -h
```

where script_name.py is to be replaced by the name of the script you would like to run.
This command will print to the screen a description of the script, a description of the necessary arguments, and a list of optional arguments.  For instance

```
python photoemission_coordinates.py -h
```

gives the output

```
usage: photoemission_coordinates.py [-h] [-n E_PER_MPART] [-t STD_T] [-R R]
                                    [-p PROFILE] [-o OUTPUT_FILE]
                                    [--fermi_energy E_FERMI]
                                    [--work_potential E_WORK]
                                    [--photon_energy E_PHOTON] [-u UNITS]
                                    [-d DELIMITER]
                                    n

Generates electron macroparticle phase coordinates according to the three setp
model. Output is t,x,y,px,py,pz.

positional arguments:
  n                     The number of macroparticles desired.

optional arguments:
  -h, --help            show this help message and exit
  -n E_PER_MPART, --number_of_electrons_per_macroparticle E_PER_MPART
                        The number of electrons simulated by each
                        macroparticle. Default is 100.
  -t STD_T, --std-t STD_T
                        The standard deviation of the gaussian of electron
                        emission in time. Default is about 21 fs.
  -R R, -r R            The radial length. For the uniform distribution, this
                        is the radius. For the Guassian, the standard
                        deviation. Default is about 96 um.
  -p PROFILE, --emission_profile PROFILE
                        Describes the profile of the electrons in the radial
                        direction. Options are Gaussian or uniform.  Default 
                        is Gaussian.
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        The path to where the output will be written. Defualt
                        is IniConditions.txt in the current working directory.
  --fermi_energy E_FERMI
                        The fermi energy of the electron source. Default is
                        5.000E-6.
  --work_potential E_WORK
                        The work potential of the electron source. Default is
                        4.450E-6.
  --photon_energy E_PHOTON
                        The energy of the lazer's photon. Default is 4.650E-6.
  -u UNITS, --units UNITS
                        Specifies the units of the momenta. Default is m/s
                        (mks). Supports mks and MeV/c (specify "MeV").
  -d DELIMITER, --delimiter DELIMITER
                        Specifies the delimiter used in the output file.
                        Default is " ".
```

To use this script to generate 100 photo-emitted 10 electron macroparticles coordinates with the output deliminated by "," and stored in the file "first_100_particles.csv", you'd issue the command

```
python photoemission_coordinates.py 100 -n 10 -d , -o first_100_particles.csv
```

This should then create (or overwrite) the file "first_100_particles.csv" in the directory where you ran the code with 100 lines of t,x,y,px,py,pz coordinates.

Similar documentation is provided for all scripts in this directory.
