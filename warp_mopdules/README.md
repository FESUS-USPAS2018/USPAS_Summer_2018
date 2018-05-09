warp_modules
Git repository of Warp simulation scripts/tools to model the UEM 
ultra-intense electron microscope experiment at 
Michigan State University. 

Professor Steven M. Lund
Physics and Astronomy Department
Facility for Rare Isotope Beams
Michigan State University
lund@frib.msu.edu
517-908-7291

Dr. Brandon Zerbe 
Physics and Astronomy Department
Michigan State University
zerbe@msu.edu 


To initialize the repository, 

To run the program, add the directory where this readme.txt file is to your 
python path:

  export PYTHONPATH=$PYTHONPATH:$PWD

This will allow python to find the modules in the import statement.

Before running these scripts, the following dependency needs to be resolved:

  Download the config2class utility and put it in your src directory (or wherever 
  you keep such things):

  git clone git@github.com:billyziege/config2class.git

and then the following directories need to be added to your python path:

  export PYTHONPATH=$PYTHONPATH:/path/to/warp_modules:/path/to/config2class/src

This directory contains the scripts:

  plot_conductors.py: Loads the section "Conductor elements" from the config file and
    plots the resulting electric field after initiallizing the grid.
  preprocess_field.py: Loads the field from the external file, ravels it for use in Fortran,
    writes it to a file, and writes the relevant statistics for loading to a neighboring 
    configuration file for easy loading.
  plot_field.py:  Loads the field from the input configuration file and outputs plots
    for the field.

To see other options for these scripts:
  % python ${script_name} -h
