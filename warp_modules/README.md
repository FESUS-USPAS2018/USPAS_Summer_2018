# Overview
The intent of the scripts in this directory is to provide additional tools to do standard
functions with warp.  As most Warp programs are written with global variables
and in-line function definitions, I have spent considerable effort in creating work-arounds
that make it clear what variables are necessary for a function while simultaneoulsy
allowing the functions to be re-usable and modularized.

# Download and "Install"

Download the warp_modules directory to your computer (clone the USPAS_Summer_2018 directory).

Add Warp modules and config2class (also in USPAS_Summer_2018 directory) to your python path with

```
export PYTHONPATH=$PYTHONPATH:/path/to/warp_modules:/path/to/config2class/src
```

That's it.  Your python program can now call these modules.  The only thing is that you will need to
add them to your PYTHONPATH every time you start a new terminal session.

# Additional details
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


To run the program, add the directory where this README file is to your 
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
