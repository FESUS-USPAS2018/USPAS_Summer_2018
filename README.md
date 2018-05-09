# Introduction
In this summer course, we will be working with 
* a Python code for photoemission process
* a FORTRAN code for N-particle simulation inside the gun
* *Warp*, a widely-used Particle-In-Cell simulation code, and python scripts using this code

Here is the instructions on software installation to be completed before coming to the course.

# Get a Linux environment (Windows Users only):
Since *Warp* only supports Linux and Mac OS, Windows users need to construct a Linux environment first. For Windows users, please check out the folder [Setup_for_Windows](https://github.com/FESUS-USPAS2018/USPAS_Summer_2018/blob/master/Setup_for_Windows) for further instructions. 

# Python preparation:
Once you have a Linux/MacOS environment, you need to install necessary libraries/packages/modules.  Directions setting up an isolated "sandbox" can be found  
[here](https://github.com/FESUS-USPAS2018/USPAS_Summer_2018/tree/master/python_sandbox) is some instructions and tips.  We will use this environment both for running *Warp* as well as running our own scripts.

# Install *Warp*
After setting up the right Python environment, we can now install *Warp*. [Here](https://github.com/FESUS-USPAS2018/USPAS_Summer_2018/tree/master/warp) are some instructions/tips.

# Python code for photoemission
One of the central codes we will discuss in the course will be a python script that produces samples of particles drawn from distributions derived from Spicer's three step model.  This script, and other other data visualization scripts, can be found [here](https://github.com/FESUS-USPAS2018/USPAS_Summer_2018/tree/master/python_scripts).

# FORTRAN code for N-particle simulation
You can find the N-partice simulation code we are going to use [here](https://github.com/FESUS-USPAS2018/USPAS_Summer_2018/tree/master/N_Particle_simulation). Since we've already installed a FORTRAN compiler during *Warp* installation, such as *gfortran*, no extra preparation is needed at this point.

# Warp supporting python modules
*Warp* is a framework that provides access to Fortran-based PIC simulations.  While working with the *Warp*, we have built tools that streamline some central functions.  We probably will not use all of these tools, but we will definitely use some of them.  We have made our scripts available to all here. 

# Contact Information
You may run into issues during installation, or you may have other concerns.  Please feel free to email us at:

Brandon Zerbe: brandon.s.zerbe@gmail.com

Xukun Xiang: xiangxuk@msu.edu
