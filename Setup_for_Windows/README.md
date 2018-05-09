# Motivation
In this course, we are going to work with 
- a python code
- a FORTRAN code
- a widely used Particle-In-Cell simulation code *WARP*

*WARP* only works with Linux and MacOS X, so the easiest way is to take the advantage of the new feature of Windows 10 *Windows Subsystem for Linux (WSL)* to create a Linux enviroment inside your windows.

If you are a Windows7 user, you can either install a Linux alongside Windows to create a dual boot, or use Cygwin in Windows7. Please let us know by email (xiangxuk@msu.edu), if you have any problems to get a Linux environment prepared.

# Install the Windows Subsystem for Linux (WSL)
[Here](https://docs.microsoft.com/en-us/windows/wsl/install-win10) is the detailed instructions from Microsoft about how to install a WSL. 

Then you need to update your Linux enviroment by
```bash
sudo apt update
```

## Install Pip
```bash
sudo apt install python-pip
```

You're done! Now you are ready for the next step with other Linux/MacOS users.
