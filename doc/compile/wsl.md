# Instructions to run Cytosim on *Windows Subsystem for Linux*

For Windows 10 or higher, Cytosim can run using the [Windows Subsystem for Linux](https://learn.microsoft.com/en-us/windows/wsl/about), 

## Installation

Please, follow the [instructions for your system](https://learn.microsoft.com/en-us/windows/wsl/install)

## Configuration

1.	Type `Turn Windows features on and off` in the search box of the taskbar, and manually enable:
	- Virtual Machine Platform
	- Windows Subsystem for Linux. 
Click OK and restart Windows.

2.	Go to Windows Microsoft Store and search for Ubuntu. Download and install the current release (Ubuntu 20.04 LTS as of 03.2022)

3.	Click Launch to start the Ubuntu Terminal. If asked to download the latest WSL2 Linux Kernel, proceed by following the instructions that will be provided. Install the update and relaunch Ubuntu. 

4.	Follow the terminal instructions to Enter new UNIX username and set a password. Once the user account is set up, launch Ubuntu from the search box of the Windows task bar. Open the home directory (aka “.”), from the command window:

        > explorer.exe .

5.	Install a WSL-compatible X-Window server, for example, [Xming](https://sourceforge.net/projects/xming/)

6.	Start the X-Window server, Xming by double-clicking on the Xming icon; select Multiple Windows. 

Hereafter, open the Ubuntu terminal and follow the instruction to download and install Cytosim on Linux.


### Authors

Carlos Lugo, Eashan Saikia, and Francois Nedelec, last updated 21.02.2024
