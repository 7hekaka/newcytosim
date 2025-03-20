# Cytosim

Cytosim is an [open-source](LICENSE.txt) cytoskeleton simulation suite designed to model large systems of flexible filaments and their associated proteins, such as molecular motors. It has been extensively used to study actin and microtubule dynamics in 1D, 2D, and 3D.

Built around an optimized C++ core engine, Cytosim runs on UNIX, macOS, GNU/Linux, and within the Windows Subsystem for Linux (WSL). Its modular and extensible design makes it highly adaptable, allowing customization for specific research needs. Additionally, a versatile set of Python programs streamlines the standard simulation workflow.

Cytosim has been used in approximately 100 research articles from the [Nedelec](doc/publications/index.md) and [other research groups](doc/publications/biblio.md).

![Cytosim](doc/data/cytosim.png)

Cytosim offers command-line tools with simulation and display capabilities. Simulations are defined using a [configuration file](doc/sim/config.md) which specifies objects, their parameters, and a sequence of operations such as advancing time, saving frames, or [generating reports](doc/sim/report.md). Below is a basic example, where parameters are expressed in [seconds, micrometers and pico-Newtons](doc/sim/units.md).

	set simul system
	{
	    time_step = 0.005
	    viscosity = 0.02
	}
	
	set space cell
	{
	    shape = sphere
	}
	
	set fiber microtubule
	{
	    rigidity = 20
	    segmentation = 0.5
	    confine = inside, 200, cell
	}
	    
	new cell
	{
	    radius = 5
	}
	
	new 5 microtubule
	{
	    length = 11
	}
	
	run 5000 system
	{
	    nb_frames = 10
	}

# Documentation

[Link to documentation](doc/index.md)

The Brownian dynamics approach was described in:  
[Collective Langevin Dynamics of Flexible Cytoskeletal Fibers](http://iopscience.iop.org/article/10.1088/1367-2630/9/11/427/meta)

A recent protocol demonstrating the use of Cytosim is:
[A Typical Workflow to Simulate Cytoskeletal Systems](https://dx.doi.org/10.3791/64125)

The documentation source files use [Markdown](https://en.wikipedia.org/wiki/Markdown) and are best viewed with [MacDown (Mac OSX only)](https://macdown.uranusjr.com) or [Typora (Cross platform)](https://typora.io) 

# Installation

Cytosim is distributed as source code and [must be compiled](doc/compile/index.md) before use. On Mac OS X and Linux this should be uncomplicated using a standard compiler. Compiling natively on Windows would require changes to the code. For Windows 10 and later, use the [Windows Subsystem for Linux](doc/compile/wsl.md). For older Windows systems, [use Cygwin](doc/compile/cygwin.md).

To download the source code, enter these commands in a terminal window:

	git clone https://gitlab.com/f.nedelec/cytosim
	cd cytosim
	
To compile using [make](https://www.gnu.org/software/make), try:
	
	make

If this fails, parameters of `makefile.inc` need to be updated.
Altermatively, it is possible to use [cmake](https://cmake.org) to configure `make` automatically:

	mkdir build
	cd build
	cmake ..
	make

For troubleshooting, please check [the compile instructions](doc/compile/index.md).
Once *cytosim* is running on your machine, proceed with the [tutorials](doc/tutorials/index.md), the page on [running simulations](doc/main/runs.md), and the examples contained in the folder `cym`. You can run them by specifying the paths:

	bin/play cym/self.cym

# Contributors

 The project was started in 1995, and received its name in 1999.
 We hope cytosim can be useful for your research. 
 Sincerely yours, The Developers of Cytosim:

*  Francois J. Nedelec     1995-
*  Dietrich Foethke        2003-2007
*  Cleopatra Kozlowski     2003-2007
*  Elizabeth Loughlin      2006-2010
*  Ludovic Brun            2008-2010
*  Beat Rupp               2008-2011
*  Jonathan Ward           2008-2014
*  Antonio Politi          2010-2012
*  Andre-Claude Clapson    2011-2013
*  Jamie-Li Rickman        2014-2019
*  Gaelle Letort           2014-2022
*  Maud Formanek           2020-2023
*  Manuel Lera-Ramirez     2017-2024
*  Serge Dmitrieff         2013-
*  Julio Belmonte          2014-

# Contact

Email: cytosim@cytosim.org


