# Cytosim's Documentation

*work in progress*

*  [**Overview**](main/overview.md)
*  [Getting started (introduction to the command line)](main/starter.md)
*  [The configuration file](sim/config.md)
*  [Simulation engine](sim/index.md)
*  [**Tutorials**](tutorials/index.md)
*  [The different executables](main/executables.md)
*  [File types](main/file_types.md)
*  [Running simulations on your computer](main/runs.md)
*  [Standard config files](main/examples.md)
*  [Graphical rendering](sim/graphics.md)
*  [Making movies](main/movies.md)
*  [Getting numbers out of Cytosim](main/report.md)
*  [Running cytosim on a cluster](main/run_slurm.md)
*  [Frequently asked questions](main/faq.md)
*  [Prior work](examples/index.md)

# Installation

Cytosim is distributed as source code and [must be compiled](compile/index.md) before use. On Mac OS X and Linux this should be straightforward if you are familiar with compilation in general. On Windows, we suggest to [compile within Cygwin](compile/cygwin.md).

To compile, enter these commands in a terminal window:

	git clone https://github.com/nedelec/cytosim.git
	cd cytosim
	make

Once cytosim is running on your machine, check the tutorials, the page on [running simulations](main/runs.md), and the examples contained in the folder `cym`. Inspect in particular the short configuration files (e.g. fiber.cym, self.cym). 

#### Troubleshooting

For more information, please check [the dedicated pages](compile/index.md).  
You may need to manually edit the makefiles depending on your platform.

# Advanced matter

*  [Code documentation](code/index.md)
*  [Doxygen documentation](doc/code/doxygen/index.html)

# Contact

cytosimATcytosimDOTorg

