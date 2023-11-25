
# Installation


Cytosim is distributed as source code and [must be compiled](doc/compile/index.md) before use. On MacOS and Linux this should be uncomplicated even if you are not familiar with program development. Compiling natively on Windows would require changes to the code, but Cytosim should [compile within Cygwin](doc/compile/cygwin.md).

To download the source code, enter these commands in a terminal window:

	git clone https://gitlab.com/f-nedelec/cytosim
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

# Dimensionality

The dimensionality is baked in at compilation time. Executables such as `sim` and `play` are thus made exclusively for 1D, 2D or 3D simulations.
