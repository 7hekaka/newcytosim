
# Dimensionality

Cytosim can perform simulations in 1D, 2D or 3D, but the dimensionality is selected at compilation time!

This means that one executable will only be able to perform simulations for one dimensionality.

### As specified in the code

The dimension by defalt is selected in `src/math/dim.h`

	#define DIM 2

If you edit this file and change `DIM`, always recompile everything:

	make clean
	make

This will replace the executables in `bin/`

You can check the dimensionality of the executable, with:

	bin/sim info
	bin/play info

### Selected at compilations

The value defined in `dim.h` can be overwritten by the compilation system (e.g. `make`).
Shortcuts are built into the makefile, and you can built directly the 2D executables:

	make bin2/sim
	make bin2/play

and similarly

	make bin3/sim
	make bin3/play

or to make all exectuables of specified dimension:

	make dim2
	make dim3

and to all everything:

	make alldim

### Recommended practice

You can compile the executables that you need and rename them accordingly:

	cp bin2/play ~/bin/play2
	cp bin3/play ~/bin/play3

It can be handy to create shortcuts or alias if you use these often.