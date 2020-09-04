# Cytosim Code Documentation

Cytosim is built around [a core C++ engine](../../src/index.md) constituted of the file in [`src/`](../../src).  
Handy [Python scripts](../../python/index.md) are located in [`python/`](../../python).

The C++ code contains documentation embedded in the comments.

A detailed HTML documentation of the C++ classes and methods can be generated with [`doxygen`](http://www.doxygen.nl).

After installing `doxygen` on your machine, in a terminal, and the project root directory, enter:

	make doc

The output is exported into [`doc/code/doxygen/`](doxygen) and accessible [from the index file](doxygen/index.html).


# Cytosim Testing

- [How to validate the executable](validation.md)

### Contact

feedbackATcytosimDOTorg


