# Cytosim Code Documentation

Cytosim is built around [a core C++ engine](../../src/index.md) constituted of the files in [`src/`](../../src) and subdirectories.  

Handy [Python scripts](../../python/index.md) are located in [`python/`](../../python).

# Code structure

- [Modularity](modularity.md)


# Code documentation

The C++ code contains documentation embedded in the comments.

This documentation can be extracted automatically with [`doxygen`](http://www.doxygen.nl) to generate HTML pages describing the C++ classes and their methods. After installing 'doxygen', from a terminal, in the project root directory, enter:

	make doc

The output will be accessible [from the index](doxygen/index.html) located in [`doc/code/doxygen`](doxygen).


# Testing Cytosim

- [How to validate the executable](validation.md)

### Contact

feedbackATcytosimDOTorg


