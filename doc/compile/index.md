# How to compile Cytosim

The core of Cytosim is written using standard [C++17](https://en.wikipedia.org/wiki/C%2B%2B17) and it is necessary to recompile the programs after each modification of the source code. Some code is platform specific, in particular optimizations relying on [SIMD](https://en.wikipedia.org/wiki/Single_instruction,_multiple_data), but this code is automatically disabled if the feature is not available. Most of the accessory tools use [Python](https://www.python.org).

Compilation requires a C++ compiler: e.g. [The GNU compiler `g++`](http://gcc.gnu.org/), [`clang`](http://clang.llvm.org) or the [`Intel compiler`](http://en.wikipedia.org/wiki/Intel_C%2B%2B_Compiler), together with a few libraries.
Compilation is started from a terminal, with a program called [`make`](http://www.gnu.org/software/make/). Optionally, [`cmake`](https://cmake.org) can be used to configure `make` on your platform. On MacOS, we recommend using [`Xcode`](https://en.wikipedia.org/wiki/Xcode), which is used for development.

### Dimensionality

The dimensionality is set during compilation. It can be querried by running `sim info`.
It is necessary to recompile to change dimensionality, [following these instructions](dimensionality.md).

### Mathematical libraries

Both `sim` or `play` use these mathematical libraries:

* [BLAS](http://netlib.org/blas)
* [LAPACK](http://netlib.org/lapack)

There exists a [public reference implementation](http://netlib.org), which can be compiled with a [FORTRAN](http://en.wikipedia.org/wiki/Fortran) compiler.

However, the precompiled library, is available on many platforms:

- [Intel Math Kernel Library](https://software.intel.com/mkl)
- [Apple's vecLib](http://developer.apple.com/hardwaredrivers/ve/vector_libraries.html)
- [OpenBLAS](https://www.openblas.net)
- Also available for many Linux distributions.

Apple's veclib is preinstalled on MacOS, and [Intel's MKL](https://en.wikipedia.org/wiki/Math_Kernel_Library) is available free of charge, but is optimal only on Intel processors.

### Interactive player
 
Cytosim's ***play*** is built on:

- [OpenGL](http://www.opengl.org/) 
- [POSIX threads](http://en.wikipedia.org/wiki/POSIX_Threads)
- [libspng](https://libspng.org/) to export images
- [GLUT](http://www.opengl.org/resources/libraries/glut/) for windowing.

Note that GLUT can be replaced by [freeGLUT](http://freeglut.sourceforge.net/).  
GLUT and OpenGL are included in MacOS.

Cytosim's `play` can export PNG images, because its code include a copy of [libspng](https://libspng.org/), but it can be linked against [libPNG]() instead. See the [compile options](options.md).

### Multi player
 
Cytosim's ***monoplay*** and ***multiplay*** are built on:

- [OpenGL](http://www.opengl.org/) 
- [POSIX threads](http://en.wikipedia.org/wiki/POSIX_Threads)
- [libspng](https://libspng.org/) to export images
- [GLFW](https://www.glfw.org/) for windowing.

GLFW must be installed.
One option is to [compile](https://www.glfw.org/docs/latest/compile_guide.html) from the [source code](https://www.glfw.org/download.html), which requires [cmake](https://cmake.org/).

# Getting Ready 

### MacOS

Install [Xcode](https://developer.apple.com/technologies/tools/), available on the Mac App Store. After installing Xcode, install the Xcode 'Command-Line Tools', an optional package providing 'make'. All necessary libraries are then installed.

We provide the Xcode project file for cytosim, which is a convenient way to access the code.


### Linux

Check the [dedicated pages](linux.md).

On Linux, you need to install the GNU compiler collection, BLAS and LAPACK. This is sufficient to compile `sim`.
Recent Linux distributions provide precompiled BLAS/LAPACK as optional installation (check your distribution).

To make `play` install the OpenGL developer libraries and FreeGLUT or OpenGLUT. 


### Windows

Native compilation on Windows can be complicated. We thus recommend installing a Linux-like environment on top of Windows, unless a dual boot system is already available. 

For Windows 10 and later, use the [Windows Subsystem for Linux](wsl.md).  
For older systems, [use Cygwin](cygwin.md).

You will need a compiler, the X window system, BLAS/LAPACK and GLUT.  


# Compilation

After installing a compiler and [gnu's make ](http://www.gnu.org/software/make/)
you are ready to compile from a terminal, with the following commands in the root directory of cytosim:

	make sim
	make play

The command `make` without arguments will build `sim` and `play`.  
If this does not work, you may need to manually edit `makefile.in` to adjust to your platform.

It is also possible to use [cmake](https://cmake.org), which finds the configuration automatically:

	mkdir build
	cd build
	cmake ..
	make sim play report

You can then check the resulting executables, normally located in subdirectory `bin`:

	bin/sim info
	bin/sim
	bin/play live

# Optimizations

To speed up the calculation, follow these steps:

On line 10 or so of  `makefile.inc`, set `MODE` to `F`:

	MODE := F

Turn off assertions by defining `NDEBUG` in `src/base/assert.h`:

	#define NDEBUG

Recompile Cytosim from scratch:

	make clean
	make

# Troubleshooting `sim`

Compilation is specified in `makefile` and `makefile.inc`, and these files may need to be adjusted manually.
Check  `makefile.inc` first and verify the [compile options](options.md).

Make attempts to automatically detect the platform:

	#---------------- MACHINE = {mac, linux, cygwin, auto}
	MACHINE := auto

But you may manually set `MACHINE` to `mac`, `linux` or `cygwin` depending on your platform,
and check the parameters set lower in the `makefile.inc`, for this platform, for example:

	ifeq ($(MACHINE),linux)
		...
	endif

To pinpoint the problem, try to build objects that have fewer depencies first:

#### Check your compiler for C++17 support and compilation switches

	make test_cxx

#### Check your BLAS/LAPACK directly:

	make test_blas
	
If you get an error at this stage, check that the compiler is linking the correct libraries.
To find the libraries on your system, this command may help:

	find /usr/lib -name libblas.*

For example, if the result is:

	/usr/lib/libblas.so

You should adjust `makefile.inc` to specify the corresponding path:

	LIBDIR := /usr/lib

Ensure that no other line changes the value of `LIBDIR`. For example, this one is commented out with the `#`:

	#LIBDIR := /usr/lib/x86_64-linux-gnu

#### Check if it can read a configuration file:

	make test_glossary

#### At this stage, you can attempt to compile `sim`:

	make sim

# Troubleshooting `play`

If you are having trouble compiling `play`, check its requirements independently:

#### Check for thread support:

	make test_thread
	
#### Check for OpenGL support:

	make test_opengl
	
#### Check for GLUT support:

	make test_glut

#### Check for GLAP (our own extension of GLUT)):

	make test_glapp

#### At this stage, you can compile `play`:

	make play

### Finally: contact us!

Please, write to `feedbackATcytosimDOTorg`.

Please, describe what fails and what you have tried.
Attach your 'makefile.inc' and tell us the platform on which you compiled.

FJN, 4.6.2022
