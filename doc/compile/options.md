# Compiler and Optimizations
 
The parameters of the compilation are set by editing the file `makefile.inc`, located in the root of the distribution. By manually editing this file, you can change:

The compiler:

	COMPILER := gcc
 
The level of optimization:

	MODE := F
 
Use the `D=debug` mode to test new code, and the 'F=fast' mode to run extensive calculations. For optimal performance, you should also disable assertions (see below).

# Assertions

Assertions (a safety mechanism used for debugging) are turned on/off in `src/base/assert_macro.h`
To make the executable faster, you can disable assertions by defining NDEBUG:

	#define NDEBUG 1


# Floating-point precision

The code is written using an alias `real` to either `float` or `double`:

- a [float](https://en.wikipedia.org/wiki/Single-precision_floating-point_format) uses 4 bytes and has 7 decimals of precision
- a [double](https://en.wikipedia.org/wiki/Double-precision_floating-point_format) uses 8 bytes and has 17 decimals of precision

The selection between float or double is done in `src/math/real.h`:
 
	#define REAL_IS_DOUBLE 1

**Using double precision is strongly advised.** 

### Using single-precision

The code might be faster in single precision, as the memory footprint is reduced by a factor 2. However, the 'solve' step **may fail** because of the reduced precision.

To use single precision, change to:

	#define REAL_IS_DOUBLE 0

After editing, recompile everything:

	make clean
	make

Always compare the results with double precision as a benchmark.

# OpenGL Window Toolkits

A cross-platform toolkit is used by all graphical programs (`play`, `multiplay`), to create a window on which to draw with OpenGL.  
Cytosim is migrating from GLUT/FreeGLUT to GLFW:

- [GLUT](http://www.opengl.org/resources/libraries/glut/)
- [GLFW](https://www.glfw.org/) for windowing.

### GLUT/FreeGLUT

GLUT and OpenGL are included in MacOS, adjust linkage to include:

	-framework GLUT -framework OpenGL 

On other platforms, install [freeGLUT](http://freeglut.sourceforge.net/).  

### GLFW

GLFW must be installed.  
One option is to [compile](https://www.glfw.org/docs/latest/compile_guide.html) from the [source code](https://www.glfw.org/download.html), which requires [cmake](https://cmake.org/).

The GLFW library and header file should be copied in folder 'ext' in the base directory:

	- ext/glfw3.a
	- ext/GLFW/glfw3.h

# PNG image support

Cytosim's source includes a copy of [libspng](https://libspng.org/), which is used to export PNG images.
To use the PNG library [http://www.libpng.org](libpng), you must install it. 
On MacOS, use [Homebrew](https://brew.sh):

	brew install libpng
 
You can then enable PNG support by editing the `makefile.inc`:

	HAS_PNG := 2

# Advanced features

* [Math Kernel Library](math_kernel_library.md)
* [SIMD optimizations](vectorization.md)
* [Parallel execution](multithreading.md)

