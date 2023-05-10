
# The Math Kernel Library (MKL)


The [Math Kernel Library](https://en.wikipedia.org/wiki/Math_Kernel_Library) is [distributed by Intel](https://software.intel.com/mkl) and free to use. It provides a version of BLAS and LAPACK that can be used with Cytosim. It was optimized for Intel processor and may give better performance on certain machines, compared to other distributions of BLAS.

This feature was implemented on 10/04/2018.


# Compilation


The code does not need to be changed. In ‘makefile.inc’ set 
    
    HAS_MKL := 2

In case of problem, check that the options of the makefile define `sequential static linking`.

Linking MKL requires multiple pass, due to cross-dependency in the library.
You must specify a 'group' of libraries:

    MKLLIB := -Wl,--start-group $(MKLDIR)/libmkl_intel_lp64.a $(MKLDIR)/libmkl_sequential.a $(MKLDIR)/libmkl_core.a -Wl,--end-group

Make sure there is no other MKLLIB defined in the makefile that would overwrite this line.

For advice on linking, check [Intel's Link Line Advisor](https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor).


# MKL on a cluster

To compile on a cluster, you may need to load the module containing the Intel MKL:

    module load imkl

Note that `MKLROOT` should be defined by the OS, when you load module `imkl`.  
The latest compiler may be required:
 
    module load foss

The rest is as usual.


On the cluster, the `imkl` module should not be needed if sim was compiled with **static linking** as specifed with `HAS_MKL := 2`.


# Verification


Check the library dependencies of cytosim with `ldd`:

    ldd bin/sim

The MacOS equivalent is
	
	otool -L bin/sim


FJN 27/04/2018
