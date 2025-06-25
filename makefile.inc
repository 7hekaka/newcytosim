# This file defines variables for make and is included by `makefile`
# Copyright 2007-2017 EMBL; 2019-2022 Sainsbury Laboratory, Cambridge University


#---------------- MACHINE = {mac, linux, cluster, cygwin, auto}
MACHINE := auto

#---------------- COMPILER = {gcc, icc, clang}
COMPILER := gcc

#---------------- MODE = {F = Fast; D = Debug; P = Profiling; C = Coverage}
MODE := F

#---------------- INTEL Math Kernel Library
# If the Intel Math Kernel library is a free library providing BLAS+LAPACK,
# If it is installed on your system, you can set:
#     HAS_MKL := 1 for dynamic linking
#     HAS_MKL := 2 for static linking
# otherwise, set HAS_MKL := 0;

HAS_MKL := 0

#---------------- PNG image export support
# `libpng` needs to be installed to save PNG images:
#     Mac OSX:
#        www.macports.org :  port install libpng 
#        Homebrew  (http://brew.sh)
#             brew install libpng
#             brew tap homebrew/dupes
#             brew install zlib
#     CentOS:  yum install libpng-devel
#     Ubuntu:  `apt-get` or `synaptic libpng-dev`
# HAS_PNG can be:
#     0 : no PNG support
#     1 : macports installation
#     2 : Homebrew installation

HAS_PNG := 0

#-------------------------------------------------------------------------------
#---------------------------  Platform Detection  ------------------------------
#-------------------------------------------------------------------------------

ifeq ($(MACHINE), auto)
   ifeq ($(shell uname), Darwin)
      MACHINE := mac
   else ifeq (login-e, $(findstring login-e, $(shell uname -n)))
      MACHINE := cluster
   else ifeq ($(shell uname), Linux)
      UNAME := $(firstword $(subst -, ,$(shell uname -n)))
      ifeq ($(UNAME), login)
          MACHINE := cluster
      else
          MACHINE := linux
      endif
   else ifeq ($(shell uname -o), Cygwin)
      MACHINE := cygwin
   else
      $(error Unknown platform)
   endif
endif

#-------------------------------------------------------------------------------
#----------------------------  Intels' OpenMP  ---------------------------------
#-------------------------------------------------------------------------------

# enable or not OpenMP (http://openmp.org):
OpenMP := -fopenmp
OpenMP :=
OMPLIB :=

# redirect OpenMP on Mac OS with Apple silicon:
ifneq ($(OpenMP),)
ifeq ($(MACHINE),mac)
ifeq ($(shell uname -m), arm64)
    # using libomp installed via Homebrew (http://brew.sh)
    OpenMP := -Xclang -fopenmp -I/opt/homebrew/include/
    OMPLIB := /opt/homebrew/lib/libomp.a
endif
endif
endif

#-------------------------------------------------------------------------------
#------------------------  Compiler Specifications  ----------------------------
#-------------------------------------------------------------------------------

# Libraries needed by almost all objects
LIBRARIES := -lpthread -lncurses

# common options for all C++ objects: use '-fopenmp' to enable OpenMP
CXXFLG := -std=c++17 -fno-rtti -Wredeclared-class-member -Wno-deprecated-declarations
CXXFLG := -std=c++17 -fno-rtti $(OpenMP) #-g -fno-inline
CXXFLG := -std=c++17 -fno-rtti

# option concerning memory alignment:
ALIGN  := -malign-double -falign-loops -falign-jumps -falign-functions
ALIGN  :=

# options concerning warnings:
WARN   := -Wundef -Wall -Wno-unknown-pragmas #-Wconversion -Wno-unused

ifeq ($(MACHINE),mac)
    LIBTOOL := /usr/bin/libtool -static -o
    CXXFLG := $(CXXFLG) -fno-pic
else
    LIBTOOL := libtool --mode=link gcc -g -O -o
    LIBTOOL := ar -rcs
endif


#-------------------------------------------------------------------------------

ifeq ($(COMPILER),gcc)

# Intel SIMD: -march=core2 (SSE3); -march=sandybridge (AVX); -march=haswell (AVX2); -march=native
# Apple Silicon: -mcpu=apple-m1
# You may use g++-10 to support Intel OMP directives
# To use GCC on OSX, you need to specify the path to header files:
# export CPATH="/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/usr/include"

    CXX := g++ -Wfatal-errors -Wno-deprecated-declarations
    FAST := -O3 -fno-trapping-math -funroll-loops -ffp-contract=fast -march=native # -mavx2 -mfma
    FAST := -O3 -fno-trapping-math -funroll-loops -ffp-contract=fast # -mcpu=apple-m1

    DEBUG := -g3 -ggdb # -Wconversion
    COVER := -fprofile-arcs -ftest-coverage
    
else ifeq ($(COMPILER),icc)

    LIBTOOL  := ar -rcs

    #CXX := icpx -vec-report0
    #CXX := icpx -bind_at_load -fvisibility=hidden -wd1224,161
    CXX := icpx -vec-report0 -wd1224,161,1478
    CXX := icpx -fargument-noalias 
    # Intel advanced instruction sets:
    # '-xHost' to optimize for host machine
    # '-msse3' for SSE3 or '-xAVX' for AVX
    # '-march=corei7-avx' for AVX
    # '-march=core-avx2' for Intel core i7 (ca. 2015)
    FAST := -O2 -march=core-avx2
    FAST := -O2 -fno-devirtualize -march=core-avx2 -mfma -Wl,-no_pie -funroll-loops -fno-math-errno -fp-model fast=1 -fp-model no-except
    WARN := -qopt-report=2 -qopt-report-phase=vec# -opt-report-file=stderr
    WARN :=

    DEBUG := -march=core-avx2 -g -wd279,383,810,869,981,1418,1419,1572,2259
    COVER :=

else ifeq ($(COMPILER),clang)

# for the standard library, use: -stdlib=libc++
# the old library: -stdlib=libstdc++ (deprecated)

    CXX := c++ -Wfatal-errors
    FAST := -O2 -fno-devirtualize -funroll-loops -mavx2 -mfma
    FAST := -O2 -fno-devirtualize -funroll-loops
    DEBUG := -g3 -ggdb
    COVER := -fprofile-arcs -ftest-coverage

endif

#-------------------------------------------------------------------------------
#------------------------------- Mac OSX ---------------------------------------
#-------------------------------------------------------------------------------

ifeq ($(MACHINE),mac)

    MKLDIR := /opt/intel/mkl/lib

    FlagsC := $(COVER)
    FlagsD := $(DEBUG) $(WARN)
    FlagsP := -O2 -fno-devirtualize $(ALIGN) -fno-inline
    FlagsF := $(FAST) $(ALIGN)

    LINK   := -framework Accelerate -framework QuartzCore $(OMPLIB) $(LIBRARIES)

    GRAFIX := -framework GLUT -framework OpenGL -framework AGL
    # Using brew's freeGLUT:
    #GRAFIX := -L/usr/local/Cellar/freeglut/3.0.0 -lglut -framework OpenGL
    
    ifeq ($(HAS_PNG), 1)

        # macports libraries:
        LIB_PNG := /opt/local/lib/libpng.a /opt/local/lib/libz.a
        INC_PNG := -I/opt/local/include

    endif
	
    ifeq ($(HAS_PNG), 2)

        # HomeBrew libraries:
        LIB_PNG := /usr/local/lib/libpng.a /usr/local/opt/zlib/lib/libz.a
        INC_PNG := -I/usr/local/include

    endif

endif

#-------------------------------------------------------------------------------
#--------------------------------- Cluster -------------------------------------
#-------------------------------------------------------------------------------

ifeq ($(MACHINE),cluster)

    # MKLROOT should be defined by 'module load imkl'
    MKLDIR := $(MKLROOT)/lib/intel64

    # Beware of AVX2 and other advanced features on heterogeneous clusters:
    FAST := -O3 -fno-devirtualize -mavx2 -mfma -ffp-contract=fast
    # Options for Intel compiler ICPX:
    FAST := -O3 -mavx2 -mfma

    FlagsD := -O0 $(DEBUG) $(WARN)
    FlagsP := $(FAST) $(ALIGN) -O0 -fno-inline -pg
    FlagsC := -O0 -fprofile-arcs -ftest-coverage
    FlagsF := $(FAST) $(ALIGN) -finline-functions
    
    ### Automatic
    LINK := -llapack -lopenblas -lgfortran $(LIBRARIES)
    ### Static linking of libC and LibStdC++
    LINK := -static-libstdc++ -static-libgcc -llapack -lopenblas -lgfortran $(LIBRARIES)

    # linkage options for programs using graphics:
    GRAFIX :=-L/usr/lib -lglut -lGL -lGLEW -lXt -lX11

endif

#-------------------------------------------------------------------------------
#---------------------------------- Linux --------------------------------------
#-------------------------------------------------------------------------------

ifeq ($(MACHINE),linux)

    ### Linux options
    FlagsD := -O0 $(DEBUG) $(WARN)
    FlagsP := -O0 -fno-inline -pg $(ALIGN) $(FAST)
    FlagsC := -O0 -fprofile-arcs -ftest-coverage
    FlagsF := $(FAST) -fno-devirtualize -finline-functions $(ALIGN) -march=native

    # specify the library path:
    LIBDIR := /usr/lib
    LIBDIR := /usr/lib/x86_64-linux-gnu

    ############## BLAS/LAPACK Libraries ###############
    ### For 64bit Centos 5.8 March 2012:
    # Centos  blas: yum install blas-devel
    # Centos  lapack: yum install lapack-devel
    ### 64bit Linux distribution Ubuntu 12.04 LTS, 2012
    # Ubuntu   blas : apt-get install libblas.dev
    # Ubuntu lapack : apt-get install liblapack.dev

    # trying to automatically detect BLAS library:
    BLASLIB := $(wildcard /usr/lib/libblas.*) $(wildcard /usr/lib/*/libblas.*)
    # extracting path for libraries:
    BLASDIR := $(word 1, $(dir $(BLASLIB)))
    ifndef BLASDIR
       $(error Error: the path for the BLAS/LAPACK library `BLASDIR` is undefined)
    else
       $(warning Using BLAS library in $(BLASDIR))
    endif

    #### For dynamic linking:
    LINK := -L$(BLASDIR) -llapack -lblas $(LIBRARIES)
    #### For static linking specify the objects directly:
    #LINK := -L$(BLASDIR) $(BLASDIR)/liblapack.so $(BLASDIR)/libblas.so $(LIBRARIES)

    ########### Install libraries needed for graphics ############
    # Centos freglut: yum install freeglut-devel
    # Centos libXi:   yum install libXi-devel
    # Centos libXmu:  yum install libXmu-devel

    # Ubuntu freglut: apt-get install freeglut3-dev
    # Ubuntu libXi:   apt-get install libXi-dev
    # Ubuntu libXmu:  apt-get install libXmu-dev

    ### linkage options for programs using graphics:
    #   check for correct library libname.so otherwise give full path
    GRAFIX := -L$(LIBDIR) -lglut -lGL -lGLU -lGLEW -lXt -lX11

    ### Specify here the PNG library
    LIB_PNG := $(LIBDIR)/libpng.a $(LIBDIR)/libz.a
    INC_PNG :=

endif


#-------------------------------------------------------------------------------
#------------------------------- Cygwin ----------------------------------------
#-------------------------------------------------------------------------------

ifeq ($(MACHINE),cygwin)

    ### Cygwin October 2016
    FlagsD := -O0 $(DEBUG) $(WARN)
    FlagsP := -O0 -fno-inline -pg $(ALIGN) $(FAST)
    FlagsC := -O0 -fprofile-arcs -ftest-coverage
    FlagsF := -O2 -fno-devirtualize -finline-functions $(ALIGN) $(FAST)

    ### Needs
    # blas:   yum install blas-devel
    # lapack: yum install lapack-devel
    LINK := -L/usr/lib -L/bin -llapack -lblas $(LIBRARIES)

    ### linkage options for programs using graphics:
    ### Needs
    ### freglut: yum install freeglut-devel
    ### libXi:   yum install libXi-devel
    ### libXmu:  yum install libXmu-devel
    ###  check for correct library yourlibname.so otherwise give full path
    GRAFIX := -L/usr/lib -lglut -lGL -lGLU -lXt -lX11 -lGLEW

    ### Specify here the PNG library
    LIB_PNG := /lib/libpng.dll.a /lib/libz.a
    INC_PNG :=

    BINEXT :=.exe

endif


#-------------------------------------------------------------------------------
#------------------------------- Native Windows --------------------------------
#-------------------------------------------------------------------------------
# THIS IS NOT MAINTAINED: it worked in 2008 but has not been tried since

ifeq ($(MACHINE),windows)

    # Use a windows-native blas/lapack compilation:
    # copying "libblas.a" and "liblapack.a" in "/usr/lib"

    FlagsD := -g  $(WARN)
    FlagsP := -pg -O2 -march=i686 $(ALIGN)
    FlagsF := -O2 -finline-functions -march=i686 $(ALIGN)
  
    LINK   := -llapack -lblas
    GRAFIX := -lglut32 -lglu32 -lopengl32
    BINEXT :=.exe

endif



#-------------------------------------------------------------------------------
#------------------ linking with Intel Math Kernel Library ---------------------
#-------------------------------------------------------------------------------

# The Intel MKL provides BLAS/LAPACK libraries optimized for Intel's CPU
# and it is free: https://software.intel.com/en-us/mkl
# For help with the linking, check:
# https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor

ifeq ($(HAS_MKL),1)

    # sequential dynamic linking:
    MKLLIB := -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
    
    # threaded dynamic linking:
    #MKLLIB := -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5

    # modify the linking command:
    LINK := -L$(MKLDIR) $(MKLLIB) $(LIBRARIES)

endif


ifeq ($(HAS_MKL),2)

    # This are options for the compiler:
    MKLOPT := -m64 -I${MKLROOT}/include

    # sequential static linking (ICC and older GCC):
    #MKLLIB := $(MKLDIR)/libmkl_intel_lp64.a $(MKLDIR)/libmkl_sequential.a $(MKLDIR)/libmkl_core.a
    
    # sequential for linux
    MKLLIB := -Wl,--start-group $(MKLDIR)/libmkl_intel_lp64.a $(MKLDIR)/libmkl_sequential.a $(MKLDIR)/libmkl_core.a -Wl,--end-group  $(LIBRARIES) -lm -ldl

    # threaded static linking:
    #MKLLIB := $(MKLDIR)/libmkl_intel_lp64.a $(MKLDIR)/libmkl_thread.a $(MKLDIR)/libmkl_core.a /opt/intel/lib/libiomp5.a

    # modify the linking command:
    LINK := $(MKLLIB) $(LIBRARIES) -lm -ldl

endif


