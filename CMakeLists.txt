cmake_minimum_required(VERSION 3.2.0)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)
set(CMAKE_BUILD_TYPE Release)

project("cytosim" C CXX)
#set(CMAKE_VERBOSE_MAKEFILE ON)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/tar)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

set(SIM_LIB_TARGET "libcytosim")
set(SIMG_LIB_TARGET "libcytosimG")
set(BASE_LIB_TARGET "libcytobase")
set(MATH_LIB_TARGET "libcytomath")
set(DISP_LIB_TARGET "libcytodisp")

if(DEFINED DIMENSION)
    add_compile_definitions(-DDIM=${DIMENSION})
    set(SIM_LIB_TARGET "cytolibD${DIMENSION}")
    set(SIMG_LIB_TARGET "cytospaceD${DIMENSION}G")
endif()

option(MAKE_SIM "build the simulation executable" ON)
option(MAKE_PLAY "build the graphical OpenGL viewer" ON)
option(MAKE_TOOLS "build all cytosim/tools executables" ON)
option(MAKE_TESTS "build all cytosim/test executables" OFF)


set(SRC_BASE_DIR "${PROJECT_SOURCE_DIR}/src/base")
set(SRC_MATH_DIR "${PROJECT_SOURCE_DIR}/src/math")
set(SRC_SIM_DIR "${PROJECT_SOURCE_DIR}/src/sim")
set(SRC_DISP_DIR "${PROJECT_SOURCE_DIR}/src/disp")
set(SRC_PLAY_DIR "${PROJECT_SOURCE_DIR}/src/play")
set(SRC_TEST_DIR "${PROJECT_SOURCE_DIR}/src/test/")
set(SRC_TOOLS_DIR "${PROJECT_SOURCE_DIR}/src/tools")
set(DOC_DIR "${PROJECT_SOURCE_DIR}/src/doc")

set(SRC_SIM_ALL_DIR
    "${SRC_SIM_DIR}"
    "${SRC_SIM_DIR}/fibers"
    "${SRC_SIM_DIR}/hands"
    "${SRC_SIM_DIR}/singles"
    "${SRC_SIM_DIR}/couples"
    "${SRC_SIM_DIR}/organizers"
    "${SRC_SIM_DIR}/spaces"
)

set(SIM_TARGET "sim")
set(PLAY_TARGET "play")

#-------------------------------- Libraries ------------------------------------

if(NOT UNIX)
    message(FATAL_ERROR "CMake File currently only setup to build on Linux and Mac")
elseif(APPLE)
    message(">>>>>> Building for APPLE")
else()
    message(">>>>>> Building for Linux")
    link_libraries("gfortran")
endif()

# static linking preferred here:
find_library(BLAS_LIB libblas.a blas REQUIRED)
if(BLAS_LIB)
    message(">>>>>> BLAS: ${BLAS_LIB}")
else()
    message(FATAL_ERROR "BLAS library not found")
endif()

find_library(LAPACK_LIB liblapack.a lapack REQUIRED)
if(LAPACK_LIB)
    message(">>>>>> LAPACK: ${LAPACK_LIB}")
else()
    message(FATAL_ERROR "LAPACK library not found")
endif()

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

#-------------------------- Graphical Libraries --------------------------------

if(MAKE_PLAY)
    if(UNIX AND NOT APPLE)
        find_package(OpenGL)
        find_package(GLEW)
        find_package(GLUT)
        find_package(X11)
        find_package(Xt)
        set(OPENGL_LIBS
            OpenGL::GL
            GLEW
            glut
            X11
            Xt
        )
    elseif(UNIX AND APPLE)
        set(OPENGL_LIBS
            "-framework Accelerate"
            "-framework GLUT"
            "-framework OpenGL"
            "-framework AGL"
        )
    endif()
endif()


if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(FAST_FLAGS -Wno-deprecated-declarations -ffast-math -funroll-loops)
    set(DEBUG_FLAGS -g3 -ggdb)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
    set(FAST_FLAGS -Wno-deprecated-declarations -ffast-math -funroll-loops)
    set(DEBUG_FLAGS -g3 -ggdb)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(RUN_FLAGS -O1 -fno-tree-vectorize -ffast-math -funroll-loops -march=native)
    set(FAST_FLAGS -O2 -fno-tree-vectorize -ffast-math -funroll-loops -march=native)
    set(DEBUG_FLAGS -g3 -ggdb -Wno-deprecated)
    set(COVER_FLAGS -fprofile-arcs -ftest-coverage)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    add_compile_options(-no-parallel -mkl=sequential -wd1224,161)
    set(FAST_FLAGS -O2 -mdynamic-no-pic -funroll-loops -fno-math-errno -mavx)
    set(DEBUG_FLAGS -g -Wall -no-pie -wd279,383,810,869,981,1418,1419,1572,2259)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    message(FATAL_ERROR "Microsoft Visual Studio not supported for compilation")
else()
    message(FATAL_ERROR "Unrecognized compiler: ${CMAKE_CXX_COMPILER_ID}")
endif()


#file(COPY ${PROJECT_SOURCE_DIR}/cym DESTINATION ${PROJECT_BINARY_DIR})

add_compile_options(-Wno-deprecated-declarations)
add_compile_options(-finline-functions)
add_compile_options(${FAST_FLAGS})
add_subdirectory("src")
