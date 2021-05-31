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
    add_definitions(-DDIM=${DIMENSION})
    set(SIM_LIB_TARGET "cytolibD${DIMENSION}")
    set(SIMG_LIB_TARGET "cytospaceD${DIMENSION}G")
endif()


option(BUILD_SIM "build the simulation executable" ON)
option(BUILD_PLAY "build the Graphical viewer" ON)
option(BUILD_REPORT "build the report executable" OFF)
option(BUILD_TOOLS "build all cytosim/tools executables" ON)
option(BUILD_TESTS "build all cytosim/test executables" ON)


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

find_library(LAPACK_LIB lapack REQUIRED)
if(LAPACK_LIB)
    message("Found LAPACK: ${LAPACK_LIB}")
endif()

if(NOT LAPACK_LIB)
    message(FATAL_ERROR "LAPACK library not found")
endif()

find_library(BLAS_LIB blas REQUIRED)
if(BLAS_LIB)
    message("Found BLAS: ${BLAS_LIB}")
endif()

if(NOT BLAS_LIB)
    message(FATAL_ERROR "BLAS library not found")
endif()

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

if(NOT UNIX)
    message(FATAL_ERROR "CMake File currently only setup to build on Linux and Mac")
endif()

if(BUILD_PLAY OR BUILD_TESTS OR BUILD_TOOLS)
    if(UNIX AND NOT APPLE)
        message("Building for Linux")
        find_package(GLUT)
        find_package(OpenGL)
        find_package(GLU)
        find_package(Xt)
        find_package(X11)
        find_package(GLEW)
        set(OPENGL_LIBS
            glut
            OpenGL::GL
            GLU
            Xt
            X11
            GLEW
        )
    elseif(UNIX AND APPLE)
        message("Building for MacOS")
        set(OPENGL_LIBS
            "-framework Accelerate"
            "-framework GLUT"
            "-framework OpenGL"
            "-framework AGL"
        )
    endif()
endif()


if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(CXX_FLAGS "-Wno-deprecated-declarations")
    set(FAST_FLAGS "-Wno-deprecated-declarations -ffast-math -funroll-loops")
    set(DEBUG_FLAGS "-g3 -ggdb")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
    set(CXX_FLAGS "-Wno-deprecated-declarations")
    set(FAST_FLAGS "-Wno-deprecated-declarations -ffast-math -funroll-loops")
    set(DEBUG_FLAGS "-g3 -ggdb")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CXX_FLAGS "-fprofile-arcs -ftest-coverage")
    set(FAST_FLAGS "-O3 -fno-tree-vectorize -ffast-math -funroll-loops -march=native -Wno-deprecated")
    set(DEBUG_FLAGS "-g3 -ggdb -Wno-deprecated")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(CXX_FLAGS "-no-parallel -mkl=sequential -wd1224,161")
    set(FAST_FLAGS "-O3 -mdynamic-no-pic -funroll-loops -fno-math-errno -mavx")
    set(DEBUG_FLAGS "-g -Wall -no-pie -wd279,383,810,869,981,1418,1419,1572,2259")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    message(FATAL_ERROR "Microsoft Visual Studio not supported for compilation")
else()
    message(FATAL_ERROR "Unrecognized compiler used: ${CMAKE_CXX_COMPILER_ID}")
endif()


file(COPY ${PROJECT_SOURCE_DIR}/cym DESTINATION ${PROJECT_BINARY_DIR})
add_definitions("-finline-functions ${FAST_FLAGS}")
add_subdirectory("src")
