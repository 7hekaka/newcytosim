set(TEST_INCLUDES
    "${SRC_BASE_DIR}"
    "${SRC_MATH_DIR}"
    "${SRC_SIM_ALL_DIR}"
    "${SRC_DISP_DIR}"
    "${SRC_PLAY_DIR}"
)

set(TEST_LIBS
    "${LAPACK_LIB}"
    "${BLAS_LIB}"
    Threads::Threads
    "${MATH_LIB_TARGET}"
    "${BASE_LIB_TARGET}"
)

set(TEST_LIST
    #"test" # Test is a reserved word in CMake 3.10.x
    "test_blas"
    "test_glossary"
    "test_random"
    "test_quaternion"
    "test_simd"
    "test_sizeof"
    "test_cxx"
    "test_math"
    "test_thread"
    "test_string"
)

foreach(TEST ${TEST_LIST})
    add_executable("${TEST}" "${SRC_TEST_DIR}/${TEST}.cc")
    set_source_files_properties(
        PROPERTIES
        EXTERNAL_OBJECT true
        GENERATED true
    )
    target_include_directories(${TEST} PUBLIC "${TEST_INCLUDES}")
    target_link_libraries(${TEST} PUBLIC "${TEST_LIBS}")
endforeach()

if(OPENGL_LIBS)

set(TEST_GL_LIBS
    rasterizerGL
    "${DISP_LIB_TARGET}"
    "${SIMG_LIB_TARGET}"
    "${MATH_LIB_TARGET}"
    "${BASE_LIB_TARGET}"
    "${LAPACK_LIB}"
    "${BLAS_LIB}"
    Threads::Threads
    "${OPENGL_LIBS}"
)

set(GL_TEST_LIST
    #"test_opengl" #Currently experiencing an issue with glext.h on Linux
    "test_glut"
    "test_glapp"
    "test_sphere"
    "test_rasterizer"
    "test_space"
    "test_grid"
    "test_solve"
)

foreach(TEST ${GL_TEST_LIST})
    add_executable(${TEST} "${SRC_TEST_DIR}/${TEST}.cc")
    target_include_directories(${TEST} PUBLIC "${TEST_INCLUDES}")
    target_link_libraries(${TEST} PUBLIC "${TEST_GL_LIBS}")
endforeach()

endif()
