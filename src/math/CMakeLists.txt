set(MATH_SOURCES
    "${SRC_MATH_DIR}/vector1.cc"
    "${SRC_MATH_DIR}/vector2.cc"
    "${SRC_MATH_DIR}/vector3.cc"
    "${SRC_MATH_DIR}/matrix11.cc"
    "${SRC_MATH_DIR}/matrix22.cc"
    "${SRC_MATH_DIR}/matrix33.cc"
    "${SRC_MATH_DIR}/matrix34.cc"
    "${SRC_MATH_DIR}/rasterizer.cc"
    "${SRC_MATH_DIR}/sparmatsym1.cc"
    "${SRC_MATH_DIR}/sparmatsym2.cc"
    "${SRC_MATH_DIR}/sparmatsymblk.cc"
    "${SRC_MATH_DIR}/sparmatsymblkdiag.cc"
    "${SRC_MATH_DIR}/sparmatblk.cc"
    "${SRC_MATH_DIR}/matfull.cc"
    "${SRC_MATH_DIR}/polygon.cc"
    "${SRC_MATH_DIR}/spherical_code.cc"
    "${SRC_MATH_DIR}/random.cc"
    "${SRC_MATH_DIR}/tesselator.cc"
    "${SRC_MATH_DIR}/random_vector.cc"
    "${SRC_MATH_DIR}/project_ellipse.cc"
    "${SRC_MATH_DIR}/quartic_solver.cc"
    "${SRC_MATH_DIR}/SFMT.c"
)

add_library(${MATH_LIB_TARGET} STATIC ${MATH_SOURCES})
target_include_directories(${MATH_LIB_TARGET} PUBLIC
    "${SRC_MATH_DIR}"
    "${SRC_BASE_DIR}"
)

if(BUILD_CYTOSIM_TESTS)
    add_library(rasterizerGL "${SRC_MATH_DIR}/rasterizer.cc")
    target_compile_definitions(rasterizerGL PRIVATE -DDISPLAY)
    target_include_directories(rasterizerGL PUBLIC
        "${SRC_MATH_DIR}"
        "${SRC_BASE_DIR}"
        "${SRC_DISP_DIR}"
    )
endif()
