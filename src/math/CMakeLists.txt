set(SOURCES_MATH
    vector1.cc
    vector2.cc
    vector3.cc
    matrix11.cc
    matrix22.cc
    matrix33.cc
    matrix34.cc
    rasterizer.cc
    sparmatsym1.cc
    sparmatsym2.cc
    sparmatsymblk.cc
    sparmatsymblkdiag.cc
    sparmatblk.cc
    matfull.cc
    polygon.cc
    spherical_code.cc
    random.cc
    tesselator.cc
    random_vector.cc
    project_ellipse.cc
    quartic_solver.cc
    SFMT.c
)

list(TRANSFORM SOURCES_MATH PREPEND "${SRC_MATH_DIR}/")

add_library(${MATH_LIB_TARGET} STATIC ${SOURCES_MATH})
target_include_directories(${MATH_LIB_TARGET} PUBLIC
    "${SRC_MATH_DIR}"
    "${SRC_BASE_DIR}"
)

if(BUILD_TESTS)
    add_library(rasterizerGL "${SRC_MATH_DIR}/rasterizer.cc")
    target_compile_definitions(rasterizerGL PRIVATE -DDISPLAY)
    target_include_directories(rasterizerGL PUBLIC
        "${SRC_MATH_DIR}"
        "${SRC_BASE_DIR}"
        "${SRC_DISP_DIR}"
    )
endif()
