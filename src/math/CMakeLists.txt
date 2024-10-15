set(SOURCES_MATH
    vector1.cc
    vector2.cc
    vector3.cc
    vector4.cc
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
    random_seed.cc
    random_pcg.cc
    modulo.cc
    tesselator.cc
    random_vector.cc
    project_ellipse.cc
    SFMT.c
)

list(TRANSFORM SOURCES_MATH PREPEND "${PROJECT_SOURCE_DIR}/src/math/")

set(MATH_INCLUDES math base)
list(TRANSFORM MATH_INCLUDES PREPEND "${PROJECT_SOURCE_DIR}/src/")

add_library(${MATH_LIBRARY} STATIC ${SOURCES_MATH})
target_include_directories(${MATH_LIBRARY} PUBLIC "${MATH_INCLUDES}" )


if(MAKE_TESTS)

    set(RASTERIZER_INCLUDES math base gym)
    list(TRANSFORM RASTERIZER_INCLUDES PREPEND "${PROJECT_SOURCE_DIR}/src/")

    add_library(rasterizerGL "${PROJECT_SOURCE_DIR}/src/math/rasterizer.cc")
    target_compile_definitions(rasterizerGL PRIVATE -DDISPLAY)
    target_include_directories(rasterizerGL PUBLIC "${RASTERIZER_INCLUDES}" )
endif()
