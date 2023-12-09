
add_subdirectory("deps")

set(GYM_DIR "${PROJECT_SOURCE_DIR}/src/gym/")

set(SOURCES_GYM
    gle.cc
    gym_matrix.cc
    gym_cap.cc
    gym_view.cc
    gym_vect.cc
    gym_draw.cc
    gym_flute.cc
    gym_check.cc
    gym_menu.cc
    gym_zoo.cc
    gym_flat.cc
    gym_image.cc
    gym_color.cc
    gym_color_list.cc
    save_image.cc
    save_image_gl.cc
    offscreen.cc
    fg_font.cc
    fg_stroke.cc
    deps/spng.c
    deps/miniz.c
)

list(TRANSFORM SOURCES_GYM PREPEND "${GYM_DIR}")

set_property(SOURCE ${GYM_DIR}/deps/spng.c PROPERTY COMPILE_DEFINITIONS SPNG_STATIC SPNG_USE_MINIZ)
set_property(SOURCE ${GYM_DIR}/save_image.cc PROPERTY INCLUDE_DIRECTORIES "${GYM_DIR}/deps")
set_property(SOURCE ${GYM_DIR}/save_image_gl.cc PROPERTY INCLUDE_DIRECTORIES "${GYM_DIR}/deps")

set(GYM_INCLUDES base math)
list(TRANSFORM GYM_INCLUDES PREPEND "${PROJECT_SOURCE_DIR}/src/")

add_library(${GYM_LIBRARY} STATIC ${SOURCES_GYM} ${SOURCES_GYM_DEPS})
target_compile_definitions(${GYM_LIBRARY} PRIVATE -DDISPLAY)
target_include_directories(${GYM_LIBRARY} PUBLIC "${GYM_INCLUDES}")
