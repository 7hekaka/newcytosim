
add_subdirectory("deps")

set(SOURCES_GYM
    gle.cc
    gym_matrix.cc
    gym_cap.cc
    gym_view.cc
    gym_vect.cc
    gym_draw.cc
    gym_flute.cc
    gym_flute_dim.cc
    gym_check.cc
    gym_menu.cc
    gym_zoo.cc
    gym_flat.cc
    gym_image.cc
    gle_color.cc
    gle_color_list.cc
    save_image.cc
    offscreen.cc
    fg_font.cc
    fg_stroke.cc
)

list(TRANSFORM SOURCES_GYM PREPEND "${PROJECT_SOURCE_DIR}/src/gym/")

set(GYM_INCLUDES base math)
list(TRANSFORM GYM_INCLUDES PREPEND "${PROJECT_SOURCE_DIR}/src/")

set(SOURCES_GYM_DEPS spng.c miniz.c)
list(TRANSFORM SOURCES_GYM_DEPS PREPEND "${PROJECT_SOURCE_DIR}/src/gym/deps/")
set_property(SOURCE spng.c PROPERTY COMPILE_DEFINITIONS SPNG_STATIC SPNG_USE_MINIZ)

add_library(${GYM_LIBRARY} STATIC ${SOURCES_GYM} ${SOURCES_GYM_DEPS})

target_compile_definitions(${GYM_LIBRARY} PRIVATE -DDISPLAY)
target_include_directories(${GYM_LIBRARY} PUBLIC "${GYM_INCLUDES}")

set_property(SOURCE save_image.cc PROPERTY INCLUDE_DIRECTORIES "${PROJECT_SOURCE_DIR}/src/gym/deps")
