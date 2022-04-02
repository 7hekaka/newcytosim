set(SOURCES_DISP
    gle.cc
    gym_text.cc
    gym_matrix.cc
    gym_flute.cc
    gym_check.cc
    gym_zoo.cc
    gle_color.cc
    gle_color_list.cc
    view.cc
    view_prop.cc
    grid_display.cc
    glapp.cc
    save_image.cc
    offscreen.cc
    fiber_disp.cc
    line_disp.cc
    point_disp.cc
    fg_font.c
    fg_stroke_roman.c
    fg_stroke_mono_roman.c
)

list(TRANSFORM SOURCES_DISP PREPEND "${SRC_DISP_DIR}/")

add_library(${DISP_LIB_TARGET} STATIC ${SOURCES_DISP})
target_compile_definitions(${DISP_LIB_TARGET} PRIVATE DISPLAY)

target_include_directories(${DISP_LIB_TARGET} PUBLIC
    "${SRC_MATH_DIR}"
    "${SRC_BASE_DIR}"
)
