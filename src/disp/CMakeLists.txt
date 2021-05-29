set(SOURCES_DISP
    "${SRC_DISP_DIR}/gle.cc"
    "${SRC_DISP_DIR}/gle_flute.cc"
    "${SRC_DISP_DIR}/gle_zoo.cc"
    "${SRC_DISP_DIR}/gle_color.cc"
    "${SRC_DISP_DIR}/gle_color_list.cc"
    "${SRC_DISP_DIR}/view.cc"
    "${SRC_DISP_DIR}/view_prop.cc"
    "${SRC_DISP_DIR}/grid_display.cc"
    "${SRC_DISP_DIR}/glapp.cc"
    "${SRC_DISP_DIR}/saveimage.cc"
    "${SRC_DISP_DIR}/offscreen.cc"
    "${SRC_DISP_DIR}/fiber_disp.cc"
    "${SRC_DISP_DIR}/line_disp.cc"
    "${SRC_DISP_DIR}/point_disp.cc"
)

add_library(${DISP_LIB_TARGET} STATIC ${SOURCES_DISP})

target_include_directories(${DISP_LIB_TARGET} PUBLIC
    "${SRC_MATH_DIR}"
    "${SRC_BASE_DIR}"
    "${SRC_SIM_ALL_DIR}"
)
