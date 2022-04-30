set(SOURCES_DISP
    view.cc
    view_prop.cc
    glapp.cc
    fiber_disp.cc
    line_disp.cc
    point_disp.cc
    grid_display.cc
    display_prop.cc
    display.cc
    display1.cc
    display2.cc
    display3.cc
)

list(TRANSFORM SOURCES_DISP PREPEND "${PROJECT_SOURCE_DIR}/src/disp/")


set(DISP_INCLUDES math base sim sim/organizers gym)
list(TRANSFORM DISP_INCLUDES PREPEND "${PROJECT_SOURCE_DIR}/src/")

add_library(${DISP_LIBRARY} STATIC ${SOURCES_DISP})

target_compile_definitions(${DISP_LIBRARY} PRIVATE -DDISPLAY)
target_include_directories(${DISP_LIBRARY} PUBLIC "${DISP_INCLUDES}" )
