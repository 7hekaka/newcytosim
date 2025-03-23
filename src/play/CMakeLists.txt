
set(SOURCES_PLAY
    player_prop.cc
    player.cc
    play.cc
)

list(TRANSFORM SOURCES_PLAY PREPEND "${PROJECT_SOURCE_DIR}/src/play/")

#-------------------------- Graphical Libraries --------------------------------

find_package(OpenGL)
message(">>>>>> OpenGL ${OPENGL_LIBRARIES}")

find_package(GLUT)
message(">>>>>> GLUT ${GLUT_LIBRARIES}")
if ( NOT GLUT_FOUND )
    message(">>>>>> Will build GLUT from source")
    add_subdirectory(deps/glut)
    add_dependencies(${PLAY_TARGET} ${GLUT_LIBRARIES})
endif()


if(UNIX AND NOT APPLE)
    find_package(GLEW)
    find_package(X11)
    find_package(Xt)
    set(MOREGL_LIBS
        ${GLEW_LIBRARIES}
        X11
        Xt
    )
endif()

#------------------------------ Play executable --------------------------------

add_executable(${PLAY_TARGET} ${SOURCES_PLAY})

target_link_libraries(${PLAY_TARGET} PUBLIC ${GLUT_LIBRARIES})
target_link_libraries(${PLAY_TARGET} PUBLIC ${OPENGL_LIBRARIES})

target_link_libraries(${PLAY_TARGET} PUBLIC
    "${SIMG_LIBRARY}"
    "${DISP_LIBRARY}"
    "${GYM_LIBRARY}"
    "${MATH_LIBRARY}"
    "${BASE_LIBRARY}"
    "${MOREGL_LIBS}"
    "${LAPACK_BLAS}"
    "${CURSES_LIBRARY}"
    Threads::Threads
)

set(PLAY_INCLUDES math base sim gym disp)
list(TRANSFORM PLAY_INCLUDES PREPEND "${PROJECT_SOURCE_DIR}/src/")

target_include_directories(${PLAY_TARGET} PUBLIC "${PLAY_INCLUDES}")
