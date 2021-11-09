set(SOURCES_PLAY
    frame_reader.cc
    display_prop.cc
    display.cc
    display1.cc
    display2.cc
    display3.cc
    sim_thread.cc
    player_prop.cc
    player.cc
    play.cc
)

list(TRANSFORM SOURCES_PLAY PREPEND "${SRC_PLAY_DIR}/")

if(MAKE_PLAY)
    add_executable(${PLAY_TARGET} ${SOURCES_PLAY})

    target_link_libraries(${PLAY_TARGET} PUBLIC
        "${SIMG_LIB_TARGET}"
        "${DISP_LIB_TARGET}"
        "${MATH_LIB_TARGET}"
        "${BASE_LIB_TARGET}"
        "${OPENGL_LIBS}"
        "${LAPACK_LIB}"
        "${BLAS_LIB}"
        Threads::Threads
    )

    target_include_directories(${PLAY_TARGET} PUBLIC
        "${SRC_MATH_DIR}"
        "${SRC_BASE_DIR}"
        "${SRC_SIM_DIR}"
        "${SRC_DISP_DIR}"
    )
endif(MAKE_PLAY)
