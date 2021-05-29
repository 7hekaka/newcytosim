set(SOURCES_PLAY
    "${SRC_PLAY_DIR}/frame_reader.cc"
    "${SRC_PLAY_DIR}/display_prop.cc"
    "${SRC_DISP_DIR}/gle_flute.cc"
    "${SRC_PLAY_DIR}/display.cc"
    "${SRC_PLAY_DIR}/display1.cc"
    "${SRC_PLAY_DIR}/display2.cc"
    "${SRC_PLAY_DIR}/display3.cc"
    "${SRC_PLAY_DIR}/sim_thread.cc"
    "${SRC_PLAY_DIR}/player_prop.cc"
    "${SRC_PLAY_DIR}/player.cc"
    "${SRC_PLAY_DIR}/play.cc"
)

if(BUILD_PLAY_EXECUTABLE)
    add_executable(${PLAY_EXE_TARGET} ${SOURCES_PLAY})

    target_link_libraries(${PLAY_EXE_TARGET} PUBLIC
        "${SIMG_LIB_TARGET}"
        "${DISP_LIB_TARGET}"
        "${MATH_LIB_TARGET}"
        "${BASE_LIB_TARGET}"
        "${OPENGL_LIBS}"
        "${LAPACK_LIB}"
        "${BLAS_LIB}"
        Threads::Threads
    )

    target_include_directories(${PLAY_EXE_TARGET} PUBLIC
        "${SRC_MATH_DIR}"
        "${SRC_BASE_DIR}"
        "${SRC_SIM_ALL_DIR}"
        "${SRC_DISP_DIR}"
    )
endif(BUILD_PLAY_EXECUTABLE)
