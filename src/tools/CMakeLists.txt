set(TOOL_INCLUDES
    "${SFMT_SRC_DIR}"
    "${MATH_SRC_DIR}"
    "${BASE_SRC_DIR}"
    "${SIM_SRC_DIR}"
    "${DISP_SRC_DIR}"
    "${PLAY_SRC_DIR}"
)

set(TOOL_LIBS
    "${SIM_LIB_TARGET}"
    "${SPACES_LIB_TARGET}"
    "${MATH_LIB_TARGET}"
    "${BASE_LIB_TARGET}"
    "${LAPACK_LIB}"
    "${BLAS_LIB}"
    Threads::Threads
)

set(TOOL_OBJS
    "${PROJECT_BINARY_DIR}/src/play/CMakeFiles/play.dir/frame_reader.cc.o"
)

set(TOOL_LIST
    "frametool"
    "sieve"
    "report"
    "reportF"
    "analyse"
	"reader"
)

# Build the Tools
foreach(TOOL_NAME ${TOOL_LIST})
    add_executable("${TOOL_NAME}" ${TOOL_OBJS} "${PROJECT_SOURCE_DIR}/src/tools/${TOOL_NAME}.cc")
    set_source_files_properties(
        ${TOOL_OBJS}
        PROPERTIES
        EXTERNAL_OBJECT true
        GENERATED true
    )
    target_link_libraries(${TOOL_NAME} PUBLIC "${TOOL_LIBS}")
    target_include_directories(${TOOL_NAME} PUBLIC ${TOOL_INCLUDES})
endforeach()

