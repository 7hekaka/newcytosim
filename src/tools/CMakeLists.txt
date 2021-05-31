set(SOURCES_TOOL
    "${SRC_PLAY_DIR}/frame_reader.cc"
)

set(TOOL_INCLUDES
    "${SRC_MATH_DIR}"
    "${SRC_BASE_DIR}"
    "${SRC_SIM_ALL_DIR}"
    "${SRC_DISP_DIR}"
    "${SRC_PLAY_DIR}"
)

set(TOOL_LIBS
    "${SIM_LIB_TARGET}"
    "${MATH_LIB_TARGET}"
    "${BASE_LIB_TARGET}"
    "${LAPACK_LIB}"
    "${BLAS_LIB}"
    Threads::Threads
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
foreach(TOOL ${TOOL_LIST})
    add_executable("${TOOL}" ${SOURCES_TOOL} "${SRC_TOOLS_DIR}/${TOOL}.cc")
    set_source_files_properties(
        PROPERTIES
        EXTERNAL_OBJECT true
        GENERATED true
    )
    target_link_libraries(${TOOL} PUBLIC "${TOOL_LIBS}")
    target_include_directories(${TOOL} PUBLIC ${TOOL_INCLUDES})
endforeach()

