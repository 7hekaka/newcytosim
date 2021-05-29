set(BASE_SOURCES
    "${SRC_BASE_DIR}/operator_new.cc"
    "${SRC_BASE_DIR}/messages.cc"
    "${SRC_BASE_DIR}/filewrapper.cc"
    "${SRC_BASE_DIR}/filepath.cc"
    "${SRC_BASE_DIR}/iowrapper.cc"
    "${SRC_BASE_DIR}/exceptions.cc"
    "${SRC_BASE_DIR}/tictoc.cc"
    "${SRC_BASE_DIR}/node.cc"
    "${SRC_BASE_DIR}/inventory.cc"
    "${SRC_BASE_DIR}/stream_func.cc"
    "${SRC_BASE_DIR}/tokenizer.cc"
    "${SRC_BASE_DIR}/glossary.cc"
    "${SRC_BASE_DIR}/property.cc"
    "${SRC_BASE_DIR}/property_list.cc"
    "${SRC_BASE_DIR}/backtrace.cc"
    "${SRC_BASE_DIR}/print_color.cc"
)

add_library(${BASE_LIB_TARGET} STATIC ${BASE_SOURCES})

set_source_files_properties(
    PROPERTIES
    EXTERNAL_OBJECT true
    GENERATED true
)

target_include_directories(${BASE_LIB_TARGET} PUBLIC
    "${SRC_BASE_DIR}"
    "${SRC_MATH_DIR}"
)
