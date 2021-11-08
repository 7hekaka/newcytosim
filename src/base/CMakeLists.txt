set(SOURCES_BASE
    operator_new.cc
    messages.cc
    filewrapper.cc
    filepath.cc
    iowrapper.cc
    exceptions.cc
    time_date.cc
    node.cc
    inventory.cc
    stream_func.cc
    tokenizer.cc
    glossary.cc
    property.cc
    property_list.cc
    backtrace.cc
    print_color.cc
)

list(TRANSFORM SOURCES_BASE PREPEND "${SRC_BASE_DIR}/")

add_library(${BASE_LIB_TARGET} STATIC ${SOURCES_BASE})

set_source_files_properties(
    PROPERTIES
    EXTERNAL_OBJECT true
    GENERATED true
)

target_include_directories(${BASE_LIB_TARGET} PUBLIC
    "${SRC_BASE_DIR}"
    "${SRC_MATH_DIR}"
)
