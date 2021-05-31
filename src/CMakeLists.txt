
add_subdirectory("base")
add_subdirectory("math")
add_subdirectory("sim")

if (BUILD_PLAY)
    add_subdirectory("disp")
    add_subdirectory("play")
endif()

if (BUILD_TOOLS)
    add_subdirectory("tools")
endif()

if (BUILD_TESTS)
    add_subdirectory("test")
endif()
