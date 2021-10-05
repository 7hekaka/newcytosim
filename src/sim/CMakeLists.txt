set(SOURCES_SPACES
    spaces/space_square.cc
    spaces/space_sphere.cc
    spaces/space_dice.cc
    spaces/space_torus.cc
    spaces/space_polygon.cc
    spaces/space_banana.cc
    spaces/space_ellipse.cc
    spaces/space_cylinder.cc
    spaces/space_cylinderZ.cc
    spaces/space_capsule.cc
    spaces/space_strip.cc
    spaces/space_periodic.cc
    spaces/space_cylinderP.cc
    spaces/space_dynamic_prop.cc
    spaces/space_dynamic_ellipse.cc
    spaces/space_dynamic_sphere.cc
    spaces/space_dynamic_disc.cc
    spaces/space_force.cc
    spaces/space_polygonZ.cc
    spaces/space_tee.cc
    spaces/space_disc.cc
    spaces/space_lid.cc
    spaces/space_ring.cc
    spaces/space_mesh.cc
)

set(SOURCES_SIM
    modulo.cc movable.cc
    mecable.cc meca.cc mecapoint.cc
    interpolation.cc interpolation4.cc
    fiber_segment.cc fiber_site.cc
    lattice.cc
    object.cc object_pool.cc object_set.cc
    sphere_prop.cc sphere.cc sphere_set.cc
    bead.cc bead_set.cc
    solid_prop.cc solid.cc solid_set.cc
    field.cc field_prop.cc field_set.cc
    event.cc event_set.cc
    chain.cc mecafil.cc
    fiber.cc fiber_prop.cc fiber_set.cc
    hand.cc hand_list.cc hand_prop.cc
    single.cc single_prop.cc single_set.cc
    couple.cc couple_prop.cc couple_set.cc
    organizer.cc organizer_set.cc
    tubule_prop.cc tubule.cc tubule_set.cc
    fiber_grid.cc point_grid.cc locus_grid.cc
    space.cc space_prop.cc space_set.cc
    simul.cc simul_prop.cc
    interface.cc parser.cc
)

set(SOURCES_HANDS
    hands/motor.cc hands/motor_prop.cc
    hands/slider.cc hands/slider_prop.cc
    hands/actor.cc hands/actor_prop.cc
    hands/nucleator.cc hands/nucleator_prop.cc
    hands/tracker.cc hands/tracker_prop.cc
    hands/cutter.cc hands/cutter_prop.cc
    hands/rescuer.cc hands/rescuer_prop.cc
    hands/mighty.cc hands/mighty_prop.cc
    hands/chewer_prop.cc hands/chewer.cc
    hands/walker_prop.cc hands/walker.cc
    hands/regulator_prop.cc hands/regulator.cc
    hands/digit_prop.cc hands/digit.cc
    hands/myosin_prop.cc hands/myosin.cc
    hands/dynein_prop.cc hands/dynein.cc
    hands/kinesin_prop.cc hands/kinesin.cc
)

set(SOURCES_FIBERS
    fibers/dynamic_fiber.cc fibers/dynamic_fiber_prop.cc
    fibers/classic_fiber.cc fibers/classic_fiber_prop.cc
    fibers/treadmilling_fiber.cc fibers/treadmilling_fiber_prop.cc
    fibers/growing_fiber.cc fibers/growing_fiber_prop.cc
)

set(SOURCES_SINGLES
    singles/picket.cc singles/picket_long.cc
    singles/wrist.cc singles/wrist_long.cc
)

set(SOURCES_COUPLES
    couples/couple_long.cc
    couples/crosslink.cc couples/crosslink_long.cc couples/crosslink_prop.cc
    couples/bridge.cc couples/bridge_prop.cc
    couples/shackle.cc couples/shackle_long.cc couples/shackle_prop.cc
    couples/fork.cc couples/fork_prop.cc
    couples/duo.cc couples/duo_prop.cc couples/duo_long.cc
)

set(SOURCES_ORGANIZERS
    organizers/aster.cc organizers/aster_prop.cc
    organizers/nucleus.cc organizers/nucleus_prop.cc
    organizers/fake.cc organizers/fake_prop.cc
    organizers/bundle.cc organizers/bundle_prop.cc
)

set(SOURCES_CYTOSIM
    "${SOURCES_SPACES}"
    "${SOURCES_SIM}"
    "${SOURCES_HANDS}"
    "${SOURCES_FIBERS}"
    "${SOURCES_SINGLES}"
    "${SOURCES_COUPLES}"
    "${SOURCES_ORGANIZERS}"
)

# objects containing graphical code under the DISPLAY keyword:
set(SOURCES_DISPLAY_SENSITIVE
    ${SOURCES_SPACES}
    space.cc
    field.cc
    fiber_grid.cc
    point_grid.cc
)

list(TRANSFORM SOURCES_CYTOSIM PREPEND "${SRC_SIM_DIR}/")
list(TRANSFORM SOURCES_DISPLAY_SENSITIVE PREPEND "${SRC_SIM_DIR}/")


if(BUILD_SIM)

    add_library(${SIM_LIB_TARGET} STATIC ${SOURCES_CYTOSIM})
    target_include_directories(${SIM_LIB_TARGET} PUBLIC
        "${SRC_MATH_DIR}"
        "${SRC_BASE_DIR}"
        "${SRC_SIM_ALL_DIR}"
    )

    add_executable(${SIM_TARGET} "${SRC_SIM_DIR}/sim.cc")
    target_link_libraries(${SIM_TARGET} PRIVATE
        "${SIM_LIB_TARGET}"
        "${MATH_LIB_TARGET}"
        "${BASE_LIB_TARGET}"
        "${LAPACK_LIB}"
        "${BLAS_LIB}"
        Threads::Threads
    )
    target_include_directories(${SIM_TARGET} PUBLIC
        "${SRC_MATH_DIR}"
        "${SRC_BASE_DIR}"
        "${SRC_SIM_DIR}"
    )
    
endif(BUILD_SIM)



if(BUILD_PLAY)

    #set(SOURCES_CYTOSIM_NO_DISPLAY ${SOURCES_CYTOSIM})
    #list(REMOVE_ITEM SOURCES_CYTOSIM_NO_DISPLAY ${SOURCES_DISPLAY_SENSITIVE})

    add_library(${SIMG_LIB_TARGET} STATIC ${SOURCES_CYTOSIM})

    target_include_directories(${SIMG_LIB_TARGET} PUBLIC
        "${SRC_MATH_DIR}"
        "${SRC_BASE_DIR}"
        "${SRC_SIM_ALL_DIR}"
        "${SRC_DISP_DIR}"
    )


    # It is only necessary to compile files in ${SOURCES_DISPLAY_SENSITIVE} with this flag
    # but I do not know how to specify this with CMake
    # Hence by setting a target_compile_definitions, all source will be compiled with -DDISPLAY

    target_compile_definitions(${SIMG_LIB_TARGET} PRIVATE DISPLAY)

    # This below does not work, as it sets the compile options for all targets including the files

    foreach(FILE ${SOURCES_DISPLAY_SENSITIVE})
        #set_source_files_properties(${FILE} TARGET_DIRECTORY ${SIMG_LIB_TARGET} PROPERTIES COMPILE_DEFINITIONS "DISPLAY")
        #message(${FILE} "will be compiled with -DDISPLAY")
    endforeach()

endif(BUILD_PLAY)
