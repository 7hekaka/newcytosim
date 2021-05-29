set(SOURCES_SPACES
    "${SRC_SIM_DIR}/space.cc"
    "${SRC_SIM_DIR}/space_prop.cc"
    "${SRC_SIM_DIR}/spaces/space_square.cc"
    "${SRC_SIM_DIR}/spaces/space_sphere.cc"
    "${SRC_SIM_DIR}/spaces/space_dice.cc"
    "${SRC_SIM_DIR}/spaces/space_torus.cc"
    "${SRC_SIM_DIR}/spaces/space_polygon.cc"
    "${SRC_SIM_DIR}/spaces/space_banana.cc"
    "${SRC_SIM_DIR}/spaces/space_ellipse.cc"
    "${SRC_SIM_DIR}/spaces/space_cylinder.cc"
    "${SRC_SIM_DIR}/spaces/space_cylinderZ.cc"
    "${SRC_SIM_DIR}/spaces/space_capsule.cc"
    "${SRC_SIM_DIR}/spaces/space_strip.cc"
    "${SRC_SIM_DIR}/spaces/space_periodic.cc"
    "${SRC_SIM_DIR}/spaces/space_cylinderP.cc"
    
    "${SRC_SIM_DIR}/spaces/space_dynamic_prop.cc"
    "${SRC_SIM_DIR}/spaces/space_dynamic_ellipse.cc"
    "${SRC_SIM_DIR}/spaces/space_dynamic_sphere.cc"
    
    "${SRC_SIM_DIR}/spaces/space_force.cc"
    "${SRC_SIM_DIR}/spaces/space_polygonZ.cc"
    "${SRC_SIM_DIR}/spaces/space_tee.cc"
    "${SRC_SIM_DIR}/spaces/space_disc.cc"
    "${SRC_SIM_DIR}/spaces/space_lid.cc"
    "${SRC_SIM_DIR}/spaces/space_ring.cc"
    "${SRC_SIM_DIR}/spaces/space_mesh.cc"
)

set(SOURCES_SIM
    "${SRC_SIM_DIR}/modulo.cc"
    "${SRC_SIM_DIR}/movable.cc"
    "${SRC_SIM_DIR}/mecable.cc"
    "${SRC_SIM_DIR}/mecapoint.cc"
    "${SRC_SIM_DIR}/interpolation.cc"
    "${SRC_SIM_DIR}/interpolation4.cc"
    "${SRC_SIM_DIR}/fiber_segment.cc"
    "${SRC_SIM_DIR}/lattice.cc"
    "${SRC_SIM_DIR}/object.cc"
    "${SRC_SIM_DIR}/object_pool.cc"
    "${SRC_SIM_DIR}/object_set.cc"
    "${SRC_SIM_DIR}/fiber_site.cc"
    
    "${SRC_SIM_DIR}/sphere_prop.cc"
    "${SRC_SIM_DIR}/sphere.cc"
    "${SRC_SIM_DIR}/sphere_set.cc"
    
    "${SRC_SIM_DIR}/bead.cc"
    "${SRC_SIM_DIR}/bead_set.cc"
    
    "${SRC_SIM_DIR}/solid_prop.cc"
    "${SRC_SIM_DIR}/solid.cc"
    "${SRC_SIM_DIR}/solid_set.cc"
    
    "${SRC_SIM_DIR}/field.cc"
    "${SRC_SIM_DIR}/field_prop.cc"
    "${SRC_SIM_DIR}/field_set.cc"
    
    "${SRC_SIM_DIR}/event.cc"
    "${SRC_SIM_DIR}/event_set.cc"
    
    "${SRC_SIM_DIR}/tubule_prop.cc"
    "${SRC_SIM_DIR}/tubule.cc"
    "${SRC_SIM_DIR}/tubule_set.cc"
    
    "${SRC_SIM_DIR}/meca.cc"
    "${SRC_SIM_DIR}/simul_prop.cc"
    "${SRC_SIM_DIR}/fiber_grid.cc"
    "${SRC_SIM_DIR}/point_grid.cc"
    "${SRC_SIM_DIR}/locus_grid.cc"
    "${SRC_SIM_DIR}/space_set.cc"
    "${SRC_SIM_DIR}/simul.cc"
    "${SRC_SIM_DIR}/interface.cc"
    "${SRC_SIM_DIR}/parser.cc"
)

set(SOURCES_HANDS
    "${SRC_SIM_DIR}/hand.cc"
    "${SRC_SIM_DIR}/hand_list.cc"
    "${SRC_SIM_DIR}/hand_prop.cc"
    "${SRC_SIM_DIR}/hands/motor.cc"
    "${SRC_SIM_DIR}/hands/motor_prop.cc"
    "${SRC_SIM_DIR}/hands/slider.cc"
    "${SRC_SIM_DIR}/hands/slider_prop.cc"
    "${SRC_SIM_DIR}/hands/actor.cc"
    "${SRC_SIM_DIR}/hands/actor_prop.cc"
    "${SRC_SIM_DIR}/hands/nucleator.cc"
    "${SRC_SIM_DIR}/hands/nucleator_prop.cc"
    "${SRC_SIM_DIR}/hands/tracker.cc"
    "${SRC_SIM_DIR}/hands/tracker_prop.cc"
    "${SRC_SIM_DIR}/hands/cutter.cc"
    "${SRC_SIM_DIR}/hands/cutter_prop.cc"
    "${SRC_SIM_DIR}/hands/rescuer.cc"
    "${SRC_SIM_DIR}/hands/rescuer_prop.cc"
    "${SRC_SIM_DIR}/hands/mighty.cc"
    "${SRC_SIM_DIR}/hands/mighty_prop.cc"
    
    "${SRC_SIM_DIR}/hands/chewer_prop.cc"
    "${SRC_SIM_DIR}/hands/chewer.cc"
    "${SRC_SIM_DIR}/hands/walker_prop.cc"
    "${SRC_SIM_DIR}/hands/walker.cc"
    "${SRC_SIM_DIR}/hands/regulator_prop.cc"
    "${SRC_SIM_DIR}/hands/regulator.cc"
    
    "${SRC_SIM_DIR}/hands/digit_prop.cc"
    "${SRC_SIM_DIR}/hands/digit.cc"
    "${SRC_SIM_DIR}/hands/myosin_prop.cc"
    "${SRC_SIM_DIR}/hands/myosin.cc"
    "${SRC_SIM_DIR}/hands/dynein_prop.cc"
    "${SRC_SIM_DIR}/hands/dynein.cc"
    "${SRC_SIM_DIR}/hands/kinesin_prop.cc"
    "${SRC_SIM_DIR}/hands/kinesin.cc"
)

set(SOURCES_FIBERS
    "${SRC_SIM_DIR}/chain.cc"
    "${SRC_SIM_DIR}/fiber.cc"
    "${SRC_SIM_DIR}/mecafil.cc"
    "${SRC_SIM_DIR}/fiber_prop.cc"
    "${SRC_SIM_DIR}/fiber_set.cc"
    "${SRC_SIM_DIR}/fibers/dynamic_fiber.cc"
    "${SRC_SIM_DIR}/fibers/dynamic_fiber_prop.cc"
    "${SRC_SIM_DIR}/fibers/classic_fiber.cc"
    "${SRC_SIM_DIR}/fibers/classic_fiber_prop.cc"
    "${SRC_SIM_DIR}/fibers/treadmilling_fiber.cc"
    "${SRC_SIM_DIR}/fibers/treadmilling_fiber_prop.cc"
    "${SRC_SIM_DIR}/fibers/growing_fiber.cc"
    "${SRC_SIM_DIR}/fibers/growing_fiber_prop.cc"
)

set(SOURCES_SINGLES
    "${SRC_SIM_DIR}/single.cc"
    "${SRC_SIM_DIR}/single_prop.cc"
    "${SRC_SIM_DIR}/single_set.cc"
    "${SRC_SIM_DIR}/singles/picket.cc"
    "${SRC_SIM_DIR}/singles/picket_long.cc"
    "${SRC_SIM_DIR}/singles/wrist.cc"
    "${SRC_SIM_DIR}/singles/wrist_long.cc"
)

set(SOURCES_COUPLES
    "${SRC_SIM_DIR}/couple.cc"
    "${SRC_SIM_DIR}/couple_prop.cc"
    "${SRC_SIM_DIR}/couple_set.cc"
    "${SRC_SIM_DIR}/couples/couple_long.cc"
    "${SRC_SIM_DIR}/couples/crosslink.cc"
    "${SRC_SIM_DIR}/couples/crosslink_long.cc"
    "${SRC_SIM_DIR}/couples/crosslink_prop.cc"
    "${SRC_SIM_DIR}/couples/bridge.cc"
    "${SRC_SIM_DIR}/couples/bridge_prop.cc"
    "${SRC_SIM_DIR}/couples/shackle.cc"
    "${SRC_SIM_DIR}/couples/shackle_long.cc"
    "${SRC_SIM_DIR}/couples/shackle_prop.cc"
    "${SRC_SIM_DIR}/couples/fork.cc"
    "${SRC_SIM_DIR}/couples/fork_prop.cc"
    "${SRC_SIM_DIR}/couples/duo.cc"
    "${SRC_SIM_DIR}/couples/duo_prop.cc"
    "${SRC_SIM_DIR}/couples/duo_long.cc"
)

set(SOURCES_ORGANIZERS
    "${SRC_SIM_DIR}/organizer.cc"
    "${SRC_SIM_DIR}/organizer_set.cc"
    "${SRC_SIM_DIR}/organizers/aster.cc"
    "${SRC_SIM_DIR}/organizers/aster_prop.cc"
    "${SRC_SIM_DIR}/organizers/nucleus.cc"
    "${SRC_SIM_DIR}/organizers/nucleus_prop.cc"
    "${SRC_SIM_DIR}/organizers/fake.cc"
    "${SRC_SIM_DIR}/organizers/fake_prop.cc"
    "${SRC_SIM_DIR}/organizers/bundle.cc"
    "${SRC_SIM_DIR}/organizers/bundle_prop.cc"
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

add_library(${SIM_LIB_TARGET} STATIC ${SOURCES_CYTOSIM})

target_include_directories(${SIM_LIB_TARGET} PUBLIC
    "${SRC_MATH_DIR}"
    "${SRC_BASE_DIR}"
    "${SRC_SIM_ALL_DIR}"
)

# objects containing graphical code under the DISPLAY keyword:
set(SOURCES_CYTOSIM_DISPLAY
    "${SOURCES_SPACES}"
    "${SRC_SIM_DIR}/field.cc"
    "${SRC_SIM_DIR}/fiber_grid.cc"
    "${SRC_SIM_DIR}/point_grid.cc"
)

#set(SOURCES_CYTOSIM_NO_DISPLAY ${SOURCES_CYTOSIM})
#list(REMOVE_ITEM SOURCES_CYTOSIM_NO_DISPLAY ${SOURCES_CYTOSIM_DISPLAY})

add_library(${SIMG_LIB_TARGET} STATIC ${SOURCES_CYTOSIM})
target_compile_definitions(${SIMG_LIB_TARGET} PRIVATE DISPLAY)
foreach(FILE ${SOURCES_CYTOSIM_DISPLAY})
    #set_source_files_properties(${FILE} TARGET_DIRECTORY ${SIMG_LIB_TARGET} PROPERTIES COMPILE_DEFINITIONS "DISPLAY")
    #message(${FILE} "will be compiled with -DDISPLAY")
endforeach()
target_include_directories(${SIMG_LIB_TARGET} PUBLIC
    "${SRC_MATH_DIR}"
    "${SRC_BASE_DIR}"
    "${SRC_SIM_ALL_DIR}"
    "${SRC_DISP_DIR}"
)


if(BUILD_SIM_EXECUTABLE)
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
        "${SRC_SIM_ALL_DIR}"
    )
endif(BUILD_SIM_EXECUTABLE)
