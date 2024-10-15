# Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University
# Copyright 2018-- Francois J. Nedelec.


OBJ_SPACE := space.o space_prop.o space_square.o space_sphere.o space_disc.o\
             space_capsule.o space_dice.o space_torus.o space_polygon.o\
             space_ellipse.o space_cylinder.o space_bicylinder.o space_ring.o space_mesh.o\
             space_cylinderZ.o space_polygonZ.o space_strip.o space_periodic.o\
             space_banana.o space_cylinderP.o space_lid.o space_force.o\
             space_dynamic_prop.o space_dynamic_sphere.o space_dynamic_ellipse.o\
             space_dynamic_disc.o space_set.o

OBJ_HAND := hand.o hand_prop.o hand_monitor.o hand_list.o\
            motor.o motor_prop.o\
            slider.o slider_prop.o\
            actor.o actor_prop.o\
            nucleator.o nucleator_prop.o\
            tracker.o tracker_prop.o\
            cutter.o cutter_prop.o\
            chewer.o chewer_prop.o\
            rescuer.o rescuer_prop.o\
            regulator.o regulator_prop.o\
            mighty.o mighty_prop.o

OBJ_DIGIT := digit.o digit_prop.o\
             walker.o walker_prop.o\
             myosin.o myosin_prop.o\
             dynein.o dynein_prop.o\
             kinesin.o kinesin_prop.o

OBJ_SINGLE := single.o single_prop.o single_set.o\
              picket.o picket_long.o wrist.o wrist_long.o

OBJ_COUPLE := couple.o couple_prop.o couple_long.o couple_set.o\
              crosslink.o crosslink_long.o crosslink_prop.o\
              bridge.o bridge_prop.o fork.o fork_prop.o\
              duo.o duo_long.o duo_prop.o\
              shackle.o shackle_long.o shackle_prop.o\
         
OBJ_FIBER := chain.o mecafil.o fiber.o fiber_prop.o\
             fiber_segment.o fiber_site.o lattice.o\
             dynamic_fiber.o dynamic_fiber_prop.o\
             classic_fiber.o classic_fiber_prop.o\
             treadmilling_fiber.o treadmilling_fiber_prop.o\
             growing_fiber.o growing_fiber_prop.o fiber_set.o

OBJ_ASTER := organizer.o aster.o aster_prop.o nucleus.o nucleus_prop.o\
             fake.o fake_prop.o bundle.o bundle_prop.o organizer_set.o

OBJ_SIM := primitives.o mecable.o object.o object_pool.o object_set.o\
           sphere_prop.o sphere.o sphere_set.o\
           solid_prop.o solid.o solid_set.o bead.o bead_set.o\
           field.o field_prop.o field_set.o event.o event_set.o\
           mecapoint.o interpolation.o interpolation4.o\
           fiber_grid.o point_grid.o locus_grid.o meca.o\
           simul_prop.o simul.o simul_part.o interface.o parser.o

# objects needed to load/run Cytosim
OBJ_READ := sim_thread.o frame_reader.o

OBJ_CYTOLIB:=$(OBJ_SPACE) $(OBJ_HAND) $(OBJ_DIGIT) $(OBJ_FIBER)\
             $(OBJ_SINGLE) $(OBJ_COUPLE) $(OBJ_ASTER) $(OBJ_SIM) $(OBJ_READ)

# objects containing graphical code under the DISPLAY keyword:
OBJ_SIMDISP:=$(OBJ_SPACE) field.o fiber_grid.o point_grid.o
OBJ_CYTOLIBG:=$(OBJ_SIMDISP:%.o=%G.o) $(filter-out $(OBJ_SIMDISP),$(OBJ_CYTOLIB))

#----------------------------targets--------------------------------------------

INC_BM:=-Isrc/base -Isrc/math

#--- normal build:

sim: sim.cc cytolib.a cytomath.a cytobase.a | bin
	$(COMPILE) $(INC_BM) $(INFO) $(OBJECTS) $(LINK) -o bin/$@
	$(DONE)
vpath sim  bin

#--- build for debuging:

simd: sim.cc cytolib.a cytomath.a cytobase.a | bin
	$(CXX) $(FlagsD) $(INC_BM) $(INFO) $(OBJECTS) $(LINK) -o bin/$@
	$(DONE)
vpath simd  bin

#--- builds with specific dimension:

bin1/sim: sim.cc cytolibD1.a cytomathD1.a cytobase.a | bin1
	$(COMPILE) $(INC_BM) $(INFO) -DDIM=1 $(OBJECTS) $(LINK) -o bin1/sim
	$(DONE)

bin2/sim: sim.cc cytolibD2.a cytomathD2.a cytobase.a | bin2
	$(COMPILE) $(INC_BM) $(INFO) -DDIM=2 $(OBJECTS) $(LINK) -o bin2/sim
	$(DONE)

bin3/sim: sim.cc cytolibD3.a cytomathD3.a cytobase.a | bin3
	$(COMPILE) $(INC_BM) $(INFO) -DDIM=3 $(OBJECTS) $(LINK) -o bin3/sim
	$(DONE)


#--- libraries:

cytolib.a: $(OBJ_CYTOLIB:%.o=build/%.o) | build lib
	$(LIBTOOL) lib/$@ $^
	$(DONE)

cytolibD1.a: $(OBJ_CYTOLIB:%.o=build/%D1.o) | build lib
	$(LIBTOOL) lib/$@ $^
	$(DONE)

cytolibD2.a: $(OBJ_CYTOLIB:%.o=build/%D2.o) | build lib
	$(LIBTOOL) lib/$@ $^
	$(DONE)

cytolibD3.a: $(OBJ_CYTOLIB:%.o=build/%D3.o) | build lib
	$(LIBTOOL) lib/$@ $^
	$(DONE)

#--- libraries with graphical code enabled

cytolibG.a: $(OBJ_CYTOLIBG:%.o=build/%.o) | build lib
	$(LIBTOOL) lib/$@ $^
	$(DONE)

cytolibGD1.a: $(OBJ_CYTOLIBG:%.o=build/%D1.o) | build lib
	$(LIBTOOL) lib/$@ $^
	$(DONE)

cytolibGD2.a: $(OBJ_CYTOLIBG:%.o=build/%D2.o) | build lib
	$(LIBTOOL) lib/$@ $^
	$(DONE)

cytolibGD3.a: $(OBJ_CYTOLIBG:%.o=build/%D3.o) | build lib
	$(LIBTOOL) lib/$@ $^
	$(DONE)

#-----------minimal include directories, DIM defined in 'dim.h'-----------------

INC_BMS=$(addprefix -Isrc/, base math sim)


$(OBJ_SPACE:%.o=build/%.o): build/%.o: %.cc %.h | build
	$(COMPILE) $(INC_BMS) -Isrc/sim/spaces -c $< -o $@

$(OBJ_HAND:%.o=build/%.o): build/%.o: %.cc %.h | build
	$(COMPILE) $(INC_BMS) -Isrc/sim/hands -c $< -o $@

$(OBJ_DIGIT:%.o=build/%.o): build/%.o: %.cc %.h | build
	$(COMPILE) $(INC_BMS) -Isrc/sim/hands -c $< -o $@

$(OBJ_FIBER:%.o=build/%.o): build/%.o: %.cc %.h | build
	$(COMPILE) $(INC_BMS) -Isrc/sim/fibers -c $< -o $@

$(OBJ_SINGLE:%.o=build/%.o): build/%.o: %.cc %.h | build
	$(COMPILE) $(INC_BMS) -Isrc/sim/singles -c $< -o $@

$(OBJ_COUPLE:%.o=build/%.o): build/%.o: %.cc %.h | build
	$(COMPILE) $(INC_BMS) -Isrc/sim/couples -c $< -o $@

$(OBJ_ASTER:%.o=build/%.o): build/%.o: %.cc %.h | build
	$(COMPILE) $(INC_BMS) -Isrc/sim/organizers -c $< -o $@

$(OBJ_SIM:%.o=build/%.o): build/%.o: %.cc %.h | build
	$(COMPILE) $(INC_BMS) -c $< -o $@

$(OBJ_READ:%.o=build/%.o): build/%.o: %.cc %.h | build
	$(COMPILE) $(INC_BMS) -c $< -o $@

#------------targets with 'DIM' defined in makefile-----------------------------

INC_SIM=$(addprefix -Isrc/, base math sim sim/fibers sim/hands sim/singles sim/couples sim/organizers sim/spaces)


$(OBJ_CYTOLIB:%.o=build/%D1.o): build/%D1.o: %.cc %.h | build
	$(COMPILE) -DDIM=1 $(INC_SIM) -c $< -o $@

$(OBJ_CYTOLIB:%.o=build/%D2.o): build/%D2.o: %.cc %.h | build
	$(COMPILE) -DDIM=2 $(INC_SIM) -c $< -o $@

$(OBJ_CYTOLIB:%.o=build/%D3.o): build/%D3.o: %.cc %.h | build
	$(COMPILE) -DDIM=3 $(INC_SIM) -c $< -o $@


#---------------------targets with DISPLAY elements-----------------------------

INC_BMSD=$(addprefix -Isrc/, base math sim sim/spaces gym disp)

$(OBJ_SIMDISP:%.o=build/%G.o): build/%G.o: %.cc %.h | build
	$(COMPILE) -DDISPLAY $(INC_BMSD) -c $< -o $@

$(OBJ_SIMDISP:%.o=build/%GD1.o): build/%GD1.o: %.cc %.h | build
	$(COMPILE) -DDIM=1 -DDISPLAY $(INC_BMSD) -c $< -o $@

$(OBJ_SIMDISP:%.o=build/%GD2.o): build/%GD2.o: %.cc %.h | build
	$(COMPILE) -DDIM=2 -DDISPLAY $(INC_BMSD) -c $< -o $@

$(OBJ_SIMDISP:%.o=build/%GD3.o): build/%GD3.o: %.cc %.h | build
	$(COMPILE) -DDIM=3 -DDISPLAY $(INC_BMSD) -c $< -o $@

