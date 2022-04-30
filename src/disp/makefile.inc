# Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#
# File src/disp/makefile.inc

OBJ_DISP := view.o view_prop.o fiber_disp.o line_disp.o point_disp.o grid_display.o glapp.o


#----------------------------targets--------------------------------------------


$(OBJ_DISP:%.o=build/%.o): build/%.o: %.cc %.h | build
	$(COMPILE) $(INCLUDE) -Isrc/math -Isrc/base -Isrc/gym -c $< -o $@

