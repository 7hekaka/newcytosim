# Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#
# File src/disp/makefile.inc

OBJ_DISP := fiber_disp.o line_disp.o point_disp.o grid_display.o tesselator.o\
            gle.o gym_text.o gle_color.o gle_color_list.o gle_zoo.o view.o view_prop.o glapp.o


#-----------------------LIB & DEF for PNG support-------------------------------

IMAGES=
IMAGES_DEF=

ifneq ($(HAS_PNG), 0)

    IMAGES_DEF+=-DHAS_PNG $(INC_PNG)
    IMAGES+=$(LIB_PNG)

endif


#----------------------------targets--------------------------------------------


cytodisp.a: $(OBJ_DISP) save_imageGL.o offscreen.o
	$(MAKELIB)
	$(DONE)


save_image.o: save_image.cc save_image.h
	$(COMPILE) $(IMAGES_DEF) -c $< -o build/$@

save_imageGL.o: save_image.cc save_image.h
	$(COMPILE) $(IMAGES_DEF) -DDISPLAY -c $< -o build/$@


offscreen.o: offscreen.cc offscreen.h offscreen_fbo.cc offscreen_glx.cc
	$(COMPILE) -c $< -o build/$@


$(OBJ_DISP): %.o: %.cc %.h | build
	$(COMPILE) -Isrc/math -Isrc/base -c $< -o build/$@
