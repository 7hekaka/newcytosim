# Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#
# File src/disp/makefile.inc

OBJ_DISP := fiber_disp.o line_disp.o point_disp.o grid_display.o tesselator.o\
            gle.o gle_color.o gle_color_list.o gym_text.o gym_zoo.o gym_check.o\
            view.o view_prop.o glapp.o

OBJ_FONT := fg_font.o fg_stroke_roman.o fg_stroke_mono_roman.o


#-----------------------LIB & DEF for PNG support-------------------------------

IMAGES=
IMAGES_DEF=

ifneq ($(HAS_PNG), 0)

    IMAGES_DEF+=-DHAS_PNG $(INC_PNG)
    IMAGES+=$(LIB_PNG)

endif


#----------------------------targets--------------------------------------------

cytodisp.a: $(OBJ_DISP) $(OBJ_FONT) save_imageGL.o offscreen.o libspng.o miniz.o
	$(MAKELIB)
	$(DONE)


# Attention: FreeGLUT font are coded in C-language (use -x c):

fg_font.o: fg_font.c fg_font.h
	$(CXX) -DNDEBUG $(FAST) -x c -c $< -o build/$@

fg_stroke_roman.o: fg_stroke_roman.c fg_font.h
	$(CXX) -DNDEBUG $(FAST) -x c -c $< -o build/$@

fg_stroke_mono_roman.o: fg_stroke_mono_roman.c fg_font.h
	$(CXX) -DNDEBUG $(FAST) -x c -c $< -o build/$@


miniz.o: miniz.c miniz.h
	$(CXX) $(FAST) -x c -c $< -o build/$@

libspng.o: spng.c spng.h
	$(CXX) -DSPNG_USE_MINIZ -DNDEBUG $(FAST) -x c -c $< -o build/$@

save_image.o: save_image.cc save_image.h
	$(COMPILE) $(IMAGES_DEF) -c $< -o build/$@

save_imageGL.o: save_image.cc save_image.h
	$(COMPILE) $(IMAGES_DEF) -DDISPLAY -c $< -o build/$@


offscreen.o: offscreen.cc offscreen.h offscreen_fbo.cc offscreen_glx.cc
	$(COMPILE) -c $< -o build/$@


$(OBJ_DISP): %.o: %.cc %.h | build
	$(COMPILE) -Isrc/math -Isrc/base -c $< -o build/$@
