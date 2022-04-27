# Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#
# File src/disp/makefile.inc

OBJ_SIMDISP := fiber_disp.o line_disp.o point_disp.o grid_display.o

OBJ_GLE := tesselator.o gle.o gle_color.o gle_color_list.o view.o view_prop.o glapp.o

OBJ_GYM := gym_menu.o gym_flute.o gym_matrix.o\
           gym_check.o gym_text.o gym_zoo.o

OBJ_FONT := fg_font.o fg_stroke_roman.o fg_stroke_mono_roman.o

OBJ_DISP := $(OBJ_GLE) $(OBJ_GYM) $(OBJ_FONT) offscreen.o spng.o miniz.o

#-----------------------LIB & DEF for PNG support-------------------------------

IMAGES=
IMAGES_DEF=

ifneq ($(HAS_PNG), 0)

    IMAGES_DEF+=-DHAS_PNG $(INC_PNG)
    IMAGES+=$(LIB_PNG)

endif


#----------------------------targets--------------------------------------------

cytodisp.a: $(OBJ_SIMDISP) $(OBJ_DISP) save_imageGL.o
	$(MAKELIB)
	$(DONE)

# Attention: For C-language objects, use '-x c':

build/miniz.o: miniz.c miniz.h | build
	$(CXX) $(FAST) -x c -c $< -o $@

build/spng.o: spng.c spng.h | build
	$(CXX) -DSPNG_USE_MINIZ -DSPNG_STATIC -DSPNG_USE_MINIZ -x c -c $< -o $@

build/save_image.o: save_image.cc save_image.h | build
	$(COMPILE) $(IMAGES_DEF) -c $< -o $@

build/save_imageGL.o: save_image.cc save_image.h | build
	$(COMPILE) $(IMAGES_DEF) -DDISPLAY -c $< -o $@

build/offscreen.o: offscreen.cc offscreen.h offscreen_fbo.cc offscreen_glx.cc | build
	$(COMPILE) -c $< -o $@


$(OBJ_SIMDISP:%.o=build/%.o): build/%.o: %.cc %.h | build
	$(COMPILE) -Isrc/math -Isrc/base -c $< -o $@

$(OBJ_GLE:%.o=build/%.o): build/%.o: %.cc %.h | build
	$(COMPILE) $(INCLUDE) -Isrc/math -Isrc/base -c $< -o $@

$(OBJ_GYM:%.o=build/%.o): build/%.o: %.cc %.h | build
	$(COMPILE) $(INCLUDE) -Isrc/math -Isrc/base -c $< -o $@

$(OBJ_FONT:%.o=build/%.o): build/%.o: %.c fg_font.h | build
	$(COMPILE) $(INCLUDE) -DNDEBUG $(FAST) -x c -c $< -o $@
