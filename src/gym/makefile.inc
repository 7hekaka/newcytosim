# cytosim/src/gym/makefile.inc   Copyright Cambridge University 2022

OBJ_GLE := tesselator.o gle.o gle_color.o gle_color_list.o

OBJ_GYM := gym_menu.o gym_flat.o gym_flute.o gym_flute_dim.o gym_matrix.o gym_view.o\
           gym_vect.o gym_draw.o gym_check.o gym_image.o gym_text.o gym_zoo.o

OBJ_FONT := fg_font.o fg_stroke.o

OBJ_GFX := $(OBJ_GLE) $(OBJ_GYM) $(OBJ_FONT) offscreen.o spng.o miniz.o

#-----------------------LIB & DEF for PNG support-------------------------------

IMAGES=
IMAGES_DEF=-DHAS_SPNG

ifneq ($(HAS_PNG), 0)

    IMAGES_DEF+=-DHAS_PNG $(INC_PNG)
    IMAGES+=$(LIB_PNG)

endif


#----------------------------targets--------------------------------------------

cytodisp.a: $(OBJ_GFX) save_imageGL.o
	$(MAKELIB)
	$(DONE)


build/save_image.o: save_image.cc save_image.h | build
	$(COMPILE) $(IMAGES_DEF) -Isrc/gym/deps -c $< -o $@

build/save_imageGL.o: save_image.cc save_image.h | build
	$(COMPILE) $(IMAGES_DEF) -Isrc/gym/deps -DDISPLAY -c $< -o $@

build/offscreen.o: offscreen.cc offscreen.h offscreen_fbo.cc offscreen_glx.cc | build
	$(COMPILE) -c $< -o $@


$(OBJ_GLE:%.o=build/%.o): build/%.o: %.cc %.h | build
	$(COMPILE) $(INCLUDE) -Isrc/math -Isrc/base -Isrc/gym -c $< -o $@

$(OBJ_GYM:%.o=build/%.o): build/%.o: %.cc %.h | build
	$(COMPILE) $(INCLUDE) -Isrc/math -Isrc/base -c $< -o $@

$(OBJ_FONT:%.o=build/%.o): build/%.o: %.cc fg_font.h | build
	$(COMPILE) $(INCLUDE) -DNDEBUG -c $< -o $@

#----------------------------C-targets------------------------------------------
# Attention: For C-language objects, use '-x c':

vpath spng.h src/gym/deps
vpath spng.c src/gym/deps
vpath miniz.h src/gym/deps
vpath miniz.c src/gym/deps

build/miniz.o: miniz.c miniz.h | build
	$(CXX) $(FAST) -x c -c $< -o $@

build/spng.o: spng.c spng.h | build
	$(CXX) $(FAST) -Isrc/gym/deps -DSPNG_STATIC -DSPNG_USE_MINIZ -x c -c $< -o $@

