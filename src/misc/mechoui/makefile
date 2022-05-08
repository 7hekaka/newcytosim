# For Cytosim/misc/mechoui. Sainsbury Laboratory; Cambridge University, 2022

#-------------------include compiler-specifications-----------------------------

include makefile.inc

#-------------------check compiler and compiler-flags---------------------------

ifeq ($CXX,)
   $(error No compiler defined for $$(MACHINE)=$(MACHINE))
endif

ifeq ($(Flags$(MODE)),)
   $(error No compiler-options defined for $$(MACHINE)=$(MACHINE) in mode $(MODE))
endif

COMPILE=$(CXX) $(CXXFLG) $(Flags$(MODE))

ifneq ($(HAS_PNG), 0)
    IMAGES_DEF = -DHAS_PNG $(INC_PNG)
    IMAGES = $(LIB_PNG)
endif

#----------------------------objects--------------------------------------------

INC_MECHOUI := $(realpath $(addprefix ../../, base math gym disp))

vpath %.h  . $(INC_MECHOUI)
vpath %.cc . $(INC_MECHOUI)
vpath %.c  $(INC_MECHOUI)
vpath %.o  build
vpath %.a  build


OBJ_MECHOUI := multi_mesh.o mechoui_param.o vector3.o glossary.o property.o\
           backtrace.o tokenizer.o stream_func.o print_color.o time_date.o filepath.o\

#----------------------------targets--------------------------------------------

INCLUDE := $(addprefix -I, $(INC_MECHOUI))

mechoui: mechoui.cc $(OBJ_MECHOUI:%=build/%) build/SFMT.o build/mechdisp.a
	$(COMPILE) $(INCLUDE) $^ $(GRAFIX) $(IMAGES) -o $@
	@printf "\x1B[35m >>>>>> made %s\x1B[0m\n" $@
vpath mechoui .

.PHONY: clean
clean:
	rm -f mechoui build/*

.PHONY: tar
tar: clean
	tar -cf mechoui.tar --dereference *

build:
	mkdir -p $@

#----------------------------rules----------------------------------------------

-include ../../gym/makefile.inc
-include ../../disp/makefile.inc

vpath spng.h ../../gym/deps
vpath spng.c ../../gym/deps
vpath miniz.h ../../gym/deps
vpath miniz.c ../../gym/deps

$(OBJ_MECHOUI:%.o=build/%.o): build/%.o: %.cc
	$(COMPILE) $(addprefix -I, $(INC_MECHOUI)) -c $< -o $@

# Attention: Mersenne-Twister is coded in C-language, hence '-x c':
build/SFMT.o: SFMT.c SFMT.h | build
	$(CXX) -DNDEBUG -DSFMT_MEXP=19937 -O2 -x c -c $< -o $@

build/mechdisp.a: $(addprefix build/, $(OBJ_GFX) view.o view_prop.o glapp.o) #save_image.o
	$(LIBTOOL) $@ $^
vpath mechdisp.a build

#----------------------------makedep--------------------------------------------

#command used to build the dependencies files automatically
MAKEDEP := g++ -std=gnu++14 -MM

build/dep: $(OBJ_MECHOUI:%.o=%.cc) | build
	(rm -f $@; for F in $^; do $(MAKEDEP) $$F >> $@; done)
	printf "> > > > > > > made %s\n" $@;
vpath dep build

-include build/dep
