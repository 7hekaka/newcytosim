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

#------------------------------files--------------------------------------------

SRC_MECHOUI := multimesh.cc mechoui_param.cc vector3.cc glossary.cc\
			   property.cc backtrace.cc tokenizer.cc stream_func.cc\
               time_date.cc filepath.cc random.cc print_color.cc

OBJ_MECHOUI := $(SRC_MECHOUI:%.cc=build/%.o) build/SFMT.o build/mechdisp.a

DIR_MECHOUI := src/mesh src/base src/math src/gym src/disp src/gym/deps

vpath %.h  src $(DIR_MECHOUI)
vpath %.cc src $(DIR_MECHOUI)
vpath %.c  $(DIR_MECHOUI)
vpath %.o  build
vpath %.a  build

#----------------------------targets--------------------------------------------

INC_MECHOUI := $(addprefix -I, $(DIR_MECHOUI))

mechoui: mechoui.cc $(OBJ_MECHOUI)
	$(COMPILE) $(INC_MECHOUI) $^ $(GRAFIX) $(IMAGES) -o $@
	@printf "\x1B[35m >>>>>> made %s\x1B[0m\n" $@
vpath mechoui .

.PHONY: clean
clean:
	rm -f mechoui build/* src/gym src/disp src/base src/math

.PHONY: tar
tar: clean
	tar -cf mechoui.tar --dereference *

build:
	mkdir -p $@

symlinks:
	if test ! -e src/gym; then ln -s ../../../gym src/gym; fi
	if test ! -e src/disp; then ln -s ../../../disp src/disp; fi
	if test ! -e src/base; then ln -s ../../../base src/base; fi
	if test ! -e src/math; then ln -s ../../../math src/math; fi
	if test ! -e src/mesh; then ln -s ../../../mesh src/mesh; fi

#----------------------------rules----------------------------------------------

-include src/gym/makefile.inc
-include src/disp/makefile.inc

$(SRC_MECHOUI:%.cc=build/%.o): build/%.o : %.cc | build symlinks
	$(COMPILE) $(INC_MECHOUI) -c $< -o $@

build/mechdisp.a: $(addprefix build/, $(OBJ_GRFX) view.o view_prop.o glapp.o)
	$(LIBTOOL) $@ $^
vpath mechdisp.a build

# Attention: Mersenne-Twister is coded in C-language, hence '-x c':
build/SFMT.o: SFMT.c SFMT.h | build symlinks
	$(CXX) -DNDEBUG -DSFMT_MEXP=19937 -O2 -x c -c $< -o $@

#----------------------------makedep--------------------------------------------

#command used to build the dependencies files automatically
MAKEDEP := g++ -std=gnu++14 -MM

build/dep: $(SRC_MECHOUI) | build
	(rm -f $@; for F in $^; do $(MAKEDEP) $$F >> $@; done)
	printf "> > > > > > > made %s\n" $@;
vpath dep build

-include build/dep
