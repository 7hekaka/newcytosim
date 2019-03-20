# Blanksim 0.1 - F. Nedelec, Copyright EMBL 2007-2015
# makefile  Francois Nedelec 2015.

#-------------------include compiler-specifications-----------------------------

include makefile.inc

#-------------------check compiler and compiler-flags---------------------------

ifeq ($CXX,)
   $(error No compiler defined for $$(MACHINE)=$(MACHINE))
endif

ifeq ($(Flags$(MODE)),)
   $(error No compiler-options defined for $$(MACHINE)=$(MACHINE) in mode $(MODE))
endif

CXXFLAGS=$(CXX) $(Flags$(MODE))

ifneq ($(HAS_PNG), 0)
    IMAGE_DEF = -DHAS_PNG $(INC_PNG)
    IMAGE_LIB = $(LIB_PNG)
endif

#----------------------------objects--------------------------------------------

INCDIR=$(addprefix ../, SFMT base math gl disp)

vpath %.h  $(INCDIR) .
vpath %.cc $(INCDIR) .
vpath %.c  $(INCDIR)
vpath %.o  build

OBJ=exceptions.o vector3.o glossary.o property.o backtrace.o offscreen.o\
    tokenizer.o stream_func.o ansi_colors.o view.o view_prop.o platonic.o\
    gle.o gle_color.o gle_color_list.o glapp.o tictoc.o filepath.o\
    mesh.o mechoui_param.o mechoui.o

#----------------------------targets--------------------------------------------

mechoui: $(OBJ) SFMT.o saveimage.o
	$(CXXFLAGS) $(addprefix build/, $(notdir $^)) $(LINKGL) $(IMAGE_LIB) -o $@
	@printf "\x1B[35m >>>>>> made %s\x1B[0m\n" $@
vpath mechoui .

.PHONY: clean
clean:
	rm -f mechoui build/*

.PHONY: tar
tar: clean
	tar -cf mechoui.tar --dereference * 

#----------------------------rules----------------------------------------------

INCLUDE=$(addprefix -I, $(INCDIR))

$(OBJ): %.o: %.cc
	$(CXXFLAGS) -c $(INCLUDE) $< -o build/$@

build/saveimage.o: saveimage.cc saveimage.h
	$(CXXFLAGS) -c $(IMAGE_DEF) $< -o $@

# Attention: Mersenne-Twister is coded in C-language,
# and we must provide this with '-x c' to the compiler:
SFMT.o: SFMT.c SFMT.h
	$(CXXFLAGS) -x c -DNDEBUG -DSFMT_MEXP=19937 -c $< -o build/$@

#----------------------------makedep--------------------------------------------

MAKEDEP=gcc -MM
SRC=$(patsubst %.o,%.cc, $(OBJ))
FILES=$(findfile $(SRC), $(INCDIR))

build/dependencies: $(SRC)
	$(MAKEDEP) $(INCLUDE) $(filter %.cc, $^) > $@

-include build/dependencies
