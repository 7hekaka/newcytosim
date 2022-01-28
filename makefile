# Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

# THIS FILE SHOULD NOT BE EDITED
# ONLY EDIT FILE makefile.inc 

# disable all implicit rules:
.SUFFIXES:

#include the compiler-specifications
include makefile.inc

#-----------------check the compilers and the compiler-flags--------------------

ifndef CXX
   $(error No compiler defined for $$(MACHINE)=$(MACHINE))
endif

ifndef Flags$(MODE)
   $(warning No compiler options defined for $$(MACHINE)=$(MACHINE) in mode $(MODE))
endif

ifndef LINK
  $(error No linkage-options defined for $$(MACHINE)=$(MACHINE))
endif

# command to invoke compiler:
COMPILE := $(CXX) $(CXXFLG) $(Flags$(MODE))

# macro to make a library:
MAKELIB = $(LIBTOOL) lib/$@ $(addprefix build/, $(notdir $^))

# macro to fix the paths of objects:
OBJECTS = $(filter %.cc, $^) $(addprefix build/, $(notdir $(filter %.o, $^))) $(addprefix lib/, $(notdir $(filter %.a, $^)))

# macro to notify that a task was completed:
DONE = printf "> > > > > > > made %s\n" $@;


SRCDIR1 := $(addprefix src/, math base sim disp play)
SRCDIR2 := $(addprefix src/sim/, spaces hands fibers singles couples organizers)
SRCDIR  := $(SRCDIR1) $(SRCDIR2)


#-----------------------GIT revision number-------------------------------------

CODE_VERSION = $(shell git rev-parse --short HEAD || echo unknown)

INFO = -D'CODE_VERSION="$(CODE_VERSION)"'

#-------------------------make's search paths-----------------------------------

vpath %.h     $(SRCDIR)
vpath %.cc    $(SRCDIR)
vpath %.o     build/
vpath %.a     lib/
vpath %.dep   dep/
vpath SFMT%   src/math/

#----------------------------targets--------------------------------------------

# calling 'make' without arguments will make sim, play and report:

.PHONY: simplay
simplay: sim report play

info:
	@echo $(MACHINE)

include src/sim/makefile.inc
include src/base/makefile.inc
include src/math/makefile.inc
include src/disp/makefile.inc
include src/play/makefile.inc

-include src/tools/makefile.inc
-include src/test/makefile.inc


# Attention: Mersenne-Twister is coded in C-language,
# and we must tell the compiler with '-x c':
SFMT.o: SFMT.c SFMT.h
	$(CXX) -DNDEBUG -DSFMT_MEXP=19937 $(FAST) -x c -c $< -o build/$@


.PHONY: all dim1 dim2 dim3 alldim allsim doc

all: sim play tools

dim1: bin1/sim bin1/report bin1/play

dim2: bin2/sim bin2/report bin2/play

dim3: bin3/sim bin3/report bin3/play

alldim: dim1 dim2 dim3

allsim: bin1/sim bin2/sim bin3/sim

doc:
	if test -d doc/code/doxygen; then rm -rf doc/code/doxygen; fi
	mkdir doc/code/doxygen;
	doxygen doc/code/doxygen.cfg > log.txt 2> /dev/null

#------------------------------- archive ---------------------------------------
.PHONY: tar tarzip tarsrc pack

tar:
	COPYFILE_DISABLE=1 tar cf cytosim.tar --exclude "*.o" --exclude "*~" --exclude xcuserdata \
	--exclude doxygen src makefile makefile.inc cym python bash doc cytosim.xcodeproj


tarzip: tar
	rm -f cytosim.tar.bz2;
	bzip2 cytosim.tar;


tarsrc:
	COPYFILE_DISABLE=1 tar cf cytosim_src.tar --exclude "*.o" --exclude ".*" \
	--exclude ".svn" --exclude "*~" src makefile makefile.inc cym python bash

pack: sterile tarsrc

#---------------------------- maintenance --------------------------------------
.PHONY: bin build clean cleaner sterile

bin:
	if ! test -d bin; then mkdir bin; fi
	
bin1:
	if ! test -d bin1; then mkdir bin1; fi
	
bin2:
	if ! test -d bin2; then mkdir bin2; fi
	
bin3:
	if ! test -d bin3; then mkdir bin3; fi

build:
	if ! test -d build; then mkdir build; fi

lib:
	if ! test -d lib; then mkdir lib; fi

clean:
	rm -f build/*.o
	rm -f lib/*.a

cleaner:
	rm -f *.cmo build/*.o lib/*.a dep/*.dep;


sterile:
	rm -rf build
	rm -rf lib
	rm -rf dep
	rm -rf bin
	rm -rf bin1
	rm -rf bin2
	rm -rf bin3
	rm -f *.cmo
	rm -f log.txt;


#---------------------------- dependencies -------------------------------------

#command used to build the dependencies files automatically
MAKEDEP := g++ -std=gnu++14 -MM $(addprefix -I, $(SRCDIR))

.PHONY: dep
dep: $(addsuffix .dep, $(addprefix dep/part, 0 1 2 3 4 5 6 7))
	$(DONE)

dep/part0.dep: src/base/*.cc
	(for F in $^; do $(MAKEDEP) $$F >> dep/part0.dep; done)

dep/part1.dep: src/math/*.cc
	(for F in $^; do $(MAKEDEP) $$F >> dep/part1.dep; done)

dep/part2.dep: src/sim/*.cc
	(for F in $^; do $(MAKEDEP) $$F >> dep/part2.dep; done)

dep/part3.dep: src/sim/*/*.cc
	(for F in $^; do $(MAKEDEP) $$F >> dep/part3.dep; done)

dep/part4.dep: src/disp/*.cc
	(for F in $^; do $(MAKEDEP) $$F >> dep/part4.dep; done)

dep/part5.dep: src/play/*.cc
	(for F in $^; do $(MAKEDEP) $$F >> dep/part5.dep; done)

dep/part6.dep: src/tools/*.cc
	(for F in $^; do $(MAKEDEP) $$F >> dep/part6.dep; done)

dep/part7.dep: src/test/*.cc
	(for F in $^; do $(MAKEDEP) $$F >> dep/part7.dep; done)

-include dep/part?.dep

