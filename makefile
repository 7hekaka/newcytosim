# This is the master makefile for Cytosim. Copyright 2020 Cambridge University

# IN MOST CASES, THIS FILE SHOULD NOT BE EDITED
# Compilation can be customized by editing `makefile.inc'.

# disable all implicit rules:
.SUFFIXES:

#include the compiler-specifications
include makefile.inc

#-----------------check the compilers and the compiler-flags--------------------

ifndef CXX
   $(error No compiler defined for MACHINE=$(MACHINE))
endif

ifndef Flags$(MODE)
   $(warning No compiler options defined for MACHINE=$(MACHINE) in mode $(MODE))
endif

ifndef LINK
  $(error No linkage-options defined for MACHINE=$(MACHINE))
endif

# command to invoke compiler:
COMPILE := $(CXX) $(CXXFLG) $(Flags$(MODE))

# macro to fix the paths of objects:
OBJECTS = $(filter %.cc, $^) $(addprefix build/, $(notdir $(filter %.o, $^))) $(addprefix lib/, $(notdir $(filter %.a, $^)))

# macro to notify that a task was completed:
DONE = printf "> > > > > > > made %s\n" $@;


SRCDIR1 := $(addprefix src/, math base sim gym disp play)
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

# local copy of Mersenne Twister:
vpath SFMT.h src/math/
vpath SFMT.c src/math/

#----------------------------targets--------------------------------------------

# calling 'make' without arguments will make sim, play and report:

.PHONY: cytosim
cytosim: sim report play

info:
	@echo $(MACHINE)

include src/sim/makefile.inc
include src/base/makefile.inc
include src/math/makefile.inc
include src/gym/makefile.inc
include src/disp/makefile.inc
include src/play/makefile.inc

-include src/tools/makefile.inc
-include src/test/makefile.inc


# Attention: Mersenne-Twister is coded in C-language,
# and we must tell the compiler with '-x c':
build/SFMT.o: SFMT.c SFMT.h
	$(CXX) -DNDEBUG -DSFMT_MEXP=19937 $(FAST) -x c -c $< -o $@


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
	--exclude doxygen src makefile makefile.inc cym python doc cytosim.xcodeproj


tarzip: tar
	rm -f cytosim.tar.bz2;
	bzip2 cytosim.tar;


tarsrc:
	COPYFILE_DISABLE=1 tar cf cytosim_src.tar --exclude "*.o" --exclude ".*" \
	--exclude ".svn" --exclude "*~" src makefile makefile.inc cym python

pack: sterile tarsrc

#---------------------------- maintenance --------------------------------------

bin:
	mkdir -p $@
	
bin1:
	mkdir -p $@
	
bin2:
	mkdir -p $@
	
bin3:
	mkdir -p $@

build:
	mkdir -p $@

lib:
	mkdir -p $@

.PHONY: clean cleaner sterile

clean:
	mkdir -p dep
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
MAKEDEP := g++ -std=c++17 -MM $(addprefix -I, $(SRCDIR))

depdir:
	mkdir -p dep

dep: depdir $(addsuffix .dep, $(addprefix dep/, 0 1 2 3 4 5 6 7 8))
	$(DONE)

dep/0.dep: src/base/*.cc
	(rm -f $@; for F in $^; do $(MAKEDEP) $$F >> $@; done)

dep/1.dep: src/math/*.cc
	(rm -f $@; for F in $^; do $(MAKEDEP) $$F >> $@; done)

dep/2.dep: src/sim/*.cc
	(rm -f $@; for F in $^; do $(MAKEDEP) $$F >> $@; done)

dep/3.dep: src/sim/*/*.cc
	(rm -f $@; for F in $^; do $(MAKEDEP) $$F >> $@; done)

dep/4.dep: src/gym/*.cc
	(rm -f $@; for F in $^; do $(MAKEDEP) $$F >> $@; done)

dep/5.dep: src/disp/*.cc
	(rm -f $@; for F in $^; do $(MAKEDEP) $$F >> $@; done)

dep/6.dep: src/play/*.cc
	(rm -f $@; for F in $^; do $(MAKEDEP) $$F >> $@; done)

dep/7.dep: src/tools/*.cc
	(rm -f $@; for F in $^; do $(MAKEDEP) $$F >> $@; done)

dep/8.dep: src/test/*.cc
	(rm -f $@; for F in $^; do $(MAKEDEP) $$F >> $@; done)

-include dep/?.dep

