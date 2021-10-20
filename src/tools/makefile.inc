# Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.


TOOLS:=frametool sieve reader report reportF analyse sweet

.PHONY: tools
tools: $(TOOLS)

vpath %.cc src/tools

#--------------------macros----------------------------------------------------

TOOL_OBJ := cytolib.a cytomath.a sparmatsymblk.o cytobase.a SFMT.o

TOOL_INC  = $(addprefix -Isrc/, math base sim sim/spaces sim/fibers sim/singles sim/organizers disp play)

TOOL_MAKE = $(COMPILE) $(TOOL_INC) $(OBJECTS) $(LINK) -o bin/$@


#-------------------targets----------------------------------------------------
 
 
frametool: frametool.cc | bin
	$(COMPILE) $^ -o bin/$@
	$(DONE)
vpath frametool bin


sieve: sieve.cc frame_reader.o $(TOOL_OBJ) | bin
	$(TOOL_MAKE)
	$(DONE)
vpath sieve bin

bin1/sieve: sieve.cc frame_reader.o cytolibD1.a cytomathD1.a cytobase.a SFMT.o | bin1
	$(COMPILE) -DDIM=1 $(TOOL_INC) $(OBJECTS) $(LINK) -o bin1/sieve
	$(DONE)
vpath bin1/sieve bin1

bin2/sieve: sieve.cc frame_reader.o cytolibD2.a cytomathD2.a cytobase.a SFMT.o | bin2
	$(COMPILE) -DDIM=2 $(TOOL_INC) $(OBJECTS) $(LINK) -o bin2/sieve
	$(DONE)
vpath bin2/sieve bin2

bin3/sieve: sieve.cc frame_reader.o cytolibD3.a cytomathD3.a cytobase.a SFMT.o | bin3
	$(COMPILE) -DDIM=3 $(TOOL_INC) $(OBJECTS) $(LINK) -o bin3/sieve
	$(DONE)
vpath bin3/sieve bin3


reader: reader.cc frame_reader.o $(TOOL_OBJ) | bin
	$(TOOL_MAKE)
	$(DONE)
vpath reader bin

cympack: cympack.cc frame_reader.o $(TOOL_OBJ) | bin
	$(TOOL_MAKE)
	$(DONE)
vpath cympack bin


report: report.cc frame_reader.o $(TOOL_OBJ) | bin
	$(COMPILE) $(TOOL_INC) $(OBJECTS) $(LINK) -o bin/$@
	$(DONE)
vpath report bin

bin1/report: report.cc frame_reader.o cytolibD1.a cytomathD1.a cytobase.a SFMT.o | bin1
	$(COMPILE) -DDIM=1 $(TOOL_INC) $(OBJECTS) $(LINK) -o bin1/report
	$(DONE)
vpath bin1/report bin1

bin2/report: report.cc frame_reader.o cytolibD2.a cytomathD2.a cytobase.a SFMT.o | bin2
	$(COMPILE) -DDIM=2 $(TOOL_INC) $(OBJECTS) $(LINK) -o bin2/report
	$(DONE)
vpath bin2/report bin2

bin3/report: report.cc frame_reader.o cytolibD3.a cytomathD3.a cytobase.a SFMT.o | bin3
	$(COMPILE) -DDIM=3 $(TOOL_INC) $(OBJECTS) $(LINK) -o bin3/report
	$(DONE)
vpath bin3/report bin3


reportF: reportF.cc frame_reader.o $(TOOL_OBJ) | bin
	$(TOOL_MAKE)
	$(DONE)
vpath reportF bin


analyse: analyse.cc frame_reader.o $(TOOL_OBJ) | bin
	$(TOOL_MAKE)
	$(DONE)
vpath analyse bin


cymart: cymart.cc frame_reader.o $(TOOL_OBJ) | bin
	$(TOOL_MAKE)
	$(DONE)
vpath cymart bin


sweet: sweet.cc frame_reader.o $(TOOL_OBJ) | bin
	$(TOOL_MAKE)
	$(DONE)
vpath sweet bin


micimage: micimage.cc save_image.cc glossary.o filepath.o tokenizer.o
micimage: stream_func.o exceptions.o random.o SFMT.o print_color.o backtrace.o
	$(COMPILE) -Isrc/math -Isrc/base -Isrc/disp -DHAS_PNG $(OBJECTS) $(LIB_PNG) -o bin/$@
	$(DONE)
vpath micimage bin


#----------------------------makedep--------------------------------------------

dep/part7.dep:
	$(MAKEDEP) $(wildcard src/tools/*.cc) > $@

