# Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.


TOOLS:=frametool sieve reader report reportF analyse

.PHONY: tools
tools: $(TOOLS)

vpath %.cc src/tools

#--------------------macros----------------------------------------------------

TOOL_LIBS := cytolib.a cytomath.a cytobase.a

TOOL_INC = $(addprefix -Isrc/, math base sim sim/spaces sim/fibers sim/singles sim/organizers gym disp play)

TOOL_MAKE = $(COMPILE) $(TOOL_INC) $(OBJECTS) $(LINK) -o bin/$@


#-------------------targets----------------------------------------------------
 
 
frametool: frametool.cc | bin
	$(COMPILE) $^ -o bin/$@
	$(DONE)
vpath frametool bin


sieve: sieve.cc $(TOOL_LIBS) | bin
	$(TOOL_MAKE)
	$(DONE)
vpath sieve bin

bin1/sieve: sieve.cc cytolibD1.a cytomathD1.a cytobase.a | bin1
	$(COMPILE) -DDIM=1 $(TOOL_INC) $(OBJECTS) $(LINK) -o bin1/sieve
	$(DONE)
vpath bin1/sieve bin1

bin2/sieve: sieve.cc cytolibD2.a cytomathD2.a cytobase.a | bin2
	$(COMPILE) -DDIM=2 $(TOOL_INC) $(OBJECTS) $(LINK) -o bin2/sieve
	$(DONE)
vpath bin2/sieve bin2

bin3/sieve: sieve.cc cytolibD3.a cytomathD3.a cytobase.a | bin3
	$(COMPILE) -DDIM=3 $(TOOL_INC) $(OBJECTS) $(LINK) -o bin3/sieve
	$(DONE)
vpath bin3/sieve bin3


reader: reader.cc $(TOOL_LIBS) | bin
	$(TOOL_MAKE)
	$(DONE)
vpath reader bin


report: report.cc $(TOOL_LIBS) | bin
	$(COMPILE) $(TOOL_INC) $(OBJECTS) $(LINK) -o bin/$@
	$(DONE)
vpath report bin

bin1/report: report.cc cytolibD1.a cytomathD1.a cytobase.a | bin1
	$(COMPILE) -DDIM=1 $(TOOL_INC) $(OBJECTS) $(LINK) -o bin1/report
	$(DONE)
vpath bin1/report bin1

bin2/report: report.cc cytolibD2.a cytomathD2.a cytobase.a  | bin2
	$(COMPILE) -DDIM=2 $(TOOL_INC) $(OBJECTS) $(LINK) -o bin2/report
	$(DONE)
vpath bin2/report bin2

bin3/report: report.cc cytolibD3.a cytomathD3.a cytobase.a | bin3
	$(COMPILE) -DDIM=3 $(TOOL_INC) $(OBJECTS) $(LINK) -o bin3/report
	$(DONE)
vpath bin3/report bin3


reportF: reportF.cc $(TOOL_LIBS) | bin
	$(TOOL_MAKE)
	$(DONE)
vpath reportF bin


analyse: analyse.cc $(TOOL_LIBS) | bin
	$(TOOL_MAKE)
	$(DONE)
vpath analyse bin


cymart: cymart.cc $(TOOL_LIBS) | bin
	$(TOOL_MAKE)
	$(DONE)
vpath cymart bin


cymwiz: cymwiz.cc | bin
	$(TOOL_MAKE)
	$(DONE)
vpath cymwiz bin


micimage: micimage.cc save_image.cc save_image_gl.cc glossary.o filepath.o tokenizer.o
micimage: stream_func.o exceptions.o random.o SFMT.o print_color.o backtrace.o
	$(COMPILE) -Isrc/math -Isrc/base -Isrc/gym -DHAS_PNG $(OBJECTS) $(LIB_PNG) -o bin/$@
	$(DONE)
vpath micimage bin

