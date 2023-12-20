# Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.


OBJ_BASE := messages.o filewrapper.o filepath.o iowrapper.o exceptions.o\
            time_date.o inventory.o stream_func.o tokenizer.o signal_handlers.o\
            glossary.o property.o property_list.o backtrace.o print_color.o\

#----------------------------rules----------------------------------------------

$(OBJ_BASE:%.o=build/%.o): build/%.o: %.cc %.h | build
	$(COMPILE) -Isrc/base -Isrc/math -c $< -o $@

build/operator_new.o: operator_new.cc | build
	$(COMPILE) -Isrc/base -Isrc/math -c $< -o $@

#----------------------------targets--------------------------------------------

cytobase.a: $(OBJ_BASE:%.o=build/%.o) build/operator_new.o | build lib
	$(LIBTOOL) lib/$@ $^
	$(DONE)

