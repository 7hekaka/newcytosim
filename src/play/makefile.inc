# Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

OBJ_PLAY := player_prop.o display_prop.o display.o display1.o display2.o display3.o

FILES_PLAYER := player.cc player.h player_disp.cc

LIB_PLAYER := cytodisp.a $(OBJ_DISP:%.o=build/%.o) cytobase.a

PLAY_INC = $(addprefix -Isrc/, math base sim gym disp sim/organizers)
PLAY_DEP = $(addprefix -Isrc/, math base sim gym disp) $(LINK) $(GRAFIX) $(IMAGES)

################################ objects with DIM defined in source code:

build/player.o: $(FILES_PLAYER) | build
	$(COMPILE) $(PLAY_INC) -c $< -o $@

$(OBJ_PLAY:%.o=build/%.o): build/%.o: %.cc %.h | build
	$(COMPILE) $(PLAY_INC) -c $< -o $@

play: play.cc build/player.o $(OBJ_PLAY:%.o=build/%.o) cytolibG.a cytomath.a $(LIB_PLAYER) | bin
	$(COMPILE) $(INFO) $(OBJECTS) $(PLAY_DEP) -o bin/play
	$(DONE)
vpath play bin

################################ objects with DIM = 1:

build/playerD1.o: $(FILES_PLAYER) | build
	$(COMPILE) -DDIM=1 $(PLAY_INC) -c $< -o $@

$(OBJ_PLAY:%.o=build/%D1.o): build/%D1.o: %.cc %.h | build
	$(COMPILE) -DDIM=1 $(PLAY_INC) -c $< -o $@

bin1/play: play.cc build/playerD1.o $(OBJ_PLAY:%.o=build/%D1.o) cytolibGD1.a cytomathD1.a $(LIB_PLAYER) | bin1
	$(COMPILE) -DDIM=1 $(INFO) $(OBJECTS) $(PLAY_DEP) -o bin1/play
	$(DONE)

################################ objects with DIM = 2:

build/playerD2.o: $(FILES_PLAYER) | build
	$(COMPILE) -DDIM=2 $(PLAY_INC) -c $< -o $@

$(OBJ_PLAY:%.o=build/%D2.o): build/%D2.o: %.cc %.h | build
	$(COMPILE) -DDIM=2 $(PLAY_INC) -c $< -o $@

bin2/play: play.cc build/playerD2.o $(OBJ_PLAY:%.o=build/%D2.o) cytolibGD2.a cytomathD2.a $(LIB_PLAYER) | bin2
	$(COMPILE) -DDIM=2 $(INFO) $(OBJECTS) $(PLAY_DEP) -o bin2/play
	$(DONE)

################################ objects with DIM = 3:

build/playerD3.o: $(FILES_PLAYER) | build
	$(COMPILE) -DDIM=3 $(PLAY_INC) -c $< -o $@

$(OBJ_PLAY:%.o=build/%D3.o): build/%D3.o: %.cc %.h | build
	$(COMPILE) -DDIM=3 $(PLAY_INC) -c $< -o $@

bin3/play: play.cc build/playerD3.o $(OBJ_PLAY:%.o=build/%D3.o) cytolibGD3.a cytomathD3.a $(LIB_PLAYER) | bin3
	$(COMPILE) -DDIM=3 $(INFO) $(OBJECTS) $(PLAY_DEP) -o bin3/play
	$(DONE)

################################ monoplay & multiplay

INC_GLFW = $(addprefix -Isrc/, math base sim gym disp) -Iext

ifeq ($(MACHINE),mac)
    LIB_GLFW := ext/libglfw3.a -framework Cocoa -framework OpenGL -framework IOKit
else
    LIB_GLFW := -lglfw
endif

monoplay: src/play/monoplay.cc $(OBJ_PLAY:%.o=build/%.o) cytolibG.a $(LIB_PLAYER) cytomath.a | bin
	$(COMPILE) $(INFO) $(INC_GLFW) $(OBJECTS) $(LIB_GLFW) $(LINK) $(GRAFIX) $(IMAGES) -o bin/monoplay
	$(DONE)
vpath monoplay bin

multiplay: src/play/multiplay.cc $(OBJ_PLAY:%.o=build/%.o) cytolibG.a $(LIB_PLAYER) cytomath.a | bin
	$(COMPILE) $(INFO) $(INC_GLFW) $(OBJECTS) $(LIB_GLFW) $(LINK) $(GRAFIX) $(IMAGES) -o bin/multiplay
	$(DONE)
vpath multiplay bin
