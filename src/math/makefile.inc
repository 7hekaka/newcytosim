# Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

OBJ_MATH := vector1.o vector2.o vector3.o vector4.o random_vector.o\
    	matrix11.o matrix22.o matrix33.o matrix.o sparmat.o sparmatsym.o\
    	sparmatsym1.o sparmatsym2.o matrix34.o matfull.o\
        project_ellipse.o polygon.o rasterizer.o\
        spherical_code.o random.o random_seed.o random_pcg.o

OBJ_MATH_DIM := modulo.o sparmatblk.o sparmatsymblk.o sparmatsymblkdiag.o

OBJ_MATH_DIM0 := $(OBJ_MATH_DIM:%.o=build/%.o)
OBJ_MATH_DIM1 := $(OBJ_MATH_DIM:%.o=build/%D1.o)
OBJ_MATH_DIM2 := $(OBJ_MATH_DIM:%.o=build/%D2.o)
OBJ_MATH_DIM3 := $(OBJ_MATH_DIM:%.o=build/%D3.o)

OBJ_MATH_NODIM := $(addprefix build/, $(OBJ_MATH) SFMT.o)

#----------------------------rules----------------------------------------------

INC_BM:=-Isrc/base -Isrc/math


$(OBJ_MATH:%.o=build/%.o): build/%.o: %.cc %.h | build
	$(COMPILE) $(INC_BM) -c $< -o $@


$(OBJ_MATH_DIM0): build/%.o: %.cc %.h | build
	$(COMPILE) $(INC_BM) -c $< -o $@

$(OBJ_MATH_DIM1): build/%D1.o: %.cc %.h | build
	$(COMPILE) $(INC_BM) -DDIM=1 -c $< -o $@

$(OBJ_MATH_DIM2): build/%D2.o: %.cc %.h | build
	$(COMPILE) $(INC_BM) -DDIM=2 -c $< -o $@

$(OBJ_MATH_DIM3): build/%D3.o: %.cc %.h | build
	$(COMPILE) $(INC_BM) -DDIM=3 -c $< -o $@


cytomath.a: $(OBJ_MATH_NODIM) $(OBJ_MATH_DIM0) | build lib
	$(LIBTOOL) lib/$@ $^
	$(DONE)

cytomathD1.a: $(OBJ_MATH_NODIM) $(OBJ_MATH_DIM1) | build lib
	$(LIBTOOL) lib/$@ $^
	$(DONE)

cytomathD2.a: $(OBJ_MATH_NODIM) $(OBJ_MATH_DIM2) | build lib
	$(LIBTOOL) lib/$@ $^
	$(DONE)

cytomathD3.a: $(OBJ_MATH_NODIM) $(OBJ_MATH_DIM3) | build lib
	$(LIBTOOL) lib/$@ $^
	$(DONE)

