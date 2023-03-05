# Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

OBJ_MATH := vector1.o vector2.o vector3.o vector4.o random_vector.o\
    	matrix11.o matrix22.o matrix33.o matrix.o sparmat.o sparmatsym.o\
    	sparmatsym1.o sparmatsym2.o matrix34.o matfull.o\
        project_ellipse.o polygon.o quartic_solver.o\
        rasterizer.o spherical_code.o random.o modulo.o

OBJ_SPARMAT := sparmatblk.o sparmatsymblk.o sparmatsymblkdiag.o

#----------------------------rules----------------------------------------------

INC_BM:=-Isrc/base -Isrc/math


$(OBJ_MATH:%.o=build/%.o): build/%.o: %.cc %.h | build
	$(COMPILE) $(INC_BM) -c $< -o $@



build/sparmatsymblk.o: sparmatsymblk.cc sparmatsymblk.h
	$(COMPILE) $(INC_BM) -c $< -o $@

build/sparmatsymblkD1.o: sparmatsymblk.cc sparmatsymblk.h
	$(COMPILE) $(INC_BM) -DDIM=1 -c $< -o $@

build/sparmatsymblkD2.o: sparmatsymblk.cc sparmatsymblk.h
	$(COMPILE) $(INC_BM) -DDIM=2 -c $< -o $@

build/sparmatsymblkD3.o: sparmatsymblk.cc sparmatsymblk.h
	$(COMPILE) $(INC_BM) -DDIM=3 -c $< -o $@


build/sparmatsymblkdiag.o: sparmatsymblkdiag.cc sparmatsymblkdiag.h
	$(COMPILE) $(INC_BM) -c $< -o $@

build/sparmatsymblkdiagD1.o: sparmatsymblkdiag.cc sparmatsymblkdiag.h
	$(COMPILE) $(INC_BM) -DDIM=1 -c $< -o $@

build/sparmatsymblkdiagD2.o: sparmatsymblkdiag.cc sparmatsymblkdiag.h
	$(COMPILE) $(INC_BM) -DDIM=2 -c $< -o $@

build/sparmatsymblkdiagD3.o: sparmatsymblkdiag.cc sparmatsymblkdiag.h
	$(COMPILE) $(INC_BM) -DDIM=3 -c $< -o $@


build/sparmatblk.o: sparmatblk.cc sparmatblk.h
	$(COMPILE) $(INC_BM) -c $< -o $@

build/sparmatblkD1.o: sparmatblk.cc sparmatblk.h
	$(COMPILE) $(INC_BM) -DDIM=1 -c $< -o $@

build/sparmatblkD2.o: sparmatblk.cc sparmatblk.h
	$(COMPILE) $(INC_BM) -DDIM=2 -c $< -o $@

build/sparmatblkD3.o: sparmatblk.cc sparmatblk.h
	$(COMPILE) $(INC_BM) -DDIM=3 -c $< -o $@


cytomath.a: $(addprefix build/, $(OBJ_MATH) $(OBJ_SPARMAT) SFMT.o ) | lib
	$(LIBTOOL) lib/$@ $^
	$(DONE)

cytomathD1.a: $(addprefix build/, $(OBJ_MATH) SFMT.o) $(OBJ_SPARMAT:%.o=build/%D1.o) | lib
	$(LIBTOOL) lib/$@ $^
	$(DONE)

cytomathD2.a: $(addprefix build/, $(OBJ_MATH) SFMT.o) $(OBJ_SPARMAT:%.o=build/%D2.o) | lib
	$(LIBTOOL) lib/$@ $^
	$(DONE)

cytomathD3.a: $(addprefix build/, $(OBJ_MATH) SFMT.o) $(OBJ_SPARMAT:%.o=build/%D3.o) | lib
	$(LIBTOOL) lib/$@ $^
	$(DONE)

