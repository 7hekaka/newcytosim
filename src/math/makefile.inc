# Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.


OBJ_MATH := vector1.o vector2.o vector3.o matrix11.o matrix22.o matrix33.o\
    	rasterizer.o platonic.o matrix.o sparmat.o sparmatsym.o\
    	sparmatsym1.o polygon.o pointsonsphere.o random.o random_vector.o\
        matrix34.o matsym.o matfull.o project_ellipse.o quartic_solver.o

OBJ_SPARMAT := sparmatblk.o sparmatsymblk.o sparmatsymblkdiag.o

#----------------------------rules----------------------------------------------

INC_BM:=-Isrc/base -Isrc/math


$(OBJ_MATH): %.o: %.cc %.h | build
	$(COMPILE) $(INC_BM) -c $< -o build/$@


sparmatsym2.o: sparmatsym2.cc sparmatsym2.h
	$(COMPILE) $(INC_BM) -c $< -o build/$@


sparmatsymblk.o: sparmatsymblk.cc sparmatsymblk.h
	$(COMPILE) $(INC_BM) -c $< -o build/$@

sparmatsymblkD1.o: sparmatsymblk.cc sparmatsymblk.h
	$(COMPILE) $(INC_BM) -DDIM=1 -c $< -o build/$@

sparmatsymblkD2.o: sparmatsymblk.cc sparmatsymblk.h
	$(COMPILE) $(INC_BM) -DDIM=2 -c $< -o build/$@

sparmatsymblkD3.o: sparmatsymblk.cc sparmatsymblk.h
	$(COMPILE) $(INC_BM) -DDIM=3 -c $< -o build/$@


sparmatsymblkdiag.o: sparmatsymblkdiag.cc sparmatsymblkdiag.h
	$(COMPILE) $(INC_BM) -c $< -o build/$@

sparmatsymblkdiagD1.o: sparmatsymblkdiag.cc sparmatsymblkdiag.h
	$(COMPILE) $(INC_BM) -DDIM=1 -c $< -o build/$@

sparmatsymblkdiagD2.o: sparmatsymblkdiag.cc sparmatsymblkdiag.h
	$(COMPILE) $(INC_BM) -DDIM=2 -c $< -o build/$@

sparmatsymblkdiagD3.o: sparmatsymblkdiag.cc sparmatsymblkdiag.h
	$(COMPILE) $(INC_BM) -DDIM=3 -c $< -o build/$@


sparmatblk.o: sparmatblk.cc sparmatblk.h
	$(COMPILE) $(INC_BM) -c $< -o build/$@

sparmatblkD1.o: sparmatblk.cc sparmatblk.h
	$(COMPILE) $(INC_BM) -DDIM=1 -c $< -o build/$@

sparmatblkD2.o: sparmatblk.cc sparmatblk.h
	$(COMPILE) $(INC_BM) -DDIM=2 -c $< -o build/$@

sparmatblkD3.o: sparmatblk.cc sparmatblk.h
	$(COMPILE) $(INC_BM) -DDIM=3 -c $< -o build/$@


cytomath.a: $(OBJ_MATH) $(OBJ_SPARMAT) | lib
	$(MAKELIB)
	$(DONE)

cytomathD1.a: $(OBJ_MATH) $(OBJ_SPARMAT:%.o=%D1.o) | lib
	$(MAKELIB)
	$(DONE)

cytomathD2.a: $(OBJ_MATH) $(OBJ_SPARMAT:%.o=%D2.o) | lib
	$(MAKELIB)
	$(DONE)

cytomathD3.a: $(OBJ_MATH) $(OBJ_SPARMAT:%.o=%D3.o) | lib
	$(MAKELIB)
	$(DONE)

