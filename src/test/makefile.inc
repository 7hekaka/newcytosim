# Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.


TESTS:=test_bits test_gillespie test_random test_math test_glossary\
       test_xpttrs test_xtbsv test_quaternion test_rigidity test_project\
       test_orthonormal test_matrix test_thread test_blas test_solver test_dot test_simd

TESTS_GL:=test_solve test_opengl test_vbo test_glut test_glapp test_tesselator\
          test_rasterizer test_space test_grid test_sphere test_pipe

.PHONY: test_gl
tests_gl: $(TESTS_GL)

.PHONY: tests
tests: $(TESTS) $(TESTS_GL)

vpath %.cc src/test


#--------------------macros----------------------------------------------------

GLTEST_MAKE = $(COMPILE) $(addprefix -Isrc/, math base sim gym disp) $(OBJECTS) $(LINK) $(GRAFIX) $(IMAGES) -o bin/$@

TEST_OBJ = $(addprefix build/, random_vector.o random.o SFMT.o backtrace.o exceptions.o)

#----------------------------Targets--------------------------------------------

vpath test_% bin

test_code: test_code.cc $(TEST_OBJ) | bin
	$(COMPILE) -Isrc/base -Isrc/math $(OBJECTS) $(LINK) -o bin/$@
	$(DONE)

test_bits: test_bits.cc $(TEST_OBJ) | bin
	$(COMPILE) -Isrc/base -Isrc/math $(OBJECTS) $(LINK) -o bin/$@
	$(DONE)

test_file: test_file.cc $(TEST_OBJ) | bin
	$(COMPILE) -Isrc/base -Isrc/math $(OBJECTS) $(LINK) -o bin/$@
	$(DONE)

test_orthonormal: test_orthonormal.cc $(TEST_OBJ) | bin
	$(COMPILE) -Isrc/base -Isrc/math $(OBJECTS) $(LINK) -o bin/$@
	$(DONE)

test_array: test_array.cc $(TEST_OBJ) | bin
	$(COMPILE) -Isrc/base -Isrc/math $(OBJECTS) $(LINK) -o bin/$@
	$(DONE)

test_blas: test_blas.cc $(TEST_OBJ) | bin
	$(COMPILE) -Isrc/base -Isrc/math $(OBJECTS) $(LINK) -o bin/$@
	$(DONE)

test_project: test_project.cc $(TEST_OBJ) | bin
	$(COMPILE) -Isrc/base -Isrc/math $(OBJECTS) $(LINK) -o bin/$@
	$(DONE)

test_xpttrs: xpttrf.h test_xpttrs.cc $(TEST_OBJ) | bin
	$(COMPILE) -Isrc/base -Isrc/math $(OBJECTS) $(LINK) -o bin/$@
	$(DONE)

test_xtbsv: xtbsv.h test_xtbsv.cc $(TEST_OBJ) | bin
	$(COMPILE) -Isrc/base -Isrc/math $(OBJECTS) $(LINK) -o bin/$@
	$(DONE)

test_rigidity: test_rigidity.cc $(TEST_OBJ) | bin
	$(COMPILE) -Isrc/base -Isrc/math $(OBJECTS) $(LINK) -o bin/$@
	$(DONE)

test_size: test_size.cc cytolib.a cytomath.a cytobase.a | bin
	$(COMPILE) -Isrc/math -Isrc/base -Isrc/sim $(OBJECTS) $(LINK) -o bin/$@
	$(DONE)

test_block: test_block.cc matrix33.o matrix34.o vector3.o vector4.o matfull.o $(TEST_OBJ) | bin
	$(COMPILE) -Isrc/base -Isrc/math $(OBJECTS) $(LINK) $(INFO) -o bin/$@
	$(DONE)

test_matrix: test_matrix.cc sparmatsymblkdiag.o sparmatsymblk.o sparmatblk.o sparmatsym.o sparmatsym1.o sparmatsym2.o matrix.o $(TEST_OBJ) operator_new.o | bin
	$(COMPILE) -Isrc/base -Isrc/math $(OBJECTS) $(OpenMP) $(OMPLIB) $(LINK) $(INFO) -o bin/$@
	$(DONE)

test_glossary: test_glossary.cc glossary.o filepath.o tokenizer.o stream_func.o exceptions.o backtrace.o print_color.o | bin
	$(COMPILE) -Isrc/base -Isrc/math $(OBJECTS) $(LINK) -o bin/$@
	$(DONE)

test_quartic: test_quartic.cc quartic_solver.o $(TEST_OBJ) | bin
	$(COMPILE) -Isrc/base -Isrc/math $(OBJECTS) $(LINK) -o bin/$@
	$(DONE)

test_quaternion: test_quaternion.cc matrix33.o $(TEST_OBJ) | bin
	$(COMPILE) -Isrc/base -Isrc/math $(OBJECTS) $(LINK) -o bin/$@
	$(DONE)

test_random: test_random.cc $(TEST_OBJ) filewrapper.o messages.o | bin
	$(COMPILE) -Isrc/base -Isrc/math $(OBJECTS) $(LINK) -o bin/$@
	$(DONE)

test_gaussian: test_gaussian.cc $(TEST_OBJ) filewrapper.o messages.o | bin
	$(COMPILE) -Isrc/base -Isrc/math $(OBJECTS) $(LINK) -o bin/$@
	$(DONE)

test_variablerate: test_variablerate.cc $(TEST_OBJ) | bin
	$(COMPILE) -Isrc/base -Isrc/math $(OBJECTS) $(LINK) -o bin/$@
	$(DONE)

test_asm: test_asm.cc | bin
	$(COMPILE) -restrict -S -c -m64 $^ -o build/test_asm.s
	$(DONE)

test_dot: test_dot.cc $(TEST_OBJ) | bin
	$(COMPILE) -Isrc/base -Isrc/math $(OBJECTS) -o bin/$@
	$(DONE)

test_simd: test_simd.cc $(TEST_OBJ) | bin
	$(COMPILE) -Isrc/base -Isrc/math $(OBJECTS) -o bin/$@
	$(DONE)

test_omp: test_omp.cc | bin
	$(COMPILE) $(OpenMP) $(OMPLIB) $^ -o bin/test_omp
	$(DONE)

test_dispatch: test_dispatch.cc | bin
	$(COMPILE) $^ -o bin/test_dispatch
	$(DONE)

test_sizeof: test_sizeof.cc | bin
	$(COMPILE) $(OBJECTS) $(LINK) -o bin/$@
	$(DONE)

test_pipe: test_pipe.cc gym_color.cc gym_color_list.cc $(TEST_OBJ) | bin
	$(COMPILE) -Isrc/base -Isrc/math -Isrc/gym $(OBJECTS) $(LINK) -o bin/$@
	$(DONE)

test_cxx: test_cxx.cc | bin
	$(COMPILE) $(OBJECTS) $(LINK) -o bin/$@
	$(DONE)

test_math: test_math.cc | bin
	$(COMPILE) -Isrc/math $(OBJECTS) $(LINK) -o bin/$@
	$(DONE)

test_thread: test_thread.cc | bin
	$(COMPILE) $(OBJECTS) $(LINK) -o bin/$@
	$(DONE)

test_string: test_string.cc | bin
	$(COMPILE) $(OBJECTS) $(LINK) -o bin/$@
	$(DONE)

test_gillespie: test_gillespie.cc $(TEST_OBJ) | bin
	$(COMPILE) -Isrc/base -Isrc/math $(OBJECTS) -o bin/$@
	$(DONE)

test_solver: test_solver.cc | bin
	$(COMPILE) -Isrc/base -Isrc/math $(OBJECTS) $(LINK) -o bin/$@
	$(DONE)

#----------------------------Targets with Graphics------------------------------

test_opengl: test_opengl.cc | bin
	$(CXX) $^ $(GRAFIX) -o bin/$@
	$(DONE)

test_vbo: test_vbo.cc | bin
	$(CXX) $^ $(GRAFIX) -o bin/$@
	$(DONE)

test_glut: test_glut.cc | bin
	$(GLTEST_MAKE)
	$(DONE)

test_solve: test_solve.cc cytomath.a cytobase.a gym_matrix.o | bin
	$(GLTEST_MAKE)
	$(DONE)

#----------------------------Targets using GLAPP--------------------------------

OBJ_GLAP := $(OBJ_DISP:%.o=build/%.o) $(OBJ_GRFX) cytobase.a cytomath.a build/SFMT.o

test_glapp: test_glapp.cc $(OBJ_GLAP) | bin
	$(GLTEST_MAKE)
	$(DONE)

test_sphere: test_sphere.cc $(OBJ_GLAP) | bin
	$(GLTEST_MAKE)
	$(DONE)

test_tesselator: test_tesselator.cc tesselator.o $(OBJ_GLAP) | bin
	$(GLTEST_MAKE)
	$(DONE)


# build the rasterizer with graphical debugging enabled:
build/rasterizerG.o: rasterizer.cc rasterizer.h | build
	$(COMPILE) -Isrc/base -Isrc/math -Isrc/gym -DDISPLAY -c $< -o $@

test_rasterizer: test_rasterizer.cc rasterizerG.o $(OBJ_GLAP) | bin
	$(GLTEST_MAKE)
	$(DONE)

test_space: test_space.cc cytolibG.a $(OBJ_GLAP) | bin
	$(GLTEST_MAKE)
	$(DONE)

test_grid: test_grid.cc $(OBJ_GLAP) grid_display.o | bin
	$(GLTEST_MAKE)
	$(DONE)
