Cytosim Todo List

# Bugs

- check diffusion of a looped filament
- Display parameter such as 'tile': command line do not overwrite
- biased movements of short filaments with sterics

# Essential

- Implement continuous attractive potential for 'steric interactions'
- Avoid discontinuity when steric attraction is on
- Calibrate end-to-end distance of semi-flexible filaments
- Calibrate rotational diffusion of filaments

# Important

- Fiber::read should allow changing 'activity'
- add units to documentation of parameters
- Export Cytosim's world to Python (pybind, Serge).
- https://github.com/pybind/pybind11 seems better than Boost::Python 
- Replace PCB and MersenneTwister by XORISHIRO: https://en.wikipedia.org/wiki/Xoroshiro128%2B

# Improvements

- remove haMonitor, using function arguments instead?
- fiber:max_length could trigger switch to a different Dynamic state
- Could add color of Hands for summary obtained with 'I'
- Implement two confinements for Mecables. That's easier than constructive geometry on Spaces
- use 'cyo' as extension for ouput file
- use Sphynx for documentation:

# Spindles

- We want the substance to accumulate at the poles:
  - saturation of chewing from the fiber tips
  - non-linear cutting from lattice?
  - saturating transport on lattice: convection stalls above a certain concentration

# v2024

- use 'hpp' for C++ header files
- use 'cpp' extension for C++ code file
- Update OpenGL code to use Shaders and eliminate immediate rendering

- Replace GLUT by GLFW:  http://www.glfw.org/index.html
- Rename gle->gym and glApp->gympad
- use Dear imGUI for rendering text and minimal GUI components (Omar Cornut)
     https://github.com/ocornut/imgui

# To reach v4

- Extend the 'properties' file to be able to store multiple time-points?
- Find a solution for 'confine = type, stiff, space' and 'confine = type, space' in Single.

- add generic classes that can be customized: OpenFiber, OpenHand, OpenCouple, OpenSingle...
- binding_range warnings are not necessary for long Couples


Reorganize the cym folder:
  -> simple examples organised as a tutorial
  -> tests + reference
  -> research
test write / read systematically on all cym files

- Review the syntax of the config file

# Dissemination

- Create a Vanilla version of Cytosim: cytosim.core,
    where Hands, Single, Couple and Fiber are final (not derived) classes,
    and Fiber implements ClassicFiber, and Single implements Picket
- Video tutorial to install and run cytosim (pombe simulation)
- Video tutorial to do something simple (self-organization / contraction?)
- Processing-style single-window GUI to start Cytosim
- IKEA/LEGO style step-by-step manual with minimal text
- maintain Manu's syntax-color highligting file for cytosim
- create 2D tutorial on how to place objects
- export to X3D format, and install viewer within a browser

# Future work:

- Biped: motor with two attachment sites, walking by a hand-over-hand mechanism
- Represents an actin filament with two strands... a microtubule with 13.
     the strands can be in register, such that vertices describe circles in 3D, 
     while the abscissaM() can be shifted appropriately

- Use placement new to create thousands of Couple / Single
	This should avoid scattering of memory, and the slow execution of thousands of tiny allocate

# Performance:

- reading Single can be pooled since the Unattached and Attached Singles are stored apart
- reduce Solid::addTriad() to add 2 points in 2D and 3 in 3D, as this is sufficient
- Solid degree of freedom in linear system reduced to 5 in 3D, as this is sufficient
- SparMatMatSym Matrices can be structured as blocks corresponding to each Mecable (1.9.2022).
	+ this would speed up block addressing: 
	+ Need to create a tuple of indices {mecaindex, pointindex}, to use as argument of add_block()
	+ Updating the matrix can be done in parallel if Mecablocks are independent
	+ column reordering by increasing index would only be needed within each Mecablock,
	+ this will also permit vecMul parallelization per diagonal, using two output vectors:
        each sub-diagonal of blocks can be treated in parallel using multiple threads
	+ todo: vecMul() takes two output vectors; SparmatSymBlock becomes rectangular
- parallel matrix-vector using SparMatMatSym with two output vector
- new Matrix33x2 two implement two consecutive blocks:
    consecutive 3x3 blocks can be handled simultaneously with SIMD vec2 operations
- Store Mecable vertices as XXXXXX...YYYYYY...ZZZZZZ... to improve vectorization?
- use VAS = Vector Array Stride; MVL = Meca Vector Length
- implement IDRS to replace BICGSTAB - convergence in matlab seems better
  http://ta.twi.tudelft.nl/nw/users/gijzen/IDR.html
- Use Pardiso matrix & solver
- multi-thread the steric engine
- try direct sparse method: HSL MA27 http://www.hsl.rl.ac.uk/index.html

- Implement Reverse CutHill-McKee reordering method in Simul::orderClustersCouple()
http://ciprian-zavoianu.blogspot.com/2009/01/project-bandwidth-reduction.html

- Implement BlockOuterProduct to speed up LongLink: (Idea 06.02.2022)
each LongLink element currently covers four 3x3 blocks = 36 scalars + diagonal elements
We shall keep the diagonal elements, but compress the off-diagonal ones:
        add_block(ii2, ii0, cc2*cc0, wT);
        add_block(ii3, ii0, cc3*cc0, wT);
        add_block(ii2, ii1, cc2*cc1, wT);
        add_block(ii3, ii1, cc3*cc1, wT);
Since wT = MatrixBlock::offsetOuterProduct(), we have wT = weight [ d.1 - A (x) A ] 
with d = diagonal scalar, A = axis vector
and with interpolation coefficients ( a, b ) the 4 block are: a.b.wT, (1-a)b.wT, a.(1-b).wT, (1-a)(1-b)wT
Each compressed element shall thus store: (a, b, d, A) since weigth can be distributed to d and A.
This reduces the loading from 36 scalars down to 6 scalars, and the caculation
        wT.V = [ d.1 - A (x) A ] V =  d.V + dot(A, V). A 
        We only need to calculate this once for A = a . PT[i] + (1-a) . PT[i+1]
        And redistribute to the two vectors: PT[j] += b.wT.A and PT[j+1] += (1-b).wT.A
This can be done on the fly in mFUL.vecMul(), if we add a list of BlockOuterProduct elements to the matrix.
This can speed up multiplication, to be tested. 
The approach can be extended to others elements used extensively: SideSlidingLinks, SideLinks, etc.


# Reporting:

Export Cytosim's simulation world to Python

https://www.boost.org/doc/libs/1_78_0/libs/python/doc/html/index.html
https://www.boost.org/doc/libs/1_76_0/libs/python/doc/html/tutorial/index.html
https://github.com/TNG/boost-python-examples


# Examples:

  Would be great to revive some old models:
- yeast_karyogamy.cym
- spindle.cym: Rose's spindle model
- pombe_bundle.cym : bundle formation
  New:
- some DNA binding/unbinding

# Input/Output:

- save all float/double to a different file 'vectors.cmo'
  objects in 'objects.cmo' will have a index to the array of vectors
  this should be much faster to read!
- reorganize the Input classes IOWrapper
- include error reporting inside Input class
   - collect errors in local string
   - clearErrors() getErrors()


# New Developments:

- Steric interaction limited to nearest neighbors, instead of a fixed cut-off distance
- Remove discontinuity in the piece-wise linear force (at steric_range + steric_radius)
- steric_range specified directly, and not 'extra' range. Test that steric_range > steric_radius
-
- A level-set based Space to import images (2D sufficient initially)
- Steric using analytical formula for potential, integrated over segment length
- a Space to represent a lamellipodium: like a capsule in Y&Z, periodic in X
- extensions of fast_diffusion:
   - Dirichlet boundary conditions: set concentration of molecules at edges of Space

- option to pull on beads with mouse
- use sidePos() for display more systematically
- Crank-Nicholson diffusion in Field (Jonathan did this)
- class Event to control the value of a parameter in time
- include simple math evaluation tool, to be able to use '10*time+5' (muParser? Lua?)

# New molecular activities:

- implement specialized motor classes: Kinesin, Dynein & Myosin

- steric between beads and fibers, but not fiber-fiber
- stabilize the steric interactions of particles: solid.cym

- new Space in 3D defined as the intersection of an arbitrary number of half-spaces
- translated/rotated Space

- Glossary could record the queries and use this information for error reporting:
  the parameter 'confined' was unused, but a query was made for 'confine'

# Graphics:

- Update Cytosim to use OpenGL ES 2.0, which can be run on METAL using Google's ANGLE
- Improve control of whether bound Hands are hidden if fiber->disp->visible==0
- SimThread could store two simulation states, such as to be able to flip back-and-forth quickly.
- stereoscopic display for 3D

# Misc:

- FlightSimulator to control Cytosim's display via pipe
- Slider to change parameter value:
   - setup by a command in config.cym: slider actin:rigidity { range = 0, 1; }
   can be done with ImGui?


FJ Nedelec, 2023
