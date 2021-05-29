Cytosim Todo List

F. Nedelec

# Bugs

- check diffusion of a looped filament
- Display parameter such as 'tile': command line do not overwrite

# Essential

- Implement continuous attractive potential for 'steric interactions'
- Avoid discontinuity when steric attraction is on
- Calibrate end-to-end distance of semi-flexible filaments
- Calibrate rotational diffusion of filaments

# Important

- Fiber::read should allow changing 'activity'
- Geometrical primitives (eg. rectangle) should use diameter and not radius
- add units to documentation of parameters

# Improvements

- fiber:max_length could trigger switch to a different Dynamic state
- Could add color of Hands for summary obtained with 'I'
- Implement two confinements for Mecables. That's easier than constructive geometry on Spaces
- use 'cyo' as extension for ouput file

# Spindles

- We want the substance to accumulate at the poles:
  - saturation of chewing from the fiber tips
  - non-linear cutting from lattice?
  - saturating transport on lattice: convection stalls above a certain concentration


# v2018

- use 'hpp' for C++ header files
- use 'cpp' extension for C++ code file
- Update OpenGL code to use Shaders and eliminate immediate rendering

- Replace GLUT by GLFW:  http://www.glfw.org/index.html
- Rename gle: gym and glApp: gymapp
- use Dear imGUI for rendering text and minimal GUI components (Omar Cornut)
     https://github.com/ocornut/imgui

- new font renderer, eg. FreeType (http://www.freetype.org/index.html)
will fix problem for off-screen: glutBitmapCharacter() requires GLUT to be initialized.

# To reach v4

- Extend 'properties.cmo' to be able to store multiple time-points?
- Find a solution for 'confine = type, stiff, space' and 'confine = type, space' in Single.

- add generic classes that can be customized: OpenFiber, OpenHand, OpenCouple, OpenSingle...
- fix SpaceBeads: list of beads not set correctly after a read state.cym
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

Build:

- build with modern 'make', eg:
    scons: http://www.scons.org instead of make
    cmake: http://www.cmake.org
    meson: https://mesonbuild.com
- Replace GLUT:
    glfw:  http://www.glfw.org/index.html
- new font renderer, eg. FreeType (http://www.freetype.org/index.html)

# Performance:

- Store Mecable vertices as XXXXXXYYYYYYZZZZZZ to improve vectorization
- use VAS = Vector Array Stride; MVL = Meca Vector Length
- implement IDRS to replace BICGSTAB - convergence in matlab seems better
  http://ta.twi.tudelft.nl/nw/users/gijzen/IDR.html
- Use Pardiso matrix & solver
- multi-thread the steric engine
- try direct sparse method: HSL MA27 http://www.hsl.rl.ac.uk/index.html

# Examples:

  Would be great to revive some old models:
- yeast_karyogamy.cym
- spindle.cym: Rose's spindle model
- pombe_bundle.cym : bundle formation
  New:
- some DNA binding/unbinding

# Input/Output:

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
   - dirichlet boundary conditions: set concentration of molecules at edges of Space

- remove variadic functions (#include <stdarg>)
- option to pull on beads with mouse
- use sidePos() for display more systematically
- Crank-Nicholson diffusion in Field (Jonathan did this)
- class Event to control the value of a parameter in time
- include simple math evaluation tool, to be able to use '10*time+5' (muParser? Lua?)

# New molecular activities:

- implement specialized motor classes: Kinesin, Dynein & Myosin
- new Hand with capping activity (Request Marko Kaksonen)

- steric between beads and fibers, but not fiber-fiber
- stabilize the steric interactions of particles: solid.cym

- new Space in 3D defined as the intersection of an arbitrary number of half-spaces
- translated/rotated Space

- Glossary could record the queries and use this information for error reporting:
  the parameter 'confined' was unused, but a query was made for 'confine'

# Graphics:

- Update OpenGL code: do not use any direct rendering method
- Improve control of whether bound Hands are hidden if fiber->disp->visible==0
- use VBO to display each fiber with style=1 and 2
- SimThread can store two simulation states, such as to be able to flip back-and-forth quickly.
- stereoscopic display for 3D

# Misc:

- Slider to change parameter value:
   - setup by a command in config.cym: slider actin:rigidity { range = 0, 1; }


