// Cytosim was created by Francois Nedelec. Copyright Cambridge University 2020
#ifndef DISPLAY_H
#define DISPLAY_H

#include "real.h"
#include "array.h"
#include "fiber.h"
#include "mecable.h"
#include "mecapoint.h"
#include "gym_color.h"
#include "display_prop.h"

class Simul;
class Mecable;
class SingleSet;
class CoupleSet;
class Couple;
class FiberSet;
class Solid;
class SolidSet;
class Organizer;
class OrganizerSet;
class Space;
class SpaceSet;
class Sphere;
class SphereSet;
class Bead;
class BeadSet;
class FieldSet;
class FiberProp;
class PropertyList;
class PointDisp;


/// defining the DISPLAY keyword enables display code in included files
#define DISPLAY 1

/**
 @brief A display element with a depth coordinate
 
 zObject is used to depth-sort and display transparent elements from back to front
 An element is a segment of a Fiber, or a sphere component from Solid/Bead/Sphere
 */
class zObject
{
    /// pointer to object
    Mecapoint point_;
    
    /// distance to the imaging plane
    real depth_;

public:
    
    /// constructor
    zObject() : depth_(0) { }

    /// constructor
    zObject(Mecable const* m, index_t i = 0) : point_(m, i), depth_(0) { }
    
    /// position
    Vector position() const { return point_.pos(); }
    
    /// query depth
    real depth() const { return depth_; }
    
    /// set depth relative to given axis
    void calculate_depth(Vector const& axis);
    
    /// display object
    void draw(Display const*) const;
};


///Base class to display Cytosim state using OpenGL
class Display
{
protected:
    
    /// array of transparent objects to be displayed last
    Array<zObject> zObjects;
    
    /// the size covered by one pixel in natural units (ie. micrometers)
    float pixelSize;
    
    /// number of pixels corresponding to one unit value of the size & width
    float unitValue;
    
    /// natural distance corresponding to one unit value of the size & width
    float sizeScale;
    
    /// direction of view
    Vector3 depthAxis;

    /// indicate which analysis should be done, that is used in the display
    size_t prep_flag;
    
    /// min and max age used to adjust color range with COLORING_AGE
    double age_scale, age_start;
    
    /// use OpenGL stencil test:
    bool stencil_;
    
    /// pointer to array of LineDisp
    LineDisp * allLineDisp;
    
    /// number  of allocated LineDisp
    size_t numLineDisp;
    
private:
    
    /// set default value of FiberProp
    void initFiberDisp(FiberProp*, PropertyList&, gym_color);
    
    /// set values of fiber's LineDisp
    void initLineDisp(Fiber const*, FiberDisp const*, LineDisp*);
    
    template < typename T >
    void initPointDisp(T * prop, PropertyList&, gym_color);
    
    /// draw translucent objects after depth-sorting
    void drawTransparentObjects(Array<zObject>&);
    
public:
    
    /// associated parameters
    DisplayProp const* prop;
    
    /// constructor
    Display(DisplayProp const*);
    
    /// virtual destructor needed, as class is base to others
    virtual ~Display();
    
    /// display opaque objects
    virtual void drawObjects(Simul const&);
    
    /// draw translucent objects after depth-sorting
    void drawTransparentObjects(Simul const&);

    /// enable/disable stencil usage
    void setStencil(bool s) { stencil_ = s; }

    /// set current pixel-size and the value of the point in pixels
    void setParameters(float pixel_size, float unit_value, Vector3 const& a);

    /// attribute a LineDisp, and set individual display values for all fibers
    void attributeLineDisp(FiberSet const&);

    /// get ready to display
    void prepareDrawing(Simul const&, PropertyList&, PropertyList&);

    /// display all simulation objects
    void drawSimul(Simul const&);
    
    /// display for periodic systems
    void drawTiled(Simul const&, int nine);

    /// scale from size / line width to natural units
    float pixscale(float w) const { return w * sizeScale; }
    
    /// scale from size into OpenGL line width (pixel) units
    float pixwidth(float w) const { return std::max(w * unitValue, 0.125f); }
    
    /// return a vector of norm `rad`, orthogonal to both `dir` and the depth axis of the current view
    Vector inViewPlane(Vector const& dir, const real rad) const
    {
#if ( DIM == 3 )
        return cross(dir, depthAxis).normalized(rad);
#elif ( DIM == 2 )
        return dir.orthogonal(rad);
#else
        return Vector(0, 1, 0);
#endif
    }
    
    /// draw primitive `obj` at given position
    static void drawObject(Vector const& pos, float rad, void (*obj)());
    
    /// draw primitive `obj` at `pos` with Z-axis oriented toward `dir`
    static void drawObject(Vector const& pos, Vector const& dir, float rad, void (*obj)());

    
    /// draw all scalar fields
    void drawFields(FieldSet const&);
    
    /// draw one Space
    void drawSpace3D(Space const*, bool back);

    /// draw all Spaces (in 3D the back side)
    void drawSpaces(SpaceSet const&);
    
    /// draw the front-side of Spaces in 3D
    void drawTransparentSpaces(SpaceSet const&);


    /// draw thin lines joining the Fiber vertices
    void drawFiberBackbone(Fiber const&, gym_color col, float width) const;

    /// draw minus end of one fiber
    virtual void drawFiberEndMinus(Fiber const&, int style, float size) const;
    
    /// draw plus end of one fiber
    virtual void drawFiberEndPlus(Fiber const&, int style, float size) const;
    
    /// draw fresh assembly near the plus ends, using white stripes
    void drawFiberGrowth(Fiber const&, float size) const;

    /// draw linear features of one fiber
    virtual void drawFiberLines(Fiber const&, int style, float width) const;
    
    /// draw one segment of one fiber (used to display transparent fibers)
    virtual void drawFiberSegmentT(Fiber const&, unsigned) const;
    
    /// Using triangles to draw a broken line with width `2*rad`
    void drawFiberWidePath(Fiber const& fib, float rad) const;

    /// draw stripes of alternating colors from segments of length `inc`, in [abs, sup]
    void drawFiberStriped2D(Fiber const&, float rad, real inc, gym_color, real onk, gym_color) const;
    
    /// draw stripes of alternating colors from segments of length `inc`, in [abs, sup]
    void drawFiberArrowed2D(Fiber const&, float rad, real inc, gym_color, real onk, gym_color) const;

    /// draw stripes of alternating colors from segments of length `inc`, in [abs, sup]
    void drawFiberStriped(Fiber const&, float rad, real inc, gym_color, real onk, gym_color) const;

    /// draw stripes of alternating colors from segments of length `inc`, in [abs, sup]
    void drawFiberStripedClip(Fiber const&, float rad, real inc, gym_color, real onk, gym_color) const;

    /// actin-like rendering using a sphere to represent each monomer
    void drawFilament(Fiber const& fib, real, gym_color, gym_color, gym_color) const;

    /// actin-like rendering using a sphere to represent each monomer
    void drawActin(Fiber const& fib, gym_color, gym_color, gym_color) const;
    
    /// microtubule-like rendering using a sphere to represent each monomer
    void drawMicrotubule(Fiber const& fib, gym_color, gym_color, gym_color) const;
    
    /// draw beads along the backbone of the fiber
    void drawFiberChevrons(Fiber const&, real rad, real gap) const;

    /// draw Fiber point-like features, eg. vertices
    virtual void drawFiberPoints(Fiber const&) const;
    
    /// draw fiduciary marks on all fibers
    virtual void drawFiberSpeckles(Fiber const&) const;
   
    /// display lattice substance using color
    virtual void drawFiberLattice1(Fiber const&, VisibleLattice const&, float rad) const;
    
    /// display lattice substance using color
    virtual void drawFiberLattice2(Fiber const&, VisibleLattice const&, float rad) const;
   
    /// display lattice cell substance and edges
    virtual void drawFiberLattice3(Fiber const&, VisibleLattice const&, float rad) const;

    /// display lattice cell edges
    virtual void drawFiberLatticeEdges(Fiber const&, VisibleLattice const&, float rad) const;

    /// display lattice cell values
    virtual void drawFiberLatticeValues(Fiber const&, VisibleLattice const&) const;

    /// display lattice bits
    virtual void drawFiberLatticeBits(Fiber const&, FiberLattice const&) const;

    /// display labels for a Fiber
    void drawFiberLabels(Fiber const&, int style) const;
    
    /// display forces acting on the fiber vertices
    void drawFiberForces(Fiber const&, real scale, float width) const;
    
    /// draw all features of one Fiber
    void drawFiber(Fiber const&);
    
    /// draw all Fibers
    void drawFibers(FiberSet const&);
    
    /// draw text associated with Fibers
    void drawFiberTexts(FiberSet const&);

    /// draw some average taken over fibers provided in list
    void drawAverageFiber(ObjectList const&, gym_color) const;
    
    /// draw some arrows computed by averaging over fibers
    void drawAverageFiber1(FiberSet const&, Property const* ) const;
    
    /// draw the average for left-pointing and right-pointing fibers
    void drawAverageFiber2(FiberSet const&, Property const* ) const;

    
    /// draw a Bead
    void drawBead(Bead const&);

    /// draw translucent elements of a Bead
    void drawBeadT(Bead const&) const;
    
    /// draw all Beads
    void drawBeads(BeadSet const&);

    
    /// draw opaque elements of a Solid
    void drawSolid(Solid const&);
    
    /// draw translucent elements of a Solid
    void drawSolidT(Solid const&, unsigned) const;

    /// draw all Solids
    void drawSolids(SolidSet const&);
    
    
    /// draw one Sphere
    void drawSphere(Sphere const&);

    /// draw translucent elements of a Sphere
    void drawSphereT(Sphere const&) const;
    
    /// draw all Spheres
    void drawSpheres(SphereSet const&);
    
    
    /// draw all free Singles
    virtual void drawSinglesF(SingleSet const&) const = 0;
    
    /// draw all attached Singles
    virtual void drawSinglesA(SingleSet const&) const = 0;
    
    /// draw all free Couples, showing Hand1
    virtual void drawCouplesF(CoupleSet const&) const = 0;

    /// draw all attached Couples
    virtual void drawCouplesA(CoupleSet const&) const = 0;
    
    /// draw one bridging Couple
    virtual void drawCoupleB(Couple const*) const { }

    /// draw all bridging Couples
    void drawCouplesB(CoupleSet const&) const;

    /// draw one Organizer
    virtual void drawOrganizer(Organizer const&) const;
    
    /// draw all Organizers
    void drawOrganizers(OrganizerSet const&);

    
    /// draw additional items
    void drawMisc(Simul const&);
    
};


#endif

