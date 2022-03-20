// Cytosim was created by Francois Nedelec. Copyright Cambridge University 2020
#ifndef DISPLAY_H
#define DISPLAY_H

#include "real.h"
#include "array.h"
#include "fiber.h"
#include "mecable.h"
#include "mecapoint.h"
#include "gle_color.h"
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
#define DISPLAY

/**
 @brief A display element with a depth coordinate
 
 zObject is used to depth-sort and display transparent objects
 */
class zObject
{
    /// pointer to object
    Mecapoint point_;
    
    /// distance to the imaging plane
    real depth_;

public:
    
    zObject() : depth_(0) { }
    
    zObject(Mecable const* m) : point_(m, 0), depth_(0) { }

    zObject(Mecable const* m, size_t i) : point_(m, i), depth_(0) { }
    
    /// position
    Vector position() const { return point_.pos(); }
    
    /// query depth
    real depth() const { return depth_; }
    
    /// set depth
    void depth(real z) { depth_ = z; }
    
    /// display object
    void draw(Display const*) const;
};


///Base class to display Cytosim state using OpenGL
class Display
{
protected:
    
    /// array of transparent objects to be displayed last
    Array<zObject> zObjects;
    
    /// the pixel size for this particular display
    float pixelSize;
    
    /// the value of the size & width unit in pixels
    float unitValue;
    
    /// the value of the size & width unit in natural units (ie. micrometers)
    float sizeScale;
    
    /// direction of view
    Vector3 depthAxis;

    /// flag used to calculate expensive analysis only once
    size_t prep_flag;
    
    /// used to calculate clusterAnalysis only once
    double prep_time;
    
    /// min and max age used to adjust color range with COLORING_AGE
    double age_scale, age_start;
    
    /// use OpenGL stencil test:
    bool stencil_;
    
private:
    
    /// set default value of FiberProp
    void prepareFiberDisp(FiberProp*, PropertyList&, gle_color);
    
    /// set values of fiber's LineDisp
    void prepareLineDisp(Fiber const*);
    
    template < typename T >
    void preparePointDisp(T * prop, PropertyList&, gle_color);
    
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
    void setPixelFactors(float pixel_size, float unit_value);
    
    /// get ready to display
    void prepareForDisplay(Simul const&, PropertyList&, Vector3 const&);

    /// display the whole simulation
    void drawSimul(Simul const&);
    
    /// display for periodic systems
    void drawTiled(Simul const&, int nine);

    /// scale from pixel size / line width to natural units
    float pixscale(float w) const { return w * sizeScale; }
    
    /// set OpenGL line width
    void lineWidth(float w) const { glLineWidth(std::max(w*unitValue, 0.25f)); }

    /// set OpenGL point size
    void pointSize(float w) const { glPointSize(std::max(w*unitValue, 0.25f)); }
    
    /// draw primitive `obj` at given position
    void drawObject(Vector const& pos, float rad, void (*obj)()) const;
    
    /// draw primitive `obj` at `pos` with Z-axis oriented toward `dir`
    void drawObject(Vector const& pos, Vector const& dir, float rad, void (*obj)()) const;
    
    /// draw primitive `obj` at `pos` with Z-axis oriented toward `dir`
    void drawFlat(Vector const& pos, float rad, void (*obj)()) const;
    
    /// draw a fine spherical object
    void drawBallT(Vector const&, real radius, gle_color const&) const;
    
    /// draw a fine spherical object
    void drawDiscT(Vector const&, real radius, gle_color const&) const;


    
    /// draw a scalar field
    void drawFields(FieldSet const&);
    
    /// draw a Space
    void drawSpace(Space const*, bool);

    /// draw all Spaces (in 3D the back side)
    void drawSpaces(SpaceSet const&);
    
    /// draw the front-side of Spaces in 3D
    void drawTransparentSpaces(SpaceSet const&);

    
    /// display string of points
    static void drawStrip(size_t cnt, real const* pts, GLenum);

    /// draw thin lines joining the Fiber vertices
    void drawFiberBackbone(Fiber const&) const;

    /// draw Fiber MINUS_END
    virtual void drawFiberMinusEnd(Fiber const&, int style, float size) const;
    
    /// draw Fiber PLUS_END
    virtual void drawFiberPlusEnd(Fiber const&, int style, float size) const;

    /// draw Fiber linear features
    virtual void drawFiberLines(Fiber const&, int style) const;
    
    /// draw one segment of a Fiber (used to display transparent fibers)
    virtual void drawFiberSegmentT(Fiber const&, size_t) const;

    /// actin-like rendering using a sphere to represent each monomer
    void         drawFilament(Fiber const& fib, gle_color const&, gle_color const&, gle_color const&) const;

    /// actin-like rendering using a sphere to represent each monomer
    void         drawActin(Fiber const& fib, gle_color const&, gle_color const&, gle_color const&) const;
    
    /// microtubule-like rendering using a sphere to represent each monomer
    void         drawMicrotubule(Fiber const& fib, gle_color const&, gle_color const&, gle_color const&) const;
    
    /// draw Fiber point-like features
    virtual void drawFiberPoints(Fiber const&) const;
    
    /// draw Fiber Speckles
    virtual void drawFiberSpeckles(Fiber const&) const;
   
    /// display lattice subtance using color
    virtual void drawFiberLattice1(Fiber const&, VisibleLattice const&, real width) const;
    
    /// display lattice subtance using color
    virtual void drawFiberLattice2(Fiber const&, VisibleLattice const&, real width) const;
   
    /// display lattice cell edges
    virtual void drawFiberLattice3(Fiber const&, VisibleLattice const&, real width) const;

    /// display lattice cell edges
    virtual void drawFiberLatticeEdges(Fiber const&, VisibleLattice const&, real width) const;

    /// display Labels for a Fiber
    void         drawFiberLabels(Fiber const&, int style) const;
    
    /// display forces acting on the fiber vertices
    void         drawFiberForces(Fiber const&, real scale) const;
    
    /// draw all features of Fiber
    virtual void drawFiber(Fiber const&);
    
    /// draw Fibers
    virtual void drawFibers(FiberSet const&);
    
    /// draw the average fiber for the pool defined by func(obj, val) == true
    void drawAverageFiber(ObjectList const&);
    
    /// draw the averaged fiber
    void drawAverageFiber1(FiberSet const&, void const*);
    
    /// draw the average for left-pointing and right-pointing fibers
    void drawAverageFiber2(FiberSet const&, void const*);

    
    /// draw a Bead
    void drawBead(Bead const&);

    /// draw translucent elements of a Bead
    void drawBeadT(Bead const&) const;
    
    /// draw the Beads
    void drawBeads(BeadSet const&);

    
    /// draw opaque elements of a Solid
    void drawSolid(Solid const&);
    
    /// draw translucent elements of a Solid
    void drawSolidT(Solid const&, size_t) const;
    
    /// draw the Solids
    void drawSolids(SolidSet const&);
    
    
    /// draw the Sphere
    void drawSphere(Sphere const&);

    /// draw translucent elements of a Sphere
    void drawSphereT(Sphere const&) const;
    
    /// draw the Spheres
    void drawSpheres(SphereSet const&);
    
    
    /// draw the free Singles
    virtual void drawSinglesF(SingleSet const&) const = 0;
    
    /// draw the attached Singles
    virtual void drawSinglesA(SingleSet const&) const = 0;
    
    /// draw the free Couples, showing Hand1
    virtual void drawCouplesF1(CoupleSet const&) const = 0;
    
    /// draw the free Couples, showing Hand1 or Hand2
    virtual void drawCouplesF2(CoupleSet const&) const = 0;
    
    /// calls drawCouplesF1 or drawCouplesF2
    void drawCouplesF(CoupleSet const&) const;
    
    /// draw the attached Couples
    virtual void drawCouplesA(CoupleSet const&) const = 0;
    
    /// draw one bridging Couple
    virtual void drawCoupleB(Couple const*) const { }

    /// draw the bridging Couples
    virtual void drawCouplesB(CoupleSet const&) const;

    /// draw Organizer
    virtual void drawOrganizer(Organizer const&) const;
    
    /// draw the Organizers
    void drawOrganizers(OrganizerSet const&);

    
    /// draw additional items
    void drawMisc(Simul const&);
    
};


#endif

