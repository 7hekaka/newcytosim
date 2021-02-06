// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef DISPLAY3_H
#define DISPLAY3_H

#include "display.h"
#include "real.h"
#include "opengl.h"
#include "vector.h"
#include "point_disp.h"


///Cytosim display class for style=3
/**
 This style is for rendering in 3D.
 It uses Lighting for better volume rendering
 */
class Display3 : public Display
{
private:
    
    /// draw a fine spherical object
    void drawTBall(Vector const&, float radius, gle_color const&) const;

    /// draw a point with a small sphere
    void drawPoint(Vector const&, PointDisp const*) const;
    
    /// draw a point with a small sphere
    void drawPoint(Vector const&, float) const;

    
    /// draw a point with a small sphere
    inline void drawHand(Vector const& p, PointDisp const* d) const
    { d->color.load_both(); drawPoint(p, d); }

    /// draw a point with a small sphere
    inline void drawHand2(Vector const& p, PointDisp const* d) const
    { if ( d->visible ) { d->color2.load_both(); drawPoint(p, d); } }
    
    /// draw Fiber model segments
    void drawFiberSegments(Fiber const&, real rad,
                           gle_color (*set_color)(Fiber const&, size_t)) const;
    
    /// draw Fiber segments not necessarily aligned with the vertices
    void drawFiberSubSegments(Fiber const&, real rad, VisibleLattice::lati_t inx, VisibleLattice::lati_t last, real abs, real inc,
                              gle_color (*set_color)(Fiber const&, long, real), real fac, real facM, real facP) const;
    
    /// display lattice subtance using specified color function
    void drawFiberLattice(Fiber const&, VisibleLattice const&, real width, gle_color (*set_color)(Fiber const&, long, real)) const;

public:
        
    ///constructor
    Display3(DisplayProp const*);
    
    ///destructor
    ~Display3() {}
    
    /// draw the given simulation state using OpenGL commands
    void drawSimul(Simul const&);
    
    /// draw Fiber MINUS_END
    void drawFiberMinusEnd(Fiber const&, int style, real size) const;
    
    /// draw Fiber PLUS_END
    void drawFiberPlusEnd(Fiber const&, int style, real size) const;
    
    /// draw Fiber linear features
    void drawFiberLines(Fiber const&) const;
    
    /// draw one segment of a Fiber
    void drawFiberSegmentT(Fiber const&, size_t) const;

    /// display lattice subtance using color
    void drawFiberLattice1(Fiber const&, VisibleLattice const&, real width) const;
    
    /// display lattice subtance using color
    void drawFiberLattice2(Fiber const&, VisibleLattice const&, real width) const;
    
    /// display lattice subtance using color
    void drawFiberLattice3(Fiber const&, VisibleLattice const&, real width) const;

    /// draw Edges of lattice
    void drawFiberLatticeEdges(Fiber const&, VisibleLattice const&, real width) const;

    /// draw Fiber point-like features
    void drawFiberPoints(Fiber const&) const;
    
    /// draw Fiber speckles
    void drawFiberSpeckles(Fiber const&) const;
    
    /// draw the Solids
    void drawSolid(Solid const&);
 
    /// draw the transparent parts of Solid
    void drawSolidT(Solid const&, size_t);

    /// draw a Bead
    void drawBead(Bead const&);
    
    /// draw transparent membrane of Bead
    void drawBeadT(Bead const&);
    
    /// draw a Sphere
    void drawSphere(Sphere const&);
    
    /// draw transparent membrane of Sphere
    void drawSphereT(Sphere const&);
    
    /// draw the free Single
    void drawSinglesF(SingleSet const&) const;

    /// draw the attached Single
    void drawSinglesA(SingleSet const&) const;
    
    /// draw the attached Single
    void drawSingleA(Single const*) const;

    /// draw free Couple
    void drawCouplesF1(CoupleSet const&) const;

    /// draw free Couple, randomizing which Hand is drawn
    void drawCouplesF2(CoupleSet const&) const;

    /// draw attached Couple
    void drawCouplesA(CoupleSet const&) const;

    /// draw one bridging Couple
    void drawCoupleB(Couple const*) const;
 
    /// draw one bridging Couple
    void drawCoupleBfast(Couple const*) const;

    /// draw an Organizer
    void drawOrganizer(Organizer const&) const;
};

#endif

