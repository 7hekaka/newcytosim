// Cytosim was created by Francois Nedelec. Copyright Cambridge University 2021

#ifndef DISPLAY3_H
#define DISPLAY3_H

#include "display.h"
#include "real.h"
#include "vector.h"
#include "gym_color.h"
#include "point_disp.h"


///Cytosim display class for style=3
/**
 This style is for rendering in 3D.
 It uses Lighting for better volume rendering
 */
class Display3 : public Display
{
private:
    
    /// draw a point with a small sphere
    void drawPoint(Vector const&, PointDisp const*) const;

    /// draw primitive `obj` at given position
    void drawObject3(Vector const&, float, void (*obj)()) const;

    /// draw a point with a small sphere
    inline void drawHand(Vector const& pos, PointDisp const* dis) const;

    /// draw a point with a small sphere
    inline void drawHandF(Vector const& pos, PointDisp const* dis) const;
    
    /// draw Fiber model segments
    void drawFiberSegmentsClip(Fiber const&, float rad,
                               gym_color (*set_color)(Fiber const&, unsigned)) const;
    
    /// draw Fiber segments not necessarily aligned with the vertices
    void drawFiberSectionsClip(Fiber const&, float rad, int inx, int last, real abs, real inc,
                               gym_color (*set_color)(Fiber const&, int, real), real fac, real facM, real facP) const;
    
    /// draw Fiber model segments
    void drawFiberSegmentsJoin(Fiber const&, float rad,
                               gym_color (*set_color)(Fiber const&, unsigned)) const;
    
    /// draw Fiber segments not necessarily aligned with the vertices
    void drawFiberSectionsJoin(Fiber const&, float rad, int inx, int last, real abs, real inc,
                               gym_color (*set_color)(Fiber const&, int, real), real fac, real facM, real facP) const;

    /// display lattice substance using specified color function
    void drawFiberLattice(Fiber const&, VisibleLattice const&, float rad, gym_color (*set_color)(Fiber const&, int, real)) const;

public:
        
    ///constructor
    Display3(DisplayProp const*);
    
    ///destructor
    ~Display3() {}
    
    /// draw the given simulation state using OpenGL commands
    void drawObjects(Simul const&);
    
    /// draw Fiber minus end
    void drawFiberEndMinus(Fiber const&, int style, float size) const;
    
    /// draw Fiber plus end
    void drawFiberEndPlus(Fiber const&, int style, float size) const;
    
    /// draw Fiber linear features
    void drawFiberLines(Fiber const&, int style, float width) const;
    
    /// draw one segment of a Fiber
    void drawFiberSegmentT(Fiber const&, unsigned) const;

    /// display lattice substance using color
    void drawFiberLattice1(Fiber const&, VisibleLattice const&, float rad) const;
    
    /// display lattice substance using color
    void drawFiberLattice2(Fiber const&, VisibleLattice const&, float rad) const;
    
    /// display lattice substance using color
    void drawFiberLattice3(Fiber const&, VisibleLattice const&, float rad) const;

    /// draw Edges of lattice
    void drawFiberLatticeEdges(Fiber const&, VisibleLattice const&, float rad) const;

    /// draw Fiber point-like features
    void drawFiberPoints(Fiber const&) const;
    
    /// draw points distributed randomly along fibers, at fixed positions
    void drawFiberSpeckles(Fiber const&) const;
    
    /// draw the free Single
    void drawSinglesF(SingleSet const&) const;

    /// draw the attached Single
    void drawSinglesA(SingleSet const&) const;
    
    /// draw an attached Single
    void drawSingleA(Single const*) const;
    
    /// draw an attached Single with a link
    void drawSingleB(Single const*) const;

    /// draw free Couple
    void drawCouplesF1(CoupleSet const&) const;

    /// draw free Couple, randomizing which Hand is drawn
    void drawCouplesF2(CoupleSet const&) const;
    
    /// draw free Couple
    void drawCouplesF(CoupleSet const&) const;

    /// draw attached Couple
    void drawCouplesA(CoupleSet const&) const;

    /// draw one bridging Couple
    void drawCoupleB(Couple const*) const;
 
    /// draw one bridging Couple
    void drawCoupleBcrude(Couple const*) const;
    
    /// draw one bridging Couple made of two identical hands
    void drawCoupleBhomo(Couple const*, PointDisp const*) const;
    
    /// draw one bridging Couple using 2 feet per hand
    void drawCoupleBwalk(Couple const*) const;

    /// draw one bridging Couple
    void drawCoupleBori(Couple const*) const;

    /// draw one bridging Couple
    void drawCoupleBside(Couple const*) const;
    
    /// draw one bridging Couple
    void drawCoupleBalt(Couple const*) const;

    /// draw one Organizer
    void drawOrganizer(Organizer const&) const;
};

#endif

