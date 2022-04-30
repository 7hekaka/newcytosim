// Cytosim was created by Francois Nedelec. Copyright Cambridge University 2021

#ifndef DISPLAY3_H
#define DISPLAY3_H

#include "display.h"
#include "real.h"
#include "vector.h"
#include "gle_color.h"
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
    void drawPoint(Vector const&, float) const;
    
    /// draw a point with a small sphere
    void drawPoint(Vector const&, PointDisp const*) const;

    /// draw primitive `obj` at given position
    void drawObject3(Vector const&, float, void (*obj)()) const;

    /// draw a point with a small sphere
    inline void drawHand(Vector const& pos, PointDisp const* dis) const;

    /// draw a point with a small sphere
    inline void drawHandF(Vector const& pos, PointDisp const* dis) const;
    
    /// draw Fiber model segments
    void drawFiberSegmentsClip(Fiber const&, real rad,
                               gle_color (*set_color)(Fiber const&, size_t)) const;
    
    /// draw Fiber segments not necessarily aligned with the vertices
    void drawFiberSectionsClip(Fiber const&, real rad, long inx, long last, real abs, real inc,
                               gle_color (*set_color)(Fiber const&, long, real), real fac, real facM, real facP) const;
    
    /// draw Fiber model segments
    void drawFiberSegmentsJoin(Fiber const&, real rad,
                               gle_color (*set_color)(Fiber const&, size_t)) const;
    
    /// draw Fiber segments not necessarily aligned with the vertices
    void drawFiberSectionsJoin(Fiber const&, real rad, long inx, long last, real abs, real inc,
                               gle_color (*set_color)(Fiber const&, long, real), real fac, real facM, real facP) const;

    /// display lattice subtance using specified color function
    void drawFiberLattice(Fiber const&, VisibleLattice const&, real width, gle_color (*set_color)(Fiber const&, long, real)) const;

public:
        
    ///constructor
    Display3(DisplayProp const*);
    
    ///destructor
    ~Display3() {}
    
    /// draw the given simulation state using OpenGL commands
    void drawObjects(Simul const&);
    
    /// draw Fiber MINUS_END
    void drawFiberMinusEnd(Fiber const&, int style, float size) const;
    
    /// draw Fiber PLUS_END
    void drawFiberPlusEnd(Fiber const&, int style, float size) const;
    
    /// draw Fiber linear features
    void drawFiberLines(Fiber const&, int style) const;
    
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
    void drawCoupleBplain(Couple const*) const;
    
    /// draw one bridging Couple
    void drawCoupleBside(Couple const*) const;
    
    /// draw one bridging Couple
    void drawCoupleBalt(Couple const*) const;

    /// draw an Organizer
    void drawOrganizer(Organizer const&) const;
};

#endif

