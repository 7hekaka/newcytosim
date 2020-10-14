// Cytosim was created by Francois Nedelec. Copyright Cambridge University 2020

#ifndef DISPLAY2_H
#define DISPLAY2_H

#include "display.h"
class PointDisp;

///Cytosim display class for style=2
/**
 This is a 2D display using Bitmap for Hands.
 It implements most of the characteristics in PointDisp and FiberDisp
 */
class Display2 : public Display
{
    /// draw a fine spherical object
    void drawBallT(Vector const&, real radius, gle_color const&) const;
    
    /// draw a point
    void drawPoint(Vector const&, PointDisp const*) const;

public:
    
    ///constructor
    Display2(DisplayProp const*);
    
    ///destructor
    ~Display2() {}
    
    
    /// draw the given simulation state using OpenGL commands
    void drawSimul(Simul const&);
    
    /// draw the Solids
    void drawSolid(Solid const&);
 
    /// draw the transparent part for the Solids
    void drawSolidT(Solid const&, size_t);
    
    /// draw a Bead
    void drawBead(Bead const&);
    
    /// draw transparent membrane of Bead
    void drawBeadT(Bead const&);
    
    /// draw a Sphere
    void drawSphere(Sphere const&);
    
    /// draw transparent membrane of Sphere
    void drawSphereT(Sphere const&);
    
    /// draw free Singles
    void drawSinglesF(SingleSet const&) const;
    
    /// draw attached Singles
    void drawSinglesA(SingleSet const&) const;

    /// draw free Couples
    void drawCouplesF1(CoupleSet const&) const;
    
    /// draw free Couples, randomizing which Hand is drawn
    void drawCouplesF2(CoupleSet const&) const;
    
    /// draw attached Couples
    void drawCouplesA(CoupleSet const&) const;
    
    /// draw bridging Couples
    void drawCoupleB(Couple const*) const;
    
    /// draw an Organizer
    void drawOrganizer(Organizer const&) const;
};

#endif

