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
public:
    
    ///constructor
    Display2(DisplayProp const*);
    
    ///destructor
    ~Display2() {}
    
    
    /// draw the given simulation state using OpenGL commands
    void drawObjects(Simul const&);

    /// draw free Singles
    void drawSinglesF(SingleSet const&) const;
    
    /// draw attached Singles
    void drawSinglesA(SingleSet const&) const;

    /// draw free Couples
    void drawCouplesF1(CoupleSet const&) const;
    
    /// draw free Couples, randomizing which Hand is drawn
    void drawCouplesF2(CoupleSet const&) const;
    
    /// draw free Couple
    void drawCouplesF(CoupleSet const&) const;

    /// draw attached Couples
    void drawCouplesA(CoupleSet const&) const;
    
    /// draw bridging Couples
    void drawCoupleB(Couple const*) const;
    
};

#endif

