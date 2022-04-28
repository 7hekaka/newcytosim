// Cytosim was created by Francois Nedelec. Copyright Cambridge University 2020

#ifndef DISPLAY1_H
#define DISPLAY1_H

#include "display.h"
class PointDisp;

///Cytosim display class for style=1
/**
 This style produces a fast 2D display.
 Some of the parameters in PointDisp are ignored.

 Point-like objects are rendered using OpenGL::Points.
 All points are displayed with the same size `point_size`.
 */
class Display1 : public Display
{
    
    /// global OpenGL line width
    float linkWidth;

    /// global OpenGL point size
    float pointSize;

public:
    
    ///constructor
    Display1(DisplayProp const*);
    
    ///destructor
    ~Display1() {}
    
    
    /// draw the given simulation state using OpenGL commands
    void drawObjects(Simul const&);
   
    /// draw Fibers with offset
    void drawFiber(Fiber const&);
    
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
    void drawCouplesB(CoupleSet const&) const;

};

#endif

