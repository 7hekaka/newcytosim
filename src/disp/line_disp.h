// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef LINE_DISP_H
#define LINE_DISP_H


#include "real.h"
#include "gym_color.h"


/// Display parameters for a Line
/**
 
 LineDisp holds the display attributes for a particular Fiber, and accordingly,
 there is one LineDisp for every Fiber. In contrast, FiberDisp holds class attributes.
 
 A user cannot set these attributes directly. Instead, the values are derived
 from the relevant FiberDisp automatically, when cytosim prepares the display.
 
 For example `end_color` will be set as one value of `FiberDisp::end_colors[]`,
 depending on the states of the ends of this fiber.

*/
class LineDisp
{
public:
    
    /// color of body
    gym_color color;
    
    /// colors of plus end and minus end
    gym_color end_color[2];
    
    /// scale to convert to color to display lines
    real color_scale;

    /// visibility flag
    int visible;

public:
    
    /// constructor
    LineDisp() { clear(); }
    
    /// destructor
    ~LineDisp() { }

    /// set to default values
    void clear();
    
};


#endif

