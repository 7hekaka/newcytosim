// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

#include "gym_cap.h"


GLboolean gym::depth_ = 0;
GLboolean gym::cull_ = 0;
GLboolean gym::blend_ = 0;


GLboolean gym::light_ = 0;
GLboolean gym::alpha_ = 0;


//-----------------------------------------------------------------------
#pragma mark - Clip Planes

void gym::enableLineStipple(short pattern)
{
#ifdef GL_VERSION_2_1
    glLineStipple(1, pattern);
    glEnable(GL_LINE_STIPPLE);
#endif
}

void gym::disableLineStipple()
{
#ifdef GL_VERSION_2_1
    glDisable(GL_LINE_STIPPLE);
#endif
}

void gym::setClipPlane(unsigned glp, double X, double Y, double Z, double S)
{
    GLdouble eq[4] = { X, Y, Z, S };
    glClipPlane(GL_CLIP_PLANE0+glp, eq);
}

