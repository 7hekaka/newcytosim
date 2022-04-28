// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

#include "gym_draw.h"

/// current color
//float gym::col_[4] = { 1, 1, 1, 1 };


GLboolean gym::depth_ = 0;
GLboolean gym::cull_ = 0;
GLboolean gym::blend_ = 0;

/**
 draw back first, and then front of object,
 CULL_FACE is temporarily enabled for this
 */
void gym::dualPass(void primitive())
{
    gym::enableCullFace(GL_FRONT);
    primitive();
    gym::switchCullFace(GL_BACK);
    primitive();
    gym::restoreCullFace();
}
