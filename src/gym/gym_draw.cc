// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

#include "gym_draw.h"
#include "gym_cap.h"

/// current color
//float gym::col_[4] = { 1, 1, 1, 1 };

/// make RGBA color current
void gym::setColor(const float col[]) { glColor4fv(col); }

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
