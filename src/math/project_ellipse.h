// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef PROJECT_ELLIPSE_H
#define PROJECT_ELLIPSE_H


#include "real.h"


/// calculate `(pX, pY)`, the projection of `(wX, wY)` on the ellipse of axes `radX, radY`
void projectEllipse(real&  pX, real&  pY,
                    real   wX, real   wY,
                    real radX, real radY);


/// calculate `p[]`, the projection of `w[]` on the ellipse of axes given in `rad[]`
void projectEllipsoid(real         p[3],
                      const real   w[3],
                      const real rad[3]);


#endif
