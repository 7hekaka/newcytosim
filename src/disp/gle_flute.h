// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#ifndef GLE_FLUTE_H
#define GLE_FLUTE_H

#include "flute.h"

/*
 Attention:
 Since the types depend on the dimensionality, this unit must be compiled
 separately for each different value of DIM.
 */

#if ( DIM == 3 )
typedef flute3 fluteV;
typedef flute6 fluteVN;
typedef flute8 fluteVC;
#else
typedef flute2 fluteV;
typedef flute4 fluteVN;
typedef flute6 fluteVC;
#endif

namespace gle
{
    fluteV* mapVertexBuffer(size_t);
    void  unmapVertexBuffer();
    void   bindVertexBuffer(size_t);

    fluteVN* mapVertexNormalBuffer(size_t);
    void  unmapVertexNormalBuffer();
    void   bindVertexNormalBuffer(size_t);

    fluteVC* mapVertexColorBuffer(size_t);
    void  unmapVertexColorBuffer();
    void   bindVertexColorBuffer(size_t);
    
    unsigned* mapIndexBuffer(size_t);
    void  unmapIndexBuffer();
    void   bindIndexBuffer(size_t);
};

#endif
