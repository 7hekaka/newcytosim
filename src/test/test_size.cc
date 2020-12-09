// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "simul.h"

#include "matrix11.h"
#include "matrix22.h"
#include "matrix33.h"
#include "matrix34.h"
#include "matrix44.h"

#include "sparmatsym.h"
#include "sparmatsym1.h"
#include "sparmatsym2.h"
#include "sparmatblk.h"
#include "sparmatsymblk.h"
#include "sparmatsymblkdiag.h"

#define PRINT(arg) printf("sizeof %16s   %lu bytes\n", #arg, sizeof(arg));

int main(int argc, char* argv[])
{
    if ( argc > 1 )
    {
        PRINT(Random);
        PRINT(Array<int>);

        PRINT(SparMatSym::Element);
        PRINT(SparMatSym1::Element);
        PRINT(SparMatSym2::Element);
        
        PRINT(SparMatBlk::Line);
        PRINT(SparMatSymBlk::Column);
        PRINT(SparMatSymBlkDiag::Column);
        
        PRINT(PointGrid);
        PRINT(FatPoint);
        PRINT(FatLocus);
        
        PRINT(PointGridF);
        PRINT(BigPoint);
        PRINT(BigLocus);
    }
    else
    {
        PRINT(Vector1);
        PRINT(Vector2);
        PRINT(Vector3);
        PRINT(Vector4);
        
        PRINT(Matrix11);
        PRINT(Matrix22);
        PRINT(Matrix33);
        PRINT(Matrix34);
        PRINT(Matrix44);
        
        PRINT(Mecapoint);
        PRINT(FiberSegment);
        PRINT(Interpolation);
        PRINT(HandMonitor);
        PRINT(FiberSite);
        PRINT(Hand);
        
        PRINT(Inventoried);
        PRINT(Movable);
        PRINT(Buddy);
        
        PRINT(Object);
        PRINT(Mecable);
        PRINT(Chain);
        PRINT(Mecafil);
        PRINT(Fiber);
        
        PRINT(Space);
        PRINT(Single);
        PRINT(Couple);
        
        PRINT(Solid);
        PRINT(Bead);
        PRINT(Sphere);
        PRINT(Simul);
    }
    return EXIT_SUCCESS;
}
