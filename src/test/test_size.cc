// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "simul.h"
#include "single_prop.h"
#include "couple_prop.h"
#include "event.h"
#include "hands/digit.h"
#include "hands/walker.h"
#include "hands/cutter.h"
#include "hands/motor.h"
#include "hands/mighty.h"
#include "single.h"
#include "singles/picket.h"
#include "singles/picket_long.h"
#include "singles/wrist.h"
#include "singles/wrist_long.h"
#include "couple.h"
#include "couples/couple_long.h"
#include "couples/crosslink.h"
#include "couples/crosslink_long.h"
#include "couples/shackle.h"
#include "couples/shackle_long.h"
#include "couples/bridge.h"
#include "couples/duo.h"
#include "couples/duo_long.h"
#include "couples/fork.h"

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
#include "rasterizer.h"

#define PRINT(arg) printf("sizeof %32s  %lu bytes\n", #arg, sizeof(arg));

int main(int argc, char* argv[])
{
    int x = ( argc>1 ? atoi(argv[1]) : 0 );
    switch ( x )
    {
        case 4:
            PRINT(Random);
            PRINT(Array<int>);
            
            PRINT(PointGrid)
            PRINT(FatVector);
            PRINT(FatPoint);
            PRINT(FatLocus);
            PRINT(PointGridCell);
            
            PRINT(BigVector);
            PRINT(BigLocus);
            PRINT(BigLocusList);
            PRINT(LocusGrid);
            
            PRINT(Rasterizer::Vertex2);
            PRINT(Rasterizer::Vertex2dZ);
            PRINT(Rasterizer::Vertex3);
            break;
        case 3:
            PRINT(Single);
            PRINT(Picket);
            PRINT(PicketLong);
            PRINT(Wrist);
            PRINT(WristLong);
            PRINT(Couple);
            PRINT(Crosslink);
            PRINT(Shackle);
            PRINT(Bridge);
            PRINT(Fork);
            PRINT(Duo);
            PRINT(CoupleLong);
            PRINT(CrosslinkLong);
            PRINT(ShackleLong);
            PRINT(DuoLong);
            break;
        case 2:
            PRINT(Mecapoint);
            PRINT(FiberSegment);
            PRINT(Interpolation);
            PRINT(HandMonitor);
            PRINT(FiberSite);
            PRINT(Hand);
            PRINT(Motor);
            PRINT(Cutter);
            PRINT(Digit);
            PRINT(Walker);
            PRINT(Mighty);
            break;
        case 1:
            PRINT(Inventoried);
            PRINT(Buddy);
            PRINT(Object);
            PRINT(Mecable);
            PRINT(Chain);
            PRINT(Mecafil);
            PRINT(Fiber);
            PRINT(Space);
            PRINT(Solid);
            PRINT(Bead);
            PRINT(Sphere);
            PRINT(Event);
            PRINT(SimulProp);
            PRINT(Simul);
            break;
        case 0:
            PRINT(Vector1);
            PRINT(Vector2);
            PRINT(Vector3);
            PRINT(Vector4);
            
            PRINT(Matrix11);
            PRINT(Matrix22);
            PRINT(Matrix33);
            PRINT(Matrix34);
            PRINT(Matrix44);
            
            PRINT(SparMatSym::Element);
            PRINT(SparMatSym1::Element);
            PRINT(SparMatSym2::Element);
            
            PRINT(SparMatBlk::Line);
            PRINT(SparMatSymBlk::Column);
            PRINT(SparMatSymBlkDiag::Block);
            //PRINT(SparMatSymBlkDiag::Pilar);
            break;
    }
}
