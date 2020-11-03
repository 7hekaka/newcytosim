// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "simul.h"

#include "matrix11.h"
#include "matrix22.h"
#include "matrix33.h"
#include "matrix34.h"
#include "matrix44.h"


#define PRINT(arg) printf("sizeof %16s   %lu bytes\n", #arg, sizeof(arg));

int main(int argc, char* argv[])
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
    
    PRINT(Random);
    PRINT(Array<int>);

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

    return EXIT_SUCCESS;
}
