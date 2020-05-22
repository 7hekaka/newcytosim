// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "simul.h"

#include "matrix11.h"
#include "matrix22.h"
#include "matrix33.h"
#include "matrix34.h"
#include "matrix44.h"


#define print_sizeof(arg) printf("sizeof %16s   %lu bytes\n", #arg, sizeof(arg));

int main(int argc, char* argv[])
{
    print_sizeof(Array<int>);
    print_sizeof(Vector1);
    print_sizeof(Vector2);
    print_sizeof(Vector3);
    print_sizeof(Vector4);
    
    print_sizeof(Matrix11);
    print_sizeof(Matrix22);
    print_sizeof(Matrix33);
    print_sizeof(Matrix34);
    print_sizeof(Matrix44);
    
    print_sizeof(Mecapoint);
    print_sizeof(FiberSegment);
    print_sizeof(Interpolation);
    print_sizeof(FiberSite);
    print_sizeof(Hand);
    print_sizeof(Single);
    print_sizeof(Couple);

    print_sizeof(Space);
    print_sizeof(Node);
    print_sizeof(Inventoried);
    print_sizeof(Movable);
    print_sizeof(Object);

    print_sizeof(Mecable);
    print_sizeof(Chain);
    print_sizeof(Mecafil);
    print_sizeof(Fiber);

    print_sizeof(Solid);
    print_sizeof(Bead);
    print_sizeof(Sphere);
    print_sizeof(Simul);

    return EXIT_SUCCESS;
}
