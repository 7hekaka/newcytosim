// Cytosim 3.0 -  Copyright Francois Nedelec et al.  EMBL 2007-2013

#include "space_lid.h"
#include "exceptions.h"
#include "mecapoint.h"
#include "meca.h"

#include "random.h"


SpaceLid::SpaceLid(const SpaceProp* p)
: Space(p), ceiling(mLength[3])
{
    if ( DIM == 1 )
        throw InvalidParameter("lid is only valid in DIM=2 or 3");
    
    force = 0;
}


void SpaceLid::setModulo(Modulo& mod) const
{
    mod.enable(0, length(0));
}

void SpaceLid::resize()
{
    checkLengths(DIM, true);
}


Vector SpaceLid::extension() const
{
    return Vector(length(0), length(1), length(2));
}


/**
 place only at upper boundary. This overrides the function in spaces
 */
Vector SpaceLid::randomPlaceNearEdge(real radius, unsigned long) const
{
    return Vector( RNG.sfloat()*length(0), ceiling, 0 );
}


//------------------------------------------------------------------------------
#pragma mark -


real SpaceLid::volume() const
{
#if ( DIM == 1 )
    return ceiling + length(0);
#elif ( DIM == 2 )
    return 2.0 * length(0) * ( ceiling + length(1) );
#else
    return 4.0 * length(0) * length(1) * ( ceiling + length(2) );
#endif
}


bool  SpaceLid::inside(Vector const& point) const
{
    return ( point[DIM-1] < ceiling );
}


bool SpaceLid::allInside(Vector const& cen, const real rad) const
{
    assert_true( rad >= 0 );
    return ( cen[DIM-1] + rad < ceiling );
}


bool SpaceLid::allOutside(Vector const& cen, const real rad) const
{
    assert_true( rad >= 0 );
    return ( cen[DIM-1] + rad > ceiling );
}


Vector SpaceLid::project(Vector const& pos) const
{
#if ( DIM > 2 )
    return Vector(pos.XX, pos.YY, ceiling);
#else
    return Vector(pos.XX, ceiling, 0);
#endif
}


//------------------------------------------------------------------------------
#pragma mark -


void SpaceLid::setInteraction(Vector const& pos, Mecapoint const& pe,
                              Meca & meca, real stiff) const
{
    Matrix::index_t inx = DIM-1 + DIM * pe.matIndex();
    
    meca.mC(inx, inx) -= stiff;
    meca.base(inx)    += stiff * ceiling;
    
#if ( DIM == 2 )
    force += stiff * ( pos.YY - ceiling );
#elif ( DIM > 2 )
    force += stiff * ( pos.ZZ - ceiling );
#endif
}


void SpaceLid::setInteraction(Vector const& pos, Mecapoint const& pe, real rad,
                              Meca & meca, real stiff) const
{
    Matrix::index_t inx = DIM-1 + DIM * pe.matIndex();
    
    meca.mC(inx, inx) -= stiff;
    meca.base(inx)    += stiff * ( ceiling - rad );
    
#if ( DIM == 2 )
    force += stiff * ( pos.YY - ceiling + rad );
#elif ( DIM > 2 )
    force += stiff * ( pos.ZZ - ceiling + rad );
#endif
}

void SpaceLid::setInteractions(Meca& meca, FiberSet const&) const
{
    force = 0;
}

void SpaceLid::step()
{
    real dc = prop->mobility_dt * force;
    
    if ( fabs(dc) < 1 )
        ceiling += dc;
    else
        std::cerr << "Error: lid displacement is too fast: " << dc << '\n';
    std::cerr << "force on lid is " << force << '\n';
    
    if ( ceiling > length(DIM-1) )
        std::cerr << "Warning: space lid has reached its maximum\n";
}


//------------------------------------------------------------------------------
//                         OPENGL  DISPLAY
//------------------------------------------------------------------------------
#pragma mark -

#ifdef DISPLAY
#include "opengl.h"
#include "gle.h"
using namespace gle;

bool SpaceLid::draw() const
{
    const real X = length(0);
    const real Y = ( DIM > 1 ) ? length(1) : 1;
    const real Z = ( DIM > 2 ) ? length(2) : 0;

#if ( DIM == 2 )
    glBegin(GL_LINES);
    gleVertex( -X, ceiling, Z );
    gleVertex(  X, ceiling, Z );
    glEnd();
#elif ( DIM > 2 )
    glBegin(GL_TRIANGLE_STRIP);
    gleVertex( -X, -Y, ceiling );
    gleVertex(  X, -Y, ceiling );
    gleVertex( -X,  Y, ceiling );
    gleVertex(  X,  Y, ceiling );
    glEnd();
#endif
    
    glLineWidth(3);
    glLineStipple(1, 0x0303);
    glEnable(GL_LINE_STIPPLE);
    glBegin(GL_LINE_LOOP);
    gleVertex( -X, -Y, Z );
    gleVertex(  X, -Y, Z );
    gleVertex(  X,  Y, Z );
    gleVertex( -X,  Y, Z );
    gleVertex( -X, -Y, Z );
    glEnd();
    
#if ( DIM > 2 )
    glBegin(GL_LINE_LOOP);
    gleVertex( -X, -Y, -Z );
    gleVertex(  X, -Y, -Z );
    gleVertex(  X,  Y, -Z );
    gleVertex( -X,  Y, -Z );
    gleVertex( -X, -Y, -Z );
    glEnd();

    glBegin(GL_LINES);
    gleVertex(  X,  Y, -Z );
    gleVertex(  X,  Y,  Z );
    gleVertex( -X,  Y, -Z );
    gleVertex( -X,  Y,  Z );
    gleVertex( -X, -Y, -Z );
    gleVertex( -X, -Y,  Z );
    gleVertex(  X, -Y, -Z );
    gleVertex(  X, -Y,  Z );
    glEnd();
#endif
  
    glDisable(GL_LINE_STIPPLE);

    return true;
}

#else

bool SpaceLid::draw() const
{
    return false;
}

#endif

