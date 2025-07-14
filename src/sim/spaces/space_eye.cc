// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "assert_macro.h"
#include "space_eye.h"
#include "exceptions.h"


SpaceEye::SpaceEye(const SpaceProp* p)
: Space(p)
{
    if ( DIM != 2 )
        throw InvalidParameter("eye is only valid in 2D at the moment");
}


void SpaceEye::resize()
{
    Space::checkLengths(2, true);
    
    mLength[2] = atan(mLength[0] / mLength[1]);
    mLength[3] = sqrt( ( mLength[0] * mLength[0] + mLength[1] * mLength[1] )) / (2 * cos(mLength[2]) );
}


Vector SpaceEye::extension() const
{
    return Vector( mLength[0], mLength[1], 0.0 );
}

#if (DIM == 1)  ///not correct in 1D: unfinished!!

real SpaceEye::volume() const
{
    throw Exception("SpaceEye is not valid in 1D");
    return 0;
}

bool  SpaceEye::inside( const real w[] ) const
{
    throw Exception("SpaceEye is not valid in 1D");
    return true;
}

void SpaceEye::project( const real w[], real p[] ) const
{
    throw Exception("SpaceEye is not valid in 1D");
}

#endif


//------------------------------------------------------------------------------

#if (DIM == 2)

real SpaceEye::volume() const
{
    return ( (M_PI - 2 * mLength[2])* mLength[3]*mLength[3] - 2*mLength[0]*(mLength[3]-mLength[1] ) );
}


bool  SpaceEye::inside( const real w[] ) const
{
    
    if ( fabs(w[0]) > mLength[0] ) return false;
    else if ( fabs(w[1]) > mLength[1] ) return false;
    else if ( w[0]*w[0] + (fabs(w[1]) - mLength[1] + mLength[3]) * (fabs(w[1]) - mLength[1] + mLength[3]) > mLength[3]*mLength[3] ) return false;
    else return true;
    
}

void SpaceEye::project( const real w[], real p[] ) const
{
    p[0] = w[0];
    p[1] = w[1];
    
    if (p[1] >= 0)
    {
        p[1] += mLength[3] - mLength[1];
        real factor = mLength[3]/sqrt(p[0]*p[0] + p[1]*p[1]);
        p[0] *= factor;
        p[1] *= factor;
        p[1] -= mLength[3] - mLength[1];
        if ( p[1] <=0 )
        {
            if (p[0] > 0.0) p[0] =  mLength[0];
            else p[0] = -1 * mLength[0];
            p[1] =  0.0;
        }
    }
    else
    {
        p[1] -= mLength[3] - mLength[1];
        real factor = mLength[3]/sqrt(p[0]*p[0] + p[1]*p[1]);
        p[0] *= factor;
        p[1] *= factor;
        p[1] += mLength[3] - mLength[1];
        if ( p[1] >=0 )
        {
            if (p[0] > 0.0) p[0] =  mLength[0];
            else p[0] = -1 * mLength[0];
            p[1] =  0.0;
        }
    }
}

#endif

//------------------------------------------------------------------------------

#if (DIM == 3)

///not correct in 3D:  unfinished!!
real SpaceEye::volume() const
{
    throw Exception("SpaceEye is not yet valid in 3D");
    return 0;
}



bool  SpaceEye::inside( const real w[] ) const
{
    throw Exception("SpaceEye is not yet valid in 3D");
    return true;
}

void SpaceEye::project( const real w[], real p[] ) const
{
    throw Exception("SpaceEye is not yet valid in 3D");
}

#endif

//------------------------------------------------------------------------------
//                         OPENGL  DISPLAY
//------------------------------------------------------------------------------

#ifdef DISPLAY
#include "opengl.h"
#include "gle.h"
using namespace gle;

bool SpaceEye::display() const
{
    // printf("%f %f %f %f\n",mLength[0],mLength[1],mLength[2],mLength[3]);
    glBegin(GL_LINE_LOOP);
    for ( real aa = 2*mLength[2]-M_PI/2; aa <= 3*M_PI/2 -2*mLength[2] ; aa += 0.05 )
        gleVertex( mLength[3]*cos(aa), mLength[3]*(sin(aa)-1)+mLength[1], 0 );
    for ( real aa = -1*(3*M_PI/2 -2*mLength[2]); aa <= -1*(2*mLength[2]-M_PI/2); aa += 0.05 )
        gleVertex( mLength[3]*cos(aa), mLength[3]*(sin(aa)+1)-mLength[1], 0 );
    glEnd();
    
    return true;
    
}

#else

bool SpaceEye::display() const
{
    return false;
}


#endif


