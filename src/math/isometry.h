// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef ISOMETRY_H
#define ISOMETRY_H

#include "dim.h"
#include "vector.h"

#if ( DIM == 1 )

   #include "matrix11.h"
   typedef Matrix11 MatrixD;

#elif ( DIM == 2 )

   #include "matrix22.h"
   typedef Matrix22 MatrixD;

#elif ( DIM == 3 )

   #include "matrix33.h"
   typedef Matrix33 MatrixD;

#else

   #include "matrix44.h"
   typedef Matrix44 MatrixD;

#endif


/// A Rotation is a matrix of dimension DIM x DIM
typedef MatrixD Rotation;


/// An affine transformation in space.
/**
 A Isometry contains a vector T and a rotation matrix M,
 and represents the affine transformation:
 
     X -> M.X + T
 
 */
class Isometry
{
public:
    
    /// rotation component
    MatrixD rot;
    
    /// translation component
    Vector  mov;

public:
    
    Isometry()
    {
        rot = MatrixD::one();
        mov.reset();
    }

    Isometry(Vector const& v)
    {
        rot = MatrixD::one();
        mov = v;
    }

    Isometry(MatrixD const& r, Vector const& v)
    {
        rot = r;
        mov = v;
    }

    bool invalid()
    {
        return ! mov.valid();
    }

    void reset()
    {
        mov.reset();
        rot = MatrixD::one();
    }

    void inverse()
    {
        rot.inverse();
        mov = -( rot * mov );
    }

    /// allow automatic conversion to a Vector
    operator Vector const& () const
    {
        return mov;
    }
    
    /// allow automatic conversion to a Rotation matrix
    operator MatrixD const& () const
    {
        return rot;
    }
    
    /// apply translation, after *this
    void translate(Vector const& v)
    {
        mov += v;
    }
    
    /// apply rotation, after *this
    void rotate(MatrixD const& mat)
    {
        rot = mat * rot;
        mov = mat * mov;
    }

    /// apply another isometry, after *this
    void combine(Isometry const& iso)
    {
        mov = iso.rot * mov + iso.mov;
        rot = iso.rot * rot;
    }
};


/// output operator
inline std::ostream& operator << (std::ostream& os, Isometry const& iso)
{
#if ( DIM > 2 )
    real angle = iso.rot.rotationAngle();
    Vector axis = iso.rot.rotationAxis();
    os << "Iso { " << iso.mov << " | " << angle << " axis " << axis << " }";
#else
    os << "Iso { " << iso.mov << " | " << iso.rot << " }";
#endif
    return os;
}


#endif
