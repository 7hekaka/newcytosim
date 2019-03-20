// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef ISOMETRY_H
#define ISOMETRY_H

#include "dim.h"


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
    Vector  vec;

public:
    
    Isometry()
    {
        vec.reset();
        rot = MatrixD::identity();
    }

    Isometry(Vector const& v)
    {
        vec = v;
        rot = MatrixD::identity();
    }

    Isometry(Vector const& v, MatrixD const& r)
    {
        vec = v;
        rot = r;
    }

    void reset()
    {
        vec.reset();
        rot = MatrixD::identity();
    }

    /// allow automatic conversion to a Vector
    operator Vector const& () const
    {
        return vec;
    }
    
    /// allow automatic conversion to a Rotation matrix
    operator MatrixD const& () const
    {
        return rot;
    }
    
    void translate(Vector const& v)
    {
        vec += v;
    }
    
    void rotate(MatrixD const& mat)
    {
        rot = mat * rot;
        vec = mat * vec;
    }

    void combine(Isometry const& iso)
    {
        vec = iso.rot * vec + iso.vec;
        rot = iso.rot * rot;
    }
};


/// output operator
inline std::ostream& operator << (std::ostream& os, Isometry const& iso)
{
    os << "Isometry { " << iso.vec << " | " << iso.rot << " }";
    return os;
}


#endif
