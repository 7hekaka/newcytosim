// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

#ifndef GYM_VECT_H
#define GYM_VECT_H

#include "gym_cap.h"
#include "gym_view.h"
#include "vector1.h"
#include "vector2.h"
#include "vector3.h"

namespace gym
{
    /// set X and Y as orthogonal vectors to Z, forming a orthonormal basis
    void orthonormal(const float Z[3], float N, float X[3], float Y[3]);

    /// translate by V and then scale by S
    inline void translate(Vector1 const& V) { translate(V.XX, 0, 0); }
    /// translate by V and then scale by S
    inline void translate(Vector2 const& V) { translate(V.XX, V.YY, 0); }
    /// translate by V and then scale by S
    inline void translate(Vector3 const& V) { translate(V.XX, V.YY, V.ZZ); }

    /// translate by V and then scale by S
    inline void transScale(Vector1 const& V, float S) { transScale(V.XX, 0, 0, S); }
    /// translate by V and then scale by S
    inline void transScale(Vector2 const& V, float S) { transScale(V.XX, V.YY, 0, S); }
    /// translate by V and then scale by S
    inline void transScale(Vector3 const& V, float S) { transScale(V.XX, V.YY, V.ZZ, S); }

    /// translate by T; rotate to align X with A and Y with B
    void transRotate(Vector2 const& T, Vector2 const& A, Vector2 const& B);

    /// translate by T; rotate to align X with A, Y with B and Z with C
    void transRotate(Vector3 const& T, Vector3 const& A, Vector3 const& B, Vector3 const& C);

    /// translate by A; rotate to align Z with AB, Z replacing X. Scale in Z to put B at (0,0,1) Scale XY plane by `rad'
    void stretchAlignZ(Vector1 const& A, Vector1 const& B, float rad);
    /// translate by A; rotate to align Z with AB, Z replacing X. Scale in Z to put B at (0,0,1) Scale XY plane by `rad'
    void stretchAlignZ(Vector2 const& A, Vector2 const& B, float rad);
    /// translate by A; rotate to align Z with AB, Z replacing X. Scale XY plane by `rad'
    void stretchAlignZ(Vector3 const& A, Vector3 const& B, float rad);
    
    /// translate by pos; rotate to align Z with dir, scale XY plane by rad
    void transAlignZ(Vector1 const& pos, float rad, Vector1 const& dir);
    void transAlignZ(Vector2 const& pos, float rad, Vector2 const& dir);
    void transAlignZ(Vector3 const& pos, float rad, Vector3 const& dir);

    /// translate by pos; rotate to align Z with dir, given norm(dir)=1, scale XY by rad and Z by fac
    void stretchAlignZ1(Vector1 const& pos, float rad, Vector1 const& dir, float fac);
    void stretchAlignZ1(Vector2 const& pos, float rad, Vector2 const& dir, float fac);
    void stretchAlignZ1(Vector3 const& pos, float rad, Vector3 const& dir, float fac);

    void rotate(Vector3 const& A, Vector3 const& B, Vector3 const& C);
  
    void rotateInverse(Vector3 const& A, Vector3 const& B, Vector3 const& C);

#pragma mark - Clip Planes

    inline void setClipPlane(unsigned glp, Vector1 const dir, Vector1 const pos)
    {
        load_ref();
        setClipPlane(glp, dir.XX, 0, 0, -dot(dir, pos));
    }
    
    inline void setClipPlane(unsigned glp, Vector2 const dir, Vector2 const pos)
    {
        load_ref();
        setClipPlane(glp, dir.XX, dir.YY, 0, -dot(dir, pos));
    }
    
    inline void setClipPlane(unsigned glp, Vector3 const dir, Vector3 const pos)
    {
        load_ref();
        setClipPlane(glp, dir.XX, dir.YY, dir.ZZ, -dot(dir, pos));
    }

}

#endif
