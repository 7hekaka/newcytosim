// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University.
// Created by F. Nedelec, Oct 2002


#ifndef QUATERNION_H
#define QUATERNION_H

#include "assert_macro.h"
#include "random.h"
#include <cmath>
#include <cstdio>
#include <iostream>

/** 
 Quaternions extend complex numbers to dimension 4.
 
 http://en.wikipedia.org/wiki/Quaternion

 The unit bases of the Quaternion space are called 1, i, j and k.
 A quaternion is therefore q[0] + i * q[1] + j * q[2] + k * q[3],
 where q[?] are four real scalars.
 q[0] is the real part, and the other parts are imaginary.
 
 While addition is standard, multiplication is anti-commutative:
 - i*i = -1,
 - i*j =  k,
 - j*i = -k,
 - etc,
 .
 
 Unit quaternions are handy to represent rotations in 3D space:
 The group of 3D rotations has three degrees of freedom,
 and is mapped directly onto the quaternions of norm 1.
 
 A rotation is represented by a symmetric matrix using 6 scalar numbers, 
 but only 4 scalars are used when a Quaternion is used.
 Thus one economize spurious scalars, and moreover quaternions are easier 
 to normalize than rotation matrices as necessary to correct for numerical errors.
 
 The rotation associated to a unit quaternion Q is:
 
     v -> Q.v.inv(Q)
 
 where the imaginary quaternion v = { 0, x, y, z } represents a 3D vector { x, y, z }.
 
 Note that it is more costly to calculate a rotated vector using this formula
 than with a 3x3 matrix-vector multiplication.
 
 The composition of two rotations thus corresponds to quaternion multiplication.
 For example, Q*P corresponds to the rotation P followed by the rotation Q.
 Thus 1/Q is the rotation that is inverse to the rotation associated with Q.
 
 The angle A of the rotation associated with the quaternion Q obeys:
 - real part of Q = std::cos(A/2),
 - norm of imaginary part of Q = std::sin(A/2).
 The rotation axis is defined by the imaginary components of Q.
 
 Quaternion<real> implements the standard mathematical operations, 
 conversions to and from 3x3 real matrices, 
 and to 4x4 transformation matrices used in OpenGL.
 */


/// a Quaternion is similar to a complex number, but in dimension four
template <typename REAL>
class Quaternion 
{

private:
    
    /// The four coordinates of a Quaternion, real part first
    /** this represents q[0] + i * q[1] + j * q[2] + k * q[3] */
    REAL q[4];
    
public:
    
    /// The default constructor does not reset any value
    Quaternion() {}
    
    /// Constructor which can be used to convert from a real
    Quaternion(REAL a, REAL b, REAL c, REAL d)
    {
        q[0] = a;
        q[1] = b;
        q[2] = c;
        q[3] = d;
    }
    
    /// Destructor (ATTENTION: non-virtual: do not derive from this class)
    ~Quaternion() {}
    
    /// setting the values from Cartesian coordinates
    void set(REAL a, REAL b, REAL c, REAL d)
    {
        q[0] = a;
        q[1] = b;
        q[2] = c;
        q[3] = d;
    }
    
    /// access to a modifiable coordinate
    REAL& operator[] (size_t n) { return q[n]; }
    
    /// access to a non-modifiable coordinate 
    REAL  operator[] (size_t n) const  { return q[n]; };
    
    /// conversion operator to a "real array"
    operator REAL*() { return q; }
    
    /// conversion to a 'real array'
    REAL *    data() { return q; }
    
    /// opposition: change sign in all coordinates 
    Quaternion operator - () const
    {
        return Quaternion(-q[0], -q[1], -q[2], -q[3]);
    }
    
    /// multiply by a real value
    Quaternion operator * (REAL f) const
    {
        return Quaternion(q[0]*f, q[1]*f, q[2]*f, q[3]*f);
    }
    
    /// divide by a real value
    Quaternion operator / (REAL f) const
    {
        return Quaternion(q[0]/f, q[1]/f, q[2]/f, q[3]/f);
    }
    
    /// add a real value in place
    void operator += (REAL f)
    {
        q[0] += f;
    }
    
    /// subtract a real value in place
    void operator -= (REAL f)
    {
        q[0] -= f;
    }
    
    /// multiply for a real value in place
    void operator *= (REAL f)
    {
        q[0] *= f;
        q[1] *= f;
        q[2] *= f;
        q[3] *= f;
    }
    
    /// divide by a real value in place
    void operator /= (REAL f)
    {
        q[0] /= f;
        q[1] /= f;
        q[2] /= f;
        q[3] /= f;
    }
    
    /// sum two quaternions
    const Quaternion  operator + (const Quaternion & X) const
    {
        return Quaternion(q[0]+X[0], q[1]+X[1], q[2]+X[2], q[3]+X[3]);
    }
    
    /// subtract two quaternions
    const Quaternion  operator - (const Quaternion & X) const
    {
        return Quaternion(q[0]-X[0], q[1]-X[1], q[2]-X[2], q[3]-X[3]);
    }
    
    /// add another quaternion in place
    void operator += (const Quaternion & X)
    {
        q[0] += X[0];
        q[1] += X[1];
        q[2] += X[2];
        q[3] += X[3];
    }
    
    /// subtract a quaternion in place
    void operator -= (const Quaternion & X)
    {
        q[0] -= X[0];
        q[1] -= X[1];
        q[2] -= X[2];
        q[3] -= X[3];
    }
    
    /// multiplication from the right side
    void operator *= (const Quaternion & X)
    {
        rightMult(X);
    }
    
    /// divide in place by another quaternion
    void operator /= (const Quaternion & X)
    {
        rightMult( X.inverted() );
    }
    
    /// multiplication between quaternions
    const Quaternion operator * (const Quaternion & X) const
    {
        Quaternion result(q[0], q[1], q[2], q[3]);
        result.rightMult(X);
        return result;
    }
    
    /// division between quaternions
    const Quaternion operator / (const Quaternion & X) const
    {
        Quaternion  result(q[0], q[1], q[2], q[3]);
        result.rightMult(X.inverted());
        return result;
    }
    
    /// extract the square of the norm, i.e. norm*norm
    REAL normSqr() const
    {
        return q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3];
    }
    
    /// extract the norm 
    REAL norm() const
    {
        return std::sqrt( q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3] );
    }
    
    /// return the normalized quaternion
    const Quaternion normalized(REAL n = 1.0) const
    {
        REAL s = n / norm();
        return Quaternion(q[0]*s, q[1]*s, q[2]*s, q[3]*s );
    }
    
    /// return the normalized Quaternion
    friend Quaternion normalize(const Quaternion & X)
    {
        REAL s = 1.0 / X.norm();
        return Quaternion(X[0]*s, X[1]*s, X[2]*s, X[3]*s );
    }

    /// scale in place to obtain norm = `n`
    void normalize(REAL n)
    {
        REAL s = n / norm();
        q[0] *= s;
        q[1] *= s;
        q[2] *= s;
        q[3] *= s;
    }
    
    /// scale to reach a norm of 1
    void normalize()
    {
        REAL s = normSqr();
        if ( s > REAL_EPSILON )
            *this *= 1.0 / std::sqrt(s);
        else
            set(1, 0, 0, 0);
    }
    
    /// conjugated quaternion +, -, -, -
    const Quaternion conjugated() const
    {
        return Quaternion(q[0], -q[1], -q[2], -q[3]);
    }
    
    /// conjugate in place +, -, -, -
    void conjugate()
    {
        q[1] = -q[1];
        q[2] = -q[2];
        q[3] = -q[3];
    }
    
    /// inversed quaternion:  1/*this
    const Quaternion inverted() const
    {
        REAL x = -normSqr();
        return Quaternion( -q[0]/x, q[1]/x, q[2]/x, q[3]/x );
    }
    
    /// inverse in place
    void inverse()
    {
        REAL x = -normSqr();
        q[0] /= -x;
        q[1] /=  x;
        q[2] /=  x;
        q[3] /=  x;
    }
    
    /// the opposed quaternion:  -*this
    const Quaternion opposed() const
    {
        return Quaternion( -q[0], -q[1], -q[2], -q[3] );
    }
    
    /// oppose in place
    void negate()
    {
        q[0] = -q[0];
        q[1] = -q[1];
        q[2] = -q[2];
        q[3] = -q[3];
    }
    
    /// this * this
    const Quaternion squared() const
    {
        return Quaternion(q[0]*q[0] - q[1]*q[1] - q[2]*q[2] - q[3]*q[3],
                          2*q[0]*q[1],
                          2*q[0]*q[2],
                          2*q[0]*q[3]);
    }
    
    /// this = this * this
    void square()
    {
        REAL a = q[0], b = q[1], c = q[2], d = q[3];
        q[0] = a*a - b*b - c*c - d*d;
        a += a;
        q[1] = a*b;
        q[2] = a*c;
        q[3] = a*d;
    }
    
    /// multiplication from the right side (this <= this * X )
    void rightMult(const Quaternion & X)
    {
        REAL q0 = q[0], q1 = q[1], q2 = q[2], q3 = q[3];
        
        q[0] = q0 * X[0] - q1 * X[1] - q2 * X[2] - q3 * X[3];
        q[1] = q0 * X[1] + q1 * X[0] + q2 * X[3] - q3 * X[2];
        q[2] = q0 * X[2] - q1 * X[3] + q2 * X[0] + q3 * X[1];
        q[3] = q0 * X[3] + q1 * X[2] - q2 * X[1] + q3 * X[0];
    }
    
    /// multiplication from the left side (this <= X * this )
    void leftMult(const Quaternion & X)
    {
        REAL q0 = q[0], q1 = q[1], q2 = q[2], q3 = q[3];
        
        q[0] = q0 * X[0] - q1 * X[1] - q2 * X[2] - q3 * X[3];
        q[1] = q0 * X[1] + q1 * X[0] - q2 * X[3] + q3 * X[2];
        q[2] = q0 * X[2] + q1 * X[3] + q2 * X[0] - q3 * X[1];
        q[3] = q0 * X[3] - q1 * X[2] + q2 * X[1] + q3 * X[0];
    }
    
    /// multiplication from the right side, different implementation
    void rightMult_fast(const Quaternion & X)
    {
        REAL E = (q[3] + q[1]) * (X[1] + X[2]);
        REAL F = (q[3] - q[1]) * (X[1] - X[2]);
        REAL G = (q[0] + q[2]) * (X[0] - X[3]);
        REAL H = (q[0] - q[2]) * (X[0] + X[3]);
        REAL A = F - E;
        REAL B = F + E;
        REAL C = (q[0] + q[1]) * (X[0] + X[1]);
        REAL D = (q[0] - q[1]) * (X[2] + X[3]);
        E = (q[3] + q[2]) * (X[0] - X[1]);
        F = (q[3] - q[2]) * (X[2] - X[3]);
        q[0] = F + (A + G + H) * 0.5;
        q[1] = C + (A - G - H) * 0.5;
        q[2] = D + (B + G - H) * 0.5;
        q[3] = E + (B - G + H) * 0.5;
    }
    
    /// multiplication from the left side, different implementation
    void leftMult_fast(const Quaternion & X)
    {
        REAL E = (X[3] + X[1])*(q[1] + q[2]);
        REAL F = (X[3] - X[1])*(q[1] - q[2]);
        REAL G = (X[0] + X[2])*(q[0] - q[3]);
        REAL H = (X[0] - X[2])*(q[0] + q[3]);
        REAL A = F - E;
        REAL B = F + E;
        REAL C = (X[0] + X[1])*(q[0] + q[1]);
        REAL D = (X[0] - X[1])*(q[2] + q[3]);
        E = (X[3] + X[2])*(q[0] - q[1]);
        F = (X[3] - X[2])*(q[2] - q[3]);
        q[0] = F + (A + G + H) * 0.5;
        q[1] = C + (A - G - H) * 0.5;
        q[2] = D + (B + G - H) * 0.5;
        q[3] = E + (B - G + H) * 0.5;
    }
    
    
    /// generate the associated 3x3 rotation matrix for unit Quaternion
    /** This assumes that norm(*this) = 1 */
    void setMatrix3(REAL mat[], int ldd) const
    {
        REAL rx, ry, rz, xx, yy, yz, xy, xz, zz, x2, y2, z2;
        
        x2 = q[1] + q[1];
        y2 = q[2] + q[2];
        z2 = q[3] + q[3];
        
        rx = q[0] * x2; ry = q[0] * y2; rz = q[0] * z2;
        xx = q[1] * x2; xy = q[1] * y2; xz = q[1] * z2;
        yy = q[2] * y2; yz = q[2] * z2; zz = q[3] * z2;
        
        mat[0      ] = 1.0 - (yy + zz);
        mat[1      ] = xy + rz;
        mat[2      ] = xz - ry;
        
        mat[0+ldd  ] = xy - rz;
        mat[1+ldd  ] = 1.0 - (xx + zz);
        mat[2+ldd  ] = yz + rx;
        
        mat[0+ldd*2] = xz + ry;
        mat[1+ldd*2] = yz - rx;
        mat[2+ldd*2] = 1.0 - (xx + yy);
    }
    
    /// generate the associated 3x3 rotation matrix for unit Quaternion
    /** This assumes that norm(*this) = 1 */
    template < typename Matrix >
    void setMatrix3(Matrix & mat) const
    {
        REAL rx, ry, rz, xx, yy, yz, xy, xz, zz, x2, y2, z2;
        
        x2 = q[1] + q[1];
        y2 = q[2] + q[2];
        z2 = q[3] + q[3];
        
        rx = q[0] * x2; ry = q[0] * y2; rz = q[0] * z2;
        xx = q[1] * x2; xy = q[1] * y2; xz = q[1] * z2;
        yy = q[2] * y2; yz = q[2] * z2; zz = q[3] * z2;
        
        mat(0,0) = 1.0 - (yy + zz);
        mat(1,0) = xy + rz;
        mat(2,0) = xz - ry;
        
        mat(0,1) = xy - rz;
        mat(1,1) = 1.0 - (xx + zz);
        mat(2,1) = yz + rx;
        
        mat(0,2) = xz + ry;
        mat(1,2) = yz - rx;
        mat(2,2) = 1.0 - (xx + yy);
    }
    
    /// Rotate a 3D vector: des = Q * src * Q.conjugated()
    /** This assumes that norm(*this) = 1 */
    void rotateVector(REAL des[3], const REAL src[3]) const
    {
        REAL two(2.0);
        
        REAL rx =  q[0]*q[1];
        REAL ry =  q[0]*q[2];
        REAL rz =  q[0]*q[3];
        REAL xx = -q[1]*q[1];
        REAL xy =  q[1]*q[2];
        REAL xz =  q[1]*q[3];
        REAL yy = -q[2]*q[2];
        REAL yz =  q[2]*q[3];
        REAL zz = -q[3]*q[3];
        
        des[0] = two * ((yy + zz)*src[0] + (xy - rz)*src[1] + (ry + xz)*src[2] ) + src[0];
        des[1] = two * ((rz + xy)*src[0] + (xx + zz)*src[1] + (yz - rx)*src[2] ) + src[1];
        des[2] = two * ((xz - ry)*src[0] + (rx + yz)*src[1] + (xx + yy)*src[2] ) + src[2];
    }
    
    
    /// set from given 3x3 rotation matrix `m`
    void setFromMatrix3(const REAL m[9])
    {
        REAL  s, trace = m[0] + m[4] + m[8];
        
        // check the diagonal
        if (trace > 0 ) {
            s = std::sqrt( trace + 1.0 );
            q[0] = s * 0.5;
            s = 0.5 / s;
            q[1] = (m[5] - m[7]) * s;
            q[2] = (m[6] - m[2]) * s;
            q[3] = (m[1] - m[3]) * s;
        }
        else {
            // trace is negative
            // find biggest coefficient on diagonal:
            int i = 0;
            if (m[1+3*1] > m[0+3*0]) i = 1;
            if (m[2+3*2] > m[i+3*i]) i = 2;
            
            s = std::sqrt( 1.0 + 2*m[i+3*i] - trace );
            q[i+1] = s * 0.5;
            if (s != 0) s = 0.5 / s;
            int j = (i+1) % 3;
            int k = (j+1) % 3;
            q[j+1] = s * ( m[j+3*i] + m[i+3*j] );
            q[k+1] = s * ( m[i+3*k] + m[k+3*i] );
            q[0]   = s * ( m[k+3*j] - m[j+3*k] );
        }
    }
    
    /// set transformation matrix: translation by T, scaling by S, rotation
    void setOpenGLMatrix(float m[16], double S, const float T[3]) const
    {
        //this code assumes that the quaternion has norm = 1,
        double x2 = S * ( q[1] + q[1] );
        double y2 = S * ( q[2] + q[2] );
        double z2 = S * ( q[3] + q[3] );
        
        double rx = q[0] * x2, ry = q[0] * y2, rz = q[0] * z2;
        double xx = q[1] * x2, xy = q[1] * y2, xz = q[1] * z2;
        double yy = q[2] * y2, yz = q[2] * z2, zz = q[3] * z2;

        m[0+4*0] = float( S - (yy + zz) );
        m[1+4*0] = float( xy + rz );
        m[2+4*0] = float( xz - ry );
        m[3+4*0] = 0.f;

        m[0+4*1] = float( xy - rz );
        m[1+4*1] = float( S - (xx + zz) );
        m[2+4*1] = float( yz + rx );
        m[3+4*1] = 0.f;

        m[0+4*2] = float( xz + ry );
        m[1+4*2] = float( yz - rx );
        m[2+4*2] = float( S - (xx + yy) );
        m[3+4*2] = 0.f;

        m[0+4*3] = T[0];
        m[1+4*3] = T[1];
        m[2+4*3] = T[2];
        m[3+4*3] = 1.f;
    }

    /// set transformation matrix: translation by T, scaling by S, rotation
    void setOpenGLMatrix(double m[16], double S, const double T[3]) const
    {
        //this code assumes that the quaternion has norm = 1,
        
        double x2 = S * ( q[1] + q[1] );
        double y2 = S * ( q[2] + q[2] );
        double z2 = S * ( q[3] + q[3] );
        
        double rx = q[0] * x2, ry = q[0] * y2, rz = q[0] * z2;
        double xx = q[1] * x2, xy = q[1] * y2, xz = q[1] * z2;
        double yy = q[2] * y2, yz = q[2] * z2, zz = q[3] * z2;
        
        m[0+4*0] = S - (yy + zz);
        m[1+4*0] = xy + rz;
        m[2+4*0] = xz - ry;
        m[3+4*0] = 0.0;

        m[0+4*1] = xy - rz;
        m[1+4*1] = S - (xx + zz);
        m[2+4*1] = yz + rx;
        m[3+4*1] = 0.0;

        m[0+4*2] = xz + ry;
        m[1+4*2] = yz - rx;
        m[2+4*2] = S - (xx + yy);
        m[3+4*2] = 0.0;

        m[0+4*3] = T[0];
        m[1+4*3] = T[1];
        m[2+4*3] = T[2];
        m[3+4*3] = 1.0;
    }
    
    /// set from polar coordinates (r, phi, theta, psi)
    void setFromPolar(const REAL v[4])
    {
        REAL a = v[0] * std::sin(v[1]);
        q[0]   = v[0] * std::cos(v[1]);   //r*std::cos(phi)
        REAL b = a * std::sin(v[2]);
        q[1]   = a * std::cos(v[2]);      //r*std::sin(phi)*std::cos(theta)
        q[2]   = b * std::cos(v[3]);      //r*std::sin(phi)*std::sin(theta)*std::cos(psi)
        q[3]   = b * std::sin(v[3]);      //r*std::sin(phi)*std::sin(theta)*std::sin(psi)
    }
    
    /// return new quaternion with polar coordinates (r, phi, theta, psi)
    static const Quaternion newFromPolar(const REAL v[4])
    {
        Quaternion result;
        result.setFromPolar(v);
        return result;
    }
    
    /// calculate the polar coordinates (r, phi, theta, psi)
    void getPolar(REAL v[4]) const
    {
        // r,  phi, theta, psi
        v[0] = norm();
        v[1] = std::acos(q[0] / v[0]);
        v[2] = std::acos(q[1] / (v[0] * std::sin(v[1])));
        v[3] = std::atan2(q[3], q[2]);
    }
    

    /// set as rotation of axis v, with angle = v.norm() in radian;
    void setFromAxis(const REAL v[3])
    {
        /** for small angles, we assume here angle ~ v.norm() */
        REAL n = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
        REAL s = std::sin(n * 0.5);
        if ( n > 0 )
            s /= n;
        q[0] = std::cos(n * 0.5);
        q[1] = v[0] * s;
        q[2] = v[1] * s;
        q[3] = v[2] * s;
    }
    
    /// set from rotation of axis v, and angle 'angle' in radian around this axis
    /** argument `v` is normalized for more security */
    void setFromAxis(const REAL v[3], REAL angle)
    {
        REAL n = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
        REAL s = std::sin(angle * 0.5) / n;
        q[0] = std::cos(angle * 0.5);
        q[1] = v[0] * s;
        q[2] = v[1] * s;
        q[3] = v[2] * s;
    }
    
    /// set as rotation of angle 'angle' and axis X, Y or Z (axis=0,1,2)
    /** along one of the unit axis specified by `axis`: ( 0: X, 1: Y, 2: Z ) */
    void setFromPrincipalAxis(int axis, REAL angle)
    {
        REAL a = angle * 0.5;
        q[0] = std::cos(a);
        q[1] = 0.0;
        q[2] = 0.0;
        q[3] = 0.0;
        q[axis+1] = std::sin(a);
    }
    
    /// return angle of the rotation
    REAL getAngle() const
    {
        REAL n = std::sqrt( q[1]*q[1] + q[2]*q[2] + q[3]*q[3] );
        return  2 * std::atan2(n, q[0]);
    }
    
    /// compute the axis and return the angle of the rotation
    REAL getAngle(REAL v[3]) const
    {
        REAL n = std::sqrt( q[1]*q[1] + q[2]*q[2] + q[3]*q[3] );
        if ( n > 0 )
        {
            REAL a = 2 * std::atan2(n, q[0]);
            n = 1.0 / n;
            v[0] = q[1] * n;
            v[1] = q[2] * n;
            v[2] = q[3] * n;
            return a;
        }
        v[0] = 0;
        v[1] = 0;
        v[2] = 1;
        return 0;
    }
    
    /// compute the axis and return the angle of the rotation
    void getAxis(REAL v[3]) const
    {
        REAL n = std::sqrt( q[1]*q[1] + q[2]*q[2] + q[3]*q[3] );
        if ( n > 0 )
        {
            n = 1.0 / n;
            v[0] = q[1] * n;
            v[1] = q[2] * n;
            v[2] = q[3] * n;
        }
        else
        {
            v[0] = 0;
            v[1] = 0;
            v[2] = 1;
        }
    }
    
    /// multiply the angle of the rotation by `s`
    const Quaternion scaledAngle(REAL s) const
    {
        REAL n = std::sqrt( q[1]*q[1] + q[2]*q[2] + q[3]*q[3] );
        if ( n > 0 )
        {
            REAL a = s * std::atan2(n, q[0]);
            n = std::sin(a) / n;
            return Quaternion(std::cos(a), n*q[1], n*q[2], n*q[3]);
        }
        return Quaternion(1, 0, 0, 0);
    }
    
    /// Linear interpolation between rotations 'this' and 'b'.
    const Quaternion slerp(const Quaternion & X, const REAL u) const
    {
        // code from Jonathan Blow
        // Calculate the cosine of the angle between the two vectors
        REAL dot = q[0]*X[0] + q[1]*X[1] + q[2]*X[2] + q[3]*X[3];
        
        // If the angle is significant, use the spherical interpolation
        if ( dot > 0.9995 ) {
            // use cheap linear interpolation
            return normalize( (*this) + (X - (*this))*u );
        }
        
        REAL tmp = std::acos(dot) * u;
        //build v2 ortogonal to *this:
        Quaternion v2 = normalize( X - (*this)*dot );
        return (*this)*std::cos(tmp) + v2*std::sin(tmp);
    }
    
    /// printf
    void print(FILE* out = stdout, bool parenthesis = false) const
    {
        if ( parenthesis )
            fprintf( out, "( %+6.3f %+6.3f %+6.3f %+6.3f )", q[0], q[1], q[2], q[3]);
        else
            fprintf( out, "  %+6.3f %+6.3f %+6.3f %+6.3f", q[0], q[1], q[2], q[3]);
    }
    
    /// printf with a new-line
    void println(FILE* out = stdout, bool parenthesis = false) const
    {
        if ( parenthesis )
            fprintf( out, "( %+6.3f %+6.3f %+6.3f %+6.3f )\n", q[0], q[1], q[2], q[3]);
        else
            fprintf( out, "  %+6.3f %+6.3f %+6.3f %+6.3f\n", q[0], q[1], q[2], q[3]);
    }
    
    /// Human friendly ouput
    void print(std::ostream& os) const
    {
        os << q[0] << " " << q[1] << " " << q[2] << " " << q[3];
    }

    
#ifdef RANDOM_H
    
    /// returns a quaternion, uniformly sampling all possible rotations
    /** James Arvo, Fast random rotation matrices. in Graphics Gems 3. */
    static const Quaternion randomRotation()
    {
        REAL u1 = RNG.preal();
        REAL u2 = M_PI*RNG.sreal();
        REAL u3 = M_PI*RNG.sreal();
        REAL s1 = std::sqrt(1-u1), s2 = std::sqrt(u1);
        
        return Quaternion<REAL>(s1*std::sin(u2), s1*std::cos(u2), s2*std::sin(u3), s2*std::cos(u3));
    }
    
#endif
    
};


/// input operator
template <typename T>
std::istream& operator >> (std::istream& is, Quaternion<T> & arg)
{
    is >> arg[0] >> arg[1] >> arg[2] >> arg[3];
    return is;
}

/// output operator
template <typename T>
std::ostream& operator << (std::ostream& os, const Quaternion<T> & arg)
{
    arg.print(os);
    return os;
}

#endif
