// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
// F. Nedelec, Strasbourg 08.06.2018

#ifndef MATRIX34
#define MATRIX34

#include "real.h"
#include "vector3.h"

#include <cstdio>
#include <iostream>

#ifdef __AVX__
#  define MATRIX34_USES_AVX REAL_IS_DOUBLE
#  include "simd.h"
#else
#  define MATRIX34_USES_AVX 0
#endif



/// 3x4 matrix class with 9 'real' elements
/*
 This matrix uses line-major storage, and is optimized for the AVX instruction
 set, with SIMD vectors that can hold 4 double-precision elements.
 It is defacto a 3x4 matrix, but can be used as a 3x3 matrix.
 The function clear_shadow() clears the last column.
 */
class alignas(32) Matrix34
{
private:
    
    /// values of the elements
    real val[12];
    
    /// access to modifiable element by index
    real& operator[](int i)       { return val[i]; }
    
    /// access element value by index
    real  operator[](int i) const { return val[i]; }

public:
    
    Matrix34() { clear_shadow(); }
    
    /// copy constructor
    Matrix34(Matrix34 const& M)
    {
        for ( int u = 0; u < 12; ++u )
            val[u] = M.val[u];
    }

    /// construct Matrix from coordinates (given by columns)
    Matrix34(real a, real b, real c,
             real d, real e, real f,
             real g, real h, real i)
    {
        val[0] = a;
        val[1] = d;
        val[2] = g;
        val[4] = b;
        val[5] = e;
        val[6] = h;
        val[8] = c;
        val[9] = f;
        val[10] = i;
        clear_shadow();
    }

    /// construct Matrix with `d` on the diagonal and other values equal to `a`
    Matrix34(real z, real d)
    {
        val[0] = d;
        val[1] = z;
        val[2] = z;
        val[4] = z;
        val[5] = d;
        val[6] = z;
        val[8] = z;
        val[9] = z;
        val[10] = d;
        clear_shadow();
    }

    ~Matrix34() {}
    
    /// human-readible identifier
    static std::string what() { return "3*4"; }

    /// set all elements to zero
    void reset()
    {
        for ( int u = 0; u < 4*3; ++u )
            val[u] = 0.0;
    }
    
    bool operator != (real zero) const
    {
        for ( int u = 0; u < 4*3; ++u )
            if ( val[u] != zero )
                return true;
        return false;
    }
    
    /// conversion to pointer of real
    operator real const*() const { return val; }

    /// conversion to array of 'real'
    real* data()             { return val; }
    real* addr(int i, int j) { return val + ( 4*i + j ); }

    /// access functions to element by line and column indices
    real& operator()(int i, int j)       { return val[4*i+j]; }
    real  operator()(int i, int j) const { return val[4*i+j]; }
    
    /// extract column vector at given index
    Vector3 column(const unsigned i) const
    {
        return Vector3(val[i], val[4+i], val[4*2+i]);
    }
    
    /// extract line vector at given index
    Vector3 line(const unsigned i) const
    {
        return Vector3(val+4*i);
    }
    
    /// extract diagonal
    Vector3 diagonal() const
    {
        return Vector3(val[0], val[5], val[10]);
    }
    
    /// sum of diagonal terms
    real trace() const
    {
        return ( val[0] + val[5] + val[10] );
    }

    /// set matrix by giving lines
    void setLines(Vector3 const& A, Vector3 const& B, Vector3 const& C)
    {
        val[0 ] = A.XX;
        val[1 ] = A.YY;
        val[2 ] = A.ZZ;
        val[4 ] = B.XX;
        val[5 ] = B.YY;
        val[6 ] = B.ZZ;
        val[8 ] = C.XX;
        val[9 ] = C.YY;
        val[10] = C.ZZ;
    }
    
    /// set matrix by giving columns
    void setColumns(Vector3 const& A, Vector3 const& B, Vector3 const& C)
    {
        val[0 ] = A.XX;
        val[1 ] = B.XX;
        val[2 ] = C.XX;
        val[4 ] = A.YY;
        val[5 ] = B.YY;
        val[6 ] = C.YY;
        val[8 ] = A.ZZ;
        val[9 ] = B.ZZ;
        val[10] = C.ZZ;
    }

    /// print matrix in human readible format
    void print(FILE * f) const
    {
        fprintf(f, " / %9.3f %+9.3f %+9.3f %+9.3f \\\n",  val[0], val[1], val[2], val[3]);
        fprintf(f, "(  %9.3f %+9.3f %+9.3f %+9.3f  )\n" , val[4], val[5], val[6], val[7]);
        fprintf(f, " \\ %9.3f %+9.3f %+9.3f %+9.3f /\n",  val[8], val[9], val[10], val[11]);
    }
    
    /// output matrix lines to std::ostream
    std::ostream& operator << (std::ostream& os)
    {
        std::streamsize w = os.width();
        os.width(1);
        os << "[";
        for ( int i = 0; i < 3; ++i )
        {
            for ( int j = 0; j < 3; ++j )
                os << " " << std::fixed << std::setw(w) << (*this)(i,j);
            if ( i < 2 )
                os << ";";
            else
                os << " ]";
        }
        os.width(w);
        return os;
    }

    /// conversion to string
    std::string to_string(int w, int p) const
    {
        std::ostringstream os;
        os.precision(p);
        os << *this;
        return os.str();
    }

    /// clear values in the 4-th column
    void clear_shadow()
    {
        val[4 ] = 0.0;
        val[8 ] = 0.0;
        val[11] = 0.0;
    }
    
    /// scale all elements
    void scale(const real alpha)
    {
        for ( int u = 0; u < 12 ++u )
            val[u] *= alpha;
    }

    /// scale matrix
    void operator *=(const real alpha)
    {
        scale(alpha);
    }
    
    /// return opposite matrix (i.e. -M)
    const Matrix34 operator -() const
    {
        Matrix34 M;
        for ( int u = 0; u < 12; ++u )
            M.val[u] = -val[u];
        return M;
    }
    
    /// scaled matrix
    const Matrix34 operator *(const real alpha) const
    {
        Matrix34 res;
        for ( int u = 0; u < 12; ++u )
            res.val[u] = val[u] * alpha;
        return res;
    }
    
    /// multiplication by scalar
    friend const Matrix34 operator *(const real alpha, Matrix34 const& mat)
    {
        return mat * alpha;
    }

    /// return sum of two matrices
    const Matrix34 operator +(Matrix34 const& M) const
    {
        Matrix34 res;
        for ( int u = 0; u < 12; ++u )
            res.val[u] = val[u] + M.val[u];
        return res;
    }

    /// return sum of two matrices
    const Matrix34 operator -(Matrix34 const& M) const
    {
        Matrix34 res;
        for ( int u = 0; u < 12; ++u )
            res.val[u] = val[u] - M.val[u];
        return res;
    }

    /// subtract given matrix
    void operator +=(Matrix34 const& M)
    {
#if MATRIX34_USES_AVX
        store4(val  , add4(load4(val  ), load4(M.val  )));
        store4(val+4, add4(load4(val+4), load4(M.val+4)));
        store4(val+8, add4(load4(val+8), load4(M.val+8)));
#else
        for ( int u = 0; u < 12; ++u )
            val[u] += M.val[u];
#endif
    }

    /// add given matrix
    void operator -=(Matrix34 const& M)
    {
#if MATRIX34_USES_AVX
        store4(val  , sub4(load4(val  ), load4(M.val  )));
        store4(val+4, sub4(load4(val+4), load4(M.val+4)));
        store4(val+8, sub4(load4(val+8), load4(M.val+8)));
#else
        for ( int u = 0; u < 12; ++u )
            val[u] -= M.val[u];
#endif
    }
    
    /// transpose matrix in place
    void transpose()
    {
        std::swap(val[1], val[4]);
        std::swap(val[2], val[8]);
        std::swap(val[6], val[9]);
    }
    
    /// return matrix where 3x3 part is transposed
    Matrix34 transposed() const
    {
        Matrix34 res;
        for ( int x = 0; x < 3; ++x )
        for ( int y = 0; y < 3; ++y )
            res[y+4*x] = val[x+4*y];
        return res;
    }
    
    /// maximum of all component's absolute values
    real norm_inf() const
    {
        real res = fabs(val[0]);
        for ( unsigned i = 1; i < 12; ++i )
            res = std::max(res, fabs(val[i]));
        return res;
    }

    /// determinant of matrix
    real determinant() const
    {
        return ( val[0]*val[5]*val[10] + val[2]*val[4]*val[9]
                +val[1]*val[6]*val[8 ] - val[2]*val[5]*val[8]
                -val[1]*val[4]*val[10] - val[0]*val[6]*val[9] );
    }
    
    /// inverse in place
    void inverse()
    {
        real det = 1.0 / determinant();
        Vector3 X = column(0);
        Vector3 Y = column(1);
        Vector3 Z = column(2);
        setLines(cross(Y,Z)*det, cross(Z,X)*det, cross(X,Y)*det);
    }

    /// return inverse matrix
    Matrix34 inverted() const
    {
        Matrix34 res;
        real det = 1.0 / determinant();
        Vector3 X = column(0);
        Vector3 Y = column(1);
        Vector3 Z = column(2);
        res.setLines(cross(Y,Z)*det, cross(Z,X)*det, cross(X,Y)*det);
        //std::clog << " mat * inverse = " << mul(res).to_string(10, 3) << "\n";
        return res;
    }
    
    /// inversion of a symmetric matrix, addressing upper triangle
    /** This methods uses a L*D*L^t factorization with:
     L = ( 1 0 0; a 1 0; b c 1 )
     D = ( u 0 0; 0 v 0; 0 0 w )
     The result is a symetric matrix
     */
    int symmetricInverse()
    {
        /*
         // solving mat =  L * D * L^t:
         val[0,0] = u;
         val[1,0] = a * u;
         val[2,0] = b * u;
         val[0,1] = a * u;
         val[1,1] = a * a * u + v;
         val[2,1] = a * b * u + c * v;
         val[0,2] = b * u;
         val[1,2] = a * b * u + c * v;
         val[2,2] = b * b * u + c * c * v + w;
         */
        real u = val[0];
        real iu = 1.0 / u;
        real a = val[1] * iu;
        real b = val[2] * iu;
        real v = val[1+4] - a * val[1];
        real iv = 1.0 / v;
        real x = val[2+4] - a * val[2];
        real c = x * iv;
        real iw = 1.0 / ( val[2+4*2] - b * val[2] - c * x );
        // inverse triangular matrix U = inverse(L^t):
        b = -b + a * c;
        a = -a;
        c = -c;
        real aiv = a * iv;
        real biw = b * iw;
        real ciw = c * iw;
        // calculate U * inverse(D) * U^t:
        val[0+4*0] = iu + a * aiv + b * biw;
        val[1+4*0] = aiv + c * biw;
        val[2+4*0] = biw;
        val[0+4*1] = val[1+4*0];
        val[1+4*1] = iv + c * ciw;
        val[2+4*1] = ciw;
        val[0+4*2] = biw;
        val[1+4*2] = ciw;
        val[2+4*2] = iw;
        return 0;
    }

    /// copy values from upper triangle to lower triangle
    void copy_half()
    {
        val[4] = val[1];
        val[8] = val[2];
        val[9] = val[6];
    }

    /// true if matrix is symmetric
    bool is_symmetric() const
    {
        return ( val[4] == val[1]
                && val[8] == val[2]
                && val[9] == val[6] );
    }

#if MATRIX34_USES_AVX
    /// multiplication by a vector: this * V
    const vec4 vecmul3(double const* V) const
    {
        vec4 vec = loadu4(V); // { x, y, z, garbage }
        vec4 s0 = mul4(load4(val  ), vec);
        vec4 s1 = mul4(load4(val+4), vec);
        vec4 s2 = mul4(load4(val+8), vec);
        vec4 s3 = setzero4();
        vec4 xy = add4(unpacklo4(s0, s1), unpackhi4(s0, s1));
        vec4 zt = add4(unpacklo4(s2, s3), unpackhi4(s2, s3));
        return add4(permute2f128(xy, zt, 0x20), permute2f128(xy, zt, 0x31));
    }
    
    /// multiplication by a vector: this * V
    const vec4 vecmul4(const vec4 vec) const
    {
        vec4 s0 = mul4(load4(val  ), vec);
        vec4 s1 = mul4(load4(val+4), vec);
        vec4 s2 = mul4(load4(val+8), vec);
        vec4 s3 = setzero4();
        vec4 xy = add4(unpacklo4(s0, s1), unpackhi4(s0, s1));
        vec4 zt = add4(unpacklo4(s2, s3), unpackhi4(s2, s3));
        return add4(permute2f128(xy, zt, 0x20), permute2f128(xy, zt, 0x31));
    }

    /// multiplication by a vector: transpose(M) * V
    const vec4 trans_vecmul3(double const* V) const
    {
        vec4 xyxy = broadcast2(V);
        vec4 xxxx = duplo4(xyxy); //broadcast1(V);
        vec4 yyyy = duphi4(xyxy); //broadcast1(V+1);
        vec4 zzzz = broadcast1(V+2);
        xxxx = mul4(load4(val), xxxx);
        yyyy = mul4(load4(val+4), yyyy);
        return fmadd4(load4(val+8), zzzz, add4(xxxx,yyyy));
    }
#endif
    
    /// multiplication by a vector: this * V
    const Vector3 vecmul0(Vector3 const& V) const
    {
        return Vector3(val[0] * V.XX + val[1] * V.YY + val[2] * V.ZZ,
                       val[4] * V.XX + val[5] * V.YY + val[6] * V.ZZ,
                       val[8] * V.XX + val[9] * V.YY + val[10] * V.ZZ);
    }
    
    /// multiplication by a vector: this * V
    const Vector3 vecmul0(real const* ptr) const
    {
        return Vector3(val[0] * ptr[0] + val[1] * ptr[1] + val[2] * ptr[2],
                       val[4] * ptr[0] + val[5] * ptr[1] + val[6] * ptr[2],
                       val[8] * ptr[0] + val[9] * ptr[1] + val[10] * ptr[2]);
    }

    /// multiplication by a vector: this * V
    inline const Vector3 vecmul(Vector3 const& vec) const
    {
#if MATRIX34_USES_AVX
        return vecmul3(vec);
#else
        return vecmul0(vec);
#endif
    }
    
    /// multiplication by a vector: this * { ptr[0], ptr[1] }
    const Vector3 vecmul(real const* ptr) const
    {
#if MATRIX34_USES_AVX
        return vecmul3(ptr);
#else
        return vecmul0(ptr);
#endif
    }

    /// multiplication with a vector: M * V
    friend Vector3 operator * (Matrix34 const& mat, Vector3 const& vec)
    {
        return mat.vecmul(vec);
    }

    /// multiplication by a vector: transpose(M) * V
    const Vector3 trans_vecmul0(Vector3 const& V) const
    {
        return Vector3(val[0] * V.XX + val[4] * V.YY + val[ 8] * V.ZZ,
                       val[1] * V.XX + val[5] * V.YY + val[ 9] * V.ZZ,
                       val[2] * V.XX + val[6] * V.YY + val[10] * V.ZZ);
    }

    /// multiplication by a vector: transpose(M) * V
    const Vector3 trans_vecmul0(real const* V) const
    {
        return Vector3(val[0] * V[0] + val[4] * V[1] + val[ 8] * V[2],
                       val[1] * V[0] + val[5] * V[1] + val[ 9] * V[2],
                       val[2] * V[0] + val[6] * V[1] + val[10] * V[2]);
    }

    /// multiplication by a vector: transpose(M) * V
    inline const Vector3 trans_vecmul(real const* V) const
    {
#if MATRIX34_USES_AVX
        return trans_vecmul3(V);
#else
        return trans_vecmul0(V);
#endif
    }

    /// multiplication by another matrix: @returns this * M
    const Matrix34 mul(Matrix34 const& M) const
    {
        Matrix34 res;
        res[0] = val[4*0] * M[0] + val[1+4*0] * M[4] + val[2+4*0] * M[8];
        res[1] = val[4*1] * M[0] + val[1+4*1] * M[4] + val[2+4*1] * M[8];
        res[2] = val[4*2] * M[0] + val[1+4*2] * M[4] + val[2+4*2] * M[8];
        
        res[0+4] = val[4*0] * M[1] + val[1+4*0] * M[5] + val[2+4*0] * M[9];
        res[1+4] = val[4*1] * M[1] + val[1+4*1] * M[5] + val[2+4*1] * M[9];
        res[2+4] = val[4*2] * M[1] + val[1+4*2] * M[5] + val[2+4*2] * M[9];
        
        res[0+4*2] = val[4*0] * M[2] + val[1+4*0] * M[6] + val[2+4*0] * M[10];
        res[1+4*2] = val[4*1] * M[2] + val[1+4*1] * M[6] + val[2+4*1] * M[10];
        res[2+4*2] = val[4*2] * M[2] + val[1+4*2] * M[6] + val[2+4*2] * M[10];
        return res;
    }
    
    /// multiplication with a matrix
    friend Matrix34 operator * (Matrix34 const& mat, Matrix34 const& mut)
    {
        return mat.mul(mut);
    }

    /// multiplication by another matrix: @returns transpose(this) * M
    const Matrix34 trans_mul(Matrix34 const& M) const
    {
        Matrix34 res;
        res[0] = val[0] * M[0] + val[0+4] * M[4] + val[0+4*2] * M[8];
        res[1] = val[1] * M[0] + val[1+4] * M[4] + val[1+4*2] * M[8];
        res[2] = val[2] * M[0] + val[2+4] * M[4] + val[2+4*2] * M[8];
        
        res[0+4] = val[0] * M[1] + val[0+4] * M[5] + val[0+4*2] * M[9];
        res[1+4] = val[1] * M[1] + val[1+4] * M[5] + val[1+4*2] * M[9];
        res[2+4] = val[2] * M[1] + val[2+4] * M[5] + val[2+4*2] * M[9];
        
        res[0+4*2] = val[0] * M[2] + val[0+4] * M[6] + val[0+4*2] * M[10];
        res[1+4*2] = val[1] * M[2] + val[1+4] * M[6] + val[1+4*2] * M[10];
        res[2+4*2] = val[2] * M[2] + val[2+4] * M[6] + val[2+4*2] * M[10];
        return res;
    }
    
    
    /// add full matrix: this <- this + M
    void add_full(Matrix34 const& M)
    {
        real const* src = M.val;
#if MATRIX34_USES_AVX
        store4(val  , add4(load4(val  ), load4(src  )));
        store4(val+4, add4(load4(val+4), load4(src+4)));
        store4(val+8, add4(load4(val+8), load4(src+8)));
#else
        for ( int u = 0; u < 12; ++u )
            val[u] += src[u];
#endif
    }
    
    /// add full matrix: this <- this + alpha * M
    void add_full(const real alpha, Matrix34 const& M)
    {
        real const* src = M.val;
#if MATRIX34_USES_AVX
        vec4 a = set4(alpha);
        store4(val  , fmadd4(a, load4(src  ), load4(val  )));
        store4(val+4, fmadd4(a, load4(src+4), load4(val+4)));
        store4(val+8, fmadd4(a, load4(src+8), load4(val+8)));
#else
        for ( int u = 0; u < 12; ++u )
            val[u] += alpha * src[u];
#endif
    }
    
    /// sub full matrix: this <- this - M
    void sub_full(Matrix34 const& M)
    {
        real const* src = M.val;
#if MATRIX34_USES_AVX
        store4(val  , sub4(load4(val  ), load4(src  )));
        store4(val+4, sub4(load4(val+4), load4(src+4)));
        store4(val+8, sub4(load4(val+8), load4(src+8)));
#else
        for ( int u = 0; u < 12; ++u )
            val[u] -= src[u];
#endif
    }

    /// subtract transposed matrix: this <- this - transposed(M)
    void sub_trans(Matrix34 const& M)
    {
        real const* src = M.val;
        for ( int x = 0; x < 3; ++x )
        for ( int y = 0; y < 3; ++y )
            val[y+4*x] -= src[x+4*y];
    }
    
    /// add transposed matrix: this <- this + alpha * transposed(M)
    void add_trans(Matrix34 const& M)
    {
        real const* src = M.val;
        for ( int x = 0; x < 3; ++x )
        for ( int y = 0; y < 3; ++y )
            val[y+4*x] += src[x+4*y];
    }
    
    /// add transposed matrix: this <- this + alpha * transposed(M)
    void add_trans(const real alpha, Matrix34 const& M)
    {
        real const* src = M.val;
        for ( int x = 0; x < 3; ++x )
        for ( int y = 0; y < 3; ++y )
            val[y+4*x] += alpha * src[x+4*y];
    }

    /// add lower triangle of matrix including diagonal: this <- this + M
    void add_half(Matrix34 const& M)
    {
        real const* src = M.val;
#if MATRIX34_USES_AVX
        store4(val  , add4(load4(val  ), load4(src  )));
        store4(val+4, add4(load4(val+4), load4(src+4)));
        store4(val+8, add4(load4(val+8), load4(src+8)));
#else
        for ( int u = 0; u < 12; ++u )
            val[u] += src[u];
#else
        for ( int i = 0; i < 3; ++i )
        for ( int j = i; j < 3; ++j )
            val[4*i+j] += src[4*i+j];
#endif
    }
    
    /// add lower triangle of matrix including diagonal: this <- this + alpha * M
    void add_half(const real alpha, Matrix34 const& M)
    {
        real const* src = M.val;
        //std::clog << "matrix alignment " << ((uintptr_t)src & 63) << "\n";
#if MATRIX34_USES_AVX
        vec4 a = set4(alpha);
        store4(val  , fmadd4(a, load4(src  ), load4(val  )));
        store4(val+4, fmadd4(a, load4(src+4), load4(val+4)));
        store4(val+8, fmadd4(a, load4(src+8), load4(val+8)));
#elif ( 1 )
        for ( int u = 0; u < 12; ++u )
            val[u] += alpha * src[u];
#else
        for ( int i = 0; i < 3; ++i )
        for ( int j = i; j < 3; ++j )
            val[4*i+j] += alpha * src[4*i+j];
#endif
    }
    
    /// subtract lower triangle of matrix including diagonal: this <- this - M
    void sub_half(Matrix34 const& M)
    {
        real const* src = M.val;
#if MATRIX34_USES_AVX
        store4(val  , sub4(load4(val  ), load4(src  )));
        store4(val+4, sub4(load4(val+4), load4(src+4)));
        store4(val+8, sub4(load4(val+8), load4(src+8)));
#elif ( 1 )
        for ( int u = 0; u < 12; ++u )
            val[u] -= src[u];
#else
        for ( int i = 0; i < 3; ++i )
        for ( int j = i; j < 3; ++j )
            val[4*i+j] -= src[4*i+j];
#endif
    }
    
    /// add alpha to diagonal
    void add_diag(real alpha)
    {
        val[0]  += alpha;
        val[5]  += alpha;
        val[10] += alpha;
    }
    
    /// add -alpha to diagonal
    void sub_diag(real alpha)
    {
        val[0]  -= alpha;
        val[5]  -= alpha;
        val[10] -= alpha;
    }

    /// add all elements of block 'S' to array 'M'
    void addto(real * M, unsigned ldd) const
    {
        M[0      ] += val[0];
        M[1      ] += val[4];
        M[2      ] += val[8];
        M[  ldd  ] += val[1];
        M[1+ldd  ] += val[5];
        M[2+ldd  ] += val[9];
        M[  ldd*2] += val[2];
        M[1+ldd*2] += val[6];
        M[2+ldd*2] += val[10];
    }
    
    /// add upper elements of this block to upper triangle of 'M'
    void addto_upper(real * M, unsigned ldd) const
    {
        M[0      ] += val[0];
        M[  ldd  ] += val[1];
        M[  ldd*2] += val[2];
        M[1+ldd  ] += val[4];
        M[1+ldd*2] += val[5];
        M[2+ldd*2] += val[10];
    }
    
    /// add upper elements of this block to both upper and lower triangles of 'M'
    void addto_symm(real * M, unsigned ldd) const
    {
        M[0      ] += val[0];
        M[1      ] += val[1];
        M[2      ] += val[2];
        M[  ldd  ] += val[1];
        M[1+ldd  ] += val[5];
        M[2+ldd  ] += val[6];
        M[  ldd*2] += val[2];
        M[1+ldd*2] += val[6];
        M[2+ldd*2] += val[10];
    }
    
    /// add all elements of this block to 'M', with transposition
    void addto_trans(real * M, unsigned ldd) const
    {
        M[0      ] += val[0];
        M[1      ] += val[1];
        M[2      ] += val[2];
        M[  ldd  ] += val[4];
        M[1+ldd  ] += val[5];
        M[2+ldd  ] += val[6];
        M[  ldd*2] += val[8];
        M[1+ldd*2] += val[9];
        M[2+ldd*2] += val[10];
    }
    
    
    /// return symmetric Matrix from coordinates (column-major, lower triangle)
    static Matrix34 symmetric(real a, real b, real c,
                              real d, real e, real f)
    {
        return Matrix34(a, b, c, b, d, e, c, e, f);
    }

    /// return diagonal Matrix from diagonal terms
    static Matrix34 diagonal(real a, real b, real c)
    {
        return Matrix34(a, 0, 0, 0, b, 0, 0, 0, c);
    }
    
    /// identity matrix
    static Matrix34 identity()
    {
        return Matrix34(0, 1);
    }

    /// return a symmetric matrix: [ dir (x) transpose(dir) ]
    static Matrix34 outerProduct(Vector3 const& dir)
    {
        return symmetric(dir.XX*dir.XX, dir.YY*dir.XX, dir.ZZ*dir.XX,
                         dir.YY*dir.YY, dir.YY*dir.ZZ,
                         dir.ZZ*dir.ZZ );
    }

    /// return a symmetric matrix: alpha * [ dir (x) transpose(dir) ]
    static Matrix34 outerProduct(Vector3 const& dir, real alpha)
    {
#if MATRIX34_USES_AVX
        Matrix34 res;
        vec4 p = permute2f128(dir, dir, 0x01);
        vec4 l = blend4(dir, p, 0b1100);
        vec4 u = blend4(dir, p, 0b0011);
        vec4 d = mul4(dir, set4(alpha));
        store4(res.val  , mul4(d, duplo4(l)));
        store4(res.val+4, mul4(d, duphi4(l)));
        store4(res.val+8, mul4(d, duplo4(u)));
        return res;
#else
        real XX = dir.XX * alpha;
        real YY = dir.YY * alpha;
        real ZZ = dir.ZZ * alpha;
        return symmetric(dir.XX*XX, dir.YY*XX, dir.ZZ*XX,
                         dir.YY*YY, dir.YY*ZZ,
                         dir.ZZ*ZZ );
#endif
    }
    
    /// return outer product: [ vec (x) transpose(dir) ]
    static Matrix34 outerProduct(Vector3 const& vec, Vector3 const& dir)
    {
#if MATRIX34_USES_AVX
        Matrix34 res;
        vec4 p = permute2f128(vec, vec, 0x01);
        vec4 l = blend4(vec, p, 0b1100);
        vec4 u = blend4(vec, p, 0b0011);
        store4(res.val  , mul4(dir, duplo4(l)));
        store4(res.val+4, mul4(dir, duphi4(l)));
        store4(res.val+8, mul4(dir, duplo4(u)));
        return res;
#else
        return Matrix34(dir.XX*vec.XX, dir.YY*vec.XX, dir.ZZ*vec.XX,
                        dir.XX*vec.YY, dir.YY*vec.YY, dir.ZZ*vec.YY,
                        dir.XX*vec.ZZ, dir.YY*vec.ZZ, dir.ZZ*vec.ZZ );
#endif
    }
    
    /// return outer product: [ vec (x) transpose(dir) ]
    static Matrix34 outerProduct(real const* vec, real const* dir)
    {
#if MATRIX34_USES_AVX
        Matrix34 res;
        vec4 d = load3(dir);
        store4(res.val  , mul4(d, broadcast1(vec  )));
        store4(res.val+4, mul4(d, broadcast1(vec+1)));
        store4(res.val+8, mul4(d, broadcast1(vec+2)));
        return res;
#else
        return Matrix34(dir[0]*vec[0], dir[1]*vec[0], dir[2]*vec[0],
                        dir[0]*vec[1], dir[1]*vec[1], dir[2]*vec[1],
                        dir[0]*vec[2], dir[1]*vec[2], dir[2]*vec[2] );
#endif
    }
    
    /// add outer product: [ vec (x) transpose(dir) ]
    void addOuterProduct(real const* vec, real const* dir)
    {
#if MATRIX34_USES_AVX
        vec4 d = load3(dir);
        store4(val  , fmadd4(d, broadcast1(vec  ), load4(val  )));
        store4(val+4, fmadd4(d, broadcast1(vec+1), load4(val+4)));
        store4(val+8, fmadd4(d, broadcast1(vec+2), load4(val+8)));
#else
        val[0    ] += dir[0]*vec[0];
        val[1    ] += dir[1]*vec[0];
        val[2    ] += dir[2]*vec[0];
        val[0+4  ] += dir[0]*vec[1];
        val[1+4  ] += dir[1]*vec[1];
        val[2+4  ] += dir[2]*vec[1],
        val[0+4*2] += dir[0]*vec[2];
        val[1+4*2] += dir[1]*vec[2];
        val[2+4*2] += dir[2]*vec[2];
#endif
    }

    /// return [ dir (x) transpose(vec) + vec (x) transpose(dir) ]
    static Matrix34 symmetricOuterProduct(Vector3 const& dir, Vector3 const& vec)
    {
        real xx = dir.XX * vec.XX;
        real yy = dir.YY * vec.YY;
        real zz = dir.ZZ * vec.ZZ;
        return symmetric(xx+xx, dir.YY*vec.XX + dir.XX*vec.YY, dir.ZZ*vec.XX + dir.XX*vec.ZZ,
                         yy+yy, dir.ZZ*vec.YY + dir.YY*vec.ZZ,
                         zz+zz);
    }
 
    /// return symmetric matrix block :  dia * I + [ dir (x) dir ] * len
    static Matrix34 offsetOuterProduct(const real dia, Vector3 const& dir, const real len)
    {
        real xl = dir.XX * len;
        real yl = dir.YY * len;
        real zl = dir.ZZ * len;
        return symmetric(xl * dir.XX + dia, yl * dir.XX, zl * dir.XX,
                         yl * dir.YY + dia, zl * dir.YY,
                         zl * dir.ZZ + dia);
    }
    
    // build the rotation matrix `M = 2 V * V' - 1` of angle 180
    static Matrix34 householder(const Vector3& axis)
    {        
        real X = axis.XX, Y = axis.YY, Z = axis.ZZ;
        real X2 = X + X, Y2 = Y + Y, Z2 = Z + Z;
        
        return symmetric(X * X2 - 1.0, X * Y2, X * Z2,
                         Y * Y2 - 1.0, Y * Z2,
                         Z * Z2 - 1.0);
    }

    /// rotation around `axis` (of norm 1) with angle set by cosinus and sinus values
    /** The values of 'c' and 's' can be scaled to obtain a matrix where the
     rotation components is also scaling. Vectors along the axis remain unchanged */
    static Matrix34 rotationAroundAxis(const Vector3& axis, const real c, const real s)
    {
        /*
         This is using Rodrigues's formula:
             I + sinus * K + ( 1 - cosinus ) * K^2
             K = -1 (x) axis
        Attention: this is correct only if norm(axis)==1
         */
        const real  X = axis.XX  ,  Y = axis.YY  ,  Z = axis.ZZ;
        const real dX = X - c * X, dY = Y - c * Y, dZ = Z - c * Z;
        const real sX = s * X    , sY = s * Y    , sZ = s * Z;

        return Matrix34(dX * X + c , dY * X + sZ, dZ * X - sY,
                        dX * Y - sZ, dY * Y + c , dZ * Y + sX,
                        dX * Z + sY, dY * Z - sX, dZ * Z + c );
    }

    /// rotation axis
    Vector3 rotationAxis() const;

    /// rotation angle
    real rotationAngle() const;

    /// calculate rotation angle and Euler angle of axis
    void getEulerAngles(real& angle, real&, real&) const;

    ///
    static Matrix34 rotationAroundAxisEuler(const real a[3]);

    /// return rotation of angle a, around axis of azimuth b and elevation c
    static Matrix34 rotationFromAngles(const real a[3]);

    
    /// return a rotation that transforms (1,0,0) into (-1,0,0)
    static Matrix34 rotation180();

    /// a rotation around the X axis of specified angle
    static Matrix34 rotationAroundX(real angle);
    
    /// a rotation around the Y axis of specified angle
    static Matrix34 rotationAroundY(real angle);
    
    /// a rotation around the Z axis of specified angle
    static Matrix34 rotationAroundZ(real angle);
    
    /// a rotation around one the axis X if `x==0`, Y if `x==1` or Z if `x==2`
    static Matrix34 rotationAroundPrincipalAxis(unsigned x, real angle);

    /// return a rotation that transforms (1,0,0) into `vec` ( norm(vec) should be > 0 )
    static Matrix34 rotationToVector(const Vector3&);
    
    /// return a random rotation that transforms (1,0,0) into `vec` ( norm(vec) should be > 0 )
    /**
     In 3D, this rotation is chosen uniformly among all the rotation transforming (1,0,0) into `vec`.
     The function will fail if ( vec == 0 ).
     */
    static Matrix34 randomRotationToVector(const Vector3&);
    
    /// a random rotation chosen uniformly
    static Matrix34 randomRotation();
    
    /// a rotation of angle 'angle' around an axis chosen randomly
    static Matrix34 randomRotation(real angle);
};

#endif

