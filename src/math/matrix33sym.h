// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University.
// FJN, Cambridge 25.02.2022

#ifndef MATRIX33SYM
#define MATRIX33SYM

#include "real.h"
#include "vector3.h"

#include <cstdio>
#include <iostream>


#if defined(__SSE3__) && REAL_IS_DOUBLE
#  include "simd.h"
#  define MATRIX33SYM_USES_SSE 1
#else
#  define MATRIX33SYM_USES_SSE 0
#endif


/// 3x3 symmetric matrix class with 6 'real' elements stored in column order
/** This class cannot represent rotations and other non-symmatric operations */
class alignas(2*sizeof(real)) Matrix33sym final
{
public:
    
    /// unsigned integer type used for indices
    typedef size_t index;

private:
    
    /// values of the elements
    real val[6];

    /// access to modifiable element by index
    real& operator[](index i)       { return val[i]; }
    
    /// access element value by index
    real  operator[](index i) const { return val[i]; }
    
public:

    /// default constructor
    Matrix33sym() { }
    
    /// copy constructor
    Matrix33sym(Matrix33sym const& M)
    {
        for ( index u = 0; u < 6; ++u )
            val[u] = M.val[u];
    }

    /// construct Matrix from coordinates (given by columns)
    Matrix33sym(real a, real b, real c,
                real d, real e, real f)
    {
        val[0] = a;
        val[1] = b;
        val[2] = c;
        val[3] = d;
        val[4] = e;
        val[5] = f;
    }

    /// construct Matrix with diagonal terms set to `d` and other terms set to `z`
    Matrix33sym(real z, real d)
    {
        val[0] = d;
        val[1] = z;
        val[2] = z;
        val[3] = d;
        val[4] = z;
        val[5] = d;
    }

    ~Matrix33sym() {}
    
#pragma mark -
    
    /// dimensionality
    static constexpr size_t dimension() { return 3; }
    
    /// human-readable identifier
#if MATRIX33SYM_USES_SSE
    static std::string what() { return "+6"; }
#else
    static std::string what() { return "6"; }
#endif
    
    /// set all elements to zero
    void reset()
    {
        for ( index u = 0; u < 6; ++u )
            val[u] = 0.0;
    }
    
    /// set diagonal to 'dia' and other elements to 'off'
    void reset(real off, real dia)
    {
        val[0] = dia;
        val[1] = off;
        val[2] = off;
        val[3] = dia;
        val[4] = off;
        val[5] = dia;
    }
    
    bool operator != (real zero) const
    {
        for ( index u = 0; u < 6; ++u )
            if ( val[u] != zero )
                return true;
        return false;
    }
    
    /// conversion to pointer of real
    //operator real const*() const { return val; }

    /// return modifiable pointer of 'real'
    real* data() { return val; }
    
#if defined(__SSE3__) && REAL_IS_DOUBLE
    vec2 data0() const { return load2(val); }
    vec2 data1() const { return load2(val+2); }
    vec2 data2() const { return load2(val+4); }
#endif

    /// unmodifiable pointer of real
    real const* data() const { return val; }

    /// address of element at line i, column j
    real* addr(const index i, const index j) { return val + ( i+((5-j)*j)/2 ); }
    /// value of element at line i, column j
    real value(const index i, const index j) const { return val[i+((5-j)*j)/2]; }

    /// access functions to element by line and column indices
    real& operator()(const index i, const index j)       { return val[i+((5-j)*j)/2]; }
    real  operator()(const index i, const index j) const { return val[i+((5-j)*j)/2]; }

    /// set elements from given array
    void load(const real ptr[])
    {
        for ( index i = 0; i < 6; ++i )
            val[i] = ptr[i];
    }

    /// copy elements to given array
    void store(real ptr[]) const
    {
        for ( index i = 0; i < 6; ++i )
            ptr[i] = val[i];
    }

    /// extract column vector at given index
    Vector3 column(const index i) const
    {
        switch( i )
        {
            case 0: return Vector3(val[0], val[1], val[2]);
            case 1: return Vector3(val[1], val[3], val[4]);
            default: return Vector3(val[2], val[4], val[5]);
        }
    }
    
    /// extract line vector at given index
    Vector3 line(const index i) const
    {
        switch( i )
        {
            case 0: return Vector3(val[0], val[1], val[2]);
            case 1: return Vector3(val[1], val[3], val[4]);
            default: return Vector3(val[2], val[4], val[5]);
        }
    }
    
    /// extract diagonal
    Vector3 diagonal() const
    {
        return Vector3(val[0], val[3], val[5]);
    }
    
    /// sum of diagonal terms
    real trace() const
    {
        return ( val[0] + val[3] + val[5] );
    }
    
#pragma mark -

    /// print matrix in human readable format
    void print(FILE * f) const
    {
        fprintf(f, " / %9.3f %+9.3f %+9.3f \\\n", val[0], val[1], val[2]);
        fprintf(f, " | %9.3f %+9.3f %+9.3f |\n" , val[1], val[3], val[4]);
        fprintf(f, " \\ %9.3f %+9.3f %+9.3f /\n", val[2], val[4], val[5]);
    }
    
    /// print [ line1; line2; line3 ]
    void print(std::ostream& os) const
    {
        const int w = (int)os.width();
        os << std::setw(1) << "[";
        for ( index i = 0; i < 3; ++i )
        {
            for ( index j = 0; j < 3; ++j )
                os << " " << std::fixed << std::setw(w) << value(i,j);
            if ( i < 2 )
                os << ";";
            else
                os << " ]";
        }
    }

    /// print [ line1; line2; line3 ]
    void print_smart(std::ostream& os) const
    {
        const int w = (int)os.width();
        os << std::setw(1) << "[";
        os << " " << std::fixed << std::setw(w) << val[0];
        os << " " << std::fixed << std::setw(w) << val[1];
        os << " " << std::fixed << std::setw(w) << val[2];
        os << ";";
        os << " " << std::fixed << std::setw(w) << "sym";
        os << " " << std::fixed << std::setw(w) << val[3];
        os << " " << std::fixed << std::setw(w) << val[4];
        os << ";";
        os << " " << std::fixed << std::setw(w) << "sym";
        os << " " << std::fixed << std::setw(w) << "sym";
        os << " " << std::fixed << std::setw(w) << val[5];
        os << " ]";
    }

    /// conversion to string
    std::string to_string(std::streamsize w, std::streamsize p) const
    {
        std::ostringstream os;
        os.precision(p);
        os.width(w);
        print(os);
        return os.str();
    }
    
    /// scale all elements
    void scale(const real alpha)
    {
        for ( index u = 0; u < 6; ++u )
            val[u] *= alpha;
    }

    /// scale matrix
    void operator *=(const real alpha)
    {
        scale(alpha);
    }
    
    /// return opposite matrix (i.e. -M)
    const Matrix33sym operator -() const
    {
        Matrix33sym M;
        for ( index u = 0; u < 6; ++u )
            M.val[u] = -val[u];
        return M;
    }
    
    /// scaled matrix
    const Matrix33sym operator *(const real alpha) const
    {
        Matrix33sym res;
        for ( index u = 0; u < 6; ++u )
            res.val[u] = val[u] * alpha;
        return res;
    }
    
    /// multiplication by scalar
    friend const Matrix33sym operator *(const real alpha, Matrix33sym const& mat)
    {
        return mat * alpha;
    }

    /// return sum of two matrices
    const Matrix33sym operator +(Matrix33sym const& M) const
    {
        Matrix33sym res;
        for ( index u = 0; u < 6; ++u )
            res.val[u] = val[u] + M.val[u];
        return res;
    }

    /// return difference of two matrices
    const Matrix33sym operator -(Matrix33sym const& M) const
    {
        Matrix33sym res;
        for ( index u = 0; u < 6; ++u )
            res.val[u] = val[u] - M.val[u];
        return res;
    }

    /// subtract given matrix
    void operator += (Matrix33sym const& M)
    {
#if MATRIX33SYM_USES_SSE
        store2(val  , add2(load2(val  ), load2(M.val  )));
        store2(val+2, add2(load2(val+2), load2(M.val+2)));
        store2(val+4, add2(load2(val+4), load2(M.val+4)));
#else
        for ( index u = 0; u < 6; ++u )
            val[u] += M.val[u];
#endif
    }

    /// add given matrix
    void operator -= (Matrix33sym const& M)
    {
#if MATRIX33SYM_USES_SSE
        store2(val  , sub2(load2(val  ), load2(M.val  )));
        store2(val+2, sub2(load2(val+2), load2(M.val+2)));
        store2(val+4, sub2(load2(val+4), load2(M.val+4)));
#else
        for ( index u = 0; u < 6; ++u )
            val[u] -= M.val[u];
#endif
    }
    
    /// transpose matrix in place
    void transpose()
    {
    }
    
    /// return transposed matrix
    Matrix33sym transposed() const
    {
        Matrix33sym res = *this;
        return res;
    }
    
    /// return scaled transposed matrix
    Matrix33sym transposed(real alpha) const
    {
        Matrix33sym res = *this * alpha;
        return res;
    }

    
    /// maximum of all component's absolute values
    real norm_inf() const
    {
        real res = abs_real(val[0]);
        for ( index i = 1; i < 6; ++i )
            res = std::max(res, abs_real(val[i]));
        return res;
    }
    
    /// inversion of a symmetric matrix, using values in lower triangle
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
        real v = val[3] - a * val[1];
        real iv = 1.0 / v;
        real x = val[4] - a * val[2];
        real c = x * iv;
        real iw = 1.0 / ( val[5] - b * val[2] - c * x );
        // inverse triangular matrix U = inverse(L^t):
        b = -b + a * c;
        a = -a;
        c = -c;
        real aiv = a * iv;
        real biw = b * iw;
        real ciw = c * iw;
        // calculate U * inverse(D) * U^t:
        val[0] = iu + a * aiv + b * biw;
        val[1] = aiv + c * biw;
        val[2] = biw;
        val[3] = iv + c * ciw;
        val[4] = ciw;
        val[5] = iw;
        return 0;
    }

    /// determinant of matrix
    real determinant() const
    {
        return ( value(0,0) * ( value(1,1)*value(2,2) - value(2,1)*value(1,2) )
                +value(0,1) * ( value(1,2)*value(2,0) - value(2,2)*value(1,0) )
                +value(0,2) * ( value(1,0)*value(2,1) - value(2,0)*value(1,1) ));
    }
    
    /// inverse in place
    void inverse()
    {
        symmetricInverse();
    }

    /// return inverse matrix
    Matrix33sym inverted() const
    {
        Matrix33sym res = *this;
        res.symmetricInverse();
        return res;
    }

    /// copy values from lower triangle to upper triangle
    void copy_lower()
    {
    }

    /// relative asymmetry of 3x3 submatrix (divided by the trace)
    real asymmetry() const
    {
        return 0;
    }
    
#pragma mark -

#if MATRIX33SYM_USES_SSE && defined(__AVX__)
    
    /// multiplication by a 3-components vector: this * V
    const vec4 vecmul3_sse(vec4 vec) const
    {
        ABORT_NOW("unfinished");
        /*
        vec4 xy = duplo2f128(vec); // xyxy
        vec4 zzzz = duplo4(duphi2f128(vec));
        vec = mul2(data0(), duplo4(xy)); // xxxx
        xy = mul2(data1(), duphi4(xy)); // yyyy
        zzzz = fmadd2(data2(), zzzz, add4(xy, vec));
        return clear4th(zzzz);
         */
    }

    /// multiplication by a vector: this * V
    const vec4 vecmul3_sse(double const* V) const
    {
        vec2 zz = loadu2(V); // xy
        vec2 xx = duplo2(zz);
        vec2 yy = duphi2(zz);
        vec2 ab = data0();
        vec2 cd = data1();
        vec2 ef = data2();
        zz = loaddup2(V+2);  // zz
        ab = fmadd2(unpackhi2(ab, cd), yy, mul2(ab, xx));
        xx = fmadd1(ef, yy, mul1(cd, xx));
        ab = fmadd2(unpacklo2(cd, ef), zz, ab);
        xx = fmadd1(unpackhi2(ef, ef), zz, xx);
        return cat22(xx, ab);
    }

#endif
    
    /// multiplication by a vector: this * V
    Vector3 vecmul_(Vector3 const& V) const
    {
        return Vector3(val[0] * V.XX + val[1] * V.YY + val[2] * V.ZZ,
                       val[1] * V.XX + val[3] * V.YY + val[4] * V.ZZ,
                       val[2] * V.XX + val[4] * V.YY + val[5] * V.ZZ);
    }
    
    /// multiplication by a vector: this * V
    Vector3 vecmul_(real const* R) const
    {
        return Vector3(val[0] * R[0] + val[1] * R[1] + val[2] * R[2],
                       val[1] * R[0] + val[3] * R[1] + val[4] * R[2],
                       val[2] * R[0] + val[4] * R[1] + val[5] * R[2]);
    }
    
    /// multiplication by a vector: transpose(M) * V
    Vector3 trans_vecmul_(Vector3 const& V) const
    {
        return Vector3(val[0] * V.XX + val[1] * V.YY + val[2] * V.ZZ,
                       val[1] * V.XX + val[3] * V.YY + val[4] * V.ZZ,
                       val[2] * V.XX + val[4] * V.YY + val[5] * V.ZZ);
    }

    /// multiplication by a vector: transpose(M) * V
    Vector3 trans_vecmul_(real const* R) const
    {
        return Vector3(val[0] * R[0] + val[1] * R[1] + val[2] * R[2],
                       val[1] * R[0] + val[3] * R[1] + val[4] * R[2],
                       val[2] * R[0] + val[4] * R[1] + val[5] * R[2]);
    }

    /// multiplication by a vector: this * V
    Vector3 vecmul(Vector3 const& vec) const
    {
#if MATRIX33SYM_USES_SSE && VECTOR3_USES_AVX
        return vecmul3_sse(vec.vec);
#else
        return vecmul_(vec);
#endif
    }
    
    /// multiplication by a vector: this * { ptr[0], ptr[1] }
    Vector3 vecmul(real const* ptr) const
    {
#if MATRIX33SYM_USES_SSE && VECTOR3_USES_AVX
        return vecmul3_sse(ptr);
#else
        return vecmul_(ptr);
#endif
    }

    /// multiplication with a vector: M * V
    friend Vector3 operator * (Matrix33sym const& mat, Vector3 const& vec)
    {
        return mat.vecmul(vec);
    }

    /// multiplication by a vector: transpose(M) * V
    Vector3 trans_vecmul(real const* V) const
    {
#if MATRIX33SYM_USES_SSE && VECTOR3_USES_AVX
        return vecmul3_sse(V);
#else
        return vecmul_(V);
#endif
    }

    /// multiplication by another matrix: @returns this * M
    const Matrix33sym mul(Matrix33sym const& M) const
    {
        ABORT_NOW("unfinished");
        Matrix33sym res;
        res(0,0) = value(0,0) * M(0,0) + value(0,1) * M(1,0) + value(0,2) * M(2,0);
        res(1,0) = value(1,0) * M(0,0) + value(1,1) * M(1,0) + value(1,2) * M(2,0);
        res(2,0) = value(2,0) * M(0,0) + value(2,1) * M(1,0) + value(2,2) * M(2,0);
        
        res(0,1) = value(0,0) * M(0,1) + value(0,1) * M(1,1) + value(0,2) * M(2,1);
        res(1,1) = value(1,0) * M(0,1) + value(1,1) * M(1,1) + value(1,2) * M(2,1);
        res(2,1) = value(2,0) * M(0,1) + value(2,1) * M(1,1) + value(2,2) * M(2,1);
        
        res(0,2) = value(0,0) * M(0,2) + value(0,1) * M(1,2) + value(0,2) * M(2,2);
        res(1,2) = value(1,0) * M(0,2) + value(1,1) * M(1,2) + value(1,2) * M(2,2);
        res(2,2) = value(2,0) * M(0,2) + value(2,1) * M(1,2) + value(2,2) * M(2,2);
        return res;
    }
    
    /// multiplication with a matrix
    friend Matrix33sym operator * (Matrix33sym const& mat, Matrix33sym const& mut)
    {
        return mat.mul(mut);
    }

    /// multiplication by another matrix: @returns transpose(this) * M
    const Matrix33sym trans_mul(Matrix33sym const& M) const
    {
        ABORT_NOW("unfinished");
        Matrix33sym res;
        res(0,0) = value(0,0) * M(0,0) + value(1,0) * M(1,0) + value(2,0) * M(2,0);
        res(1,0) = value(0,1) * M(0,0) + value(1,1) * M(1,0) + value(2,1) * M(2,0);
        res(2,0) = value(0,2) * M(0,0) + value(1,2) * M(1,0) + value(2,2) * M(2,0);
        
        res(0,1) = value(0,0) * M(0,1) + value(1,0) * M(1,1) + value(2,0) * M(2,1);
        res(1,1) = value(0,1) * M(0,1) + value(1,1) * M(1,1) + value(2,1) * M(2,1);
        res(2,1) = value(0,2) * M(0,1) + value(1,2) * M(1,1) + value(2,2) * M(2,1);
        
        res(0,2) = value(0,0) * M(0,2) + value(1,0) * M(1,2) + value(2,0) * M(2,2);
        res(1,2) = value(0,1) * M(0,2) + value(1,1) * M(1,2) + value(2,1) * M(2,2);
        res(2,2) = value(0,2) * M(0,2) + value(1,2) * M(1,2) + value(2,2) * M(2,2);
        return res;
    }
    
    
    /// add full matrix: this <- this + M
    void add_full(Matrix33sym const& M)
    {
        real const* src = M.val;
#if MATRIX33SYM_USES_SSE
        store2(val  , add2(load2(src  ), load2(val  )));
        store2(val+2, add2(load2(src+2), load2(val+2)));
        store2(val+4, add2(load2(src+4), load2(val+4)));
#else
        for ( index u = 0; u < 6; ++u )
            val[u] += src[u];
#endif
    }
    
    /// add full matrix: this <- this + alpha * M
    void add_full(const real alpha, Matrix33sym const& M)
    {
        real const* src = M.val;
#if MATRIX33SYM_USES_SSE
        vec2 a = set2(alpha);
        store2(val  , fmadd2(a, load2(src  ), load2(val  )));
        store2(val+2, fmadd2(a, load2(src+2), load2(val+2)));
        store2(val+4, fmadd2(a, load2(src+4), load2(val+4)));
#else
        for ( index u = 0; u < 6; ++u )
            val[u] += alpha * src[u];
#endif
    }
    
    /// sub full matrix: this <- this - M
    void sub_full(Matrix33sym const& M)
    {
        real const* src = M.val;
#if MATRIX33SYM_USES_SSE
        store2(val  , sub2(load2(val  ), load2(src  )));
        store2(val+2, sub2(load2(val+2), load2(src+2)));
        store2(val+4, sub2(load2(val+4), load2(src+4)));
#else
        for ( index u = 0; u < 6; ++u )
            val[u] -= src[u];
#endif
    }
    
    /// subtract full matrix: this <- this - alpha * M
    void sub_full(const real alpha, Matrix33sym const& M)
    {
        real const* src = M.val;
#if MATRIX33SYM_USES_SSE
        vec2 a = set2(alpha);
        store2(val  , fnmadd2(a, load2(src  ), load2(val  )));
        store2(val+2, fnmadd2(a, load2(src+2), load2(val+2)));
        store2(val+4, fnmadd2(a, load2(src+4), load2(val+4)));
#else
        for ( index u = 0; u < 6; ++u )
            val[u] -= alpha * src[u];
#endif
    }

    /// subtract transposed matrix: this <- this - transposed(M)
    void sub_trans(Matrix33sym const& M)
    {
        real const* src = M.val;
        for ( index u = 0; u < 6; ++u )
            val[u] -= src[u];
    }
    
    /// add transposed matrix: this <- this + alpha * transposed(M)
    void add_trans(Matrix33sym const& M)
    {
        real const* src = M.val;
        for ( index u = 0; u < 6; ++u )
            val[u] += src[u];
    }
    
    /// add transposed matrix: this <- this + alpha * transposed(M)
    void add_trans(const real alpha, Matrix33sym const& M)
    {
        real const* src = M.val;
        for ( index u = 0; u < 6; ++u )
            val[u] += alpha * src[u];
    }

    /// add lower triangle of matrix including diagonal: this <- this + M
    void add_half(Matrix33sym const& M)
    {
        real const* src = M.val;
#if MATRIX33SYM_USES_SSE
        store2(val  , add2(load2(src  ), load2(val  )));
        store2(val+2, add2(load2(src+2), load2(val+2)));
        store2(val+4, add2(load2(src+4), load2(val+4)));
#else
        for ( index u = 0; u < 6; ++u )
            val[u] += src[u];
#endif
    }
    
    /// add lower triangle of matrix including diagonal: this <- this + alpha * M
    void add_half(const real alpha, Matrix33sym const& M)
    {
        real const* src = M.val;
        //std::clog << "matrix alignment " << ((uintptr_t)src & 63) << "\n";
#if MATRIX33SYM_USES_SSE
        vec2 a = set2(alpha);
        store2(val  , fmadd2(a, load2(src  ), load2(val  )));
        store2(val+2, fmadd2(a, load2(src+2), load2(val+2)));
        store2(val+4, fmadd2(a, load2(src+4), load2(val+4)));
#else
        for ( index u = 0; u < 6; ++u )
            val[u] += alpha * src[u];
#endif
    }
    
    /// subtract lower triangle of matrix including diagonal: this <- this - M
    void sub_half(Matrix33sym const& M)
    {
        real const* src = M.val;
#if MATRIX33SYM_USES_SSE
        store2(val  , sub2(load2(val  ), load2(src  )));
        store2(val+2, sub2(load2(val+2), load2(src+2)));
        store2(val+4, sub2(load2(val+4), load2(src+4)));
#else
        for ( index u = 0; u < 6; ++u )
            val[u] -= src[u];
#endif
    }
    
    /// add alpha to diagonal
    void add_diag(real alpha)
    {
        val[0] += alpha;
        val[3] += alpha;
        val[5] += alpha;
    }
    
    /// add -alpha to diagonal
    void sub_diag(real alpha)
    {
        val[0] -= alpha;
        val[3] -= alpha;
        val[5] -= alpha;
    }

    /// add all elements of block 'S' to array 'M'
    void addto(real * M, size_t ldd) const
    {
        M[0      ] += val[0];
        M[1      ] += val[1];
        M[2      ] += val[2];
        M[  ldd  ] += val[1];
        M[1+ldd  ] += val[3];
        M[2+ldd  ] += val[4];
        M[  ldd*2] += val[2];
        M[1+ldd*2] += val[4];
        M[2+ldd*2] += val[5];
    }
    
    /// add scaled elements of block 'S' to array 'M'
    void addto(real * M, size_t ldd, real alpha) const
    {
        M[0      ] += alpha * val[0];
        M[1      ] += alpha * val[1];
        M[2      ] += alpha * val[2];
        M[  ldd  ] += alpha * val[1];
        M[1+ldd  ] += alpha * val[3];
        M[2+ldd  ] += alpha * val[4];
        M[  ldd*2] += alpha * val[2];
        M[1+ldd*2] += alpha * val[4];
        M[2+ldd*2] += alpha * val[5];
    }
    
    /// add lower elements of this block to lower triangle of 'M'
    void addto_lower(real * M, size_t ldd) const
    {
        M[0      ] += val[0];
        M[1      ] += val[1];
        M[2      ] += val[2];
        M[1+ldd  ] += val[3];
        M[2+ldd  ] += val[4];
        M[2+ldd*2] += val[5];
    }
    
    /// add scaled lower elements of this block to lower triangle of 'M'
    void addto_lower(real * M, size_t ldd, real alpha) const
    {
        M[0      ] += alpha * val[0];
        M[1      ] += alpha * val[1];
        M[2      ] += alpha * val[2];
        M[1+ldd  ] += alpha * val[3];
        M[2+ldd  ] += alpha * val[4];
        M[2+ldd*2] += alpha * val[5];
    }
    
    /// add lower elements of this block to both upper and lower triangles of 'M'
    void addto_symm(real * M, size_t ldd) const
    {
        M[0      ] += val[0];
        M[1      ] += val[1];
        M[2      ] += val[2];
        M[  ldd  ] += val[1];
        M[1+ldd  ] += val[3];
        M[2+ldd  ] += val[4];
        M[  ldd*2] += val[2];
        M[1+ldd*2] += val[4];
        M[2+ldd*2] += val[5];
    }
    
    /// add all elements of this block to 'M', with transposition
    void addto_trans(real * M, size_t ldd) const
    {
        M[0      ] += val[0];
        M[1      ] += val[1];
        M[2      ] += val[2];
        M[  ldd  ] += val[1];
        M[1+ldd  ] += val[3];
        M[2+ldd  ] += val[4];
        M[  ldd*2] += val[2];
        M[1+ldd*2] += val[4];
        M[2+ldd*2] += val[5];
    }
    
#pragma mark -

    /// return diagonal Matrix from diagonal terms
    static Matrix33sym diagonal(real a, real b, real c)
    {
        return Matrix33sym(a, 0, 0, b, 0, c);
    }
    
    /// identity matrix
    static Matrix33sym identity()
    {
        return Matrix33sym(0, 1);
    }

    /// return a symmetric matrix: [ dir (x) transpose(dir) ]
    static Matrix33sym outerProduct(Vector3 const& dir)
    {
        return Matrix33sym(dir.XX*dir.XX, dir.YY*dir.XX, dir.ZZ*dir.XX,
                           dir.YY*dir.YY, dir.YY*dir.ZZ,
                           dir.ZZ*dir.ZZ );
    }

    /// return a symmetric matrix: alpha * [ dir (x) transpose(dir) ]
    static Matrix33sym outerProduct(Vector3 const& dir, real alpha)
    {
        real XX = dir.XX * alpha;
        real YY = dir.YY * alpha;
        real ZZ = dir.ZZ * alpha;
        return Matrix33sym(dir.XX*XX, dir.YY*XX, dir.ZZ*XX,
                           dir.YY*YY, dir.YY*ZZ,
                           dir.ZZ*ZZ );
    }
    
    /// return outer product: [ dir (x) transpose(vec) ]
    static Matrix33sym outerProduct(Vector3 const& dir, Vector3 const& vec)
    {
        ABORT_NOW("non symmetric");
        return Matrix33sym(0, 0);
    }
    
    /// return outer product: [ dir (x) transpose(vec) ]
    static Matrix33sym outerProduct(real const* dir, real const* vec)
    {
        ABORT_NOW("non symmetric");
        return Matrix33sym(0, 0);
    }
    
    /// add outer product: [ dir (x) transpose(vec) ]
    void addOuterProduct(real const* dir, real const* vec)
    {
        ABORT_NOW("non symmetric");
    }

    /// return [ dir (x) transpose(vec) + vec (x) transpose(dir) ]
    static Matrix33sym symmetricOuterProduct(Vector3 const& dir, Vector3 const& vec)
    {
        real X = dir.XX * vec.XX;
        real Y = dir.YY * vec.YY;
        real Z = dir.ZZ * vec.ZZ;
        return Matrix33sym(X+X, dir.YY*vec.XX + dir.XX*vec.YY, dir.ZZ*vec.XX + dir.XX*vec.ZZ,
                           Y+Y, dir.ZZ*vec.YY + dir.YY*vec.ZZ,
                           Z+Z);
    }
 
    /// return symmetric matrix block :  dia * Id + [ dir (x) dir ] * len
    static Matrix33sym offsetOuterProduct(const real dia, Vector3 const& dir, const real len)
    {
        real X = dir.XX * len;
        real Y = dir.YY * len;
        real Z = dir.ZZ * len;
        return Matrix33sym(X * dir.XX + dia, Y * dir.XX, Z * dir.XX,
                           Y * dir.YY + dia, Z * dir.YY,
                           Z * dir.ZZ + dia);
    }
    
    /// build the rotation matrix `M = 2 * dir (x) dir - Id` of axis `dir` and angle 180
    static Matrix33sym householder(const Vector3& axis)
    {        
        real X = axis.XX, Y = axis.YY, Z = axis.ZZ;
        real X2 = X + X, Y2 = Y + Y, Z2 = Z + Z;
        
        return Matrix33sym(X * X2 - 1.0, X * Y2, X * Z2,
                           Y * Y2 - 1.0, Y * Z2,
                           Z * Z2 - 1.0);
    }

    
    /// build the matrix `dia * Id + vec (x) Id`
    /**
     Thus applying M to V results in `dia * V + vec (x) V`
     */
    static Matrix33sym vectorProduct(const real dia, const Vector3& vec)
    {
        ABORT_NOW("non symmetric");
        return Matrix33sym(dia, 0);
    }
    
    /// 2D Rotation around `axis` (of norm 1) with angle set by cosine and sine values
    /**
     This is not a 3D rotation and the images of 'axis' is zero!
     This is equivalent to rotationAroundAxis(axis, c, s) - outerProduct(axis)
     Attention: This is meant to be called with `norm(axis)==1` and `c*c + s*s == 1`
     but the values of 'c' and 's' can be tweaked to scale the resulting matrix.
     */
    static Matrix33sym planarRotation(const Vector3& axis, const real c, const real s)
    {
        ABORT_NOW("non symmetric");
        return Matrix33sym(0, 0);
    }
    
    static Matrix33sym rotationAroundAxis(const Vector3& axis, const real c, const real s)
    {
        ABORT_NOW("non symmetric");
        return Matrix33sym(0, 0);
    }
};


/// output operator
inline std::ostream& operator << (std::ostream& os, Matrix33sym const& arg)
{
    std::ios::fmtflags fgs = os.flags();
    arg.print_smart(os);
    os.setf(fgs);
    return os;
}

#endif

