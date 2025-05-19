// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University.
// FJN, Strasbourg 08.06.2018

#ifndef MATRIX33
#define MATRIX33

#include "real.h"
#include "vector3.h"

#include <cstdio>
#include <iostream>
#include <charconv>

/// BLD is the leading dimension of the matrix
/**
 The code works with BLD = 3 or 4, and typically memory storage is less with 3,
 but performance might be better with 4, as SIMD-AVX calls handle doubles by 4.
 */
#include "simd.h"

#if defined(__AVX__) && REAL_IS_DOUBLE
#  define MATRIX33_USES_AVX 1
#  define BLD 3
#elif defined(__SSE3__) && !REAL_IS_DOUBLE
#  include "simd_float.h"
#  define MATRIX33_USES_AVX 0
#  define BLD 3
#else
#  define MATRIX33_USES_AVX 0
#  define BLD 3
#endif


/// 3x3 matrix class with 9 'real' elements stored in column order
#if ( BLD == 4 )
class alignas(4*sizeof(real)) Matrix33 final
#else
class Matrix33 final
#endif
{
private:
    
    /// values of the elements
    real val[BLD*3];

    /// access to modifiable element by index
    real& operator[](index_t i)       { return val[i]; }
    
    /// access element value by index
    real  operator[](index_t i) const { return val[i]; }
    
public:
    
    /// clear values that do not represent matrix elements
    void clear_shadow()
    {
#if ( BLD == 4 )
        val[3      ] = 0.0;
        val[3+BLD  ] = 0.0;
        val[3+BLD*2] = 0.0;
#endif
    }

    /// default constructor
    Matrix33() { clear_shadow(); }
    
    /// copy constructor
    Matrix33(Matrix33 const& M)
    {
        for ( index_t u = 0; u < BLD*3; ++u )
            val[u] = M.val[u];
    }

    /// construct Matrix from coordinates (given by columns)
    Matrix33(real a, real b, real c,
             real d, real e, real f,
             real g, real h, real i)
    {
        val[0] = a;
        val[1] = b;
        val[2] = c;
        val[0+BLD] = d;
        val[1+BLD] = e;
        val[2+BLD] = f;
        val[0+BLD*2] = g;
        val[1+BLD*2] = h;
        val[2+BLD*2] = i;
        clear_shadow();
    }

    /// construct Matrix with diagonal terms set to `d` and other terms set to `z`
    Matrix33(real z, real d)
    {
        val[0] = d;
        val[1] = z;
        val[2] = z;
        val[0+BLD] = z;
        val[1+BLD] = d;
        val[2+BLD] = z;
        val[0+BLD*2] = z;
        val[1+BLD*2] = z;
        val[2+BLD*2] = d;
        clear_shadow();
    }

    ~Matrix33() {}
    
#pragma mark -
    
    /// dimensionality
    static constexpr size_t dimension() { return 3; }
    
    /// human-readable identifier
#if ( BLD == 3 )
    static std::string what() { return "9"; }
#else
    static std::string what() { return "+9"; }
#endif
    
    /// set all elements to zero
    void reset()
    {
        for ( index_t u = 0; u < BLD*3; ++u )
            val[u] = 0.0;
    }
    
    /// set diagonal to 'dia' and other elements to 'off'
    void reset1(real off, real dia)
    {
        for ( index_t u = 0; u < BLD*3; ++u )
            val[u] = off;
        for ( index_t u = 0; u < 3; ++u )
            val[u+BLD*u] = dia;
    }
    
    bool operator != (real zero) const
    {
        for ( index_t u = 0; u < BLD*3; ++u )
            if ( val[u] != zero )
                return true;
        return false;
    }
    
    /// conversion to pointer of real
    //operator real const*() const { return val; }

    /// return modifiable pointer of 'real'
    real* data() { return val; }
    
#if defined(__AVX__) && REAL_IS_DOUBLE
#  if ( BLD == 4 )
    vec4 data0() const { return streamload4(val  ); }
    vec4 data1() const { return streamload4(val+4); }
    vec4 data2() const { return streamload4(val+8); }
#  else
    vec4 data0() const { return loadu4(val); }
    vec4 data1() const { return loadu4(val+BLD); }
    vec4 data2() const { return load3Z(val+BLD*2); } // last value may be garbage
#  endif
#elif USE_SIMD && !REAL_IS_DOUBLE
#  if ( BLD == 4 )
    vec4f data0() const { return streamload4f(val  ); }
    vec4f data1() const { return streamload4f(val+4); }
    vec4f data2() const { return streamload4f(val+8); }
#  else
    vec4f data0() const { return loadu4f(val); }
    vec4f data1() const { return loadu4f(val+BLD); }
    vec4f data2() const { return load3Zf(val+BLD*2); } // last value may be garbage
#  endif
#endif

    /// unmodifiable pointer of real
    real const* data() const { return val; }

    /// address of element at line i, column j
    real* addr(const index_t i, const index_t j) { return val + ( i + BLD*j ); }
    /// value of element at line i, column j
    real value(const index_t i, const index_t j) const { return val[i+BLD*j]; }

    /// access functions to element by line and column indices
    real& operator()(const index_t i, const index_t j)       { return val[i+BLD*j]; }
    real  operator()(const index_t i, const index_t j) const { return val[i+BLD*j]; }

    /// set elements from given array
    void load(const real ptr[])
    {
        for ( index_t i = 0; i < 3; ++i )
        for ( index_t j = 0; j < 3; ++j )
            val[i+BLD*j] = ptr[i+3*j];
    }

    /// copy elements to given array
    void store(real ptr[]) const
    {
        for ( index_t i = 0; i < 3; ++i )
        for ( index_t j = 0; j < 3; ++j )
            ptr[i+3*j] = val[i+BLD*j];
    }

    /// extract column vector at given index
    Vector3 column(const index_t i) const
    {
        return Vector3(val+BLD*i);
    }
    
    /// extract line vector at given index
    Vector3 line(const index_t i) const
    {
        return Vector3(val[i], val[BLD+i], val[BLD*2+i]);
    }
    
    /// extract diagonal
    Vector3 diagonal() const
    {
        return Vector3(val[0], val[BLD+1], val[BLD*2+2]);
    }
    
    /// sum of diagonal terms
    real trace() const
    {
        return ( val[0] + val[BLD+1] + val[BLD*2+2] );
    }
    
#pragma mark -
    
    /// set matrix by giving lines
    void setLines(Vector3 const& A, Vector3 const& B, Vector3 const& C)
    {
        val[0      ] = A.XX;
        val[1      ] = B.XX;
        val[2      ] = C.XX;
        val[0+BLD  ] = A.YY;
        val[1+BLD  ] = B.YY;
        val[2+BLD  ] = C.YY;
        val[0+BLD*2] = A.ZZ;
        val[1+BLD*2] = B.ZZ;
        val[2+BLD*2] = C.ZZ;
    }
    
    /// set matrix by giving columns
    void setColumns(Vector3 const& A, Vector3 const& B, Vector3 const& C)
    {
        val[0      ] = A.XX;
        val[1      ] = A.YY;
        val[2      ] = A.ZZ;
        val[0+BLD  ] = B.XX;
        val[1+BLD  ] = B.YY;
        val[2+BLD  ] = B.ZZ;
        val[0+BLD*2] = C.XX;
        val[1+BLD*2] = C.YY;
        val[2+BLD*2] = C.ZZ;
    }

    /// print matrix in human readable format
    void print(FILE * f) const
    {
        fprintf(f, " / %9.3f %+9.3f %+9.3f \\\n", val[0], val[0+BLD], val[0+BLD*2]);
        fprintf(f, " | %9.3f %+9.3f %+9.3f |\n" , val[1], val[1+BLD], val[1+BLD*2]);
        fprintf(f, " \\ %9.3f %+9.3f %+9.3f /\n", val[2], val[2+BLD], val[2+BLD*2]);
    }
    
    /// print [ line1; line2; line3 ]
    void print(std::ostream& os) const
    {
        const int w = (int)os.width();
        os << std::setw(1) << "[";
        for ( index_t i = 0; i < 3; ++i )
        {
            for ( index_t j = 0; j < 3; ++j )
                os << " " << std::fixed << std::setw(w) << value(i,j);
            if ( i < 2 )
                os << ";";
            else
                os << " ]";
        }
    }

    /// value of element at line i, column j
    std::string format_value(const index_t i, const index_t j) const
    {
        char tmp[32];
        real x = val[i+BLD*j];
        snprintf(tmp, sizeof(tmp), " %+6.3f", x);
        return std::string(tmp);
    }

    /// print [ line1; line2; line3 ]
    void print_smart(std::ostream& os) const
    {
        bool s0 = ( abs_real(value(0,1)-value(1,0)) < REAL_EPSILON );
        bool s1 = ( abs_real(value(0,2)-value(2,0)) < REAL_EPSILON );
        bool s2 = ( abs_real(value(1,2)-value(2,1)) < REAL_EPSILON );
        bool sym = ( s0 && s1 && s2 );
        os << "[" << format_value(0,0);
        os << format_value(1,0);
        os << format_value(2,0);
        os << ";";
        if ( sym )
        {
            os << " sym";
            os << format_value(1,1);
            os << format_value(2,1);
            os << ";";
            os << " sym";
            os << " sym";
        }
        else
        {
            os << format_value(0,1);
            os << format_value(1,1);
            os << format_value(2,1);
            os << ";";
            os << format_value(0,2);
            os << format_value(1,2);
        }
        os << format_value(2,2) << " ]";
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
        for ( index_t u = 0; u < BLD*3; ++u )
            val[u] *= alpha;
    }

    /// scale matrix
    void operator *=(const real alpha)
    {
        scale(alpha);
    }
    
    /// return opposite matrix (i.e. -M)
    const Matrix33 operator -() const
    {
        Matrix33 M;
        for ( index_t u = 0; u < BLD*3; ++u )
            M.val[u] = -val[u];
        return M;
    }
    
    /// scaled matrix
    const Matrix33 operator *(const real alpha) const
    {
        Matrix33 res;
        for ( index_t u = 0; u < BLD*3; ++u )
            res.val[u] = val[u] * alpha;
        return res;
    }
    
    /// multiplication by scalar
    friend const Matrix33 operator *(const real alpha, Matrix33 const& mat)
    {
        return mat * alpha;
    }

    /// return sum of two matrices
    const Matrix33 operator +(Matrix33 const& M) const
    {
        Matrix33 res;
        for ( index_t u = 0; u < BLD*3; ++u )
            res.val[u] = val[u] + M.val[u];
        return res;
    }

    /// return difference of two matrices
    const Matrix33 operator -(Matrix33 const& M) const
    {
        Matrix33 res;
        for ( index_t u = 0; u < BLD*3; ++u )
            res.val[u] = val[u] - M.val[u];
        return res;
    }

    /// subtract given matrix
    void operator += (Matrix33 const& M)
    {
#if MATRIX33_USES_AVX && ( BLD == 4 )
        store4(val  , add4(load4(val  ), load4(M.val  )));
        store4(val+4, add4(load4(val+4), load4(M.val+4)));
        store4(val+8, add4(load4(val+8), load4(M.val+8)));
#else
        for ( index_t u = 0; u < BLD*3; ++u )
            val[u] += M.val[u];
#endif
    }

    /// add given matrix
    void operator -= (Matrix33 const& M)
    {
#if MATRIX33_USES_AVX && ( BLD == 4 )
        store4(val  , sub4(load4(val  ), load4(M.val  )));
        store4(val+4, sub4(load4(val+4), load4(M.val+4)));
        store4(val+8, sub4(load4(val+8), load4(M.val+8)));
#else
        for ( index_t u = 0; u < BLD*3; ++u )
            val[u] -= M.val[u];
#endif
    }
    
    /// transpose matrix in place
    void transpose()
    {
        std::swap(val[1], val[0+BLD]);
        std::swap(val[2], val[0+BLD*2]);
        std::swap(val[2+BLD], val[1+BLD*2]);
    }
    
    /// return transposed matrix
    Matrix33 transposed() const
    {
        Matrix33 res;
#if MATRIX33_USES_AVX && ( BLD == 4 )
        vec4 m012 = load4(val  );
        vec4 m345 = load4(val+4);
        vec4 m678 = load4(val+8);
        vec4 z = shuffle4(m012, m345, 0b0011);
        vec4 u = catshift2(z, m678);
        vec4 t = shuffle4(m012, m345, 0b1000);
        store4(res.val  , blend0010(t, u));
        store4(res.val+4, blend22(z, shuffle4(u, m345, 0b1100)));
        store4(res.val+8, blend22(u, m678));
#else
        for ( index_t x = 0; x < 3; ++x )
        for ( index_t y = 0; y < 3; ++y )
            res[y+BLD*x] = val[x+BLD*y];
#endif
        return res;
    }
    
    /// return scaled transposed matrix
    Matrix33 transposed(real alpha) const
    {
        Matrix33 res;
#if MATRIX33_USES_AVX && ( BLD == 4 )
        vec4 a = set4(alpha);
        vec4 m012 = mul4(a, load4(val  ));
        vec4 m345 = mul4(a, load4(val+4));
        vec4 m678 = mul4(a, load4(val+8));
        vec4 z = shuffle4(m012, m345, 0b0011);
        vec4 u = catshift2(z, m678);
        vec4 t = shuffle4(m012, m345, 0b1000);
        store4(res.val  , blend0010(t, u));
        store4(res.val+4, blend22(z, shuffle4(u, m345, 0b1100)));
        store4(res.val+8, blend22(u, m678));
#else
        for ( index_t x = 0; x < 3; ++x )
        for ( index_t y = 0; y < 3; ++y )
            res[y+BLD*x] = alpha * val[x+BLD*y];
#endif
        return res;
    }

    
    /// maximum of all component's absolute values
    real norm_inf() const
    {
        real res = abs_real(val[0]);
        for ( index_t i = 1; i < 3*BLD; ++i )
            res = max_real(res, abs_real(val[i]));
        return res;
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
        real det = 1.0 / determinant();
        Vector3 X = column(0);
        Vector3 Y = column(1);
        Vector3 Z = column(2);
        setLines(cross(Y,Z)*det, cross(Z,X)*det, cross(X,Y)*det);
    }

    /// return inverse matrix
    Matrix33 inverted() const
    {
        Matrix33 res;
        real det = 1.0 / determinant();
        Vector3 X = column(0);
        Vector3 Y = column(1);
        Vector3 Z = column(2);
        res.setLines(cross(Y,Z)*det, cross(Z,X)*det, cross(X,Y)*det);
        //std::clog << " mat * inverse = " << mul(res).to_string(10, 3) << "\n";
        return res;
    }
    
    /// inversion of a symmetric matrix, using values in lower triangle
    /** This methods uses a L*D*L^t factorization with:
     L = ( 1 0 0; a 1 0; b c 1 )
     D = ( u 0 0; 0 v 0; 0 0 w )
     The result is a symmetric matrix
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
        real v = val[1+BLD] - a * val[1];
        real iv = 1.0 / v;
        real x = val[2+BLD] - a * val[2];
        real c = x * iv;
        real iw = 1.0 / ( val[2+BLD*2] - b * val[2] - c * x );
        // inverse triangular matrix U = inverse(L^t):
        b = -b + a * c;
        a = -a;
        c = -c;
        real aiv = a * iv;
        real biw = b * iw;
        real ciw = c * iw;
        // calculate U * inverse(D) * U^t:
        val[0+BLD*0] = iu + a * aiv + b * biw;
        val[1+BLD*0] = aiv + c * biw;
        val[2+BLD*0] = biw;
        val[0+BLD*1] = val[1+BLD*0];
        val[1+BLD*1] = iv + c * ciw;
        val[2+BLD*1] = ciw;
        val[0+BLD*2] = biw;
        val[1+BLD*2] = ciw;
        val[2+BLD*2] = iw;
        return 0;
    }

    /// copy values from lower triangle to upper triangle
    void copy_lower()
    {
        val[0+BLD  ] = val[1];
        val[0+BLD*2] = val[2];
        val[1+BLD*2] = val[2+BLD];
    }

    /// relative asymmetry of 3x3 submatrix (divided by the trace)
    real asymmetry() const
    {
        real t = abs_real(value(0,0)) + abs_real(value(1,1)) + abs_real(value(2,2));
        real a = abs_real(value(0,1) - value(1,0));
        real b = abs_real(value(0,2) - value(2,0));
        real c = abs_real(value(2,1) - value(1,2));
        return ( a + b + c ) / t;
    }
    
#pragma mark -

#if MATRIX33_USES_AVX
    
    /// multiplication by a 3-components vector: this * V
    const vec4 vecmul3_avx(vec4 vec) const
    {
        vec4 xy = duplo2f128(vec); // xyxy
        vec4 zzzz = duplo4(duphi2f128(vec));
        vec = mul4(data0(), duplo4(xy)); // xxxx
        xy = mul4(data1(), duphi4(xy)); // yyyy
        zzzz = fmadd4(data2(), zzzz, add4(xy, vec));
        return clear4th(zzzz);
    }

    /// multiplication by a vector: this * V
    const vec4 vecmul3_avx(double const* V) const
    {
        vec4 zzzz = broadcast2(V); // xyxy
        vec4 xxxx = duplo4(zzzz); //broadcast1(V);
        vec4 yyyy = duphi4(zzzz); //broadcast1(V+1);
        zzzz = broadcast1(V+2);
        xxxx = mul4(data0(), xxxx);
        yyyy = mul4(data1(), yyyy);
        zzzz = fmadd4(data2(), zzzz, add4(xxxx, yyyy));
        return clear4th(zzzz);
    }

    /// transpose-multiplication by a vector: transpose(M) * V
    const vec4 trans_vecmul3_avx(double const* V) const
    {
        vec4 zz = setzero4();
        vec4 vec = blend31(loadu4(V), zz); // { x, y, z, 0 }
        vec4 s0 = mul4(data0(), vec);
        vec4 s1 = mul4(data1(), vec);
        vec4 s2 = mul4(data2(), vec);
        s0 = add4(unpacklo4(s0, s1), unpackhi4(s0, s1));
        s2 = add4(unpacklo4(s2, zz), unpackhi4(s2, zz));
        return add4(catshift2(s0, s2), blend22(s0, s2));
    }
#endif
    
    /// multiplication by a vector: this * V
    Vector3 vecmul_(Vector3 const& V) const
    {
        return Vector3(val[0] * V.XX + val[  BLD] * V.YY + val[  BLD*2] * V.ZZ,
                       val[1] * V.XX + val[1+BLD] * V.YY + val[1+BLD*2] * V.ZZ,
                       val[2] * V.XX + val[2+BLD] * V.YY + val[2+BLD*2] * V.ZZ);
    }
    
    /// multiplication by a vector: this * V
    Vector3 vecmul_(real const* R) const
    {
        return Vector3(val[0] * R[0] + val[  BLD] * R[1] + val[  BLD*2] * R[2],
                       val[1] * R[0] + val[1+BLD] * R[1] + val[1+BLD*2] * R[2],
                       val[2] * R[0] + val[2+BLD] * R[1] + val[2+BLD*2] * R[2]);
    }
    
    /// transpose-multiplication by a vector: transpose(M) * V
    Vector3 trans_vecmul_(Vector3 const& V) const
    {
        return Vector3(val[0    ] * V.XX + val[1      ] * V.YY + val[2      ] * V.ZZ,
                       val[BLD  ] * V.XX + val[1+BLD  ] * V.YY + val[2+BLD  ] * V.ZZ,
                       val[BLD*2] * V.XX + val[1+BLD*2] * V.YY + val[2+BLD*2] * V.ZZ);
    }

    /// transpose-multiplication by a vector: transpose(M) * V
    Vector3 trans_vecmul_(real const* R) const
    {
        return Vector3(val[0    ] * R[0] + val[1      ] * R[1] + val[2      ] * R[2],
                       val[BLD  ] * R[0] + val[1+BLD  ] * R[1] + val[2+BLD  ] * R[2],
                       val[BLD*2] * R[0] + val[1+BLD*2] * R[1] + val[2+BLD*2] * R[2]);
    }

    /// multiplication by a vector: this * V
    Vector3 vecmul(Vector3 const& vec) const
    {
#if MATRIX33_USES_AVX && VECTOR3_USES_AVX
        return vecmul3_avx(vec.xyz);
#elif MATRIX33_USES_AVX
        return vecmul3_avx(vec.data());
#else
        return vecmul_(vec);
#endif
    }
    
    /// multiplication by a vector: this * { ptr[0], ptr[1] }
    Vector3 vecmul(real const* ptr) const
    {
#if MATRIX33_USES_AVX
        return vecmul3_avx(ptr);
#else
        return vecmul_(ptr);
#endif
    }

    /// multiplication with a vector: M * V
    friend Vector3 operator * (Matrix33 const& mat, Vector3 const& vec)
    {
        return mat.vecmul(vec);
    }

    /// transpose-multiplication by a vector: transpose(M) * V
    Vector3 trans_vecmul(real const* V) const
    {
#if MATRIX33_USES_AVX
        return trans_vecmul3_avx(V);
#else
        return trans_vecmul_(V);
#endif
    }

    /// multiplication by another matrix: @returns this * M
    const Matrix33 mul(Matrix33 const& M) const
    {
        Matrix33 res;
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
    friend Matrix33 operator * (Matrix33 const& mat, Matrix33 const& mut)
    {
        return mat.mul(mut);
    }

    /// multiplication by another matrix: @returns transpose(this) * M
    const Matrix33 trans_mul(Matrix33 const& M) const
    {
        Matrix33 res;
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
    void add_full(Matrix33 const& M)
    {
        real const* src = M.val;
#if MATRIX33_USES_AVX && ( BLD == 4 )
        store4(val  , add4(load4(val  ), load4(src  )));
        store4(val+4, add4(load4(val+4), load4(src+4)));
        store4(val+8, add4(load4(val+8), load4(src+8)));
#elif MATRIX33_USES_AVX
        storeu4(val  , add4(loadu4(val  ), loadu4(src  )));
        storeu4(val+4, add4(loadu4(val+4), loadu4(src+4)));
        store1(val+8, add1(load1(val+8), load1(src+8)));
#else
        for ( index_t u = 0; u < BLD*3; ++u )
            val[u] += src[u];
#endif
    }
    
    /// add full matrix: this <- this + alpha * M
    void add_full(const real alpha, Matrix33 const& M)
    {
        real const* src = M.val;
#if MATRIX33_USES_AVX && ( BLD == 4 )
        vec4 a = set4(alpha);
        store4(val  , fmadd4(a, load4(src  ), load4(val  )));
        store4(val+4, fmadd4(a, load4(src+4), load4(val+4)));
        store4(val+8, fmadd4(a, load4(src+8), load4(val+8)));
#else
        for ( index_t u = 0; u < BLD*3; ++u )
            val[u] += alpha * src[u];
#endif
    }
    
    /// sub full matrix: this <- this - M
    void sub_full(Matrix33 const& M)
    {
        real const* src = M.val;
#if MATRIX33_USES_AVX && ( BLD == 4 )
        store4(val  , sub4(load4(val  ), load4(src  )));
        store4(val+4, sub4(load4(val+4), load4(src+4)));
        store4(val+8, sub4(load4(val+8), load4(src+8)));
#else
        for ( index_t u = 0; u < BLD*3; ++u )
            val[u] -= src[u];
#endif
    }
    
    /// subtract full matrix: this <- this - alpha * M
    void sub_full(const real alpha, Matrix33 const& M)
    {
        real const* src = M.val;
#if MATRIX33_USES_AVX && ( BLD == 4 )
        vec4 a = set4(alpha);
        store4(val  , fnmadd4(a, load4(src  ), load4(val  )));
        store4(val+4, fnmadd4(a, load4(src+4), load4(val+4)));
        store4(val+8, fnmadd4(a, load4(src+8), load4(val+8)));
#else
        for ( index_t u = 0; u < BLD*3; ++u )
            val[u] -= alpha * src[u];
#endif
    }

    /// subtract transposed matrix: this <- this - transposed(M)
    void sub_trans(Matrix33 const& M)
    {
        real const* src = M.val;
        for ( int x = 0; x < 3; ++x )
        for ( int y = 0; y < 3; ++y )
            val[y+BLD*x] -= src[x+BLD*y];
    }
    
    /// add transposed matrix: this <- this + alpha * transposed(M)
    void add_trans(Matrix33 const& M)
    {
        real const* src = M.val;
        for ( int x = 0; x < 3; ++x )
        for ( int y = 0; y < 3; ++y )
            val[y+BLD*x] += src[x+BLD*y];
    }
    
    /// add transposed matrix: this <- this + alpha * transposed(M)
    void add_trans(const real alpha, Matrix33 const& M)
    {
        real const* src = M.val;
        for ( int x = 0; x < 3; ++x )
        for ( int y = 0; y < 3; ++y )
            val[y+BLD*x] += alpha * src[x+BLD*y];
    }

    /// add lower triangle of matrix including diagonal: this <- this + M
    void add_half(Matrix33 const& M)
    {
        real const* src = M.val;
#if MATRIX33_USES_AVX && ( BLD == 4 )
        store4(val  , add4(load4(val  ), load4(src  )));
        store4(val+4, add4(load4(val+4), load4(src+4)));
        store4(val+8, add4(load4(val+8), load4(src+8)));
#elif ( 1 )
        for ( index_t u = 0; u < BLD*3; ++u )
            val[u] += src[u];
#else
        for ( int x = 0; x < 3; ++x )
        for ( int y = x; y < 3; ++y )
            val[y+BLD*x] += src[y+BLD*x];
#endif
    }
    
    /// add lower triangle of matrix including diagonal: this <- this + alpha * M
    void add_half(const real alpha, Matrix33 const& M)
    {
        real const* src = M.val;
        //std::clog << "matrix alignment " << ((uintptr_t)src & 63) << "\n";
#if MATRIX33_USES_AVX && ( BLD == 4 )
        vec4 a = set4(alpha);
        store4(val  , fmadd4(a, load4(src  ), load4(val  )));
        store4(val+4, fmadd4(a, load4(src+4), load4(val+4)));
        store4(val+8, fmadd4(a, load4(src+8), load4(val+8)));
#elif ( 1 )
        val[0] += alpha * src[0];
        val[1] += alpha * src[1];
        val[2] += alpha * src[2];
        val[1+BLD] += alpha * src[1+BLD];
        val[2+BLD] += alpha * src[2+BLD];
        val[2+BLD*2] += alpha * src[2+BLD*2];
#elif ( 0 )
        for ( index_t u = 0; u < BLD*3; ++u )
            val[u] += alpha * src[u];
#else
        for ( int x = 0; x < 3; ++x )
        for ( int y = x; y < 3; ++y )
            val[y+BLD*x] += alpha * src[y+BLD*x];
#endif
    }
    
    /// add lower triangle of matrix including diagonal: this <- this + alpha * M
    void add_half(const real alpha, Matrix33 const& M, const real diag)
    {
        real const* src = M.val;
        //std::clog << "matrix alignment " << ((uintptr_t)src & 63) << "\n";
#if ( 1 )
        val[0] += alpha * ( src[0] + diag );
        val[1] += alpha * src[1];
        val[2] += alpha * src[2];
        val[1+BLD] += alpha * ( src[1+BLD] + diag );
        val[2+BLD] += alpha * src[2+BLD];
        val[2+BLD*2] += alpha * ( src[2+BLD*2] + diag );
#else
        for ( int x = 0; x < 3; ++x )
        {
            val[y+BLD*x] += alpha * ( src[y+BLD*x] + dia );
            for ( int y = x+1; y < 3; ++y )
                val[y+BLD*x] += alpha * src[y+BLD*x];
        }
#endif
    }

    /// subtract lower triangle of matrix including diagonal: this <- this - M
    void sub_half(Matrix33 const& M)
    {
        real const* src = M.val;
#if MATRIX33_USES_AVX && ( BLD == 4 )
        store4(val  , sub4(load4(val  ), load4(src  )));
        store4(val+4, sub4(load4(val+4), load4(src+4)));
        store4(val+8, sub4(load4(val+8), load4(src+8)));
#elif ( 1 )
        for ( index_t u = 0; u < BLD*3; ++u )
            val[u] -= src[u];
#else
        for ( int x = 0; x < 3; ++x )
        for ( int y = x; y < 3; ++y )
            val[y+BLD*x] -= src[y+BLD*x];
#endif
    }
    
    /// add alpha to diagonal
    void add_diag(real alpha)
    {
        val[0]       += alpha;
        val[1+BLD]   += alpha;
        val[2+BLD*2] += alpha;
    }
    
    /// return copy of *this, with `alpha` added to the diagonal
    Matrix33 plus_diagonal(real alpha) const
    {
        Matrix33 res;
        res.val[0+BLD*0] = val[0+BLD*0] + alpha;
        res.val[1+BLD*0] = val[1+BLD*0];
        res.val[2+BLD*0] = val[2+BLD*0];
        res.val[0+BLD*1] = val[0+BLD*1];
        res.val[1+BLD*1] = val[1+BLD*1] + alpha;
        res.val[2+BLD*1] = val[2+BLD*1];
        res.val[0+BLD*2] = val[0+BLD*2];
        res.val[1+BLD*2] = val[1+BLD*2];
        res.val[2+BLD*2] = val[2+BLD*2] + alpha;
        return res;
    }
    
    /// add all elements of block 'S' to array 'M'
    void addto(real * M, size_t ldd) const
    {
        M[0      ] += val[0];
        M[1      ] += val[1];
        M[2      ] += val[2];
        M[  ldd  ] += val[0+BLD];
        M[1+ldd  ] += val[1+BLD];
        M[2+ldd  ] += val[2+BLD];
        M[  ldd*2] += val[0+BLD*2];
        M[1+ldd*2] += val[1+BLD*2];
        M[2+ldd*2] += val[2+BLD*2];
    }
    
    /// add scaled elements of block 'S' to array 'M'
    void addto(real * M, size_t ldd, real alpha) const
    {
        M[0      ] += alpha * val[0];
        M[1      ] += alpha * val[1];
        M[2      ] += alpha * val[2];
        M[  ldd  ] += alpha * val[0+BLD];
        M[1+ldd  ] += alpha * val[1+BLD];
        M[2+ldd  ] += alpha * val[2+BLD];
        M[  ldd*2] += alpha * val[0+BLD*2];
        M[1+ldd*2] += alpha * val[1+BLD*2];
        M[2+ldd*2] += alpha * val[2+BLD*2];
    }
    
    /// add lower elements of this block to lower triangle of 'M'
    void addto_lower(real * M, size_t ldd) const
    {
        M[0      ] += val[0];
        M[1      ] += val[1];
        M[2      ] += val[2];
        M[1+ldd  ] += val[1+BLD];
        M[2+ldd  ] += val[2+BLD];
        M[2+ldd*2] += val[2+BLD*2];
    }
    
    /// add scaled lower elements of this block to lower triangle of 'M'
    void addto_lower(real * M, size_t ldd, real alpha) const
    {
        M[0      ] += alpha * val[0];
        M[1      ] += alpha * val[1];
        M[2      ] += alpha * val[2];
        M[1+ldd  ] += alpha * val[1+BLD];
        M[2+ldd  ] += alpha * val[2+BLD];
        M[2+ldd*2] += alpha * val[2+BLD*2];
    }
    
    /// add lower elements of this block to both upper and lower triangles of 'M'
    void addto_symm(real * M, size_t ldd) const
    {
        M[0      ] += val[0];
        M[1      ] += val[1];
        M[2      ] += val[2];
        M[  ldd  ] += val[1];
        M[1+ldd  ] += val[1+BLD];
        M[2+ldd  ] += val[2+BLD];
        M[  ldd*2] += val[2];
        M[1+ldd*2] += val[2+BLD];
        M[2+ldd*2] += val[2+BLD*2];
    }
    
    /// add all elements of this block to 'M', with transposition
    void addto_trans(real * M, size_t ldd) const
    {
        M[0      ] += val[0];
        M[1      ] += val[0+BLD];
        M[2      ] += val[0+BLD*2];
        M[  ldd  ] += val[1];
        M[1+ldd  ] += val[1+BLD];
        M[2+ldd  ] += val[1+BLD*2];
        M[  ldd*2] += val[2];
        M[1+ldd*2] += val[2+BLD];
        M[2+ldd*2] += val[2+BLD*2];
    }
    
    
    /// return symmetric Matrix from coordinates (column-major, lower triangle)
    static Matrix33 symmetric(real a, real b, real c,
                              real d, real e, real f)
    {
        return Matrix33(a, b, c, b, d, e, c, e, f);
    }
    
#pragma mark -

    /// return diagonal Matrix from diagonal terms
    static Matrix33 diagonal(real a, real b, real c)
    {
        return Matrix33(a, 0, 0, 0, b, 0, 0, 0, c);
    }
    
    /// identity matrix
    static Matrix33 one()
    {
        return Matrix33(0, 1);
    }

    /// return a symmetric matrix: [ dir (x) transpose(dir) ]
    static Matrix33 outerProduct(Vector3 const& dir)
    {
        return symmetric(dir.XX*dir.XX, dir.YY*dir.XX, dir.ZZ*dir.XX,
                         dir.YY*dir.YY, dir.YY*dir.ZZ,
                         dir.ZZ*dir.ZZ );
    }

    /// return a symmetric matrix: alpha * [ dir (x) transpose(dir) ]
    static Matrix33 outerProduct(Vector3 const& dir, real alpha)
    {
#if MATRIX33_USES_AVX && ( BLD == 4 )
        Matrix33 res;
        vec4 d = dir;
        vec4 p = swap2f128(d);
        vec4 l = blend22(d, p);
        vec4 u = blend22(p, d);
        vec4 a = mul4(d, set4(alpha));
        store4(res.val  , mul4(a, duplo4(l)));
        store4(res.val+4, mul4(a, duphi4(l)));
        store4(res.val+8, mul4(a, duplo4(u)));
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
    
    /// return outer product: [ dir (x) transpose(vec) ]
    static Matrix33 outerProduct(Vector3 const& dir, Vector3 const& vec)
    {
#if MATRIX33_USES_AVX && ( BLD == 4 )
        Matrix33 res;
        vec4 v = vec;
        vec4 p = swap2f128(v);
        vec4 l = blend22(v, p);
        vec4 u = blend22(p, v);
        vec4 d = dir;
        store4(res.val  , mul4(d, duplo4(l)));
        store4(res.val+4, mul4(d, duphi4(l)));
        store4(res.val+8, mul4(d, duplo4(u)));
        return res;
#else
        return Matrix33(dir.XX*vec.XX, dir.YY*vec.XX, dir.ZZ*vec.XX,
                        dir.XX*vec.YY, dir.YY*vec.YY, dir.ZZ*vec.YY,
                        dir.XX*vec.ZZ, dir.YY*vec.ZZ, dir.ZZ*vec.ZZ );
#endif
    }
    
    /// return outer product: [ dir (x) transpose(vec) ]
    static Matrix33 outerProduct(real const* dir, real const* vec)
    {
#if MATRIX33_USES_AVX && ( BLD == 4 )
        Matrix33 res;
        vec4 d = load3(dir);
        store4(res.val  , mul4(d, broadcast1(vec  )));
        store4(res.val+4, mul4(d, broadcast1(vec+1)));
        store4(res.val+8, mul4(d, broadcast1(vec+2)));
        return res;
#else
        return Matrix33(dir[0]*vec[0], dir[1]*vec[0], dir[2]*vec[0],
                        dir[0]*vec[1], dir[1]*vec[1], dir[2]*vec[1],
                        dir[0]*vec[2], dir[1]*vec[2], dir[2]*vec[2] );
#endif
    }
    
    /// add outer product: [ dir (x) transpose(vec) ]
    void addOuterProduct(real const* dir, real const* vec)
    {
#if MATRIX33_USES_AVX && ( BLD == 4 )
        vec4 d = load3(dir);
        store4(val  , fmadd4(d, broadcast1(vec  ), load4(val  )));
        store4(val+4, fmadd4(d, broadcast1(vec+1), load4(val+4)));
        store4(val+8, fmadd4(d, broadcast1(vec+2), load4(val+8)));
#else
        val[0      ] += dir[0]*vec[0];
        val[1      ] += dir[1]*vec[0];
        val[2      ] += dir[2]*vec[0];
        val[0+BLD  ] += dir[0]*vec[1];
        val[1+BLD  ] += dir[1]*vec[1];
        val[2+BLD  ] += dir[2]*vec[1];
        val[0+BLD*2] += dir[0]*vec[2];
        val[1+BLD*2] += dir[1]*vec[2];
        val[2+BLD*2] += dir[2]*vec[2];
#endif
    }

    /// return [ dir (x) transpose(vec) + vec (x) transpose(dir) ]
    static Matrix33 symmetricOuterProduct(Vector3 const& dir, Vector3 const& vec)
    {
        real X = dir.XX * vec.XX;
        real Y = dir.YY * vec.YY;
        real Z = dir.ZZ * vec.ZZ;
        return symmetric(X+X, dir.YY*vec.XX + dir.XX*vec.YY, dir.ZZ*vec.XX + dir.XX*vec.ZZ,
                         Y+Y, dir.ZZ*vec.YY + dir.YY*vec.ZZ,
                         Z+Z);
    }
    
    /// return symmetric matrix block :  -dir^2 * Id + [ dir (x) dir ]
    static Matrix33 offsetOuterProduct(Vector3 const& dir)
    {
        real X = dir.XX;
        real Y = dir.YY;
        real Z = dir.ZZ;
        return symmetric(-Y*Y - Z*Z, Y*X, Z*X,
                         -X*X - Z*Z, Z*Y,
                         -X*X - Y*Y);
    }

    /// return symmetric matrix block :  dia * Id + [ dir (x) dir ] * len
    static Matrix33 offsetOuterProduct(const real dia, Vector3 const& dir, const real len)
    {
        real X = dir.XX * len;
        real Y = dir.YY * len;
        real Z = dir.ZZ * len;
        return symmetric(X * dir.XX + dia, Y * dir.XX, Z * dir.XX,
                         Y * dir.YY + dia, Z * dir.YY,
                         Z * dir.ZZ + dia);
    }
    
    /// build the rotation matrix `M = 2 * dir (x) dir - Id` of axis `dir` and angle 180
    static Matrix33 householder(const Vector3& axis)
    {        
        real X = axis.XX, Y = axis.YY, Z = axis.ZZ;
        real X2 = X + X, Y2 = Y + Y, Z2 = Z + Z;
        
        return symmetric(X * X2 - 1.0, X * Y2, X * Z2,
                         Y * Y2 - 1.0, Y * Z2,
                         Z * Z2 - 1.0);
    }

    
    /// build the matrix `dia * Id + vec (x) Id`
    /**
     Thus applying M to V results in `dia * V + vec (x) V`
     */
    static Matrix33 vectorProduct(const real dia, const Vector3& vec)
    {
        return Matrix33(    dia,  vec.ZZ, -vec.YY,
                        -vec.ZZ,     dia,  vec.XX,
                         vec.YY, -vec.XX,     dia);
    }
    
    /// 2D Rotation around `axis` (of norm 1) with angle set by cosine and sine values
    /**
     This is not a 3D rotation and the images of 'axis' is zero!
     This is equivalent to rotationAroundAxis(axis, c, s) - outerProduct(axis)
     Attention: This is meant to be called with `norm(axis)==1` and `c*c + s*s == 1`
     but the values of 'c' and 's' can be tweaked to scale the resulting matrix.
     */
    static Matrix33 planarRotation(const Vector3& axis, const real c, const real s)
    {
        const real  X = axis.XX,  Y = axis.YY,  Z = axis.ZZ;
        const real cX = -c * X,  cY = -c * Y,  cZ = -c * Z;
        const real sX =  s * X,  sY =  s * Y,  sZ =  s * Z;

        return Matrix33(cX * X + c , cY * X + sZ, cZ * X - sY,
                        cX * Y - sZ, cY * Y + c , cZ * Y + sX,
                        cX * Z + sY, cY * Z - sX, cZ * Z + c );
    }

#pragma mark - Rotations
    
    /// rotation around `axis` (of norm 1) with angle set by cosine and sine values
    /**
     Attention: the result is a rotation only if `norm(axis)==1` and `c*c + s*s == 1`
     but the values of 'c' and 's' can be scaled to obtain a matrix where the
     rotated components are also scaled. Vectors along the axis remain unchanged
     */
    static Matrix33 rotationAroundAxis(const Vector3& axis, const real c, const real s)
    {
        /*
         This is using Rodrigues's formula:
             Id + sine * K + ( 1 - cosine ) * K^2
             K = -Id (x) axis       (-K is the cross product matrix)
             K^2 = axis (x) axis - Id
         
         https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
         */
        const real  X = axis.XX  ,  Y = axis.YY  ,  Z = axis.ZZ;
        const real dX = X - c * X, dY = Y - c * Y, dZ = Z - c * Z;
        const real sX = s * X    , sY = s * Y    , sZ = s * Z;

        return Matrix33(dX * X + c , dY * X + sZ, dZ * X - sY,
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
    static Matrix33 rotationAroundAxisEuler(const real a[3]);

    /// return rotation of angle a, around axis of azimuth b and elevation c
    static Matrix33 rotationFromAngles(const real a[3]);

    
    /// return a rotation that transforms (1,0,0) into (-1,0,0)
    static Matrix33 rotation180();
    
    /// mirror image: X -> -X
    static Matrix33 flipX();

    /// a rotation that brings (1, 1, 1) into (1, 0, 0)
    static Matrix33 align111();

    /// return rotation of axis Z with angle defined by cosine, sine values
    static Matrix33 rotationAroundZ(const real C, const real S)
    {
        return Matrix33(C, S, 0, -S, C, 0, 0, 0, 1);
    }

    /// a rotation around the X axis of specified angle
    static Matrix33 rotationAroundX(real angle);
    
    /// a rotation around the Y axis of specified angle
    static Matrix33 rotationAroundY(real angle);
    
    /// a rotation around the Z axis of specified angle
    static Matrix33 rotationAroundZ(real angle);
    
    /// a rotation around one the axis: X if `i=0`, Y if `i=1` or Z if `i=2`
    static Matrix33 rotationAroundPrincipalAxis(index_t i, real angle);

    /// return a rotation that transforms (1,0,0) into `vec` ( norm(vec) should be > 0 )
    static Matrix33 rotationToVector(const Vector3&);
    
    /// return a matrix that transform (1,0,0) into X, and (0,0,1) into Z
    static Matrix33 rotationToVectors(Vector3 X, Vector3 Z);

    /// return a random rotation that transforms (1,0,0) into `vec` ( norm(vec) should be > 0 )
    static Matrix33 randomRotationToVector(const Vector3&);
    
    /// a random rotation chosen uniformly
    static Matrix33 randomRotation();
    
    /// a rotation of angle 'angle' around an axis chosen randomly
    static Matrix33 randomRotation(real angle);
};


/// output operator
inline std::ostream& operator << (std::ostream& os, Matrix33 const& arg)
{
    std::ios::fmtflags fgs = os.flags();
    arg.print_smart(os);
    os.setf(fgs);
    return os;
}

#endif

