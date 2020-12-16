// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
// F. Nedelec, Strasbourg 08.06.2018

#ifndef MATRIX22
#define MATRIX22

#include "real.h"
#include "vector2.h"
#include <cstdio>
#include <iostream>

/*
 This matrix can use AVX instructions if 'real == double'
 */
#ifdef __AVX__
#  define MATRIX22_USES_AVX REAL_IS_DOUBLE
#  include "simd.h"
#else
#  define MATRIX22_USES_AVX 0
#endif


/// 2x2 matrix class with 4 'real' elements stored in column order
class alignas(4*sizeof(real)) Matrix22 final
{
public:
    
    /// unsigned integer type used for indices
    typedef size_t index;
    
private:
    
    union {
        real val[4];
#if MATRIX22_USES_AVX
        vec4 mat;
#endif
    };

    /// access to modifiable element by index
    real& operator[](index i)       { return val[i]; }
    
    /// access element value by index
    real  operator[](index i) const { return val[i]; }

public:
    
    Matrix22() {}
    
    /// copy constructor
    Matrix22(Matrix22 const& M)
    {
#if MATRIX22_USES_AVX
        mat = M.mat;
#else
        for ( index u = 0; u < 4; ++u )
            val[u] = M.val[u];
#endif
    }
    
    /// copy constructor with scaling
    Matrix22(Matrix22 const& M, real alpha)
    {
#if MATRIX22_USES_AVX
        mat = mul4(M.mat, set4(alpha));
#else
        for ( index u = 0; u < 4; ++u )
            val[u] = alpha * M.val[u];
#endif
    }

    /// construct Matrix from coordinates (column-major)
    Matrix22(real a, real b, real c, real d)
    {
        val[0] = a;
        val[1] = b;
        val[2] = c;
        val[3] = d;
    }

    /// construct Matrix with `d` on the diagonal and other values equal to `z`
    Matrix22(real z, real d)
    {
        val[0] = d;
        val[1] = z;
        val[2] = z;
        val[3] = d;
    }

    /// constructor from array
    Matrix22(real const M[])
    {
#if MATRIX22_USES_AVX
        mat = loadu4(M);
#else
        val[0] = M[0];
        val[1] = M[1];
        val[2] = M[2];
        val[3] = M[3];
#endif
    }
    
#if MATRIX22_USES_AVX
    /// constructor from SIMD vector
    Matrix22(vec4 const& M)
    {
        mat = M;
    }
#endif

    /// destructor
    ~Matrix22() {}
    
    /// dimensionality
    static constexpr size_t dimension() { return 2; }
    
    /// human-readable identifier
#if MATRIX22_USES_AVX
    static std::string what() { return "+4"; }
#else
    static std::string what() { return "4"; }
#endif
    
    /// set all elements to zero
    void reset()
    {
#if MATRIX22_USES_AVX
        mat = setzero4();
#else
        for ( index u = 0; u < 4; ++u )
            val[u] = 0.;
#endif
    }

    /// set diagonal to 'dia' and other elements to 'off'
    void reset(real dia, real off)
    {
#if MATRIX22_USES_AVX
        mat = set4(dia, off, off, dia);
#else
        val[0] = dia;
        val[1] = off;
        val[2] = off;
        val[3] = dia;
#endif
    }
    
    /// true if any element is different than 'zero'
    bool operator != (real zero) const
    {
        for ( index u = 0; u < 4; ++u )
            if ( val[u] != zero )
                return true;
        return false;
    }

    /// copy values from lower triangle to upper triangle
    void copy_lower()
    {
        val[2] = val[1];
    }
    
    /// conversion to pointer of real
    operator real const*() const { return val; }

    /// return modifiable pointer of 'real'
    real* data() { return val; }

    /// return unmodifiable pointer of real
    real const* data() const { return val; }

    /// return address of element at (i, j)
    real* addr(const index i, const index j) { return val + ( i + 2*j ); }
    
    /// access functions to element by line and column indices
    real& operator()(const index i, const index j)       { return val[i+2*j]; }
    real  operator()(const index i, const index j) const { return val[i+2*j]; }
    
    /// set elements from given array
    void load(const real ptr[])
    {
        for ( index u = 0; u < 4; ++u )
            val[u] = ptr[u];
    }

    /// copy elements to given array
    void store(real ptr[]) const
    {
        for ( index u = 0; u < 4; ++u )
            ptr[u] = val[u];
    }

    /// extract column vector at given index
    Vector2 column(const index i) const
    {
        return Vector2(val+2*i);
    }
    
    /// extract line vector at given index
    Vector2 line(const index i) const
    {
        return Vector2(val[i], val[2+i]);
    }

    /// extract diagonal
    Vector2 diagonal() const
    {
        return Vector2(val[0], val[3]);
    }

    /// sum of diagonal terms
    real trace() const
    {
        return ( val[0] + val[3] );
    }

    /// human-friendly ouput
    void print(FILE * f) const
    {
        fprintf(f, "/ %9.3f %9.3f \\\n", val[0], val[2]);
        fprintf(f, "\\ %9.3f %9.3f /\n", val[1], val[3]);
    }
    
    /// print on one line
    void print(std::ostream& os) const
    {
        const int w = (int)os.width();
        os << std::setw(2) << "[ ";
        os << std::setw(w) << (*this)(0,0) << " ";
        os << std::setw(w) << (*this)(0,1) << "; ";
        os << std::setw(w) << (*this)(1,0) << " ";
        os << std::setw(w) << (*this)(1,1) << " ]";
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

    /// relative asymmetry = abs( difference of off-diagonal elements ) / trace
    real asymmetry() const
    {
        real t = abs_real(val[0]) + abs_real(val[3]);
        return abs_real(val[2]-val[1]) / t;
    }
    
    /// scale all elements
    void scale(const real alpha)
    {
#if MATRIX22_USES_AVX
        mat = mul4(mat, set4(alpha));
#else
        for ( index u = 0; u < 4; ++u )
            val[u] *= alpha;
#endif
    }
    
    /// scale matrix
    void operator *=(const real alpha)
    {
        scale(alpha);
    }
    
    /// return opposite matrix (i.e. -M)
    const Matrix22 operator -() const
    {
        Matrix22 M;
        for ( index u = 0; u < 4; ++u )
            M.val[u] = -val[u];
        return M;
    }

    /// scaled matrix
    const Matrix22 operator *(const real alpha) const
    {
#if MATRIX22_USES_AVX
        return Matrix22(mul4(mat, set4(alpha)));
#else
        Matrix22 M;
        for ( index u = 0; u < 4; ++u )
            M.val[u] = val[u] * alpha;
        return M;
#endif
    }
    
    /// multiplication by scalar
    friend const Matrix22 operator *(const real alpha, Matrix22 const& mat)
    {
        return mat * alpha;
    }
    
    /// return sum of two matrices
    const Matrix22 operator +(Matrix22 const& M) const
    {
#if MATRIX22_USES_AVX
        return Matrix22(add4(mat, M.mat));
#else
        Matrix22 res;
        for ( index u = 0; u < 4; ++u )
            res.val[u] = val[u] + M.val[u];
        return res;
#endif
    }

    /// return difference of two matrices
    const Matrix22 operator -(Matrix22 const& M) const
    {
#if MATRIX22_USES_AVX
        return Matrix22(sub4(mat, M.mat));
#else
        Matrix22 res;
        for ( index u = 0; u < 4; ++u )
            res.val[u] = val[u] - M.val[u];
        return res;
#endif
    }

    /// subtract given matrix
    void operator +=(Matrix22 const& M)
    {
#if MATRIX22_USES_AVX
        mat = add4(mat, M.mat);
#else
        for ( index u = 0; u < 4; ++u )
            val[u] += M.val[u];
#endif
    }

    /// add given matrix
    void operator -=(Matrix22 const& M)
    {
#if MATRIX22_USES_AVX
        mat = sub4(mat, M.mat);
#else
        for ( index u = 0; u < 4; ++u )
            val[u] -= M.val[u];
#endif
    }

    /// transpose matrix in place
    void transpose()
    {
        real t = val[1];
        val[1] = val[2];
        val[2] = t;
    }
    
    /// return transposed matrix
    Matrix22 transposed() const
    {
#if MATRIX22_USES_AVX && defined __AVX2__
        return Matrix22(permute4x64(mat, 0xD8));
#else
        return Matrix22(val[0], val[2], val[1], val[3]);
#endif
    }
       
    /// return scaled transposed matrix
    Matrix22 transposed(real alpha) const
    {
#if MATRIX22_USES_AVX && defined __AVX2__
        return Matrix22(mul4(set4(alpha), permute4x64(mat, 0xD8)));
#else
        return Matrix22(alpha*val[0], alpha*val[2], alpha*val[1], alpha*val[3]);
#endif
    }
 
    /// maximum of all component's absolute values
    real norm_inf() const
    {
        real a = std::max(abs_real(val[0]), abs_real(val[1]));
        real b = std::max(abs_real(val[2]), abs_real(val[3]));
        return std::max(a, b);
    }

    /// determinant
    real determinant() const
    {
        return ( val[0] * val[3] - val[1] * val[2] );
    }
    
    /// inverse in place
    void inverse()
    {
        real s = 1.0 / determinant();
        real x =  s * val[0];
        val[0] =  s * val[3];
        val[1] = -s * val[1];
        val[2] = -s * val[2];
        val[3] = x;
    }

    /// return inverse matrix
    Matrix22 inverted() const
    {
        real s = 1.0 / determinant();
        return Matrix22(s * val[3], -s * val[1], -s * val[2], s * val[0]);
    }

#if MATRIX22_USES_AVX
    static const vec4 transposed(vec4 const& mat)
    {
#ifdef __AVX2__
        return permute4x64(mat, 0xD8);
#else
        return blend4(mat, permute4(swap2f128(mat),0b1100), 0b0110);
#endif
    }

    /// multiplication by a vector: this * V
    static const vec2 vecmul2(vec4 const& mat, vec2 const& vec)
    {
#if 0
        vec2 res;
        res[0] = mat[0] * vec[0] + mat[2] * vec[1];
        res[1] = mat[1] * vec[0] + mat[3] * vec[1];
        return res;
#elif defined(__FMA__)
        vec2 u = mul2(getlo(mat), unpacklo2(vec,vec));
        return fmadd2(gethi(mat), unpackhi2(vec,vec), u);
#else
        vec2 u = mul2(getlo(mat), unpacklo2(vec,vec));
        vec2 w = mul2(gethi(mat), unpackhi2(vec,vec));
        return add2(u, w);
#endif
    }
    
    /// multiplication by a vector: this * V
    static const vec2 vecmul2(vec4 const& mat, real const* V)
    {
#if 0
        vec2 res;
        res[0] = mat[0] * V[0] + mat[2] * V[1];
        res[1] = mat[1] * V[0] + mat[3] * V[1];
        return res;
#elif defined(__FMA__)
        vec2 u = mul2(getlo(mat), loaddup2(V));
        return fmadd2(gethi(mat), loaddup2(V+1), u);
#else
        vec2 u = mul2(getlo(mat), loaddup2(V));
        vec2 w = mul2(gethi(mat), loaddup2(V+1));
        return add2(u, w);
#endif
    }

    /// multiplication by a vector: transpose(M) * V
    static const vec2 trans_vecmul2(vec4 const& mat, real const* V)
    {
#if 0
        vec2 res;
        res[0] = mat[0] * V[0] + mat[1] * V[1];
        res[1] = mat[2] * V[0] + mat[3] * V[1];
        return res;
#elif defined __AVX2__
        vec4 m = mul4(mat, broadcast2(V));
        vec4 l = permute4x64(m, 0x88);
        vec4 h = permute4x64(m, 0xDD);
        return add2(getlo(l), getlo(h));
#else
        vec4 m = mul4(mat, broadcast2(V));
        vec2 h = gethi(m);
        return add2(unpacklo2(getlo(m), h), unpackhi2(getlo(m), h));
#endif
    }
    
    /// multiplication by another matrix: @returns val * mat
    static const vec4 mul(vec4 const& val, vec4 const& mat)
    {
#ifdef __FMA__
        vec4 s = mul4(permute2f128(val,val,0x20), duplo4(mat));
        return fmadd4(permute2f128(val,val,0x31), duphi4(mat), s);
#else
        vec4 s = mul4(permute2f128(val,val,0x20), duplo4(mat));
        vec4 t = mul4(permute2f128(val,val,0x31), duphi4(mat));
        return add4(s, t);
#endif
    }
    
    /// multiplication by another matrix: @returns val * mat
    static const vec4 mul(real const* val, vec4 const& mat)
    {
        vec4 s = mul4(broadcast2(val), duplo4(mat));
        return fmadd4(broadcast2(val+2), duphi4(mat), s);
    }

    /// multiplication by another matrix: @returns transpose(this) * mat
    static const vec4 trans_mul(vec4 const& val, vec4 const& mat)
    {
#ifdef __FMA__
        vec4 s = mul4(permute4x64(val, 0x88), duplo4(mat));
        return fmadd4(permute4x64(val, 0xDD), duphi4(mat), s);
#elif defined __AVX2__
        vec4 s = mul4(permute4x64(val, 0x88), duplo4(mat));
        vec4 t = mul4(permute4x64(val, 0xDD), duphi4(mat));
        return add4(s, t);
#else
        vec4 a = transposed(val);
        vec4 s = mul4(permute2f128(a,a,0x20), duplo4(mat));
        return fmadd4(permute2f128(a,a,0x31), duphi4(mat), s);
#endif
    }

#endif

    /// multiplication by a vector: this * V
    Vector2 vecmul0(Vector2 const& V) const
    {
        return Vector2(val[0] * V.XX + val[2] * V.YY,
                       val[1] * V.XX + val[3] * V.YY);
    }

    /// multiplication by a vector: this * V
    Vector2 vecmul0(real const* ptr) const
    {
        return Vector2(val[0] * ptr[0] + val[2] * ptr[1],
                       val[1] * ptr[0] + val[3] * ptr[1]);
    }

    /// multiplication by a vector: this * V
    Vector2 vecmul(Vector2 const& V) const
    {
#if MATRIX22_USES_AVX
        return Vector2(vecmul2(mat, V.vec));
#else
        return vecmul0(V);
#endif
    }

    /// multiplication by a vector: this * { ptr[0], ptr[1] }
    Vector2 vecmul(real const* ptr) const
    {
#if MATRIX22_USES_AVX
        return Vector2(vecmul2(mat, ptr));
#else
        return vecmul0(ptr);
#endif
    }
    
    friend Vector2 operator * (Matrix22 const& mat, Vector2 const& vec)
    {
        return mat.vecmul(vec);
    }

    /// multiplication by a vector: this * V
    Vector2 trans_vecmul0(Vector2 const& V) const
    {
        return Vector2(val[0] * V.XX + val[1] * V.YY,
                       val[2] * V.XX + val[3] * V.YY);
    }

    /// multiplication by a vector: transpose(M) * V
    Vector2 trans_vecmul0(real const* ptr) const
    {
        return Vector2(val[0] * ptr[0] + val[1] * ptr[1],
                       val[2] * ptr[0] + val[3] * ptr[1]);
    }

    /// multiplication by a vector: transpose(M) * V
    Vector2 trans_vecmul(real const* ptr) const
    {
#if MATRIX22_USES_AVX
        return Vector2(trans_vecmul2(mat, ptr));
#else
        return trans_vecmul0(ptr);
#endif
    }

    /// multiplication by another matrix: @returns this * M
    const Matrix22 mul(Matrix22 const& M) const
    {
#if MATRIX22_USES_AVX
        return mul(val, M.mat);
#else
        Matrix22 res;
        res.val[0] = val[0] * M[0] + val[2] * M[1];
        res.val[1] = val[1] * M[0] + val[3] * M[1];
        res.val[2] = val[0] * M[2] + val[2] * M[3];
        res.val[3] = val[1] * M[2] + val[3] * M[3];
        return res;
#endif
    }
    
    /// matrix-matrix multiplication
    friend Matrix22 operator * (Matrix22 const& mat, Matrix22 const& mut)
    {
        return mat.mul(mut);
    }

    /// multiplication by another matrix: @returns transpose(this) * M
    const Matrix22 trans_mul(Matrix22 const& M) const
    {
#if MATRIX22_USES_AVX
        return trans_mul(mat, M.mat);
#else
        Matrix22 res;
        res.val[0] = val[0] * M[0] + val[1] * M[1];
        res.val[1] = val[2] * M[0] + val[3] * M[1];
        res.val[2] = val[0] * M[2] + val[1] * M[3];
        res.val[3] = val[2] * M[2] + val[3] * M[3];
        return res;
#endif
    }

    /// add full matrix: this <- this + M
    void add_full(Matrix22 const& M)
    {
#if MATRIX22_USES_AVX
        mat = add4(mat, M.mat);
#else
        for ( index u = 0; u < 4; ++u )
            val[u] += M.val[u];
#endif
    }

    /// add full matrix: this <- this + alpha * M
    void add_full(const real alpha, Matrix22 const& M)
    {
#if MATRIX22_USES_AVX
        mat = fmadd4(set4(alpha), M.mat, mat);
#else
        for ( index u = 0; u < 4; ++u )
            val[u] += alpha * M.val[u];
#endif
    }
    
    /// sub full matrix: this <- this - M
    void sub_full(Matrix22 const& M)
    {
#if MATRIX22_USES_AVX
        mat = sub4(mat, M.mat);
#else
        for ( index u = 0; u < 4; ++u )
            val[u] -= M.val[u];
#endif
    }

    /// subtract full matrix: this <- this - alpha * M
    void sub_full(const real alpha, Matrix22 const& M)
    {
#if MATRIX22_USES_AVX
        mat = fnmadd4(set4(alpha), M.mat, mat);
#else
        for ( index u = 0; u < 4; ++u )
            val[u] -= alpha * M.val[u];
#endif
    }

    /// subtract transposed matrix: this <- this - transposed(M)
    void sub_trans(Matrix22 const& M)
    {
        val[0] -= M.val[0];
        val[1] -= M.val[2];
        val[2] -= M.val[1];
        val[3] -= M.val[3];
    }

    /// add transposed matrix: this <- this + alpha * transposed(M)
    void add_trans(Matrix22 const& M)
    {
        val[0] += M.val[0];
        val[1] += M.val[2];
        val[2] += M.val[1];
        val[3] += M.val[3];
    }

    /// add transposed matrix: this <- this + alpha * transposed(M)
    void add_trans(const real alpha, Matrix22 const& M)
    {
        val[0] += alpha * M.val[0];
        val[1] += alpha * M.val[2];
        val[2] += alpha * M.val[1];
        val[3] += alpha * M.val[3];
    }

    /// add lower triangle of matrix including diagonal: this <- this + M
   void add_half(Matrix22 const& M)
    {
#if MATRIX22_USES_AVX
        mat = add4(mat, M.mat);
#elif ( 1 )
        for ( index u = 0; u < 4; ++u )
            val[u] += M.val[u];
#else
        val[0] += M.val[0];
        val[1] += M.val[1];
        val[3] += M.val[3];
#endif
    }
    
    /// add lower triangle of matrix including diagonal: this <- this + alpha * M
    void add_half(const real alpha, Matrix22 const& M)
    {
#if MATRIX22_USES_AVX
        mat = fmadd4(set4(alpha), M.mat, mat);
#elif ( 1 )
        for ( index u = 0; u < 4; ++u )
            val[u] += alpha * M.val[u];
#else
        val[0] += alpha * M.val[0];
        val[1] += alpha * M.val[1];
        val[3] += alpha * M.val[3];
#endif
    }
    
    /// subtract lower triangle of matrix including diagonal: this <- this - M
    void sub_half(Matrix22 const& M)
    {
#if MATRIX22_USES_AVX
        mat = sub4(mat, M.mat);
#elif ( 1 )
        for ( index u = 0; u < 4; ++u )
            val[u] -= M.val[u];
#else
        val[0] -= M.val[0];
        val[1] -= M.val[1];
        val[3] -= M.val[3];
#endif
    }
    
    /// add alpha to diagonal
    void add_diag(real alpha)
    {
        val[0] += alpha;
        val[3] += alpha;
    }
    
    /// add -alpha to diagonal
    void sub_diag(real alpha)
    {
        val[0] -= alpha;
        val[3] -= alpha;
    }

    
    /// add all elements of block 'S' to array 'M'
    void addto(real * M, size_t ldd) const
    {
        M[0    ] += val[0];
        M[1    ] += val[1];
        M[  ldd] += val[2];
        M[1+ldd] += val[3];
    }
    
    /// add lower elements of this block to upper triangle of 'M'
    void addto_lower(real * M, size_t ldd) const
    {
        M[0    ] += val[0];
        M[1    ] += val[1];
        M[1+ldd] += val[3];
        assert_true( val[2] == 0 );
    }
    
    /// add lower elements of this block to both upper and lower triangles of 'M'
    void addto_symm(real * M, size_t ldd) const
    {
        M[0    ] += val[0];
        M[1    ] += val[1];
        M[  ldd] += val[1];
        M[1+ldd] += val[3];
        //assert_true( val[2] == 0 );
    }
    
    /// add all elements of this block to 'M', with transposition
    void addto_trans(real * M, size_t ldd) const
    {
        M[0    ] += val[0];
        M[1    ] += val[2];
        M[  ldd] += val[1];
        M[1+ldd] += val[3];
    }
    
#pragma mark -

    /// return diagonal Matrix from diagonal terms
    static Matrix22 diagonal(real a, real b)
    {
        return Matrix22(a, 0, 0, b);
    }

    /// identity matrix
    static Matrix22 identity()
    {
        return Matrix22(0, 1);
    }


    /// set a symmetric matrix: [ dir (x) transpose(dir) ]
    static Matrix22 outerProduct(Vector2 const& dir)
    {
        
        real xy = dir.XX * dir.YY;
        return Matrix22(dir.XX * dir.XX, xy,
                        xy, dir.YY * dir.YY);
    }

    /// set a symmetric matrix: alpha * [ dir (x) transpose(dir) ]
    static Matrix22 outerProduct(Vector2 const& dir, real alpha)
    {
        real XX = dir.XX * dir.XX;
        real XY = dir.XX * dir.YY;
        real YY = dir.YY * dir.YY;
        return Matrix22(XX, XY, XY, YY) * alpha;
    }
    
    /// return outer product: [ dir (x) transpose(vec) ]
    static Matrix22 outerProduct(Vector2 const& dir, Vector2 const& vec)
    {
        return Matrix22(dir.XX * vec.XX, dir.YY * vec.XX,
                        dir.XX * vec.YY, dir.YY * vec.YY);
    }
    
    /// return [ dir (x) transpose(vec) + vec (x) transpose(dir) ]
    static Matrix22 symmetricOuterProduct(Vector2 const& dir, Vector2 const& vec)
    {
        real xx = dir.XX * vec.XX;
        real yy = dir.YY * vec.YY;
        real xy = dir.YY * vec.XX + dir.XX * vec.YY;
        return Matrix22(xx+xx, xy, xy, yy+yy);
    }

    /// return symmetric matrix block :  dia * I + [ dir (x) dir ] * len
    static Matrix22 offsetOuterProduct(const real dia, Vector2 const& dir, const real len)
    {
        real xl = dir.XX * len;
        real yl = dir.YY * len;
        real xy = dir.XX * yl;
        return Matrix22(xl * dir.XX + dia, xy, xy, yl * dir.YY + dia);
    }
    
    /// build the matrix `dia * Id + vec (x) Id`
    static Matrix22 vectorProduct(const real dia, const real off)
    {
        return Matrix22(dia, off, -off, dia);
    }

    /// return rotation matrix of angle defined by cosinus and sinus
    static Matrix22 rotation(const real c, const real s)
    {
        return Matrix22(c, s, -s, c);
    }
    
    /// return rotation matrix of given angle
    static Matrix22 rotation(const real ang)
    {
        real c = std::cos(ang), s = std::sin(ang);
        return Matrix22(c, s, -s, c);
    }

    /// angle of rotation
    real rotationAngle() const;
    
    /// return a rotation that transforms (1,0,0) into (-1,0,0)
    static Matrix22 rotation180();

    /// return a rotation that transforms (1,0,0) into `vec` ( norm(vec) should be > 0 )
    static Matrix22 rotationToVector(const Vector2& vec);
    
    /// return a random rotation that transforms (1,0,0) into `vec` ( norm(vec) should be > 0 )
    static Matrix22 randomRotationToVector(const Vector2& vec)
    {
        return rotationToVector(vec);
    }
    
    /// a random rotation chosen uniformly
    static Matrix22 randomRotation();
    
    /// a rotation of angle '+/- angle'
    static Matrix22 randomRotation(real angle);

};

/// output operator
inline std::ostream& operator << (std::ostream& os, Matrix22 const& arg)
{
    std::ios::fmtflags fgs = os.flags();
    arg.print(os);
    os.setf(fgs);
    return os;
}

#endif

