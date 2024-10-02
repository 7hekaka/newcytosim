// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University.
// FJN, Strasbourg 08.06.2018

#ifndef MATRIX11
#define MATRIX11

#include "real.h"
#include "vector1.h"
#include <cstdio>
#include <iostream>

/// 1x1 matrix class with 1 'real' elements
/**
 This matrix can use AVX instructions if 'real == double'
 */
class Matrix11 final
{
private:
    
    real val_;
    
    /// access to modifiable element by index
    real& operator[](index_t i)       { return val_; }
    
    /// access element value by index
    real  operator[](index_t i) const { return val_; }

public:
    
    Matrix11() {}
    
    /// construct Matrix with value `v`
    Matrix11(real v)
    {
        val_ = v;
    }

    /// construct Matrix with term set to `d`
    Matrix11(real, real d)
    {
        val_ = d;
    }

    ~Matrix11() {}
    
#pragma mark -
    
    /// dimensionality
    static constexpr size_t dimension() { return 1; }
    
    /// human-readable identifier
    static std::string what() { return "1"; }
    
    /// set all elements to zero
    void reset()
    {
        val_ = 0.;
    }
    
    /// set diagonal to 'dia' and other elements to 'off'
    void reset(real, real dia)
    {
        val_ = dia;
    }

    /// true if element is different from 'zero'
    bool operator != (real zero) const
    {
        return ( val_ != zero );
    }
    
    /// conversion to real
    //operator real() const { return val_; }

    /// copy values from lower triangle to upper triangle
    void copy_lower() { }

    /// direct access to the scalar
    real& value()       { return val_; }

    /// direct access to the scalar
    real  value() const { return val_; }
    
    /// modifiable pointer of 'real'
    real* data() { return &val_; }
    
    /// unmodifiable pointer of 'real'
    real const* data() const { return &val_; }

    /// address of element at line i, column j
    real* addr(const index_t i, const index_t j) { return &val_; }
    /// value of element at line i, column j
    real value(const index_t i, const index_t j) const { return val_; }

    /// access functions to element by line and column indices
    real& operator()(const index_t i, const index_t j)       { return val_; }
    real  operator()(const index_t i, const index_t j) const { return val_; }
    
    /// set elements from given array
    void load(const real ptr[])
    {
        val_ = ptr[0];
    }

    /// copy elements to given array
    void store(real ptr[]) const
    {
        ptr[0] = val_;
    }

    /// extract column vector at given index
    Vector1 column(const index_t) const
    {
        return Vector1(val_);
    }
    
    /// extract line vector at given index
    Vector1 line(const index_t) const
    {
        return Vector1(val_);
    }

    /// extract diagonal
    Vector1 diagonal() const
    {
        return Vector1(val_);
    }

    /// sum of diagonal terms
    real trace() const
    {
        return val_;
    }

    /// human-friendly output
    void print(FILE * f) const
    {
        fprintf(f, "[ %9.3f ]\n", val_);
    }

    /// conversion to string
    std::string to_string(int w, std::streamsize p) const
    {
        std::ostringstream os;
        os.precision(p);
        os << "[ " << std::setw(w) << value() << " ]";
        return os.str();
    }

    /// always zero
    real asymmetry() const
    {
        return 0;
    }
    
    /// scale all elements
    void scale(const real alpha)
    {
        val_ *= alpha;
    }
    
    /// scale matrix
    void operator *=(const real alpha)
    {
        scale(alpha);
    }
    
    /// return opposite matrix (i.e. -M)
    const Matrix11 operator -() const
    {
        return Matrix11(-val_);
    }
    
    /// scaled matrix
    const Matrix11 operator *(const real alpha) const
    {
        return Matrix11( val_ * alpha );
    }

    /// return sum of two matrices
    const Matrix11 operator +(Matrix11 const& M) const
    {
        return Matrix11( val_ + M.val_ );
    }

    /// return difference of two matrices
    const Matrix11 operator -(Matrix11 const& M) const
    {
        return Matrix11( val_ - M.val_ );
    }

    /// subtract given matrix
    void operator +=(Matrix11 const& M)
    {
        val_ += M.val_;
    }

    /// add given matrix
    void operator -=(Matrix11 const& M)
    {
        val_ -= M.val_;
    }
    
    /// transpose matrix in place
    void transpose()
    {
    }
    
    /// return transposed matrix
    Matrix11 transposed() const
    {
        return Matrix11(val_);
    }
        
    /// return scaled transposed matrix
    Matrix11 transposed(real alpha) const
    {
        return Matrix11(alpha*val_);
    }

    /// maximum of all component's absolute values
    real norm_inf() const
    {
        return abs_real(val_);
    }
    
    /// inverse in place
    void inverse()
    {
        val_ = 1.0 / val_;
    }

#pragma mark -

    /// multiplication by a vector: this * V
    Vector1 vecmul(Vector1 const& V) const
    {
        return Vector1(val_ * V.XX);
    }
    
    /// multiplication by a vector: this * V
    Vector1 vecmul(real const* ptr) const
    {
        return Vector1(val_ * ptr[0]);
    }

    /// matrix-vector multiplication
    friend Vector1 operator * (Matrix11 const& mat, Vector1 const& vec)
    {
        return mat.vecmul(vec);
    }

    /// multiplication by a vector: transpose(M) * V
    Vector1 trans_vecmul(real const* V) const
    {
        return Vector1(val_ * V[0]);
    }

    /// multiplication by another matrix: @returns this * B
    const Matrix11 mul(Matrix11 const& B) const
    {
        return Matrix11(val_ * B.val_);
    }
    
    /// matrix-matrix multiplication
    friend Matrix11 operator * (Matrix11 const& mat, Matrix11 const& mut)
    {
        return mat.mul(mut);
    }

    /// multiplication by another matrix: @returns transpose(this) * B
    const Matrix11 trans_mul(Matrix11 const& B) const
    {
        return Matrix11(val_ * B.val_);
    }

    /// add full matrix: this <- this + M
    void add_full(Matrix11 const& M)
    {
        val_ += M.val_;
    }
    
    /// add full matrix: this <- this + alpha * M
    void add_full(const real alpha, Matrix11 const& M)
    {
        val_ += alpha * M.val_;
    }
    
    /// subtract full matrix: this <- this + M
    void sub_full(Matrix11 const& M)
    {
        val_ -= M.val_;
    }
    
    /// subtract full matrix: this <- this - alpha * M
    void sub_full(const real alpha, Matrix11 const& M)
    {
        val_ -= alpha * M.val_;
    }

    /// add lower triangle of matrix including diagonal: this <- this + M
    void add_half(Matrix11 const& M)
    {
        val_ += M.val_;
    }
    
    /// add lower triangle of matrix including diagonal: this <- this + alpha * M
    void add_half(const real alpha, Matrix11 const& M)
    {
        val_ += alpha * M.val_;
    }
    
    /// add lower triangle of matrix including diagonal: this <- this + alpha * M
    void add_half(const real alpha, Matrix11 const& M, const real dia)
    {
        val_ += alpha * ( M.val_ + dia );
    }

    /// subtract lower triangle of matrix including diagonal: this <- this - M
    void sub_half(Matrix11 const& M)
    {
        val_ -= M.val_;
    }
    
    /// add alpha to diagonal
    void add_diag(real alpha)
    {
        val_ += alpha;
    }
    
    /// return copy of *this, with `alpha` added to the diagonal
    Matrix11 plus_diagonal(real alpha) const
    {
        return Matrix11(val_+alpha);
    }
    
    /// add all elements of block 'S' to array 'M'
    void addto(real * M, size_t ldd) const
    {
        M[0] += val_;
    }
    
    /// add scaled elements of block 'S' to array 'M'
    void addto(real * M, size_t ldd, real alpha) const
    {
        M[0] += alpha * val_;
    }

    /// add element of this block to 'M'
    void addto_lower(real * M, size_t ldd) const
    {
        M[0] += val_;
    }
    
    /// add scaled element of this block to 'M'
    void addto_lower(real * M, size_t ldd, real alpha) const
    {
        M[0] += alpha * val_;
    }

    /// add element of this block to 'M'
    void addto_symm(real * M, size_t ldd) const
    {
        M[0] += val_;
    }
    
    /// add all elements of this block to 'M', with transposition
    void addto_trans(real * M, size_t ldd) const
    {
        M[0] += val_;
    }
    
#pragma mark -

    /// return `a * Identity`
    static Matrix11 diagonal(real a)
    {
        return Matrix11(a);
    }
    
    /// identity matrix
    static Matrix11 one()
    {
        return Matrix11(1);
    }

    /// set a symmetric matrix: [ dir (x) transpose(dir) ]
    static Matrix11 outerProduct(Vector1 const& dir)
    {
        return Matrix11(dir.XX * dir.XX);
    }

    /// set a symmetric matrix: alpha * [ dir (x) transpose(dir) ]
    static Matrix11 outerProduct(Vector1 const& dir, real alpha)
    {
        return Matrix11(dir.XX * dir.XX * alpha);
    }
    
    /// return outer product: [ dir (x) transpose(vec) ]
    static Matrix11 outerProduct(Vector1 const& dir, Vector1 const& vec)
    {
        return Matrix11(dir.XX * vec.XX);
    }
    
    /// return null matrix
    static Matrix11 offsetOuterProduct(Vector1 const& dir)
    {
        return Matrix11(0);
    }

    /// return symmetric matrix block :  dia * I + [ dir (x) dir ] * len
    static Matrix11 offsetOuterProduct(const real dia, Vector1 const& dir, const real len)
    {
        return Matrix11(len * dir.XX * dir.XX + dia);
    }
    
    /// build the matrix `dia * Id`
    static Matrix11 vectorProduct(const real dia, const real)
    {
        return Matrix11(dia);
    }
    
    /// meaningless function
    static Matrix11 planarRotation(const real axis, const real c, const real s)
    {
        return Matrix11(sign_real(axis));
    }

    /// return rotation matrix of angle defined by c = +/- 1
    static Matrix11 rotation(const real c, const real s)
    {
        return Matrix11(sign_real(c));
    }
    
    /// return rotation matrix of given angle
    static Matrix11 rotation(const real ang)
    {
        return Matrix11(sign_real(std::cos(ang)));
    }

    /// rotation angle
    real rotationAngle() const;

    /// return a rotation that transforms (1,0,0) into (-1,0,0)
    static Matrix11 rotation180();
    
    /// mirror image: X -> -X
    static Matrix11 flipX();

    /// the identity
    static Matrix11 align111() { return Matrix11(1); }

    /// return a rotation that transforms (1,0,0) into `vec` ( norm(vec) should be > 0 )
    static Matrix11 rotationToVector(const Vector1& vec);
    
    /// return a random rotation that transforms (1,0,0) into `vec` ( norm(vec) should be > 0 )
    static Matrix11 randomRotationToVector(const Vector1& vec)
    {
        return rotationToVector(vec);
    }
    
    /// a random rotation chosen uniformly
    static Matrix11 randomRotation();
    
    /// a random rotation of given angle
    static Matrix11 randomRotation(real);

};


/// output operator
inline std::ostream& operator << (std::ostream& os, Matrix11 const& mat)
{
    os << "[ " << mat.value() << " ]";
    return os;
}

#endif

