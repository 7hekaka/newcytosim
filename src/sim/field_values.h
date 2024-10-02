// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef FIELD_VALUES_H
#define FIELD_VALUES_H

#include "real.h"
#include "iowrapper.h"
#include <iostream>


/// Scalar type (real) for a Field
/**
 Example of type that can be used in Field
 */
class FieldScalar
{
public:
    
    /// single scalar value
    real val;
    
public:
    
    /// constructor
    FieldScalar()                     { val = 0; }
    
    /// implicit conversion from real
    FieldScalar(real a)               { val = a; }
    
    /// implicit conversion to real
    operator real&()                  { return val; }

    /// set to zero
    void clear()                      { val = 0; }
    
    /// comparison
    bool negative() const { return val < 0; }

    /// write
    void write(Outputter& out) const  { out.writeFloat(val); }
    
    /// read
    void read(Inputter& in)           { val = in.readFloat(); }

};


/// Vector of N scalar, suitable for a Field
/**
 Example of type that can be used in Field
 */
template < int N >
class FieldVector
{
    /// vector of values
    real val[N];
    
public:
    
    /// the dimensionality (given as template argument)
    static const int dimension_ = N;
    
    /// constructor
    FieldVector<N>() { for (int i=0; i<N; ++i) val[i] = 0; }

    /// constructor
    FieldVector<N>(real x) { for (int i=0; i<N; ++i) val[i] = x; }

    /// set to zero
    void clear() { for (int i=0; i<N; ++i) val[i] = 0; }

    /// access to the vector components
    real const& operator[](int i) const { assert_true(i<N); return val[i]; }

    /// access to modifiable vector components
    real& operator[](int i)        { assert_true(i<N); return val[i]; }

    /// multiplication
    FieldVector<N> operator * (real x) const { FieldVector R; for (int i=0; i<N; ++i) R[i] = val[i] * x; return R; }
    
    /// addition
    FieldVector<N> operator + (FieldVector<N> const& X) { FieldVector R; for (int i=0; i<N; ++i) R[i] = val[i] + X[i]; return R; }
    
    /// substractions
    FieldVector<N> operator - (FieldVector<N> const& X) { FieldVector R; for (int i=0; i<N; ++i) R[i] = val[i] - X[i]; return R; }

    /// comparison
    bool negative() const { for (int i=1; i<N; ++i) if ( val[i] < 0 ) return true; return false; }

    /// assignment
    void operator = (real x) { for (int i=0; i<N; ++i) val[i] = x; }
    
    /// accumulation
    void operator += (FieldVector<N> const& X) { for (int i=0; i<N; ++i) val[i] += X[i]; }
    

    /// minimization
    void e_min(FieldVector<N> const& X) { for (int i=0; i<N; ++i) val[i] = std::min(val[i], X[i]); }
    
    /// maximization
    void e_max(FieldVector<N> const& X) { for (int i=0; i<N; ++i) val[i] = std::max(val[i], X[i]); }

    /// read
    void read(Inputter& in)          { for (int i=0; i<N; ++i) val[i] = in.readFloat(); }
    
    /// write
    void write(Outputter& out) const { for (int i=0; i<N; ++i) out.writeFloat(val[i]); }

};

/// output
template < int N >
std::ostream& operator << (std::ostream& os, FieldVector<N> const& FV) { for (size_t i=0; i<N; ++i) os << " " << FV[i]; return os; }

/// input
template < int N >
std::istream& operator >> (std::istream& is, FieldVector<N>& FV) { for (size_t i=0; i<N; ++i) is >> FV[i]; return is; }


#endif
