// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University

#ifndef INTERPOLATION4_H
#define INTERPOLATION4_H

#include "vector.h"
#include "iowrapper.h"
#include "mecapoint.h"

class Interpolation;
class Simul;
class Meca;

/// This represents an interpolation over up to 4 vertices of a Mecable
/** FJN 17.9.2018 */
class Interpolation4
{
public:
    
    /// Mecable from which points are interpolated
    Mecable const* mec_;
    
    /// interpolation coefficients for points [ref, ref+1, ref+2, ref+3]
    /** The sum of these 4 coefficients is equal to one */
    real coef_[4];

    /// index of first interpolated point
    index_t prime_;
    
    /// number of interpolated points (order), expected to be DIM+1 at most
    index_t rank_;
    
    /// set coefficients
    void set_coef(real a, real b, real c) { coef_[1]=a; coef_[2]=b; coef_[3]=c; polish(); }

    /// calculate rank by considering which coefficients are not null
    void polish();
    
public:
    
    /// set pointers to null
    void clear();
    
    /// constructor
    Interpolation4() : mec_(nullptr), prime_(0), rank_(0) { set_coef(0, 0, 0); }

    /// set as pointing to vertex `p` of `mec`
    void set(Mecable const* mec, index_t p);
    
    /// set as interpolated between vertices `p` and `q` of `mec`
    void set(Mecable const* mec, index_t p, index_t q, real coef);

    /// set as interpolated over 4 vertices, defined by position 'vec'
    void set(Mecable const*, index_t, Vector const& vec);

    /// attachment mecable
    Mecable const* mecable() const { return mec_; }
    
    /// position in space calculated from interpolation
    Vector pos() const { return mec_->interpolatePoints(prime_, coef_, rank_); }
    
    /// direction normal to the sphere, at the point of interpolation
    Vector normal() const;

    /// number of points beeing interpolated
    index_t rank() const { return rank_; }

    /// first attachement point
    Mecapoint vertex0() const { return Mecapoint(mec_, prime_); }

    /// interpolation coefficients
    const real* coef() const { return coef_; }
    
    /// call Meca::addLink to given Mecapoint
    void addLink(Meca&, Mecapoint const&, real weight) const;

    /// call Meca::addLink with given Interpolation and *this
    void addLink(Meca&, Interpolation const&, real weight) const;
    
    /// check validity
    int  invalid() const;
    
    /// output
    void write(Outputter& out) const;
    
    /// input
    void read(Inputter& in, Simul&);
    
    /// printout
    void print(std::ostream& out) const;
};

#endif
