// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University

#ifndef INTERPOLATION_H
#define INTERPOLATION_H

#include "real.h"
#include "vector.h"
#include "mecable.h"
#include "mecapoint.h"
#include "chain.h"

class FiberSegment;

/// Indicates an intermediate position between two points of a Mecable
/**
 This class defines an interpolation between two points on the same Mecable.
 The Mecable is designated by a pointer, and two vertices by their indices.
 An interpolation coefficient in [0,1] then defines an intermediate position
 between the two vertices:
 
     pos = (1-coef) * pt1_ + coef * pt2_
 
 or equivalently:
 
     pos = pt1_ + coef * ( pt2_ - pt1_ )

 This class used to be called 'PointInterpolated' before 12.2017
 */
class Interpolation final
{
private:
    
    /// Mecable from which the points are interpolated 
    Mecable const * mec_;
    
    /// interpolation coefficient: pos = (1-coef) * pt1_ + coef * pt2_
    real coef_;

    /// 1.0 / ( distance between points )
    real ilen_;
    
    /// index of vertex 1 in mec_
    index_t pt1_;

    /// index of vertex 2 in mec_
    index_t pt2_;
    
public:
    
    /// reset member variables
    Interpolation() : mec_(nullptr), coef_(0), ilen_(0), pt1_(0), pt2_(0) { }
    
    /// set to interpolate P and P+1 on `m`, with coefficient `c`
    Interpolation(const Chain * m, real c, index_t P)
    : mec_(m), coef_(c), ilen_(m->segmentationInv()), pt1_(P), pt2_(P+1) { assert_true(P+1 < m->nbPoints() ); }
    
    /// set to interpolate P and Q on `m`, with coefficient `c`
    Interpolation(const Mecable * m, real c, index_t P, index_t Q)
    : mec_(m), coef_(c), ilen_(0), pt1_(P), pt2_(Q) { assert_true(Q < m->nbPoints() ); }

    /// disabled old-style constructor
    //Interpolation(const Mecable*, unsigned, unsigned, real) = delete;

    /// disabled constructor avoid conversions from int to real
    //Interpolation(const Mecable*, unsigned, real) = delete;

    /// set to interpolate given fiber segment, at abscissa `abs` in [0,1]
    Interpolation(FiberSegment const&, real abs);
    
    
    /// Reset member variables
    void clear()
    {
        mec_ = nullptr; pt1_ = 0; pt2_ = 0; coef_ = 0;
    }
    
    /// Set to interpolate p1 and p2 on ps, with coefficient c, on the same Mecable
    void setPoints(unsigned p, unsigned q, const real c)
    {
        pt1_ = p; pt2_ = q; coef_ = c;
    }
    
    /// Index of point 1 in the matrix of dynamics (Meca::mISO)
    index_t matIndex1() const { return mec_->matIndex() + pt1_; }
    
    /// Index of point 2 in the matrix of dynamics (Meca::mISO)
    index_t matIndex2() const { return mec_->matIndex() + pt2_; }
    
    /// true if the pointer seems to be valid.
    bool valid() const { return (mec_!=nullptr) && (pt1_<mec_->nbPoints()) && (pt2_<mec_->nbPoints()); }

    /// Constant pointer to the Mecable
    Mecable const* mecable() const { return mec_; }

    /// Mecapoint corresponding to first point
    Mecapoint vertex1() const { return Mecapoint(mec_, pt1_); }
    
    /// Mecapoint corresponding to second point
    Mecapoint vertex2() const { return Mecapoint(mec_, pt2_); }
    
    /// Index of point 1 in object
    index_t point1() const { return pt1_; }
  
    /// Index of point 2 in object
    index_t point2() const { return pt2_; }
    
    /// Index of point with the smallest weight (ie. point1 or point2)
    index_t lightest_point() const { return ( coef_ > 0.5 ) ? pt1_ : pt2_; }

    /// interpolation coefficient on first point
    real coef0() const { return 1 - coef_; }

    /// interpolation coefficient on second point
    real coef1() const { return coef_; }

    /// interpolation coefficient on first point (historical function)
    //real coef2()  const { return 1 - coef_; }

    /// Set interpolation coefficient
    void coef(real c) { coef_ = c; }
    
    /// Interpolated position in space
    Vector pos()  const { return mec_->interpolatePoints(pt1_, pt2_, coef_); }
    
    /// position of first point
    Vector pos1() const { return mec_->posP(pt1_); }
    
    /// position of second point 
    Vector pos2() const { return mec_->posP(pt2_); }
    
    /// that is pos2() - pos1()
    Vector diff() const { return mec_->diffPoints(pt1_, pt2_); }
    
    /// distance between point1 and point2
    real len()    const { return diff().norm(); }

    /// squared distance between point1 and point2
    real lenSqr() const { return diff().normSqr(); }

    /// normalize(pos2() - pos1())
    Vector dir()  const { return normalize(diff()); }
    
    /// distance between point1 and point2
    real lenInv() const { assert_true(ilen_>0); return ilen_; }

    /// true if the coefficient is in [0, 1]
    bool inside() const { return ( 0 <= coef_ ) && ( coef_ <= 1.0 ); }

    /// test if `this` has a common point with argument
    bool overlapping(Mecapoint const&) const;
    
    /// test if `this` has a common point with argument
    bool overlapping(Interpolation const&) const;

    /// Human friendly ouput
    void print(std::ostream&) const;
    
    /// check validity
    int invalid() const;

};

/// output operator for debugging purpose
std::ostream& operator << (std::ostream&, Interpolation const&);

#endif
