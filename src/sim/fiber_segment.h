// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef FIBER_SEGMENT_H
#define FIBER_SEGMENT_H


#include "real.h"
#include "vector.h"
#include "mecapoint.h"
#include "fiber.h"


/// Indicates the segment of a Fiber located between two consecutive vertices
/** 
 FiberSegment is used to refer to the section of a Fiber located between two vertices.
 The i-th FiberSegment covers the section [i, i+1[

 This is used to calculate the distance to this segment,
 the intersection of the segment with a plane.
*/
class FiberSegment final
{
private:
    
    /// Fiber to which the segment belongs to
    Fiber const * fib_;
    
    /// index of segment's first point
    index_t sgi_;
    
public:
    
    /// construct without initialization
    FiberSegment() : sgi_(0) {}
    
    /// constructor
    FiberSegment(Fiber const* f, unsigned p) : fib_(f), sgi_(p) {}
    
    /// setter
    void set(Fiber const* f, unsigned p) { fib_ = f; sgi_ = p; }

    /// the Fiber
    Fiber const* fiber() const { return fib_; }
    
    /// the Fiber's identity
    ObjectID identity() const { if ( fib_ ) return fib_->identity(); return 0; }
    
    /// index of segment
    index_t point() const { return sgi_; }
    
    /// set segment index
    void point(unsigned p) { sgi_ = p; }

    /// Index of segment's first vertex in the isotropic matrix (Meca::mISO)
    index_t matIndex0()  const { return fib_->matIndex() + sgi_; }

    /// abscissa at start of segment (i.e. corresponding to point())
    real abscissa1()    const { return fib_->abscissaPoint(sgi_); }
    
    /// abscissa of second point
    real abscissa2()    const { return fib_->abscissaPoint(sgi_+1); }

    /// the length of the segment
    real len()          const { return fib_->segmentation(); }
    
    /// should return 1.0 / len()
    real lenInv()       const { return fib_->segmentationInv(); }
    
    /// true if abscissa 'a', counted from point 0 is within the segment
    bool within(real a) const { return ( 0 <= a ) & ( a <= fib_->segmentation() ); }
    
    /// position of first point
    Vector pos1()       const { return fib_->posP(sgi_); }
    
    /// position of second point
    Vector pos2()       const { return fib_->posP(sgi_+1); }

    /// interpolated position, where c is in [0, 1]
    Vector midPoint(real c) const { return fib_->midPoint(sgi_, c); }
    
    /// that is [ pos2() + pos1() ] / 2
    Vector middle()     const { return fib_->midPoint(sgi_); }

    /// that is pos2() - pos1()
    Vector diff()       const { return fib_->diffPoints(sgi_); }

    /// normalized tangent to Fiber
    Vector dir()        const { return fib_->dirSegment(sgi_); }
    
    /// Mecapoint corresponding to first point
    Mecapoint vertex1() const { return Mecapoint(fib_, sgi_); }
    
    /// Mecapoint corresponding to second point
    Mecapoint vertex2() const { return Mecapoint(fib_, sgi_+1); }
    
    /// true if the segment is the first of the Fiber
    bool isFirst()  const { return ( sgi_ == 0U ); }

    /// true if the segment is not the first of the Fiber
    bool notFirst() const { return ( sgi_ > 0U ); }
    
    /// true if the segment is the last of the fiber
    bool isLast()   const { return ( sgi_+2U == fib_->nbPoints() ); }
    
    /// true if the segment is not the last of the fiber
    bool notLast()  const { return ( sgi_+2U < fib_->nbPoints() ); }

    
    /// return abscissa of the projection of `w` on the line supporting the segment, and set distance
    real projectPoint0(Vector w, real& dis2) const;

    /// return abscissa in [0, 1] of the closest point to `w` on the segment, and set distance
    real projectPoint(Vector w, real& dis2) const;

    /// old hand-crafted projectionPoint(), but incompatible with periodic boundary conditions
    real projectPointF(const real[], real& dis2) const;
    
    /// calculates the closest distance between two infinite lines and set corresponding abscissa
    real shortestDistanceSqr(FiberSegment const&, real& a, real& b) const;

    /// Human friendly ouput
    void print(std::ostream&) const;

    /// Human friendly ouput
    std::string to_string() const;
};

/// print for debug purpose
std::ostream& operator << (std::ostream&, FiberSegment const&);

#endif

