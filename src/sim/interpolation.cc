// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "fiber_segment.h"
#include "interpolation.h"
#include "mecapoint.h"


/** This converts the dimensionfull abs into an interpolation coefficient in [0,1] */
Interpolation::Interpolation(FiberSegment const& loc, real abs)
{
    mec_ = loc.fiber();
    pt1_ = loc.point();
    pt2_ = loc.point() + 1;
    // convert abscissa to coefficient:
    coef_ = abs * loc.lenInv();
    assert_true( 0 <= coef_ && coef_ <= 1 );
    //coef_ = max_real(0, min_real(1, coef_));
}


bool Interpolation::overlapping(Mecapoint const& p) const
{
    return ( mec_==p.mecable() && ( pt1_==p.point() || pt2_==p.point() ));
}


bool Interpolation::overlapping(Interpolation const& p) const
{
    return ( mec_==p.mec_ &&
            ( pt1_==p.pt1_ || pt1_==p.pt2_ || pt2_==p.pt1_ || pt2_==p.pt2_ ));
}


void Interpolation::print(std::ostream& os) const
{
    if ( mec_ )
        os << "(" << mec_->reference() << " " << pt1_ << " " << pt2_ << " " << coef_ << ")";
    else
        os << "(null)";
}


std::ostream& operator << (std::ostream& os, Interpolation const& arg)
{
    arg.print(os);
    return os;
}


int Interpolation::invalid() const
{
    if ( !mec_ )
        return 1;

    if ( pt1_ >= mec_->nbPoints() )
        return 2;

    if ( pt2_ > mec_->nbPoints() )
        return 3;

    if ( coef_ < 0 )
        return 4;
    
    if ( coef_ > 1 )
        return 5;
    
    return 0;
}

