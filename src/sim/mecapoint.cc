// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "mecapoint.h"
#include "assert_macro.h"
#include "interpolation.h"
#include "iowrapper.h"
#include "simul.h"


void Mecapoint::read(Inputter& in, Simul& sim)
{
    ObjectTag g;
    Object * obj = sim.readReference(in, g);
    mec_ = Simul::toMecable(obj);
    pti_ = 0;
    if ( obj )
    {
        pti_ = in.readUInt16();
        if ( ! mec_ )
            throw InvalidIO("invalid pointer while reading Mecapoint");
    }
}


void Mecapoint::write(Outputter& out) const
{
    Object::writeReference(out, mec_);
    if ( mec_ ) out.writeUInt16(pti_);
}


bool Mecapoint::overlapping(Mecapoint const& p) const
{
    return ( mec_ == p.mec_  &&  pti_ == p.pti_ );
}


bool Mecapoint::adjacent(Mecapoint const& p) const
{
    return ( mec_ == p.mec_  &&
            ( pti_ == p.pti_ || pti_ == p.pti_+1 || pti_+1 == p.pti_ ));
}


void Mecapoint::print(std::ostream& os) const
{
    if ( mec_ )
        os << "(" << mec_->reference() << "  " << point() << ")";
    else
        os << "(void)";
}


std::ostream& operator << (std::ostream& os, Mecapoint const& arg)
{
    arg.print(os);
    return os;
}


int Mecapoint::invalid() const
{
    if ( !mec_ )
        return 1;

    if ( point() >= mec_->nbPoints() )
        return 2;
    
    return 0;
}

