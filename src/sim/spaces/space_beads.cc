// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "dim.h"
#include "space_beads.h"
#include "bead_prop.h"
#include "object_set.h"
#include "glossary.h"
#include "simul.h"


/**
 The parameters BeadSet and BeadProp define the inner volume for this Space
*/
SpaceBeads::SpaceBeads(SpaceProp const* p)
: Space(p), mBeadProp(nullptr)
{
    mBeadName = "undefined";
    for ( int d = 0; d < 3; ++d )
    {
        bbMin[d] = 0;
        bbMax[d] = 0;
    }
}


void SpaceBeads::setBeads()
{
    Simul const& sim = simul();
    mBeadProp = sim.findProperty<BeadProp>("bead", mBeadName);

    mBeads.clear();
    for ( Bead * b = sim.beads.first(); b; b=b->next() )
    {
        if ( b->property() == mBeadProp )
            mBeads.push_back(b);
    }
    std::clog << "SpaceBeads has " << mBeads.size() << " beads\n";
}


/**
 refresh the list of Beads 
 */
void SpaceBeads::resize(Glossary& opt)
{
    opt.set(mBeadName, "bead") || opt.set(mBeadName, "beads");
}

/**
 Calculate a box in which the box are entirely inside
 */
void SpaceBeads::setBoundaries()
{
    for ( unsigned d = 0; d < DIM; ++d )
    {
        bbMin[d] =  INFINITY;
        bbMax[d] = -INFINITY;
    }
    
    for ( Bead * i : mBeads )
    {
        real rad = i->radius();
        Vector pos = i->position();
        for ( unsigned d = 0; d < DIM; ++d )
        {
            bbMin[d] = std::min(bbMin[d], pos[d]-rad);
            bbMax[d] = std::max(bbMax[d], pos[d]+rad);
        }
    }
    //std::clog << "SpaceBeads::boundaries " << bbMin[0] << "  " << bbMax[0] << '\n';
}


void SpaceBeads::step()
{
    if ( !mBeadProp ) setBeads();
    setBoundaries();
}


void SpaceBeads::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(bbMin[0], bbMin[1], bbMin[2]);
    sup.set(bbMax[0], bbMax[1], bbMax[2]);
}


/**
 Calculates the sum of the bead's volumes,
 thus assuming that the beads do no overlap
 */
real SpaceBeads::volume() const
{
    real res = 0;
    for ( Bead * i : mBeads )
        res += i->volume();
    return res;
}


bool SpaceBeads::inside(Vector const& pos) const
{
    for ( unsigned d = 0; d < DIM; ++d )
    {
        if ( pos[d] < bbMin[d] ) return false;
        if ( bbMax[d] < pos[d] ) return false;
    }
    
    for ( Bead * i : mBeads )
        if ( distanceSqr(pos, i->position()) < i->radiusSqr() )
            return true;
    
    return false;
}


Vector SpaceBeads::project(Vector const& pos) const
{
    Bead * bid = nullptr;
    real dmin = INFINITY;

    for ( Bead * i : mBeads )
    {
        real d = abs_real(distance(pos, i->position()) - i->radius());
        if ( d < dmin )
        {
            dmin = d;
            bid  = i;
        }
    }
    
    Vector vec = ( pos - bid->position() ).normalized(bid->radius());
    return bid->position() + vec;
}


//------------------------------------------------------------------------------

void SpaceBeads::setInteraction(Vector const& pos, Mecapoint const& pe, Meca& meca, real stiff) const
{
    ABORT_NOW("SpaceBeads is incomplete");
}


void SpaceBeads::setInteraction(Vector const& pos, Mecapoint const& pe, real rad, Meca& meca, real stiff) const
{
    ABORT_NOW("SpaceBeads is incomplete");
}

