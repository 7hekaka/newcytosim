// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "dim.h"
#include "space_beads.h"
#include "bead_prop.h"
#include "object_set.h"
#include "simul.h"



/**
 The parameters BeadSet and BeadProp define the inner volume for this Space
*/
SpaceBeads::SpaceBeads(const SpaceProp* p)
: Space(p)
{
    for ( int d = 0; d < 3; ++d )
    {
        bbMin[d] = 0;
        bbMax[d] = 0;
    }
}

/**
 refresh the list of Beads 
 */
void SpaceBeads::resize(Glossary& opt)
{
    if ( objset() )
    {
        std::string name;
        opt.set(name, "bead");
        Simul const& sim = simul();
        Property * bip = sim.properties.find("bead", name);

        mBeads.clear();
        if ( !bip )
        {
            for ( Bead * bd = sim.beads.first(); bd; bd=bd->next() )
                mBeads.push_back(bd);
        }
        else
        {
            for ( Bead * bd = sim.beads.first(); bd; bd=bd->next() )
            {
                if ( bd->property() == bip )
                    mBeads.push_back(bd);
            }
        }
#if ( 1 )
        static size_t nb = 0;
        if ( nb != mBeads.size() )
        {
            nb = mBeads.size();
            std::clog << "SpaceBeads has " << nb << " beads\n";
        }
#endif
    }
}

/**
 Calculate a box in which the box are entirely inside
 */
void SpaceBeads::setBoundaries()
{
    for ( int d = 0; d < DIM; ++d )
    {
        bbMin[d] =  INFINITY;
        bbMax[d] = -INFINITY;
    }
    
    for ( Bead * i : mBeads )
    {
        real rad = i->radius();
        Vector pos = i->position();
        for ( int d = 0; d < DIM; ++d )
        {
            bbMin[d] = std::min(bbMin[d], pos[d]-rad);
            bbMax[d] = std::max(bbMax[d], pos[d]+rad);
        }
    }
    //std::clog << "SpaceBeads::bounding " << bbMin[0] << "  " << bbMax[0] << std::endl;
}


void SpaceBeads::step()
{
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


bool SpaceBeads::inside(Vector const& point) const
{
    for ( int d = 0; d < DIM; ++d )
    {
        if ( point[d] < bbMin[d] ) return false;
        if ( bbMax[d] < point[d] ) return false;
    }
    
    Vector pos(point);
    
    for ( Bead * i : mBeads )
        if ( pos.distanceSqr(i->position()) < i->radiusSqr() )
            return true;
    
    return false;
}


Vector SpaceBeads::project(Vector const& pos) const
{
    Bead * bid = nullptr;
    real dmin = INFINITY;

    for ( Bead * i : mBeads )
    {
        real d = fabs(pos.distance(i->position()) - i->radius());
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

void SpaceBeads::setInteraction(Vector const& pos, Mecapoint const& pe, Meca & meca, real stiff) const
{
    ABORT_NOW("SpaceBeads is incomplete");
}


void SpaceBeads::setInteraction(Vector const& pos, Mecapoint const& pe, real rad, Meca & meca, real stiff) const
{
    ABORT_NOW("SpaceBeads is incomplete");
}

//------------------------------------------------------------------------------
//                         OPENGL  DISPLAY
//------------------------------------------------------------------------------

bool SpaceBeads::draw() const
{
    return true;
}


