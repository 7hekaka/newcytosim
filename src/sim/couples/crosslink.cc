// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "crosslink.h"
#include "crosslink_prop.h"
#include "object_set.h"
#include "exceptions.h"
#include "random.h"
#include "modulo.h"
#include "space.h"
#include "meca.h"

//------------------------------------------------------------------------------
Crosslink::Crosslink(CrosslinkProp const* p, Vector const& w)
: Couple(p, w)
{
}


Crosslink::~Crosslink()
{
}


//------------------------------------------------------------------------------
#pragma mark -

/**
 Simulates:
 - diffusive motion
 - attachment
 .
 */
void Crosslink::stepFF()
{
    diffuse();
    
    // activity:
    if ( RNG.flip() )
        cHand1->stepUnattached(simul(), cPos);
    else
        cHand2->stepUnattached(simul(), cPos);
}


void Crosslink::setInteractions(Meca& meca) const
{
    assert_true( attached1() && attached2() );
    
    meca.addLink(cHand1->interpolation(), cHand2->interpolation(), prop()->stiffness);
}

