// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "bridge.h"
#include "bridge_prop.h"
#include "exceptions.h"
#include "modulo.h"
#include "meca.h"
#include "random.h"

extern Modulo const* modulo;

//------------------------------------------------------------------------------

Bridge::Bridge(BridgeProp const* p, Vector const& w)
: Couple(p, w), prop(p)
{
}


Bridge::~Bridge()
{
    prop = nullptr;
}


void Bridge::stepAA()
{
    Vector f = Bridge::force();
    real fn = f.norm();
    
    if ( cHand1->checkKramersDetachment(fn) )
        cHand1->detach();
    else
        cHand1->stepLoaded( f);
    
    if ( cHand2->checkKramersDetachment(fn) )
        cHand2->detach();
    else
        cHand2->stepLoaded(-f);
}

//------------------------------------------------------------------------------

/**
 Calculates the force for the interLongLink()
 */
Vector Bridge::force() const
{
    Vector d = cHand2->pos() - cHand1->pos();
    
    //correct for periodic space:
    if ( modulo )
        modulo->fold(d);
    
    real dn = d.norm();
    // the norm of the force vector is always known, even if 'd' is nearly zero:
    real f = prop->stiffness * ( dn - prop->length );
    
    if ( dn > REAL_EPSILON )
        return d * ( f / dn );
    
    return Vector(f,0,0);
}


/**
 This uses interLongLink().
 */
void Bridge::setInteractions(Meca& meca) const
{
    meca.addLongLink(cHand1->interpolation(), cHand2->interpolation(), prop->length, prop->stiffness);
}
