// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "shackle.h"
#include "shackle_prop.h"
#include "meca.h"


Shackle::Shackle(ShackleProp const* p, Vector const& w)
: Couple(p, w)
{
}


Shackle::~Shackle()
{
}


/**
 The interaction is slipery on hand1
 */
void Shackle::setInteractions(Meca& meca) const
{
    Interpolation const& pt1 = cHand1->interpolation();
    Interpolation const& pt2 = cHand2->interpolation();

    meca.addSlidingLink(pt1, pt2, prop()->stiffness);
}


void Shackle::stepAA()
{
    real dis = INFINITY;
    
    // project the position of cHand2 to set abscissa of cHand1
    real abs = cHand1->fiber()->projectPoint(cHand2->pos(), dis);
    
    //std::clog << "Shackle delta-abs = " << abs - cHand1->abscissa() << '\n';
    cHand1->moveTo(abs);
    // need to update to calculate force:
    cHand1->reinterpolate();
    
    Vector f = Couple::force();
    real mag = f.norm();
    
    if ( cHand1->checkKramersDetachment(mag) )
        cHand1->detach();
    else
        cHand1->stepLoaded( f);

    if ( cHand2->checkKramersDetachment(mag) )
        cHand2->detach();
    else
        cHand2->stepLoaded(-f);
}

