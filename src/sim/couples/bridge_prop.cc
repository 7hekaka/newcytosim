// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "dim.h"
#include "messages.h"
#include "exceptions.h"
#include "glossary.h"
#include "simul_prop.h"
#include "bridge_prop.h"
#include "bridge.h"


Couple * BridgeProp::newCouple() const
{
    //std::clog << "BridgeProp::newCouple" << '\n';
    return new Bridge(this);
}


void BridgeProp::clear()
{
    CoupleProp::clear();
}


void BridgeProp::read(Glossary& glos)
{
    CoupleProp::read(glos);
}


void BridgeProp::complete(Simul const& sim)
{
    CoupleProp::complete(sim);
}


void BridgeProp::write_values(std::ostream& os) const
{
    CoupleProp::write_values(os);
}

