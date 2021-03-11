// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "actor.h"
#include "actor_prop.h"
#include "glossary.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "simul.h"

//------------------------------------------------------------------------------

Actor::Actor(ActorProp const* p, HandMonitor* h)
: Hand(p,h), prop(p)
{
    throw InvalidParameter("the actor class in unfinished");
}


//------------------------------------------------------------------------------

void Actor::stepUnloaded()
{
    assert_true( attached() );
}


void Actor::stepLoaded(Vector const& force)
{
    assert_true( attached() );

    // do something:
}

