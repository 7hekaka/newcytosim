// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "event.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "simul.h"


//------------------------------------------------------------------------------
Event::~Event()
{
    //Cytosim::log("destroying Event %p\n", this);
}



//------------------------------------------------------------------------------
void Event::write(Outputter& out) const
{
}


void Event::read(Inputter & in, Simul& sim, ObjectTag tag)
{
}
