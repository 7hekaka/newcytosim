// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "organizer.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "simul.h"


Organizer::~Organizer()
{
    //Cytosim::log("destroying Organizer ", this, "\n");
}


void Organizer::grasp(Mecable * m)
{
    Buddy::connect(m);
    mObjects.push_back(m);
}


void Organizer::grasp(Mecable * m, size_t ix)
{
    assert_true( ix < mObjects.size() );

    if ( m != mObjects[ix] )
    {
        Buddy::disconnect(mObjects[ix]);
        Buddy::connect(m);
    }
    
    mObjects[ix] = m;
}


bool Organizer::check(Mecable const* m) const
{
    MecableList::const_iterator i = std::find(mObjects.begin(), mObjects.end(), m);
    return ( i != mObjects.end() );
}


void Organizer::goodbye(Buddy const* b)
{
    //std::clog << reference() << " lost " << b << "\n";
    MecableList::iterator i = std::find(mObjects.begin(), mObjects.end(), b);
    if ( i != mObjects.end() )
        *i = nullptr;
}

/* This is not necessary anymore, and object should all be linked */
void Organizer::checkOrganized() const
{
    for ( Mecable * i : mObjects )
    {
        if ( i && ! i->linked() )
        {
            std::clog << " Organized " << i->reference() << " is not linked\n";
        }
    }
}

//------------------------------------------------------------------------------
#pragma mark - Placement

/**
 \return The centroid from all the object positions
 */
Vector Organizer::position() const
{
    Vector res(0,0,0);
    for ( Mecable const* i : mObjects )
        res += i->position();
    return res / (real)mObjects.size();
}


Vector Organizer::positionP(index_t ix) const
{
    Vector res(0,0,0);
    for ( Mecable const* i : mObjects )
        res += i->posPoint(ix);
    return res / (real)mObjects.size();
}

/*
void Organizer::translate(Vector const& T)
{
    for ( Mecable * mec : mObjects )
    {
        if ( mec && mec->mobile() & 1 )
        {
            mec->translate(T);
            mec->flag(0);
        }
    }
}


void Organizer::rotate(Rotation const& T)
{
    for ( Mecable * mec : mObjects )
    {
        if ( mec && mec->mobile() & 2 )
        {
            mec->rotate(T);
            mec->flag(0);
        }
    }
}
*/

real Organizer::sumDragCoefficient() const
{
    real res = 0;
    for ( Mecable const* i : mObjects )
        res += i->dragCoefficient();
    return res;
}


//------------------------------------------------------------------------------
#pragma mark - I/O

void Organizer::writeOrganized(Outputter& out) const
{
    out.writeUInt16(mObjects.size());
    for ( Mecable const* i : mObjects )
    {
        Mecable::writeReference(out, i);
    }
}


void Organizer::readOrganized(Inputter& in, Simul& sim, size_t nbo)
{
    nbOrganized(nbo);
    
    //std::clog << " Organizer::read with " << nb << " objects" << '\n';
    ObjectTag g;
    for ( size_t i = 0; i < nbo; ++i )
    {
        Object * w = sim.readReference(in, g);
        if ( w )
        {
            //std::clog << " Organized(" << i << ") is " << w->reference() << '\n';
            grasp(Simul::toMecable(w), i);
        }
        else
            grasp(nullptr, i);
    }
}


void Organizer::read(Inputter& in, Simul& sim, ObjectTag)
{
    readOrganized(in, sim, in.readUInt16());
}

