// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "dim.h"
#include "space_rotate.h"
#include "exceptions.h"


SpaceRotate::SpaceRotate(Space * spc)
: Space(spc->prop)
{
    if ( DIM <= 1 )
        throw InvalidParameter("space:rotate is only valid in DIM=2 or 3");
    if ( !spc )
        throw InvalidParameter("space:rotate invoked with void argument");

    mSpace = spc;
}


SpaceRotate::~SpaceRotate()
{
    delete(mSpace);
    mSpace = nullptr;
}


Vector SpaceRotate::forward(Vector const& pos) const
{
    return Vector(-pos[2], pos[1], pos[0]);
}


Vector SpaceRotate::backward(Vector const& pos) const
{
    return Vector(pos[2], pos[1], -pos[0]);
}


void SpaceRotate::boundaries(Vector& inf, Vector& sup) const
{
    mSpace->boundaries(inf, sup);
    inf = forward(inf);
    sup = forward(sup);
}


bool SpaceRotate::inside(Vector const& pos) const
{
    return mSpace->inside(backward(pos));
}


Vector SpaceRotate::project(Vector const& pos) const
{
    return forward(mSpace->project(backward(pos)));
}


/* It would be necessary to use a temporary Meca, to swap the indices between X and Z */
void SpaceRotate::setInteraction(Vector const& pos, Mecapoint const&, Meca&, real stiff) const
{
    ABORT_NOW("unfinished SpaceRotate");
}


/* It would be necessary to use a temporary Meca, to swap the indices between X and Z */
void SpaceRotate::setInteraction(Vector const& pos, Mecapoint const&, real rad, Meca&, real stiff) const
{
    ABORT_NOW("unfinished SpaceRotate");
}

//------------------------------------------------------------------------------
//                         OPENGL  DISPLAY
//------------------------------------------------------------------------------

#ifdef DISPLAY
#include "opengl.h"

bool SpaceRotate::draw() const
{
    glRotated(90, 0,  1, 0);
    mSpace->draw();
    glRotated(90, 0, -1, 0);
    return true;
}

#else

bool SpaceRotate::draw() const
{
    return false;
}

#endif
