// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include "assert_macro.h"
#include "bead.h"
#include "bead_prop.h"
#include "exceptions.h"
#include "single.h"
#include "hand_prop.h"
#include "iowrapper.h"
#include "glossary.h"
#include "meca.h"
#include "simul.h"
#include "space.h"
#include "modulo.h"


//------------------------------------------------------------------------------

Bead::Bead(BeadProp const* p, Vector pos, real rad)
: paRadius(rad), paDrag(0), prop(p)
{
    assert_true(rad >= 0);
    setNbPoints(1);
    setPoint(0, pos);
}


Bead::~Bead()
{
    if ( objset() )
        simul().singles.deleteWrists(this);
    prop = nullptr;
}


//------------------------------------------------------------------------------

void Bead::setInteractions(Meca& meca) const
{
#if NEW_SOLID_CLAMP
    if ( prop->clamp_stiff > 0 )
        meca.addPointClamp(Mecapoint(this,0), prop->clamp_pos, prop->clamp_stiff);
#endif

    if ( prop->confine != CONFINE_OFF )
    {
        Space const* spc = prop->confine_space_ptr;
        
        switch ( prop->confine )
        {
            case CONFINE_INSIDE:
            {
                // Confine only the center
                Vector cen(pPos);
                if ( ! spc->inside(cen) )
                    spc->setConfinement(cen, Mecapoint(this, 0), meca, prop->confine_stiffness);
            } break;
                
            case CONFINE_OUTSIDE:
            {
                // confine the center outside
                Vector cen(pPos);
                if ( spc->inside(cen) )
                    spc->setConfinement(cen, Mecapoint(this, 0), meca, prop->confine_stiffness);
            } break;
                
            case CONFINE_ALL_INSIDE:
            {
                // Confine the entire bead
                Vector cen(pPos);
                if ( ! spc->allInside(cen, paRadius) )
                    spc->setConfinement(cen, Mecapoint(this, 0), paRadius, meca, prop->confine_stiffness);
            } break;
                
            case CONFINE_ON:
                spc->setConfinement(position(), Mecapoint(this, 0), meca, prop->confine_stiffness);
                break;
                
            default:
                throw InvalidParameter("Invalid bead::confine");
        }
    }
}


real Bead::addBrownianForces(real const* rnd, real alpha, real* rhs) const
{
    // Brownian amplitude:
    real b = std::sqrt( alpha * paDrag );

    for ( size_t d = 0; d < DIM; ++d )
        rhs[d] += b * rnd[d];
    
    //the amplitude is needed in Meca
    return b / paDrag;
}


/**
 If `drag` is not specified, its value is calculated using Stokes' law:

       drag = 6 * M_PI * viscosity * radius;

*/
void Bead::setDragCoefficient()
{
    paDrag = ( 6 * M_PI ) * ( prop->viscosity * paRadius );
    if ( prop->drag > 0 )
    {
        //std::clog << "drag set for `" << prop->name() << "' bypassing Stokes' law\n";
        paDrag = prop->drag;
    }
    
#if ( 0 )
    static bool virgin = true;
    if ( paRadius > 0  &&  virgin )
    {
        std::clog << "Bead `" << prop->name() << "' (radius " << paRadius << ") has drag " << paDrag << '\n';
        virgin = false;
    }
#endif
}


/**
 The projection here just scales by the mobility
 */
void Bead::projectForces(const real* X, real* Y) const
{
    assert_true( paDrag > 0 );
    real s = 1.0 / paDrag;
    for ( size_t d = 0; d < DIM; ++d )
        Y[d] = s * X[d];
}


//------------------------------------------------------------------------------

void Bead::write(Outputter& out) const
{
    writeMarker(out, TAG);
    out.writeFloats(position(), DIM, '\n');
    out.writeSoftSpace();
    out.writeFloat(paRadius);
}


void Bead::read(Inputter& in, Simul&, ObjectTag)
{
    Vector pos;
    in.readFloats(pos, DIM);
    setPoint(0, pos);
    real r = in.readFloat();
    resize(r);
}

