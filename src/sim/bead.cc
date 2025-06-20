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
    setDragCoefficient();
#if NEW_SOLID_CLAMP
    clamp_place.reset();
    clamp_stiff = 0;
#endif
}


Bead::~Bead()
{
    if ( objset() )
        simul().singles.deleteWrists(this);
    prop = nullptr;
}


void Bead::allocateMecable(const index_t nbp)
{
    assert_true( nbp == 1 );
    //Mecable::allocateMemory(nbp, 0);
    /* Directly setting the Vertex pointer to the member data,
     avoids the minuscule allocation (3 reals) that would occur
     for each Bead in the system.
     This will work as long as Vector can directly match a real[]
     */
    pPos = paCenter.data();
}

//------------------------------------------------------------------------------

void Bead::setInteractions(Meca& meca) const
{
#if NEW_SOLID_CLAMP
    if ( clamp_stiff > 0 )
        meca.addPointClamp(Mecapoint(this,0), clamp_place, clamp_stiff);
#endif

    if ( prop->confine != CONFINE_OFF )
    {
        Space const* spc = prop->confine_space;
        assert_true(spc);
        
        switch ( prop->confine )
        {
            case CONFINE_INSIDE:
            {
                // Confine only the center
                if ( ! spc->inside(paCenter) )
                    spc->setConfinement(paCenter, Mecapoint(this, 0), meca, prop->confine_stiff);
            } break;
                
            case CONFINE_OUTSIDE:
            {
                // confine the center outside
                if ( spc->inside(paCenter) )
                    spc->setConfinement(paCenter, Mecapoint(this, 0), meca, prop->confine_stiff);
            } break;
                
            case CONFINE_ALL_INSIDE:
            {
                // Confine the entire bead
                if ( ! spc->allInside(paCenter, paRadius) )
                    spc->setConfinement(paCenter, Mecapoint(this, 0), paRadius, meca, prop->confine_stiff);
            } break;
                
            case CONFINE_ON:
                spc->setConfinement(paCenter, Mecapoint(this, 0), meca, prop->confine_stiff);
                break;
                
            default:
                throw InvalidParameter("Invalid bead::confine");
        }
    }
}


real Bead::addBrownianForces(real const* fce, real alpha, real* rhs) const
{
    // Brownian amplitude:
    real b = std::sqrt( alpha * paDrag );

    for ( size_t d = 0; d < DIM; ++d )
        rhs[d] = fce[d] + b * rhs[d];
    
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
    if ( paRadius > 0 )
    {
        static std::string msg;
        std::ostringstream oss;
        oss << "Bead `" << prop->name() << "' (radius " << paRadius << ") has drag " << paDrag << '\n';
        if ( msg != oss.str() )
        {
            msg = oss.str();
            std::clog << msg;
        }
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
    out.writeFloats(paCenter, DIM, '\n');
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

