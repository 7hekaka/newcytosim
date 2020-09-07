// Cytosim 3.0 -  Copyright Francois Nedelec et al.  EMBL 2007-2013

#include "space_lid.h"
#include "exceptions.h"
#include "mecapoint.h"
#include "iowrapper.h"
#include "glossary.h"
#include "meca.h"

#include "random.h"


SpaceLid::SpaceLid(SpaceDynamicProp const* p)
: Space(p), prop(p)
{
    if ( DIM == 1 )
        throw InvalidParameter("lid is only valid in DIM=2 or 3");
    
    for ( int d = 0; d < 3; ++d )
        halflength_[d] = 0;
    top_ = 0;
    force_ = 0;
}


void SpaceLid::resize(Glossary& opt)
{
    for ( unsigned d = 0; d < DIM; ++d )
    {
        real len = halflength_[d];
        if ( opt.set(len, "length", d) )
            len *= 0.5;
        if ( len <= 0 )
            throw InvalidParameter("lid:length_[] must be > 0");
        halflength_[d] = len;
    }
    
    opt.set(top_, "top");
    if ( top_ < 0 )
        throw InvalidParameter("lid:top must be >= 0");
    
    update();
}


void SpaceLid::update()
{
    modulo_.reset();
    modulo_.enable(0, 2*halflength_[0]);
}


void SpaceLid::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(-halflength_[0],-halflength_[1],-halflength_[2]);
    sup.set( halflength_[0], halflength_[1], top_);
}


void SpaceLid::bounce(Vector& pos) const
{
    if ( !SpaceLid::inside(pos) )
        bounceOnEdges(pos);
    
    // periodic in all except the last dimension:
#if ( DIM > 1 )
    pos.XX = fold_real(pos.XX, modulo_.period_[0]);
#endif
#if ( DIM > 2 )
    pos.YY = fold_real(pos.YY, modulo_.period_[1]);
#endif
}


/**
 place only at upper boundary. This overrides the function in spaces
 */
Vector SpaceLid::randomPlaceOnEdge(real) const
{
    return Vector( RNG.sfloat()*halflength_[0], top_, 0 );
}


//------------------------------------------------------------------------------
#pragma mark -


real SpaceLid::volume() const
{
#if ( DIM == 1 )
    return top_ + halflength_[0];
#elif ( DIM == 2 )
    return 2.0 * halflength_[0] * ( top_ + halflength_[1] );
#else
    return 4.0 * halflength_[0] * halflength_[1] * ( top_ + halflength_[2] );
#endif
}


bool SpaceLid::inside(Vector const& pos) const
{
    return ( pos[DIM-1] < top_ );
}


bool SpaceLid::allInside(Vector const& cen, const real rad) const
{
    assert_true( rad >= 0 );
    return ( cen[DIM-1] + rad < top_ );
}


bool SpaceLid::allOutside(Vector const& cen, const real rad) const
{
    assert_true( rad >= 0 );
    return ( cen[DIM-1] + rad > top_ );
}


Vector SpaceLid::project(Vector const& pos) const
{
#if ( DIM > 2 )
    return Vector(pos.XX, pos.YY, top_);
#else
    return Vector(pos.XX, top_, 0);
#endif
}


//------------------------------------------------------------------------------
#pragma mark -


void SpaceLid::setInteraction(Vector const& pos, Mecapoint const& pe,
                              Meca& meca, real stiff) const
{
#if ( DIM == 2 )
    meca.addPlaneClampY(pe, top_, stiff);
#elif ( DIM > 2 )
    meca.addPlaneClampZ(pe, top_, stiff);
#endif
    
#if ( DIM == 2 )
    force_ += stiff * ( pos.YY - top_ );
#elif ( DIM > 2 )
    force_ += stiff * ( pos.ZZ - top_ );
#endif
}


void SpaceLid::setInteraction(Vector const& pos, Mecapoint const& pe, real rad,
                              Meca& meca, real stiff) const
{
#if ( DIM == 2 )
    meca.addPlaneClampY(pe, top_-rad, stiff);
#elif ( DIM > 2 )
    meca.addPlaneClampZ(pe, top_-rad, stiff);
#endif
    
#if ( DIM == 2 )
    force_ += stiff * ( pos.YY - top_ + rad );
#elif ( DIM > 2 )
    force_ += stiff * ( pos.ZZ - top_ + rad );
#endif
}

void SpaceLid::setInteractions(Meca& meca) const
{
    force_ = 0;
}

void SpaceLid::step()
{
    real dc = prop->mobility_dt * force_;
    
    if ( abs_real(dc) < 1 )
        top_ += dc;
    else
        std::cerr << "Error: lid displacement is too fast: " << dc << '\n';
    std::cerr << "force on lid is " << force_ << '\n';
    
    if ( top_ > halflength_[DIM-1] )
        std::cerr << "Warning: space lid has reached its maximum\n";
}


//------------------------------------------------------------------------------

void SpaceLid::write(Outputter& out) const
{
    out.put_characters("lid", 16);
    out.writeUInt16(6);
    out.writeFloat(halflength_[0]);
    out.writeFloat(halflength_[1]);
    out.writeFloat(halflength_[2]);
    out.writeFloat(top_);
    out.writeFloat(force_);
    out.writeFloat(0.f);
}


void SpaceLid::setLengths(const real len[])
{
    halflength_[0] = len[0];
    halflength_[1] = len[1];
    halflength_[2] = len[2];
    top_       = len[3];
    force_     = len[4];
    update();
}

void SpaceLid::read(Inputter& in, Simul&, ObjectTag)
{
    real len[8] = { 0 };
    read_data(in, len, "lid");
    setLengths(len);
}

//------------------------------------------------------------------------------
//                         OPENGL  DISPLAY
//------------------------------------------------------------------------------
#pragma mark -

#ifdef DISPLAY
#include "opengl.h"
#include "gle.h"
using namespace gle;

bool SpaceLid::draw() const
{
    const real X = halflength_[0];
    const real Y = ( DIM > 1 ) ? halflength_[1] : 1;
    const real Z = ( DIM > 2 ) ? halflength_[2] : 0;

#if ( DIM == 2 )
    glBegin(GL_LINES);
    gleVertex( -X, top_, Z );
    gleVertex(  X, top_, Z );
    glEnd();
#elif ( DIM > 2 )
    glBegin(GL_TRIANGLE_STRIP);
    gleVertex( -X, -Y, top_ );
    gleVertex(  X, -Y, top_ );
    gleVertex( -X,  Y, top_ );
    gleVertex(  X,  Y, top_ );
    glEnd();
#endif
    
    glLineWidth(3);
    glLineStipple(1, 0x0303);
    glEnable(GL_LINE_STIPPLE);
    glBegin(GL_LINE_LOOP);
    gleVertex( -X, -Y, Z );
    gleVertex(  X, -Y, Z );
    gleVertex(  X,  Y, Z );
    gleVertex( -X,  Y, Z );
    gleVertex( -X, -Y, Z );
    glEnd();
    
#if ( DIM > 2 )
    glBegin(GL_LINE_LOOP);
    gleVertex( -X, -Y, -Z );
    gleVertex(  X, -Y, -Z );
    gleVertex(  X,  Y, -Z );
    gleVertex( -X,  Y, -Z );
    gleVertex( -X, -Y, -Z );
    glEnd();

    glBegin(GL_LINES);
    gleVertex(  X,  Y, -Z );
    gleVertex(  X,  Y,  Z );
    gleVertex( -X,  Y, -Z );
    gleVertex( -X,  Y,  Z );
    gleVertex( -X, -Y, -Z );
    gleVertex( -X, -Y,  Z );
    gleVertex(  X, -Y, -Z );
    gleVertex(  X, -Y,  Z );
    glEnd();
#endif
  
    glDisable(GL_LINE_STIPPLE);

    return true;
}

#else

bool SpaceLid::draw() const
{
    return false;
}

#endif

