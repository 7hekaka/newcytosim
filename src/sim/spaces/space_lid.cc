// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University.

#include "dim.h"
#include "space_lid.h"
#include "exceptions.h"
#include "mecapoint.h"
#include "iowrapper.h"
#include "glossary.h"
#include "random.h"
#include "meca.h"


SpaceLid::SpaceLid(SpaceDynamicProp const* p)
: Space(p), prop(p)
{
    if ( DIM == 1 )
        throw InvalidParameter("lid  is not usable in 1D");
    halflength_[0] = 0;
    halflength_[1] = 0;
    bot_ = 0;
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
    
    real bot = bot_, top = top_;
    if ( opt.set(top, "length", DIM-1) )
    {
        bot = -0.5 * top;
        top =  0.5 * top;
    }
    else
    {
        opt.set(bot, "bottom");
        opt.set(top, "top");
    }

    if ( top < bot )
        throw InvalidParameter("lid:top must be >= lid:bottom");
    
    bot_ = bot;
    top_ = top;
    
    update();
}


void SpaceLid::update()
{
    modulo_.reset();
    for ( unsigned d = 0; d < DIM-1; ++d )
        modulo_.enable(d, 2*halflength_[d]);
}


void SpaceLid::boundaries(Vector& inf, Vector& sup) const
{
#if ( DIM >= 3 )
    inf.set(-halflength_[0],-halflength_[1], bot_);
    sup.set( halflength_[0], halflength_[1], top_);
#else
    inf.set(-halflength_[0], bot_, 0);
    sup.set( halflength_[0], top_, 0);
#endif
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
 place only at upper boundary. This overrides the function in Space
 */
Vector SpaceLid::randomPlaceOnEdge(real) const
{
    return Vector(RNG.sfloat()*halflength_[0], top_, 0);
}


//------------------------------------------------------------------------------
#pragma mark -


real SpaceLid::volume() const
{
#if ( DIM == 1 )
    return ( top_ - bot_ );
#elif ( DIM == 2 )
    return 2.0 * halflength_[0] * ( top_ - bot_ );
#else
    return 4.0 * halflength_[0] * halflength_[1] * ( top_ - bot_ );
#endif
}


bool SpaceLid::inside(Vector const& pos) const
{
#if ( DIM == 1 )
    return (( bot_ <= pos.XX ) & ( pos.XX <= top_ ));
#elif ( DIM == 2 )
    return (( bot_ <= pos.YY ) & ( pos.YY <= top_ ));
#else
    return (( bot_ <= pos.ZZ ) & ( pos.ZZ <= top_ ));
#endif
}


bool SpaceLid::allInside(Vector const& pos, const real rad) const
{
    assert_true( rad >= 0 );
#if ( DIM == 1 )
    return (( bot_+rad <= pos.XX ) & ( pos.XX+rad <= top_ ));
#elif ( DIM == 2 )
    return (( bot_+rad <= pos.YY ) & ( pos.YY+rad <= top_ ));
#else
    return (( bot_+rad <= pos.ZZ ) & ( pos.ZZ+rad <= top_ ));
#endif
}


bool SpaceLid::allOutside(Vector const& pos, const real rad) const
{
    assert_true( rad >= 0 );
#if ( DIM == 1 )
    return (( bot_ > pos.XX+rad ) | ( pos.XX > top_+rad ));
#elif ( DIM == 2 )
    return (( bot_ > pos.YY+rad ) | ( pos.YY > top_+rad ));
#else
    return (( bot_ > pos.ZZ+rad ) | ( pos.ZZ > top_+rad ));
#endif
}


Vector SpaceLid::project(Vector const& pos) const
{
#if ( DIM == 1 )
    real X = sign_select(2 * pos.XX - bot_ - top_, bot_, top_);
    return Vector(X);
#elif ( DIM == 2 )
    real Y = sign_select(2 * pos.YY - bot_ - top_, bot_, top_);
    return Vector(pos.XX, Y);
#else
    real Z = sign_select(2 * pos.ZZ - bot_ - top_, bot_, top_);
    return Vector(pos.XX, pos.YY, Z);
#endif
}


//------------------------------------------------------------------------------
#pragma mark - setInteraction


void SpaceLid::setInteraction(Vector const& pos, Mecapoint const& pe,
                              Meca& meca, real stiff) const
{
#if ( DIM == 2 )
    real S = 2 * pos.YY - bot_ - top_;
#elif ( DIM > 2 )
    real S = 2 * pos.ZZ - bot_ - top_;
#endif

    // record force only on top edge:
#if ( DIM == 2 )
    meca.addPlaneClampY(pe, sign_select(S, bot_, top_), stiff);
    force_ += ( S > 0 ) * stiff * ( pos.YY - top_ );
#elif ( DIM > 2 )
    meca.addPlaneClampZ(pe, sign_select(S, bot_, top_), stiff);
    force_ += ( S > 0 ) * stiff * ( pos.ZZ - top_ );
#endif
}


void SpaceLid::setInteraction(Vector const& pos, Mecapoint const& pe, real rad,
                              Meca& meca, real stiff) const
{
#if ( DIM == 2 )
    real S = 2 * pos.YY - bot_ - top_;
#elif ( DIM > 2 )
    real S = 2 * pos.ZZ - bot_ - top_;
#endif

    // record force only on top edge:
#if ( DIM == 2 )
    real Y = sign_select(S, bot_+rad, top_-rad);
    meca.addPlaneClampY(pe, Y, stiff);
    force_ += ( S > 0 ) * stiff * ( pos.YY - Y );
#elif ( DIM > 2 )
    real Z = sign_select(S, bot_+rad, top_-rad);
    meca.addPlaneClampZ(pe, Z, stiff);
    force_ += ( S > 0 ) * stiff * ( pos.ZZ - Z );
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
    out.writeFloat(bot_);
    out.writeFloat(top_);
    out.writeFloat(0.f);
    out.writeFloat(force_);
}


void SpaceLid::setLengths(const real len[])
{
    halflength_[0] = len[0];
    halflength_[1] = len[1];
    bot_   = len[2];
    top_   = len[3];
    force_ = len[5];
    update();
}


void SpaceLid::read(Inputter& in, Simul&, ObjectTag)
{
    real len[8] = { 0 };
    read_data(in, 8, len, "lid");
    setLengths(len);
}


//------------------------------------------------------------------------------
#pragma mark - Display

#ifdef DISPLAY
#include "opengl.h"
#include "gle.h"
using namespace gle;

bool SpaceLid::draw() const
{
    const real X = halflength_[0];
    const real T = top_;
    const real B = bot_;
    
#if ( DIM >= 3 )
    const real Y = halflength_[1];
    // draw faces:
    glBegin(GL_TRIANGLE_STRIP);
    glNormal3f(0, 0, 1);
    gleVertex(-X,  Y, B);
    gleVertex( X,  Y, B);
    gleVertex(-X, -Y, B);
    gleVertex( X, -Y, B);
    glEnd();
    glBegin(GL_TRIANGLE_STRIP);
    glNormal3f(0, 0, -1);
    gleVertex(-X,  Y, T);
    gleVertex(-X, -Y, T);
    gleVertex( X,  Y, T);
    gleVertex( X, -Y, T);
    glEnd();
    // draw outline:
    glBegin(GL_LINE_STRIP);
    gleVertex(-X,  Y, B);
    gleVertex( X,  Y, B);
    gleVertex( X, -Y, B);
    gleVertex(-X, -Y, B);
    gleVertex(-X,  Y, B);
    glEnd();
    glBegin(GL_LINE_STRIP);
    gleVertex(-X,  Y, T);
    gleVertex(-X, -Y, T);
    gleVertex( X, -Y, T);
    gleVertex( X,  Y, T);
    gleVertex(-X,  Y, T);
    glEnd();
#else
    glBegin(GL_LINES);
    gleVertex(-X, T, 0);
    gleVertex( X, T, 0);
    gleVertex( X, B, 0);
    gleVertex(-X, B, 0);
    glEnd();
#endif
    
    // draw periodic boundaries:
    glLineStipple(1, 0x000F);
    glEnable(GL_LINE_STIPPLE);
    glBegin(GL_LINES);
#if ( DIM >= 3 )
    gleVertex( X,  Y, T);
    gleVertex( X,  Y, B);
    gleVertex( X, -Y, T);
    gleVertex( X, -Y, B);
    gleVertex(-X,  Y, T);
    gleVertex(-X,  Y, B);
    gleVertex(-X, -Y, T);
    gleVertex(-X, -Y, B);
#else
    gleVertex( X, T, 0);
    gleVertex( X, B, 0);
    gleVertex(-X, T, 0);
    gleVertex(-X, B, 0);
#endif
    glEnd();
    glDisable(GL_LINE_STIPPLE);
    return true;
}

#else

bool SpaceLid::draw() const
{
    return false;
}

#endif

