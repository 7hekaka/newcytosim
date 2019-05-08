// Cytosim 3.0 -  Copyright Francois Nedelec et al.  EMBL 2007-2013

#include "space_lid.h"
#include "exceptions.h"
#include "mecapoint.h"
#include "iowrapper.h"
#include "glossary.h"
#include "meca.h"

#include "random.h"


SpaceLid::SpaceLid(SpaceProp const* p)
: Space(p)
{
    if ( DIM == 1 )
        throw InvalidParameter("lid is only valid in DIM=2 or 3");
    
    for ( int d = 0; d < 3; ++d )
        length_[d] = 0;
    top_ = 0;
    force_ = 0;
}


void SpaceLid::resize(Glossary& opt)
{
    for ( int d = 0; d < DIM; ++d )
    {
        real len = length_[d];
        if ( opt.set(len, "length", d) )
            len *= 0.5;
        if ( len <= 0 )
            throw InvalidParameter("lid:length_[] must be > 0");
        length_[d] = len;
    }
    
    opt.set(top_, "ceiling");
    if ( top_ < 0 )
        throw InvalidParameter("lid:ceiling must be >= 0");
}


void SpaceLid::setModulo(Modulo& mod) const
{
    mod.enable(0, length_[0]);
}

Vector SpaceLid::extension() const
{
    return Vector(length_[0], length_[1], length_[2]);
}


/**
 place only at upper boundary. This overrides the function in spaces
 */
Vector SpaceLid::randomPlaceNearEdge(real radius, unsigned long) const
{
    return Vector( RNG.sfloat()*length_[0], top_, 0 );
}


//------------------------------------------------------------------------------
#pragma mark -


real SpaceLid::volume() const
{
#if ( DIM == 1 )
    return top_ + length_[0];
#elif ( DIM == 2 )
    return 2.0 * length_[0] * ( top_ + length_[1] );
#else
    return 4.0 * length_[0] * length_[1] * ( top_ + length_[2] );
#endif
}


bool  SpaceLid::inside(Vector const& point) const
{
    return ( point[DIM-1] < top_ );
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
                              Meca & meca, real stiff) const
{
    Matrix::index_t inx = DIM-1 + DIM * pe.matIndex();
    
    meca.mC(inx, inx) -= stiff;
    meca.base(inx)    += stiff * top_;
    
#if ( DIM == 2 )
    force_ += stiff * ( pos.YY - top_ );
#elif ( DIM > 2 )
    force_ += stiff * ( pos.ZZ - top_ );
#endif
}


void SpaceLid::setInteraction(Vector const& pos, Mecapoint const& pe, real rad,
                              Meca & meca, real stiff) const
{
    Matrix::index_t inx = DIM-1 + DIM * pe.matIndex();
    
    meca.mC(inx, inx) -= stiff;
    meca.base(inx)    += stiff * ( top_ - rad );
    
#if ( DIM == 2 )
    force_ += stiff * ( pos.YY - top_ + rad );
#elif ( DIM > 2 )
    force_ += stiff * ( pos.ZZ - top_ + rad );
#endif
}

void SpaceLid::setInteractions(Meca& meca, FiberSet const&) const
{
    force_ = 0;
}

void SpaceLid::step()
{
    real dc = prop->mobility_dt * force_;
    
    if ( fabs(dc) < 1 )
        top_ += dc;
    else
        std::cerr << "Error: lid displacement is too fast: " << dc << '\n';
    std::cerr << "force on lid is " << force_ << '\n';
    
    if ( top_ > length_[DIM-1] )
        std::cerr << "Warning: space lid has reached its maximum\n";
}


//------------------------------------------------------------------------------

void SpaceLid::write(Outputter& out) const
{
    out.put_line(" "+prop->shape+" ");
    out.writeUInt16(3);
    out.writeFloat(length_[0]);
    out.writeFloat(length_[1]);
    out.writeFloat(length_[2]);
    out.writeFloat(top_);
    out.writeFloat(force_);
}


void SpaceLid::setLengths(const real len[])
{
    length_[0] = len[0];
    length_[1] = len[1];
    length_[2] = len[2];
    top_       = len[3];
    force_     = len[4];
}

void SpaceLid::read(Inputter& in, Simul&, ObjectTag)
{
    real len[8] = { 0 };
    read_data(in, len);
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
    const real X = length_[0];
    const real Y = ( DIM > 1 ) ? length_[1] : 1;
    const real Z = ( DIM > 2 ) ? length_[2] : 0;

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

