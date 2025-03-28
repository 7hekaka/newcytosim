// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.
#include "space_dice.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "glossary.h"
#include "meca.h"


SpaceDice::SpaceDice(SpaceProp const* p)
: Space(p)
{
    if ( DIM == 1 )
        throw InvalidParameter("dice is not usable in 1D");

    for ( int d = 0; d < 4; ++d )
        half_[d] = 0;
    edge_ = 0;
}


void SpaceDice::resize(Glossary& opt)
{
    real edg = edge_;
    
    opt.set(edg, "edge");
    if ( edg < 0 )
        throw InvalidParameter("dice:edge must be >= 0");
    if ( edg == 0 )
        std::clog << "Warning: dice:edge is unset or null";

    for ( unsigned d = 0; d < DIM; ++d )
    {
        real len = half_[d];
        if ( opt.set(len, "length", d) )
            len *= 0.5;
        if ( len < edg )
            throw InvalidParameter("dice:length[] must be >= 2 * edge");
        half_[d] = len;
    }
    
    edge_ = edg;
    update();
}


/**
 The `dice` is included in the rectangle
 */
void SpaceDice::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(-half_[0],-half_[1],-half_[2]);
    sup.set( half_[0], half_[1], half_[2]);
}


/**
 If `radius==0`, the volume should be the volume of a rectangle
 */
real SpaceDice::volume() const
{
#if ( DIM == 1 )
    return 2 * half_[0];
#elif ( DIM == 2 )
    return 4 * half_[0]*half_[1] + (M_PI-4.0) * square(edge_);
#else
    return 8 * half_[0]*half_[1]*half_[2]
    + (2.0*M_PI-8.0) * ( half_[0] + half_[1] + half_[2] - 3 * edge_ ) * square(edge_)
    + (4.0*M_PI/3.0-8.0) * cube(edge_);
#endif
}


real SpaceDice::surface() const
{
#if ( DIM == 1 )
    return 2;
#elif ( DIM == 2 )
    return 4 * ( half_[0]+half_[1] ) + (2.0*M_PI-8.0) * edge_;
#else
    return 8 * ( half_[0]*half_[1] + half_[0]*half_[2] + half_[1]*half_[2] )
    + (4.0*M_PI-16.0) * ( half_[0] + half_[1] + half_[2] - 3 * edge_ ) * edge_
    + (4.0*M_PI-24.0) * square(edge_);
#endif
}

//------------------------------------------------------------------------------

bool SpaceDice::inside(Vector const& W) const
{
    real dis = 0;
    for ( unsigned d = 0; d < DIM; ++d )
        dis += square(max_real(0, abs_real(W[d]) - half_[d] + edge_));
    return ( dis <= edgeSqr_ );
}


bool SpaceDice::allInside(Vector const& W, real rad) const
{
    assert_true( rad >= 0 );
    real E = max_real(0, edge_-rad);  // remaining edge
    real R = E + rad;  // size reduction

    real dis = 0;
    for ( unsigned d = 0; d < DIM; ++d )
        dis += square(max_real(0, abs_real(W[d]) - max_real(0, half_[d]-R)));
    return ( dis <= E * E );
}

//------------------------------------------------------------------------------

#if ( DIM == 1 )

Vector SpaceDice::project(Vector const& W) const
{
    return Vector(std::copysign(half_[0], W.XX), 0, 0);
}

#else

Vector SpaceDice::project(Vector const& W) const
{
    Vector P(W);
    bool in = true;

    real X = half_[0] - abs_real(W.XX);
    if ( X < edge_ ) { P.XX = std::copysign(half_[0]-edge_, W.XX); in=false; }
    
    real Y = half_[1] - abs_real(W.YY);
    if ( Y < edge_ ) { P.YY = std::copysign(half_[1]-edge_, W.YY); in=false; }

#if ( DIM > 2 )
    real Z = half_[2] - abs_real(W.ZZ);
    if ( Z < edge_ ) { P.ZZ = std::copysign(half_[2]-edge_, W.ZZ); in=false; }
#endif
    
    if ( in )
    {
        // find the dimensionality corresponding to the closest face
#if ( DIM > 2 )
        if ( Z < Y )
        {
            if ( Z < X )
                P.ZZ = std::copysign(half_[2], W.ZZ);
            else
                P.XX = std::copysign(half_[0], W.XX);
        }
        else
#endif
        {
            if ( Y < X )
                P.YY = std::copysign(half_[1], W.YY);
            else
                P.XX = std::copysign(half_[0], W.XX);
        }
        return P;
    }

    //normalize to radius(), and add to p to get the real projection
    real dis = edge_ / norm(W-P);
    return dis * ( W - P ) + P;
}
#endif

//------------------------------------------------------------------------------

/**
 Here 'dim' is the reduced dimension: half_[d] - edge_
 */
void SpaceDice::setConfinement(Vector const& w, Mecapoint const& mp,
                               Meca& meca, real stiff, const real dim[], real edg)
{
#if ( DIM == 1 )
    meca.addPlaneClampX(mp, std::copysign(dim[0], w.XX), stiff);
#else
    real dX = dim[0] - abs_real(w.XX);
    real dY = dim[1] - abs_real(w.YY);
    bool inX = ( dX > edg );
    bool inY = ( dY > edg );
#if ( DIM > 2 )
    real dZ = dim[2] - abs_real(w.ZZ);
    bool inZ = ( dZ > edg );
    if ( inX && inY && inZ )
#else
    if ( inX && inY )
#endif
    {
        // find the dimensionality corresponding to the closest face
#if ( DIM > 2 )
        if ( dZ < dY )
        {
            if ( dZ < dX )
                meca.addPlaneClampZ(mp, std::copysign(dim[2], w.ZZ), stiff);
            else
                meca.addPlaneClampX(mp, std::copysign(dim[0], w.XX), stiff);
        }
        else
#endif
        {
            if ( dY < dX )
                meca.addPlaneClampY(mp, std::copysign(dim[1], w.YY), stiff);
            else
                meca.addPlaneClampX(mp, std::copysign(dim[0], w.XX), stiff);
        }
    }
#if ( DIM > 2 )
    else if ( inY && inZ )
    {
        meca.addPlaneClampX(mp, std::copysign(dim[0], w.XX), stiff);
    }
    else if ( inX && inZ )
    {
        meca.addPlaneClampY(mp, std::copysign(dim[1], w.YY), stiff);
    }
    else if ( inX && inY )
    {
        meca.addPlaneClampZ(mp, std::copysign(dim[2], w.ZZ), stiff);
    }
#endif
    else if ( inX )
    {
#if ( DIM > 2 )
        real cY = std::copysign(dim[1]-edg, w.YY);
        real cZ = std::copysign(dim[2]-edg, w.ZZ);
        meca.addCylinderClamp(mp, Vector(1, 0, 0), Vector(0, cY, cZ), edg, stiff);
#else
        meca.addPlaneClampY(mp, std::copysign(dim[1], w.YY), stiff);
#endif
    }
    else if ( inY )
    {
#if ( DIM > 2 )
        real cX = std::copysign(dim[0]-edg, w.XX);
        real cZ = std::copysign(dim[2]-edg, w.ZZ);
        meca.addCylinderClamp(mp, Vector(0, 1, 0), Vector(cX, 0, cZ), edg, stiff);
#else
        meca.addPlaneClampX(mp, std::copysign(dim[0], w.XX), stiff);
#endif
    }
#if ( DIM > 2 )
    else if ( inZ )
    {
        real cX = std::copysign(dim[0]-edg, w.XX);
        real cY = std::copysign(dim[1]-edg, w.YY);
        meca.addCylinderClamp(mp, Vector(0, 0, 1), Vector(cX, cY, 0), edg, stiff);
    }
#endif
    else
    {
        real cX = std::copysign(dim[0]-edg, w.XX);
        real cY = std::copysign(dim[1]-edg, w.YY);
#if ( DIM > 2 )
        real cZ = std::copysign(dim[2]-edg, w.ZZ);
#else
        real cZ = 0;
#endif
        meca.addSphereClamp(mp, Vector(cX, cY, cZ), edg, stiff);
    }
#endif
}


void SpaceDice::setConfinement(Vector const& pos, Mecapoint const& mp, Meca& meca, real stiff) const
{
    setConfinement(pos, mp, meca, stiff, half_, edge_);
}


void SpaceDice::setConfinement(Vector const& pos, Mecapoint const& mp, real rad, Meca& meca, real stiff) const
{
    real E = max_real(0, edge_-rad);  // remaining edge

    real dim[DIM];
    for ( unsigned d = 0; d < DIM; ++d )
        dim[d] = max_real(0, half_[d] - rad);

    setConfinement(pos, mp, meca, stiff, dim, E);
}

//------------------------------------------------------------------------------

void SpaceDice::write(Outputter& out) const
{
    writeMarker(out, TAG);
    writeShape(out, "LLLE");
    out.writeUInt16(4);
    out.writeFloat(half_[0]);
    out.writeFloat(half_[1]);
    out.writeFloat(half_[2]);
    out.writeFloat(edge_);
}


void SpaceDice::setLengths(const real len[8])
{
    half_[0] = len[0];
    half_[1] = len[1];
    half_[2] = len[2];
    edge_ = len[3];
    update();
}

void SpaceDice::read(Inputter& in, Simul&, ObjectTag)
{
    real len[8] = { 0 };
    readShape(in, 8, len, "LLLE");
    setLengths(len);
}

//------------------------------------------------------------------------------
#pragma mark - OpenGL display

#ifdef DISPLAY

#include "gle.h"
#include "gym_draw.h"
#include "gym_view.h"
#include "gym_flute.h"
#include "tesselator.h"

void SpaceDice::draw2D(float width) const
{
    const unsigned CNT = 128;
    drawSection(2, 0, CNT, width);
}

void SpaceDice::draw3D() const
{
    const float X(half_[0] - edge_);
    const float Y(half_[1] - edge_);
    const float Z(half_[2] - edge_);
    
#if 1
    Tesselator mesh;
    mesh.buildDice(X, Y, Z, edge_, 2*gle::finesse, 1);
    unsigned cnt = mesh.num_vertices();
    flute3 * fl3 = gym::mapBufferV3(cnt);
    mesh.store_vertices((float*)fl3);
    gym::unmapBufferV3N0();
    //gym::drawPoints(4, 0, cnt);
    // do not copy the last 12 triangles, corresponding to the faces drawn below
    unsigned tri = 3 * ( mesh.num_faces() - 12 );
    unsigned short* inx = gym::mapIndexBuffer(tri);
    memcpy(inx, mesh.face_data(), tri*sizeof(Tesselator::INDEX));
    gym::unmapIndexBuffer();

    static_assert(std::is_same<Tesselator::INDEX, GLushort>::value, "Index type mismatch");
    glDrawElements(GL_TRIANGLES, tri, GL_UNSIGNED_SHORT, nullptr);
    gym::unbind2();
#endif

#if 1
    // this gets the normals right on the planar surfaces
    const float XR(half_[0]);
    const float YR(half_[1]);
    const float ZR(half_[2]);
    flute6* flu = gym::mapBufferV3N3(36);
    gle::setExplodedCube(flu, X, Y, Z, XR, YR, ZR);
    gym::unmapBufferV3N3();
    gym::drawTriangles(0, 36);
#endif
    gym::cleanupVN();
}

#else

void SpaceDice::draw2D(float) const {}
void SpaceDice::draw3D() const {}

#endif


