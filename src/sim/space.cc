// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include "dim.h"
#include "space.h"
#include "space_prop.h"
#include "exceptions.h"
#include "interpolation.h"
#include "messages.h"
#include "iowrapper.h"
#include "meca.h"


Space::Space(SpaceProp const* p) 
: prop(p)
{
    assert_true(prop);
}


Space::~Space()
{
    //std::clog << "~Space(" << prop->name() << ")\n";
    prop = nullptr;
}

//------------------------------------------------------------------------------
#pragma mark - Random Places

/**
 Provide a uniform random distribution in the volume by Monte-Carlo.
 
 Method: throws a point in the rectangular volume provided by boundaries()
 until inside() returns true. bailing out after max_trials attempts
*/
Vector Space::place() const
{
    size_t ouf = 0, max_trials = 1<<14;
    Vector res, inf, sup;
    boundaries(inf, sup);
    Vector dif = sup - inf;
    
    do {
        res = inf + dif.e_mul(Vector::randP());
        if ( ++ouf > max_trials )
        {
            throw InvalidParameter("random placement failed for space `"+prop->name()+"'");
            //Cytosim::warn("random placement failed for space `",prop->name(),"'\n");
            return Vector(0,0,0);
        }
    } while ( ! inside(res) );
    
    return res;
}


/**
 Return a `point` for which:
 - inside(point) = true
 - inside(point, radius) = false
 */
Vector Space::placeNearEdge(real rad, size_t max_trials) const
{
    size_t ouf = 0;
    Vector res;
    do {
        res = place();
        assert_true( inside(res) );
        if ( ++ouf > max_trials )
            throw InvalidParameter("edge placement failed for space `"+prop->name()+"'");
    } while ( allInside(res, rad) );
    return res;
}


/**
 Return a random point on the edge, using this method:
 - toss a random point `pos` within the range extended by `rad`.
 - project `pos` on the edge
 - return projection if the distance to `pos` is less than `rad`
 .
 This does not give a uniform distribution on the surface, and accumulation can
 occur in regions of high curvature, since more volume projects on the surface.
 */
Vector Space::onSurface(real rad, size_t max_trials) const
{
    size_t ouf = 0;
    real D = abs_real(rad);
    real RR = rad * rad;
    Vector pos, res, inf, dif;
    
    boundaries(inf, dif);
    inf -= Vector(D, D, D);
    dif += Vector(D, D, D) - inf;
    
    do {
        pos = inf + dif.e_mul(Vector::randP());
        res = project(pos);
        D = ( pos - res ).normSqr();
        if ( ++ouf > max_trials )
            throw InvalidParameter("surface placement failed for `"+prop->name()+"'");
    } while ( D > RR );
    
    return res;
}

//------------------------------------------------------------------------------
#pragma mark - Inside/Outside

/**
 A bead is entirely inside if:
 - its center is inside,
 - the minimal distance (center-to-edge) is greater than the radius
 .
 */
bool Space::allInside(Vector const& cen, const real rad) const
{
    if ( ! inside(cen) )
        return false;

    return ( distanceToEdgeSqr(cen) >= rad * rad );
}

/**
 A bead is entirely outside if:
 - its center is outside,
 - the minimal distance (center-to-edge) is greater than the radius
 .
 
 Attention: this is not equivalent to !allInside(center, radius)
 */
bool Space::allOutside(Vector const& cen, const real rad) const
{
    if ( inside(cen) )
        return false;
        
    return ( distanceToEdgeSqr(cen) >= rad * rad );
}

//------------------------------------------------------------------------------
#pragma mark - Project

/**
this code is equivalent to SpaceInflate::project(), with a negative radius
 */
Vector Space::projectDeflated(Vector const& pos, const real rad) const
{
    if ( rad < 0 )
        ABORT_NOW("radius should not be negative");

    Vector prj = project(pos);
    
    ///\todo problem in project() with radius if point is exactly on the box (n==0)
    //if (n==0) we do not know the orthogonal direction to follow. We should
    //take another point near by, and project from there.
    
    Vector dif = pos - prj;
    real n = dif.normSqr();
    
    if ( n > 0 )
        n = ( inside(pos) ? +rad : -rad ) / std::sqrt(n);
    else {
        throw Exception("in project(..., radius): the point is on the edge");
        //printf("point % .3f % .3f % .3f :", pos[0], pos[1], pos[2]);
        //printf("inside = %i :", inside(point));
        //printf("proj  % .3f % .3f % .3f\n", prj[0], prj[1], prj[2]);
    }
    
    return prj + n * dif;
}

//------------------------------------------------------------------------------
#pragma mark - Misc


real Space::maxExtension() const
{
    Vector inf, sup;
    boundaries(inf, sup);
    return std::max(inf.norm_inf(), sup.norm_inf());
}

/**
 The volume is estimated with a simple Monte-Carlo approach:
 - throw points in the rectangular volume provided by boundaries()
 - count how many are inside the volume with inside()
 .
 Then
 
     volume ~ ( number-of-points-inside / number-of-point ) * volume-of-rectangle
 
 */
real Space::estimateVolume(size_t cnt) const
{
    Vector inf, sup;
    boundaries(inf, sup);
    Vector dif = sup - inf;
    
    size_t in = 0;
    for ( size_t i = 0; i < cnt; ++i )
    {
        Vector pos = inf + dif.e_mul(Vector::randP());
        in += inside(pos);
    }
    
    real vol = real(in) / real(cnt);
    for ( size_t d = 0; d < DIM; ++d )
        vol *= dif[d];

    return vol;
}


/**
 This uses `Space::project()` to reflect `w` on the edge of the Space,
 until the result eventually falls inside.
 
 In most geometries, this works well, but if the distance from the point
 to the edge is very large compared to the width of the space, the number
 of iterations may be large.
 
 In general, this does not satisfies detailed balance, and may lead to artificial
 accumulations in certain part of space, such as the center of a sphere.
*/
Vector Space::bounceOnEdges(Vector const& pos) const
{
    Vector P = pos;
    // bounce on the edge, and return if inside
    int cnt = 0;
    do {
        P = 2*project(P) - P;
        if ( inside(P) )
            return P;
    } while ( ++cnt < 8 );
    
    Vector Q, R;
    do {
        Q = project(P);
        P = 2*Q - P;
        if ( inside(P) )
            return P;
        R = project(P);
        P = 2*R - P;
        if ( inside(P) )
            return P;
        if ( distanceSqr(Q, R) < REAL_EPSILON )
            return ( Q + R ) * 0.5;
    } while ( ++cnt < 16 );

    static size_t msg = 0;
    if ( ++msg < 16 )
    {
        std::cerr << "Warning: "+prop->name()+":bounce fails?\n";
        do {
            Q = project(P);
            P = 2*Q - P;
            std::cerr << pos << " :  " << cnt << "  " << P << " --proj--> " << Q << '\n';
            if ( inside(P) )
                return P;
        } while ( ++cnt < 24 );
    }
    
    // if space is convex, the midpoint should be inside
    return ( Q + R ) * 0.5;
}


Vector Space::bounce(Vector const& pos) const
{
    if ( !inside(pos) )
        return bounceOnEdges(pos);
    return pos;
}


/** 
 `normalToEdge(Vector const&)` uses an iterative method to find
 the normal to the edge, using Space::project().
 
 If you know for certain that `point[]` is far from the edge,
 the normal can be more directly obtained from the projection:

     project(point, proj);
     normal = normalize( proj - point )
 
 */
Vector Space::normalToEdge(Vector const& pos) const
{
    Vector dir;
    Vector prj = project(pos);

    if ( distanceSqr(pos, prj) > 0.1 )
        dir = pos - prj;
    else
    {
        const real goal = square(1000*REAL_EPSILON);
        Vector P, M;
        real H = 0.5;
        for ( unsigned i = 0; i < 8; ++i )
        {
            H /= 2;
            for ( unsigned j = 0; j < 256; ++j )
            {
                dir = Vector::randU(H);
                if ( inside(prj+dir) != inside(prj-dir) )
                    break;
            }
            for ( unsigned j = 0; j < 32; ++j )
            {
                P = project(prj+dir);
                M = project(prj-dir);
                
                // refine the estimate:
                Vector ref = 0.5 * ( M - P );
                dir = ( dir + ref ).normalized(H);
                
                // check convergence:
                if ( ref.normSqr() < goal )
                    goto finish;
            }
        }
        printf("warning: Space::normalToEdge(%9.3f %9.3f %9.3f) failed\n", pos.x(), pos.y(), pos.z());
        printf("         error = %e at height = %8.4f\n", distance(P, M), H);
    }

finish:
    if ( inside(prj+dir) )
        return -normalize(dir);
    else
        return normalize(dir);
}


//------------------------------------------------------------------------------
real  Space::distanceToEdgeSqr(Vector const& pos) const
{
    return distanceSqr(project(pos), pos);
}


real  Space::signedDistanceToEdge(Vector const& pos) const
{
    if ( inside(pos) )
        return -distanceToEdge(pos);
    else
        return +distanceToEdge(pos);
}


//------------------------------------------------------------------------------
#pragma mark - Interactions

/**
 Call the appropriate interaction from `meca`, to force `mp` to be on the edge of the Space.
 
 This implementation uses `pos` to find the local normal to the edge of the Space.
 and then calls Meca::addPlaneClamp, with the approprimate aguments.
 This generates a friction-less potential centered on the edge.
 */

void Space::setConfinement(Vector const& pos, Mecapoint const& mp, Meca& meca, real stiff) const
{
    Vector prj = project(pos);
    assert_true(prj.valid());
    Vector dir = pos - prj;
    real n = dir.normSqr();
    if ( n > REAL_EPSILON )
        meca.addPlaneClamp(mp, prj, dir, stiff/n);
}

/**
 Call the appropriate interaction from `meca`, to confine `mp`, which is at position `pos`.
 
 The default implementation projects `pos`,
 to calculate the direction of the normal to the edge of the Space,
 and then calls Meca::addPlaneClamp, with the approprimate aguments.
 This generates a friction-less potential centered on the edge.
 */

void Space::setConfinement(Vector const& pos, Mecapoint const& mp, real rad, Meca& meca, real stiff) const
{
    Vector prj = projectDeflated(pos, rad);
    Vector dir = pos - prj;
    real n = dir.normSqr();
    if ( n > 0 )
        meca.addPlaneClamp(mp, prj, dir, stiff/n);
}

#if ( 0 )

/**
 This calls Space::setConfinement(pos, Mecapoint, meca, stiff) twice,
 to generate a force on `pi` (which is at position `pos`) toward the surface.
 */
void Space::setConfinement(Vector const& pos, Interpolation const& pi, Meca& meca, real stiff) const
{
    setConfinement(pos, pi.vertex1(), meca, pi.coef0()*stiff);
    setConfinement(pos, pi.vertex2(), meca, pi.coef1()*stiff);
}


/**
 This will add a force component if:
 - ( conf == inside ) && ! Space::inside(pos)
 - ( conf == outside ) &&  Space::inside(pos)
 - ( conf == surface )
 .
 */
void Space::setConfinement(Interpolation const& pi, Meca& meca, real stiff, Confinement conf) const
{
    if ( conf == CONFINE_ON )
    {
        setConfinement(pi.pos(), pi, meca, stiff);
    }
    else if ( conf == CONFINE_INSIDE )
    {
        Vector pos = pi.pos();
        if ( ! inside(pos) )
            setConfinement(pos, pi, meca, stiff);
    }
    else if ( conf == CONFINE_OUTSIDE )
    {
        Vector pos = pi.pos();
        if ( inside(pos) )
            setConfinement(pos, pi, meca, stiff);
    }
}

#endif


//------------------------------------------------------------------------------
#pragma mark - IO

void Space::writeShape(Outputter& out, std::string const& arg)
{
    out.put_characters(arg, 16);
}


void Space::readShape(Inputter& in, std::string const& expected)
{
    std::string str;
#if BACKWARD_COMPATIBILITY < 52
    if ( in.formatID() < 52 )
        str = in.get_word(); // stored as a space-terminated string
    else
#endif
        str = in.get_characters(16); // stored as 16 characters
    
    // compare with expected layout:
    if ( str != expected && in.formatID() > 56 )
    {
        std::cerr << "Notice: space layout mismatch: found `" << str;
        std::cerr << "' but `" << expected << "' was expected\n";
    }
}


size_t Space::readLengths(Inputter& in, size_t n_len, real len[])
{
    size_t n = 0;
#if BACKWARD_COMPATIBILITY < 43
    if ( in.formatID() < 43 )
        n = in.readUInt8();
    else
#endif
        n = in.readUInt16();
    
    size_t d = 0;
    for ( ; d < std::min(n_len,n); ++d )
        len[d] = in.readFloat();
    for ( ; d < n; ++d )
        in.readFloat();
    return n;
}


size_t Space::readShape(Inputter& in, size_t n_len, real len[], std::string const& expected)
{
    readShape(in, expected);
    return readLengths(in, n_len, len);
}

void Space::write(Outputter& out) const
{
    writeMarker(out, TAG);
    writeShape(out, prop->shape);
    out.writeUInt16(0);
}


void Space::read(Inputter& in, Simul&, ObjectTag)
{
    real len[8] = { 0 };
    readShape(in, 8, len, prop->shape);
    setLengths(len);
}

//------------------------------------------------------------------------------
#pragma mark - Display


#ifdef DISPLAY

#include "gym_draw.h"
#include "gym_flute.h"
#include "gym_flute_dim.h"


void Space::drawSection(int dim, real pos, size_t cnt, float width) const
{
    Vector inf, sup;
    boundaries(inf, sup);
    Vector ppp(pos, pos, pos);
    int xxx = ( dim + 1 ) % DIM;
    int yyy = ( xxx + 1 ) % DIM;
    real ix = inf[xxx], sx = sup[xxx];
    real iy = inf[yyy], sy = sup[yyy];
    real dx = ( sx - ix ) / cnt;
    real dy = ( sy - iy ) / cnt;

    size_t i = 0;
    fluteD* pts = gym::mapBufferVD(4*cnt+4);
    ppp[xxx] = sx;
    for ( real y = iy; y <= sy; y += dy )
    {
        ppp[yyy] = y;
        pts[i++] = project(ppp);
    };
    ppp[yyy] = sy;
    for ( real x = sx; x >= ix; x -= dx )
    {
        ppp[xxx] = x;
        pts[i++] = project(ppp);
    }
    ppp[xxx] = ix;
    for ( real y = sy; y >= iy; y -= dy )
    {
        ppp[yyy] = y;
        pts[i++] = project(ppp);
    };
    ppp[yyy] = iy;
    for ( real x = ix; x <= sx; x += dx )
    {
        ppp[xxx] = x;
        pts[i++] = project(ppp);
    }
    gym::unmapBufferVD();
    gym::drawLineStrip(width, 0, i);
}

#else

void Space::drawSection(int, real, size_t, float width) const
{
    //you will get this output if objects for play were not compiled properly:
    //DISPLAY should be defined on the compiler command, with: -DDISPLAY
    printf("dummy Space::drawSection()");
}

#endif
