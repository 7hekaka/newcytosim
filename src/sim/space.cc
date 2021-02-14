// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

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
 
 Algorithm: throw a point in the rectangular volume provided by boundaries()
 until inside() returns true.
*/
Vector Space::randomPlace() const
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
            //Cytosim::warn << "random placement failed for space `"+prop->name()+"'\n";
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
Vector Space::randomPlaceNearEdge(real rad, size_t nb_trials) const
{
    size_t ouf = 0;
    Vector res;
    do {
        res = randomPlace();
        assert_true( inside(res) );
        if ( ++ouf > nb_trials )
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
 */
Vector Space::randomPlaceOnEdge(real rad, size_t max_trials) const
{
    size_t ouf = 0;
    real D = abs_real(rad), RR = rad * rad;
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


real Space::max_extension() const
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
*/
void Space::bounceOnEdges(Vector& pos) const
{
    Vector P;
    // bounce on the edge, and return if inside
    int cnt = 0;
    do {
        P = project(pos);
        pos = 2*P - pos;
        if ( inside(pos) )
            return;
    } while ( ++cnt < 8 );
    
    Vector Q;
    do {
        Q = project(pos);
        pos = 2*Q - pos;
        if ( inside(pos) )
            return;
        if ( distanceSqr(P, Q) < REAL_EPSILON )
        {
            pos = ( P + Q ) * 0.5;
            return;
        }
        P = project(pos);
        pos = 2*P - pos;
        if ( inside(pos) )
            return;
    } while ( ++cnt < 16 );

    static size_t msg = 0;
    if ( ++msg < 16 )
    {
        std::cerr << "Warning: "+prop->name()+":bounce fails?\n";
        do {
            P = project(pos);
            pos = 2*P - pos;
            std::cerr << cnt << "  " << pos << " proj " << P << '\n';
            if ( inside(pos) )
                return;
        } while ( ++cnt < 24 );
    }
    
    // if space is convex, the midpoint should be inside
    pos = ( P + Q ) * 0.5;
}


void Space::bounce(Vector& pos) const
{
    if ( !inside(pos) )
        bounceOnEdges(pos);
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
    const real goal = square(100*REAL_EPSILON);
    
    Vector P, M, res;
    Vector prj = project(pos);

    real H = 1;
    for ( unsigned i = 0; i < 12; ++i )
    {
        H /= 2;
        for ( unsigned j = 0; j < 16; ++j )
        {
            //start from a random vector:
            res = Vector::randU(H);
            
            for ( unsigned n = 0; n < 32; ++n )
            {
                P = project(prj+res);
                M = project(prj-res);
                
                // refine the estimate:
                Vector ref = 0.5 * ( M - P );
                res += ref;
                
                // check convergence:
                if ( ref.normSqr() < goal )
                {
                    if ( 2 * res.norm() < H )
                        res.normalize(H);
                    else
                    {
                        if ( inside(prj+res) )
                            return -normalize(res);
                        else
                            return normalize(res);
                    }
                }
            }
        }
    }
    
    printf("warning: normalToEdge() failed\n");
    printf("         error = %e at height = %e\n", distance(P, prj), H);
    if ( inside(prj+res) )
        return -normalize(res);
    else
        return normalize(res);
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
 Call the appropriate interaction from `meca`, to force `pe` to be on the edge of the Space.
 
 This implementation uses `pos` to find the local normal to the edge of the Space.
 and then calls Meca::addPlaneClamp, with the approprimate aguments.
 This generates a friction-less potential centered on the edge.
 */

void Space::setInteraction(Vector const& pos, Mecapoint const& pe, Meca& meca, real stiff) const
{
    Vector prj = project(pos);
    assert_true(prj.valid());
    Vector dir = pos - prj;
    real n = dir.normSqr();
    if ( n > REAL_EPSILON )
        meca.addPlaneClamp(pe, prj, dir, stiff/n);
}

/**
 Call the appropriate interaction from `meca`, to confine `pe`, which is at position `pos`.
 
 The default implementation projects `pos`,
 to calculate the direction of the normal to the edge of the Space,
 and then calls Meca::addPlaneClamp, with the approprimate aguments.
 This generates a friction-less potential centered on the edge.
 */

void Space::setInteraction(Vector const& pos, Mecapoint const& pe, real rad, Meca& meca, real stiff) const
{
    Vector prj = projectDeflated(pos, rad);
    Vector dir = pos - prj;
    real n = dir.normSqr();
    if ( n > 0 )
        meca.addPlaneClamp(pe, prj, dir, stiff/n);
}

#if ( 0 )

/**
 This calls Space::setInteraction(pos, Mecapoint, meca, stiff) twice,
 to generate a force on `pi` (which is at position `pos`) toward the surface.
 */
void Space::setInteraction(Vector const& pos, Interpolation const& pi, Meca& meca, real stiff) const
{
    setInteraction(pos, pi.exact1(), meca, pi.coef0()*stiff);
    setInteraction(pos, pi.exact2(), meca, pi.coef1()*stiff);
}


/**
 This will add a force component if:
 - ( conf == inside ) && ! Space::inside(pos)
 - ( conf == outside ) &&  Space::inside(pos)
 - ( conf == surface )
 .
 */
void Space::setInteraction(Interpolation const& pi, Meca& meca, real stiff, Confinement conf) const
{
    if ( conf == CONFINE_ON )
    {
        setInteraction(pi.pos(), pi, meca, stiff);
    }
    else if ( conf == CONFINE_INSIDE )
    {
        Vector pos = pi.pos();
        if ( ! inside(pos) )
            setInteraction(pos, pi, meca, stiff);
    }
    else if ( conf == CONFINE_OUTSIDE )
    {
        Vector pos = pi.pos();
        if ( inside(pos) )
            setInteraction(pos, pi, meca, stiff);
    }
}

#endif


//------------------------------------------------------------------------------
#pragma mark - IO

void Space::writeShape(Outputter& out, std::string const& arg)
{
    out.put_characters(arg, 16);
}


void Space::readShape(Inputter& in, size_t n_len, real len[], std::string const& expected)
{
#ifdef BACKWARD_COMPATIBILITY
    if ( in.formatID() < 35 )
        return;
    
    if ( in.formatID() < 36 )
    {
        for ( unsigned d = 0; d < 3; ++d )
            len[d] = in.readFloat();
        return;
    }
    
    if ( in.formatID() < 41 )
    {
        unsigned n = in.readUInt8();
        for ( unsigned d = 0; d < n; ++d )
            len[d] = in.readFloat();
        return;
    }
#endif
    
    std::string str;
#ifdef BACKWARD_COMPATIBILITY
    if ( in.formatID() < 52 )
        str = in.get_word(); // stored as a space-terminated string
    else
#endif
        str = in.get_characters(16); // stored as 16 characters
    
    // compare with expected shape:
    if ( str != expected )
    {
#if 0
        InvalidIO e("shape mismatch");
        e << "found `" << str << "' in file when `" << expected << "' was expected";
        throw e;
#else
        std::cerr << "Warning: shape mismatch: ";
        std::cerr << "found `" << str << "' in file when `" << expected << "' was expected\n";
#endif
    }
    
    // read the dimensions:
    size_t n = 0;
#ifdef BACKWARD_COMPATIBILITY
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
}


void Space::write(Outputter& out) const
{
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
#include "gle.h"


void Space::drawSection(int dim, real pos, size_t cnt) const
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
    fluteD* pts = gle::mapVertexBuffer(4*cnt+4);
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
    gle::unmapVertexBuffer();
    glDrawArrays(GL_LINE_STRIP, 0, i);
}

#else

void Space::drawSection(int, real, size_t) const
{
    //you will get this output if objects for play was not compiled properly:
    //DISPLAY should be defined on the compiler command, with: -DDISPLAY
    printf("dummy Space::drawSection()");
}

#endif
