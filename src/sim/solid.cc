// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University.

#include "assert_macro.h"
#include "primitives.h"
#include "random_vector.h"
#include "solid.h"
#include "solid_prop.h"
#include "exceptions.h"
#include "messages.h"
#include "hand_prop.h"
#include "iowrapper.h"
#include "tokenizer.h"
#include "glossary.h"
#include "meca.h"
#include "simul.h"
#include "space.h"
#include "singles/wrist.h"
#include "lapack.h"
#include "quaternion.h"
#include "matrix33.h"

//------------------------------------------------------------------------------
#pragma mark - Step


void Solid::step()
{
}


void Solid::setInteractions(Meca& meca) const
{
#if NEW_RADIAL_FLOW
    LOG_ONCE("NEW_RADIAL_FLOW enabled: Solids converge to the same point\n");
    /// Special code for Maria Burdyniuk
    double now = simul().time();
    if ( prop->flow_time[0] > now )
    {
        Mecapoint pt(this,0);
        Vector dir = prop->flow_center - pt.pos();
        real s = dragCoefficient() / ( prop->flow_time[1] - now );
        meca.addForce(pt, dir * s);
    }
#endif
#if NEW_SOLID_CLAMP
    if ( clamp_stiff > 0 )
    {
        // this attaches to first point of Solid:
        meca.addPointClamp(Mecapoint(this,0), clamp_place, clamp_stiff);
    }
#endif
    if ( prop->confine != CONFINE_OFF )
    {
        Space const* spc = prop->confine_space;
        assert_true(spc);

        switch ( prop->confine )
        {
            case CONFINE_INSIDE:
                for ( index_t i = 0; i < nPoints; ++i )
                {
                    // confine all massive points:
                    if ( radius(i) > 0 )
                    {
                        Vector pos = posP(i);
                        if ( ! spc->inside(pos) )
                            spc->setConfinement(pos, Mecapoint(this, i), meca, prop->confine_stiff);
                    }
                }
                break;
                
            case CONFINE_OUTSIDE:
                for ( index_t i = 0; i < nPoints; ++i )
                {
                    // confine all massive points:
                    if ( radius(i) > 0 )
                    {
                        Vector pos = posP(i);
                        if ( spc->inside(pos) )
                            spc->setConfinement(pos, Mecapoint(this, i), meca, prop->confine_stiff);
                    }
                }
                break;
                
            case CONFINE_ALL_INSIDE:
                for ( index_t i = 0; i < nPoints; ++i )
                {
                    const real rad = radius(i);
                    // confine all massive points:
                    if ( rad > 0 )
                    {
                        Vector pos = posP(i);
                        if ( ! spc->allInside(pos, rad) )
                            spc->setConfinement(pos, Mecapoint(this, i), rad, meca, prop->confine_stiff);
                    }
                }
                break;
                
            case CONFINE_ON:
                for ( index_t i = 0; i < nPoints; ++i )
                {
                    // only confine massive points:
                    if ( radius(i) > 0 )
                        spc->setConfinement(posP(i), Mecapoint(this, i), meca, prop->confine_stiff);
                }
                break;
                
            case CONFINE_POINT:
                spc->setConfinement(posP(0), Mecapoint(this, 0), meca, prop->confine_stiff);
                break;
                
            case CONFINE_POINT_INSIDE:
            {
                Vector pos = posP(0);
                if ( !spc->inside(pos) )
                    spc->setConfinement(pos, Mecapoint(this, 0), meca, prop->confine_stiff);
            } break;

            default:
                throw InvalidParameter("Invalid solid::confine");
        }
    }
}

//------------------------------------------------------------------------------
#pragma mark - Construction

void Solid::reset()
{
    soRadius = nullptr;
    soShape  = nullptr;
    soAmount = 0;
    soDrag   = 0;
#if ( DIM > 2 )
    soMomentum = Matrix33(0, 1);
#endif
#if NEW_SOLID_CLAMP
    clamp_place.reset();
    clamp_stiff = 0;
#endif
}


Solid::Solid()
{
    reset();
    prop = nullptr;
}


Solid::Solid(SolidProp const* p)
: prop(p)
{
    reset();
}


Solid::Solid(const Solid & o)
: Mecable(o)
{
    reset();
    prop = o.prop;
    allocateMecable(o.nPoints);
    for ( index_t p = 0; p < nPoints; ++p )
        soRadius[p] = o.soRadius[p];
    fixShape();
}


Solid & Solid::operator = (const Solid & o)
{
    reset();
    prop = o.prop;
    allocateMecable(o.nPoints);
    Mecable::operator = (o);
    for ( index_t p = 0; p < nPoints; ++p )
        soRadius[p] = o.soRadius[p];
    fixShape();
    return *this;
}


Solid::~Solid()
{
    if ( objset() )
        simul().singles.deleteWrists(this);
    release();
    prop = nullptr;
}


/**
 This extends Mecable::allocateMecable().
 */
void Solid::allocateMecable(const index_t nbp)
{
    real * ptr = Mecable::allocateMemory(nbp, 1+DIM);
    if ( ptr )
    {
        index_t all = allocated();
        //std::clog << "Solid::allocateMecable " << all << '\n';
        soRadius = ptr;
        soShape = ptr + all;
        //print(std::clog, 0);
    }
}


void Solid::release()
{
    soRadius = nullptr;
    soShape = nullptr;
}


//------------------------------------------------------------------------------
#pragma mark - Build


index_t Solid::makePoint(ObjectList& objs, Glossary& opt, std::string const& var, Simul& sim)
{
    std::string str;
    index_t inx = 0;
    unsigned nbp = 1;
    // optionally specify a number of points
    if ( opt.set_positive_integer(nbp, var) )
        ++inx;
    
    if ( nbp > 0 )
    {
        // get sphere radius:
        real rad = 0;
        opt.set(rad, var, inx+1);
        
        index_t fip = nPoints;
        str = opt.value(var, inx);
        // add 'nbp' points:
        for ( index_t n = 0; n < nbp; ++n )
        {
            Vector vec = Cytosim::findPosition(str, nullptr);
            addSphere(vec, rad);
        }
        
        // attach Single to this set of points:
        ++inx;
        while ( opt.set(str, var, ++inx) )
            sim.singles.makeWrists(objs, this, fip, nbp, str);
    }
    return nbp;
}

/**
 single? = SINGLE, POINT1, POINT2, COEFFICIENT
 */
void Solid::makeWrist(ObjectList& objs, Glossary& opt, std::string const& var, Simul& sim)
{
    index_t num = 1;
    std::string str;
    index_t a = 0, b = 0;
    real c = 0.0;
    
    opt.set(str, var, 0);
    Tokenizer::split_integer(num, str);
    SingleProp const* sip = sim.findProperty<SingleProp>("single", str);

    // get index of point A
    opt.set(str, var, 1);
    b = a = point_index(str);

    // get index of point B
    if ( opt.set(str, var, 2) )
    {
        b = point_index(str);
        
        // get interpolation coefficient
        opt.set(c, var, 3);
        if ( c < 0 || 1 < c )
            throw InvalidParameter("interpolation coefficient must be in [0, 1]");
    }
    
    for ( index_t i = 0; i < num; ++i )
    {
        // add a Wrist anchored between 'a' and 'b':
        Wrist * w = sip->newWrist(this, a);
        if ( b != a )
            w->rebase(this, a, b, c);
        objs.push_back(w);
    }
}


/*
 add Wrists anchored on the local coordinate system of a sphere at index 'ref':
 using unit vectors here since the Triad is build already with a scale 'rad'
 */
void Solid::addWrists(ObjectList& objs, index_t num, SingleProp const* sip, index_t ref, std::string const& spec)
{
    for ( index_t i = 0; i < num; ++i )
    {
        Vector vec(0,0,0);
        try {
            vec = Cytosim::findPosition(spec, nullptr);
        } catch( InvalidParameter& e ) {
            //e << " in `" << spec << "'\n";
            throw;
        }
        Wrist * w = sip->newWrist(this, ref, vec);
        objs.push_back(w);
    }
}


index_t Solid::makeBall(ObjectList& objs, Glossary& opt, std::string const& var, Simul& sim)
{
    std::string str;
    // get sphere radius:
    real rad = 0;
    if ( opt.set_block(str, '[', var, 1) ) // some code specified
    {
        // get a range of radius
        float a = 0, b = 0;
        if ( 2 != sscanf(str.c_str(), "%f, %f", &a, &b) )
            throw InvalidParameter("expected range ([REAL, REAL]) in Solid's radius");
        rad = RNG.real_uniform(a, b);
    }
    else
        opt.set(rad, var, 1);
    if ( rad <= 0 )
        throw InvalidParameter("radius of solid:ball must be > 0");

    // get position of center:
    Vector cen = Cytosim::findPosition(opt.value(var, 0), nullptr);
    // add a bead with a local coordinate system
    index_t ref = addSphere(cen, rad);
    addTriad(rad);
    rad = abs_real(rad);

#if ( DIM > 2 )
    real sep = 1.0;
    if ( opt.set(sep, "separation") )
    {
        // count number of singles to be attached:
        index_t nbs = 0;
        index_t inx = 2;
        while ( opt.peek(str, var, inx++) )
        {
            index_t num = 1;
            Tokenizer::split_integer(num, str);
            nbs += num;
        }
        std::vector<Vector> pts(nbs, Vector(0,0,0));
        real sup = 2 * std::pow(M_PI/(3*M_SQRT2*nbs), 1.0/3);
        index_t nbp = tossPointsBall(pts, std::min(sep/rad, sup), 128);
        if ( nbp < nbs )
        {
            Cytosim::warn("Could only fit ", nbp, "/", nbs, " points in Sphere at separation ", sep, "\n");
            Cytosim::warn(" Target separation for ", nbs, " points is ", sup * rad, "\n");
        }
        nbs = 0;
        inx = 2;
        while ( opt.set(str, var, inx++) )
        {
            index_t num = 1;
            Tokenizer::split_integer(num, str);
            SingleProp const* sip = sim.findProperty<SingleProp>("single", str);
            for ( index_t u = 0; u < num; ++u )
            {
                Wrist * w = sip->newWrist(this, ref, pts[nbs]);
                nbs = ( nbs + 1 ) % nbp;
                objs.push_back(w);
            }
        }
    }
    else
#endif
    {
        index_t inx = 2;
        while ( opt.set(str, var, inx++) )
        {
            index_t num = 1;
            Tokenizer::split_integer(num, str);
            std::string nam = Tokenizer::split_symbol(str);
            if ( nam.empty() )
                throw InvalidParameter("the name of a single should be specified in `"+str+"'");
            SingleProp const* sip = sim.findProperty<SingleProp>("single", nam);
            if ( str.size() )
                addWrists(objs, num, sip, ref, str);
            else
            {
                for ( index_t i = 0; i < num; ++i )
                {
                    Wrist * w = sip->newWrist(this, ref, Vector::randB());
                    objs.push_back(w);
                }
            }
        }
    }
    return ref;
}


index_t Solid::makeSphere(ObjectList& objs, Glossary& opt, std::string const& var, Simul& sim)
{
    std::string str;
    real rad = 0;
    if ( opt.set_block(str, '[', var, 1) ) // some code specified
    {
        // get a range of radius
        float a = 0, b = 0;
        if ( 2 != sscanf(str.c_str(), "%f, %f", &a, &b) )
            throw InvalidParameter("expected range ([REAL, REAL]) in Solid's radius");
        rad = RNG.real_uniform(a, b);
    }
    else
        opt.set(rad, var, 1);

    // get position of center:
    Vector cen = Cytosim::findPosition(opt.value(var, 0), nullptr);
    // add a bead with a local coordinate system
    index_t ref = addSphere(cen, rad);
    addTriad(rad);
    rad = abs_real(rad);

#if ( DIM > 2 )
    real sep = 1.0, dev = 0.0;
    if ( opt.set(sep, "separation") && opt.set(dev, "deviation") )
    {
        if ( dev > rad )
            throw InvalidParameter("solid:deviation should be <= radius");
        // attach Single on the surface of this sphere:
        size_t nbs = opt.num_values(var) - 2;
        std::vector<Vector> pts(nbs, Vector(0,0,0));
        distributePointsSphere(pts, sep/rad, 128);
        index_t inx = 2;
        while ( opt.set(str, var, inx++) )
        {
            // get a number and the name of a class:
            index_t num = 1;
            Tokenizer::split_integer(num, str);
            SingleProp const* sip = sim.findProperty<SingleProp>("single", str);
            for ( index_t u = 0; u < num; ++u )
            {
                Vector pos = pts[inx-3];
                Vector vec = normalize(pos+pos.randOrthoB(dev/rad));
                Wrist * w = sip->newWrist(this, ref, vec);
                objs.push_back(w);
            }
        }
    }
    else if ( opt.set(sep, "separation") )
    {
        // count number of singles to be attached:
        index_t nbs = 0;
        index_t inx = 2;
        while ( opt.peek(str, var, inx++) )
        {
            index_t num = 1;
            Tokenizer::split_integer(num, str);
            nbs += num;
        }
        // generate a set of unit vectors:
        std::vector<Vector> pts(nbs, Vector(0,0,0));
        distributePointsSphere(pts, sep/rad, 128);
        nbs = 0;
        inx = 2;
        while ( opt.set(str, var, inx++) )
        {
            index_t num = 1;
            Tokenizer::split_integer(num, str);
            SingleProp const* sip = sim.findProperty<SingleProp>("single", str);
            for ( index_t u = 0; u < num; ++u )
            {
                Wrist * w = sip->newWrist(this, ref, pts[nbs++]);
                objs.push_back(w);
            }
        }
    }
    else
#endif
    {
        index_t inx = 2;
        while ( opt.set(str, var, inx++) )
        {
            index_t num = 1;
            Tokenizer::split_integer(num, str);
            std::string nam = Tokenizer::split_symbol(str);
            //std::clog << num << " [" << nam << "] |" << str << "|\n";
            if ( nam.empty() )
                throw InvalidParameter("the name of a single should be specified in `"+str+"'");
            SingleProp const* sip = sim.findProperty<SingleProp>("single", nam);
            if ( str.size() )
                addWrists(objs, num, sip, ref, str);
            else
            {
                for ( index_t i = 0; i < num; ++i )
                {
                    Wrist * w = sip->newWrist(this, ref, Vector::randU());
                    objs.push_back(w);
                }
            }
        }
    }
    return ref;
}


Fiber* Solid::makeFiber(ObjectList& objs, Glossary& opt, std::string const& var, Simul& sim)
{
    index_t ref = 0;
    const real rad = radius(ref);
    Vector A, B;
    std::string str;
    
    opt.set(str, var, 0);
    FiberProp const* fip = sim.findProperty<FiberProp>("fiber", str);
    
    opt.set(str, var, 2);
    SingleProp const* sip = sim.findProperty<SingleProp>("single", str);

    // get position of anchoring points
    if ( !opt.set(A, var, 3) || !opt.set(B, var, 4) )
        throw InvalidParameter("points must be specified as fiber?[1] and fiber?[2]");
    
    Wrist * w1 = sip->newWrist(this, ref, A/rad);
    objs.push_back(w1);
    
    Wrist * w2 = sip->newWrist(this, ref, B/rad);
    objs.push_back(w2);

    real len = distance(A, B);
    if ( len < REAL_EPSILON )
        throw InvalidParameter("anchor points must be distinct");
    
    if ( !opt.set(str, var, 1) )
        throw InvalidParameter("fiber specs must be specified as fiber?[1]");
    
    ObjectList list;
    Fiber * F = sim.fibers.newFiber(list, fip, str);
    objs.append(list);

    Vector dir = ( B - A ) / len;
    if ( opt.set(str, var, 5) && str=="plus_end" )
    {
        ObjectSet::rotateObjects(list, Rotation::rotationToVector(-dir));
        ObjectSet::translateObjects(list, A - F->posEndP());
        w1->attach(FiberSite(F, F->abscissaP()));
        w2->attach(FiberSite(F, F->abscissaP()-len));
    }
    else
    {
        ObjectSet::rotateObjects(list, Rotation::rotationToVector(dir));
        ObjectSet::translateObjects(list, A - F->posEndM());
        w1->attach(FiberSite(F, F->abscissaM()));
        w2->attach(FiberSite(F, F->abscissaM()+len));
    }
    return F;
}


/**
 @ingroup NewObject
 
 There are different ways to specify the number and positions of points in a Solid:
 
     new solid NAME
     {
       point1 = [INTEGER,] POSITION, RADIUS [, SINGLE]
       point2 = [INTEGER,] POSITION, RADIUS [, SINGLE]
       point3 = [INTEGER,] POSITION, RADIUS [, SINGLE]
       etc.
     }
 
 each `point#` specifies a number of points to be added.
 The first parameter (`INTEGER`) specifies the number of points.
 The second argument (`POSITION`) specifies their position with respect to the center.
 The keywords are the same as for other position in cytosim (see examples below).
 The last argument (`RADIUS`) specifies the radius of the bead attached at this point,
 and it can be zero.
 
 Examples:
 
     new solid blob
     {
       point1 = center, 1.0
       point2 = 10, sphere 1, 0, grafted
       ...
     }
 
 `POSITION` can be a `VECTOR`, or the usual keywords:
 - `center`
 - `ball RADIUS`
 - `sphere RADIUS`
 - `equator RADIUS`
 .
 
 Another way to specify points of a Solid:

     new solid NAME
     {
       sphere1 = POSITION, RADIUS [, SINGLE]
       sphere2 = POSITION, RADIUS [, SINGLE]
       etc.
     }
 
 each `sphere#` specifies one sphere to be added.
 The first argument (`POSITION`) specifies the position with respect to the center.
 The keywords are the same as for other position in cytosim (see examples below).
 The second argument (`RADIUS`) specifies the radius of the bead attached at this point,
 and it should not be zero.
 

 <h3> Add Singles to a Solid </h3>
 
 Extra parameters can be used to add Single to the points of a Solid:
 
     new solid NAME
     {
       point1  = ... , SINGLE
       sphere1 = ... , SINGLE
       etc.
       anchor  = SINGLE [, SINGLE] ...
       etc.
     }
 
 Where `SINGLE` is string containing at most 3 words: `[INTEGER] NAME [each]`,
 where the `INTEGER` specifies the number of Singles, `NAME` specifies their name,
 and the optional word `each` species that the command applies to every point.
 
 The command `attach` applies to all the points of the Solid, while `attach1`,
 `attach2`, etc. apply to the points specified by `point1`, `point2`, etc. only.
 With `attach`, the Singles are distributed randomly on all the points,
 and if `each` is specified, the specification is repeated for each point.
 
 For example if `grafted` is the name of a Single, one can use:

     new solid NAME
     {
        point1 = center, 1, grafted
        sphere1 = 1 0 0, 1, 7 grafted
     }
 */

ObjectList Solid::build(Glossary& opt, Simul& sim)
{
    ObjectList objs;
    std::string str;
    
#if NEW_SOLID_CLAMP
    if ( opt.set(str, "clamp") )
    {
        if ( str == "position" )
            opt.set(clamp_place, "position");
        else
            clamp_place = Cytosim::readPosition(str, prop->confine_space);
        //std::clog << reference() << " clamped at " << clamp_place << "\n";
    }
    opt.set(clamp_stiff, "clamp", 1) || opt.set(clamp_stiff, "clamp_stiffness");
    if ( clamp_stiff < 0 )
        throw InvalidParameter("clamp[0] (stiffness) should be >= 0");
#endif
    
    if ( opt.has_key("point0") || opt.has_key("sphere0") )
        throw InvalidParameter("point indices start at 1 (use `point1`, `point2`, etc.)");
    
    bool more = true;
    std::string var = "";
    for ( index_t inp = 1; more; ++inp)
    {
        more = false;
        // parameter 'sphere###' will add a sphere, with Singles on their surface
        var = "sphere" + std::to_string(inp);
        if ( opt.has_key(var) )
        {
            index_t ref = makeSphere(objs, opt, var, sim);
            if ( opt.has_key("twin") )
                rotateTriad(ref, Rotation::align111());
            more = true;
        }
        
        // parameter 'ball###' will add a sphere, with Singles in their volume
        var = "ball" + std::to_string(inp);
        if ( opt.has_key(var) )
        {
            makeBall(objs, opt, var, sim);
            more = true;
        }
        
        // parameter 'point###' will add a point, with Single attached
        var = "point" + std::to_string(inp);
        if ( opt.has_key(var) )
        {
            makePoint(objs, opt, var, sim);
            more = true;
        }
    }
         
#if BACKWARD_COMPATIBILITY < 100
    /* distribute Singles over all points. Deprecated, since 03.2017 */
    size_t inp = 0;
    while ( opt.set(str, "anchor", inp++) )
        sim.singles.makeWrists(objs, this, 0, nPoints, str);
#endif
#if BACKWARD_COMPATIBILITY < 100
    /* Anchor Single interpolated between two vertices, deprecated 24.11.2022 */
    inp = 1;
    var = "anchor1";
    while ( opt.has_key(var) )
    {
        std::cerr << "please use `single1` instead of removed `anchor1`\n";
        var = "anchor" + std::to_string(++inp);
    }
#endif
    /* Anchor Single to intermediate positions between two vertices (24.11.2022) */
    inp = 1;
    var = "single1";
    while ( opt.has_key(var) )
    {
        makeWrist(objs, opt, var, sim);
        var = "single" + std::to_string(++inp);
    }
    /* Attach fibers using two interpolated anchors (24.11.2022) */
    /* Except that if 'fibers' is specified, this is building an Aster */
    if ( !opt.has_key("fibers") )
    {
        inp = 1;
        var = "fiber1";
        while ( opt.has_key(var) )
        {
            makeFiber(objs, opt, var, sim);
            var = "fiber" + std::to_string(++inp);
        }
    }
    // verify the number of points:
    inp = 0;
    if ( opt.set(inp, "nb_points")  &&  inp != nPoints )
    {
        throw InvalidParameter("could not find the number of points specified in solid:nb_points");
    }

    objs.push_back(this);
    
    return objs;
}


//------------------------------------------------------------------------------
#pragma mark - Shape

index_t Solid::addSphere(Vector const& vec, real rad)
{
    index_t inx = addPoint(vec);
    soRadius[inx] = rad;
    //std::clog << "addSphere(" << vec << ", " << rad << ") for " << reference() << " index " << inx << "\n";
    return inx;
}

// the coordinate system can be inverted by specifying negative 'arm'
index_t Solid::addTriad(real arm)
{
    if ( nPoints < 1 )
        throw InvalidParameter("cannot add Triad to a Solid with no point");
    
    index_t inx = lastPoint();
    real ARM = abs_real(arm);
    
    //std::clog << "Solid::addTriad(" << arm << ") at index " << inx << "\n";
    Vector vec = posP(inx);
    
    if ( DIM > 0 ) addPoint(vec+Vector(arm,0,0));
    if ( DIM > 1 ) addPoint(vec+Vector(0,ARM,0));
    if ( DIM > 2 ) addPoint(vec+Vector(0,0,ARM));
    
    return inx;
}


void Solid::rotateTriad(index_t ref, Rotation const& rot)
{
    if ( nPoints < ref + DIM )
        throw InvalidParameter("cannot identify Triad at this index");

    Vector vec = posP(ref);
    for ( index_t i = 1; i <= DIM; ++i )
    {
        Vector off = posPoint(ref+i) - vec;
        setPoint(ref+i, vec + rot*off);
    }
}


/** will return the size of the base vector, negated if triad is indirect */
real Solid::hasTriad(index_t inx, real epsilon) const
{
    if ( inx + DIM >= nbPoints() )
        return 0;
    Vector P = posPoint(inx), Y;
    Vector X = posPoint(inx+1) - P;
    real res = X.normSqr();
    if ( DIM >= 2 )
    {
        Y = posPoint(inx+2) - P;
        if ( abs_real(Y.normSqr() - res) > epsilon )
            return 0;
        if ( abs_real(dot(X, Y)) > epsilon )
            return 0;
    }
    if ( DIM >= 3 )
    {
        Vector Z = posPoint(inx+3) - P;
        if ( abs_real(Z.normSqr() - res) > epsilon )
            return 0;
        if ( abs_real(dot(X, Z)) > epsilon )
            return 0;
        if ( abs_real(dot(Y, Z)) > epsilon )
            return 0;
    }
    return res;
}


void Solid::setRadius(const index_t indx, const real rad)
{
    assert_true( indx < nPoints );
    if ( rad < 0 )
        throw InvalidParameter("solid:radius must be positive");
    soRadius[indx] = rad;
}


real Solid::sumRadius() const
{
    real res = 0;
    
    for ( index_t i = 0; i < nPoints; ++i )
        res += radius(i);
    
    return res;
}


/**
 find the 3 nearest spheres to sphere at `inx` on the same Solid
 The result should be independent of the actual position, since the shape
 of a Solid is fixed, and it could be precalculated
 */
index_t Solid::closestSpheres(const index_t inx, index_t& i1, index_t& i2, index_t& i3) const
{
    assert_true( inx < nbPoints() );
    const index_t sup = nbPoints();
    Vector X = posP(inx);
    i1 = sup;
    i2 = sup;
    i3 = sup;
    real d1 = INFINITY;
    real d2 = INFINITY;
    real d3 = INFINITY;
    for ( index_t i = 0; i < sup; ++i )
    {
        if ( i != inx && radius(i) > 0 )
        {
            real d = distanceSqr(X, posP(i));
            if ( d < d1 )
            {
                d3 = d2;
                i3 = i2;
                d2 = d1;
                i2 = i1;
                d1 = d;
                i1 = i;
            }
            else if ( d < d2 )
            {
                d3 = d2;
                i3 = i2;
                d2 = d;
                i2 = i;
            }
            else if ( d < d3 )
            {
                d3 = d;
                i3 = i;
            }
        }
    }
    return ( i1 < sup ) + ( i2 < sup ) + ( i3 < sup );
}


Vector Solid::centroid() const
{
    if ( nPoints == 0 )
        throw InvalidParameter("cannot calculate centroid of a Solid without point");
    
    if ( nPoints == 1 )
        return posP(0);
    
    real sum = 0;
    Vector res(0,0,0);
    for ( index_t i = 0; i < nPoints; ++i )
    {
        real R = radius(i);
        if ( R > 0 )
        {
            res += R * posP(i);
            sum += R;
        }
    }
    if ( sum < REAL_EPSILON )
        throw InvalidParameter("cannot calculate centroid of a Solid without drag sphere");
    
    return res / sum;
}


/**
  
 @return a vector aligned with the principal component of the spheres in the Solid

 */
Vector Solid::orientation() const
{
    real S = 0, M[9] = { 0 };
    
    for ( index_t i = 0; i < nPoints; ++i )
    {
        const real w = radius(i);
        Vector p = posP(i);
            
        M[0] += w * ( DIM * p.XX * p.XX - 1 );
#if ( DIM > 1 )
        M[1] += w * ( DIM * p.YY * p.XX );
        M[4] += w * ( DIM * p.YY * p.YY - 1 );
#endif
#if ( DIM > 2 )
        M[2] += w * ( DIM * p.ZZ * p.XX );
        M[5] += w * ( DIM * p.ZZ * p.YY );
        M[8] += w * ( DIM * p.ZZ * p.ZZ - 1 );
#endif
        S += w;
    }
    
    if ( nPoints < 2 )
        throw InvalidParameter("cannot calculate orientation of a Solid with less than 2 spheres");
    
    int nbv = 0;
    real vec[9] = { 0 };

    if ( S > REAL_EPSILON )
    {
        // rescale matrix:
        for ( index_t d = 0; d < 9; ++d )
            M[d] = M[d] / S;
        
        real val[3] = { 0 };
        real work[32];
        int iwork[16];
        int ifail[4];
        int info = 0;
        
        // find vector with the largest eigenvalue:
        lapack::xsyevx('V','I','L', DIM, M, 3, 0, 0, DIM, DIM, REAL_EPSILON,
                       &nbv, val, vec, 3, work, 32, iwork, ifail, &info);
    }
    if ( nbv > 0 )
        return Vector(vec[0], vec[1], vec[2]);
    return Vector(0, 0, 0);
}

//------------------------------------------------------------------------------
#pragma mark - Shape


/**
 fixShape() copies the current shape in the array soShape[],
 and calculates the moment of inertia of the ensemble of points.
 The reference soShape[] is used by 'reshape()'.
 */
void Solid::fixShape()
{
    if ( nPoints == 0 )
        throw InvalidParameter("Cannot fix Solid with no points!");
    
    //std::clog << "Fixing Solid " << reference() << " with " << nPoints << " points\n";
    
    Vector avg, dev;
    calculateMomentum(avg, dev);
    
    // store momentum of the current shape:
    soVariance = dev.e_sum();
    
    // store the current points:
    soAmount = nPoints;
    // set reference to current shape translated to be centered:
    for ( index_t p = 0; p < nPoints; ++p )
    {
        ( Vector(pPos+DIM*p) - avg ).store(soShape+DIM*p);
        //for ( index_t d = 0; d < DIM; ++d )
        //    soShape[DIM*p+d] = pPos[DIM*p+d] - avg[d];
    }
}


/**
 The function rescale the reference shape `soShape[]`, calculated when fixShape() was called.
 The i-th dimension is scaled according to `mag[i]`.
 Note that `reshape()` should be called to apply the new reference to the current shape.
 */
void Solid::scaleShape(const real mag[DIM])
{
    //recalculate the momentum needed in rescale():
    soVariance = 0;
    // scale the specified dimension
    for ( index_t p = 0; p < soAmount; ++p )
    {
        for ( index_t d = 0; d < DIM; ++d )
        {
            soShape[DIM*p+d] *= mag[d];
            soVariance += square(soShape[DIM*p+d]);
        }
    }
}


/**
 Rescale the current cloud of points around its center of gravity,
 to recover the same 'size' as stored in `soVariance`.
 Size is measured as sum( ( x - g )^2 ).
 */
void Solid::rescale()
{
    Vector avg, dev;
    calculateMomentum(avg, dev);
    
    // calculate the momentum of the current shape:
    real M = dev.e_sum();
    
    if ( M > 0 )
    {
        // calculate the scaling factor to restore the size to 'soVariance':
        real mag = std::sqrt( soVariance / M );
    
        // scale the shape around the center of gravity:
        for ( index_t p = 0; p < nPoints; ++p )
        {
            real * ptr = pPos + DIM * p;
            (avg+mag*(Vector(ptr)-avg)).store(ptr);
            //for ( index_t d = 0; d < DIM; ++d )
            //    ptr[d] = mag * ( ptr[d] - avg[d] ) + avg[d];
        }
    }
}


/**
 reshape() finds the best isometric transformation = rotation + translation
 to bring `shap[]` onto `Mecable::pPos[]`, and replaces pPos[] by the transformed shap[].
 This restores the shape of the cloud of points given in shap[], into the current
 position and orientation of the object.
 The best translation is the ones that conserves the center of gravity,
 The best rotation is obtained differently in 2D and 3D, but it is unique,
 except for special cases.

 @todo: store the rotation and translation calculated by reshape()
*/

void Solid::reshape1D(real const* shap)
{
    assert_true( DIM == 1 );
    assert_true( soAmount == nPoints );
         
    real X = 0, C = 0;
    for ( index_t i = 0; i < nPoints; ++i )
    {
        C += pPos[i] * shap[i];
        X += pPos[i];
    }
    
    X /= real( nPoints );
    C /= abs_real(C);
    
    for ( index_t i = 0; i < nPoints; ++i )
        pPos[i] = C * shap[i] + X;
}


/**
 In 2D, the best rotation is obtained by simple math on the cross products
 and vector products of shap[] and pPos[]: (see it on paper)
*/
void Solid::reshape2D(real const* shap)
{
    assert_true( DIM == 2 );
    assert_true( soAmount == nPoints );
    
    real X = 0, Y = 0;
    real C = 0, S = 0;
    
    for ( index_t i = 0; i < nPoints; ++i )
    {
        X += pPos[0+DIM*i];
        Y += pPos[1+DIM*i];
        C += pPos[DIM*i] * shap[0+DIM*i] + pPos[1+DIM*i] * shap[1+DIM*i];
        S += shap[DIM*i] * pPos[1+DIM*i] - shap[1+DIM*i] * pPos[0+DIM*i];
    }
    X /= nPoints;
    Y /= nPoints;
    
    real n = std::sqrt( C*C + S*S );
    // cosine and sine of the rotation:
    if ( n > REAL_EPSILON ) {
        C = C / n;
        S = S / n;
    }
    else {
        C = 1;
        S = 0;
    }
    //printf(" n %8.3f, c %8.3f, s %8.3f norm = %8.3f\n", n, C, S, C*S + S*S);
    
    // apply transformation = rotation + translation:
    for ( index_t i = 0; i < nPoints; ++i )
    {
        pPos[DIM*i  ] = C * shap[DIM*i] - S * shap[DIM*i+1] + X;
        pPos[DIM*i+1] = S * shap[DIM*i] + C * shap[DIM*i+1] + Y;
    }
}


/**
 In 3D, we follow the procedure described by Berthold K.P. Horn in
 "Closed-form solution of absolute orientation using unit quaternions"
 Journal of the optical society of America A, Vol 4, Page 629, April 1987
*/
void Solid::reshape3D(real const* shap)
{
    assert_true( DIM == 3 );
    assert_true( soAmount == nPoints );
    Vector3 avg(Mecable::position());
    
    Matrix33 S(0,0);
    for ( index_t i = 0; i < nPoints; ++i )
        S.addOuterProduct(shap+DIM*i, pPos+DIM*i);
    
    // rescale matrix to keep its magnitude in range
    real alpha = S.diagonal().abs().e_sum();
    S.scale(1.0/alpha);

    real N[4*4] = { 0 };
    
    // set lower triangle of the 4x4 traceless matrix:
    N[0*4+0] = S(0,0) + S(1,1) + S(2,2);
    N[0*4+1] = S(1,2) - S(2,1);
    N[0*4+2] = S(2,0) - S(0,2);
    N[0*4+3] = S(0,1) - S(1,0);
    
    N[1*4+1] = S(0,0) - S(1,1) - S(2,2);
    N[1*4+2] = S(0,1) + S(1,0);
    N[1*4+3] = S(2,0) + S(0,2);
    
    N[2*4+2] = S(1,1) - S(0,0) - S(2,2);
    N[2*4+3] = S(1,2) + S(2,1);
    
    N[3*4+3] = S(2,2) - S(1,1) - S(0,0);

    //VecPrint::full(4, 4, N, 4, 3);

    /* 
     Use LApack to find the largest Eigenvalue (unused), and associated Eigenvector,
     which is the quaternion corresponding to the best rotation.
     */
    
    int nbv;
    real val[4];
    Quaternion<real> quat;
    const int lwork = 35*4;
    real work[lwork];
    int iwork[5*4];
    int ifail[4];
    int info = 0;
    
    lapack::xsyevx('V','I','L', 4, N, 4, 0, 0, 4, 4, REAL_EPSILON,
                   &nbv, val, quat.data(), 4, work, lwork, iwork, ifail, &info);
    
    //printf("optimal LWORK = %f\n", work[0]);
    if ( info == 0 )
    {
        //get the rotation matrix corresponding to the quaternion:
        quat.setMatrix3(S);
        
        if ( 0 )
        {
            // measure remaining error:
            real dev = 0;
            for ( index_t i = 0; i < nPoints; ++i )
                dev += distanceSqr(avg+S*Vector3(shap+DIM*i), pPos+DIM*i);
            printf("%s error %6.4f :", reference().c_str(), dev);
            printf(" eigenvalue %6.2f ", val[0]); quat.println();
        }

        // apply rotation + translation:
        for ( index_t i = 0; i < nPoints; ++i )
            (avg+S*Vector3(shap+DIM*i)).store(pPos+DIM*i);
    }
    else
    {
        // apply translation:
        for ( index_t i = 0; i < nPoints; ++i )
            (avg+Vector3(shap+DIM*i)).store(pPos+DIM*i);
        
        printf("Solid::reshape(): lapack::xsyevx() failed with code %i\n", info);
    }
}

/**
 
 Solid::getPoints() calls rescale() often and reshape() occasionally, because:
 - reshape() corrects for all kind of numerical drift but is computationally expensive
 - rescale() corrects for one kind of numerical drift that is dominant.
 .
 
 The calls for different solids are shifted by using the identity() of each Solid.
 */
void Solid::getPoints(real const* ptr)
{
    Mecable::getPoints(ptr);
    
    // for one point, nothing should be done
    if ( nPoints < 2 )
        return;
#if ( DIM >= 3 )
    reshape3D(soShape);
#elif ( DIM == 2 )
    reshape2D(soShape);
#else
    reshape1D(soShape);
#endif
}


//------------------------------------------------------------------------------
#pragma mark - Drag

/**
The mobility is that of a set of spheres in an infinite fluid (Stokes law):

Translation:
    dposition/dtime = mu_T * force
    mu_T = 6 * PI * viscosity * radius

Rotation:
    dangle/dtime = mu_R * torque
    mu_R = 8 * PI * viscosity * radius^3

 This sums up the drag of the sphere, ignoring any hydrodynamic interaction.
 Thus drag will generally be overestimated.
 */
void Solid::setDragCoefficient()
{
    real sumR = 0;
    real sumR3 = 0;
    Vector cen(0,0,0);
#if ( DIM < 3 )
    real roti = 0;     //in 2D, the total rotational inertia
#endif
    
    for ( index_t i = 0; i < nPoints; ++i )
    {
        real R = radius(i);
        if ( R > 0 )
        {
            sumR   += R;
            sumR3  += R * R * R;
            cen    += R * posP(i);
#if ( DIM < 3 )
            roti   += R * posP(i).normSqr();
#endif
        }
    }
    
    soCenter = cen / sumR;
    soDrag   = sumR * ( 6 * M_PI );
    
#if ( DIM > 2 )
    soDragRot = sumR3 * ( 8 * M_PI );
#else
    // in 2D, reduce to centroid:
    soDragRot = sumR3 * ( 8 * M_PI ) + roti * ( 6 * M_PI ) - soDrag * soCenter.normSqr();
#endif
    
    // sanity check:
    if ( soDrag < REAL_EPSILON )
        throw InvalidParameter("solid "+reference()+" has no drag (no bead with radius > 0)");

    if ( soDragRot < REAL_EPSILON )
        throw InvalidParameter("ill-formed Solid "+reference()+" has zero rotational drag");

#if ( 0 )
    std::clog << "Solid " << reference() << " (viscosity " << prop->viscosity << ") has:\n";
    std::clog << "     translation drag " << soDrag * prop->viscosity << "\n";
    std::clog << "     rotation drag    " << soDragRot * prop->viscosity << "\n";
#endif
}


void Solid::prepareMecable()
{
    setDragCoefficient();
    makeProjection();
}


real Solid::addBrownianForces(real const* fce, real alpha, real* rhs) const
{
    // mobility of entire Solid
    const real drag = prop->viscosity * soDrag;

    if ( std::isinf(drag) )
        return INFINITY;
    
    // amplitude of Brownian motion, scaled as it is added to every points:
    const real b = std::sqrt( alpha * drag / (real)nPoints );

    // noise is equally distributed over all points:
    for ( index_t j = 0; j < DIM*nPoints; ++j )
        rhs[j] = b * rhs[j] + fce[j];
    
    //std::clog << reference() << " drag " << b << " " << drag << "   " << std::sqrt(nPoints)*b/drag << "\n";
    // typical displacement of one point due to Brownian motion:
    return std::sqrt( alpha / drag );
}


//------------------------------------------------------------------------------
#pragma mark - Projection


#if ( DIM == 1 )

/**
 The projection in 1D is just summing all the forces,
 and distributing equally to all the points:
*/
void Solid::makeProjection()
{
}

void Solid::projectForces(const real* X, real* Y) const
{
    real T = 0;
    for ( index_t p = 0; p < nPoints; ++p )
        T += X[p];
    
    T /= prop->viscosity * soDrag;
    
    for ( index_t p = 0; p < nPoints; ++p )
        Y[p] = T;
}

#elif ( DIM == 2 )

/**
 Recalculate soCenter, and the rotational moment of inertia.
 */
void Solid::makeProjection()
{
    soCenter = centroid();
    
#if ( 0 )
    /*
     In 2D the rotational moment of inertia is a scalar that is invariant
     by rotation, and it is not normally necessary to recalculate it here
     
     The code below checks that the value has not changed:
     */
    Vector cen(0,0,0);
    real roti = 0;
    real sumR = 0;
    real sumR3 = 0;
    
    for ( index_t i = 0; i < nPoints; ++i )
    {
        real R = radius(i);
        if ( R > 0 )
        {
            sumR  += R;
            sumR3 += R * R * R;
            cen   += R * posP(i);
            roti  += R * posP(i).normSqr();
        }
    }
    
    cen /= sumR;
    real m = sumR3 * ( 8 * M_PI ) + roti * ( 6 * M_PI ) - soDrag * cen.normSqr();
    std::clog << "Solid2D::error(rotational_drag) " << abs_real(soDragRot-m) << " for " << reference() << "\n";
#endif
}


void Solid::projectForces(const real* X, real* Y) const
{
    Vector T(0.0,0.0);  // Translation
    real R = 0;         // Infinitesimal Rotation (a vector in Z)
    
    for ( index_t p = 0; p < nPoints; ++p )
    {
        real const* pos = pPos + DIM * p;
        real const* xxx = X + DIM * p;
        
        T.XX += xxx[0];
        T.YY += xxx[1];
        
        R += pos[0] * xxx[1] - pos[1] * xxx[0];
    }
    
    const real A = 1.0 / ( prop->viscosity * soDragRot );
    const real B = 1.0 / ( prop->viscosity * soDrag );

    R = A * ( R + cross(T,soCenter) );
    T = B * T + cross(soCenter,R);
    
    for ( index_t p = 0; p < nPoints; ++p )
    {
        real const* pos = pPos + DIM * p;
        real * yyy = Y + DIM * p;
        
        yyy[0] = T.XX - R * pos[1];
        yyy[1] = T.YY + R * pos[0];
    }
}

#elif ( DIM >= 3 )

/**
 To project in 3D, we calculate the resulting tensor by summing all
 the forces on all points, reducing it at the center of gravity.
 From this, we can deduce the forces compatible with solid motion,
 which is a combination of translation and rotation.
 */
void Solid::makeProjection()
{
    ///\todo: from reshape, we know the rotation matrix from the stored shape
    //to the current shape. We could use it to transform the inertia matrix
    real sum = 0;
    index_t cnt = 0;
    real cX = 0, cY = 0, cZ = 0;
    real XX=0, XY=0, XZ=0, YY=0, YZ=0, ZZ=0;
    
    for ( index_t i = 0; i < nPoints; ++i )
    {
        const real R = radius(i);
        if ( R > 0 )
        {
            ++cnt;
            sum += R;
            const Vector pos = posP(i);
            const Vector vec = R * pos;
            cX += vec.XX;
            cY += vec.YY;
            cZ += vec.ZZ;
            XX += vec.XX * pos.XX;
            XY += vec.XX * pos.YY;
            XZ += vec.XX * pos.ZZ;
            YY += vec.YY * pos.YY;
            YZ += vec.YY * pos.ZZ;
            ZZ += vec.ZZ * pos.ZZ;
        }
    }
    cX /= sum;
    cY /= sum;
    cZ /= sum;
    soCenter.set(cX, cY, cZ);
    if ( cnt == 1 )
    {
        // in this case there is a single bead, and the matrix is isotropic
        soMomentum = Matrix33(0, 1.0/soDragRot);
        return;
    }
    
    // scale to get the correct mobility:
    const real A = 6 * M_PI;
    const real B = soDrag;
    const real D = soDragRot - soDrag * ( cX*cX + cY*cY + cZ*cZ );
    
    // finally set the matrix in front of R in projectForces()
    // the 3x3 matrix is symmetric, and we only set its lower half:
    soMomentum(0,0) = D + A * (YY+ZZ) + B * cX * cX;
    soMomentum(1,0) =   - A *  XY     + B * cX * cY;
    soMomentum(2,0) =   - A *  XZ     + B * cX * cZ;
    soMomentum(1,1) = D + A * (XX+ZZ) + B * cY * cY;
    soMomentum(2,1) =   - A *  YZ     + B * cY * cZ;
    soMomentum(2,2) = D + A * (XX+YY) + B * cZ * cZ;
    //Matrix33 mat = soMomentum;
    soMomentum.symmetricInverse();

#if ( 0 )
    // checking here that the calculated inverse is the inverse!
    mat.copy_lower();
    const real dif = ( mat * soMomentum - Matrix33(0,1) ).norm_inf();
    //std::clog << "Solid " << reference() << ' ' << std::fixed << dif << '\n';
    if ( dif > 0.001 )
    {
        std::clog << "Solid " << reference() << ' ' << cnt << '\n';
        std::clog << std::setw(10) << soMomentum << '\n';
        std::clog << std::setw(10) << mat << '\n';
    }
#endif
}


/**
 This calculated Y <- P * X, where
 P is the projection associated with the constraints of motion without
 deformation (solid object)
 
 This calculates the total force and momentum in the center of mobility,
 scale to get speed, and distribute according to solid motion mechanics.
*/
void Solid::projectForces0(const real* X, real* Y) const
{
    Vector T(0,0,0);    //Translation
    Vector R(0,0,0);    //Rotation
    
    for ( index_t p = 0; p < nPoints; ++p )
    {
        real * pos = pPos + DIM * p;
        real const* xxx = X + DIM * p;
        
        // T = T + xxx
        T.XX += xxx[0];
        T.YY += xxx[1];
        T.ZZ += xxx[2];
        
        // R = R + cross(pos, xxx)
        R.XX += pos[1] * xxx[2] - pos[2] * xxx[1];
        R.YY += pos[2] * xxx[0] - pos[0] * xxx[2];
        R.ZZ += pos[0] * xxx[1] - pos[1] * xxx[0];
    }
    
    Vector V = R + cross(T, soCenter);
    
    const real visc = prop->viscosity;
    R = soMomentum.vecmul(V) / visc;
    // reduce Torque to center of mobility:
    T = cross(soCenter, R) + T / ( visc * soDrag );

    for ( index_t p = 0; p < nPoints; ++p )
    {
        real const* pos = pPos + DIM * p;
        real * yyy = Y + DIM * p;
        
        yyy[0] = T.XX + R.YY * pos[2] - R.ZZ * pos[1];
        yyy[1] = T.YY + R.ZZ * pos[0] - R.XX * pos[2];
        yyy[2] = T.ZZ + R.XX * pos[1] - R.YY * pos[0];
    }
}


void Solid::projectForces(const real* X, real* Y) const
{
    Vector T(0,0,0);   //Translation
    Vector R(0,0,0);   //Rotation
    
    for ( index_t p = 0; p < nPoints; ++p )
    {
        Vector F(X+DIM*p);
        T += F;
        R += cross(Vector(pPos+DIM*p), F);
    }
    
    Vector V = R + cross(T, soCenter);

    const real visc = prop->viscosity;
    R = soMomentum.vecmul(V) / visc;
    // reduce Torque to center of mobility:
    T = cross(soCenter, R) + T / ( visc * soDrag );
    
    for ( index_t p = 0; p < nPoints; ++p )
    {
        // Y <- T + cross(R, X)
        (T+cross(R, Vector(pPos+DIM*p))).store(Y+DIM*p);
    }
}

#endif

//------------------------------------------------------------------------------
#pragma mark - I/O & print


void Solid::write(Outputter& out) const
{
    writeMarker(out, TAG);
    out.writeUInt16(nPoints);
    for ( index_t p = 0; p < nPoints ; ++p )
    {
        out.writeFloats(pPos + DIM * p, DIM, '\n');
        out.writeSoftSpace();
        out.writeFloat(soRadius[p]);
    }
#if NEW_SOLID_CLAMP
    writeMarker(out, CLAMP_TAG);
    out.writeFloat(clamp_stiff);
    out.writeFloats(clamp_place, DIM, '\n');
#endif
}


void Solid::read(Inputter& in, Simul& sim, ObjectTag tag)
{
    if ( tag == TAG )
    {
        try
        {
            index_t nbp = in.readUInt16();
            setNbPoints(nbp);
            for ( index_t i = 0; i < nbp ; ++i )
            {
                in.readFloats(pPos+DIM*i, DIM);
                soRadius[i] = in.readFloat();
            }
        }
        catch( Exception & e )
        {
            clearPoints();
            throw;
        }
        fixShape();
    }
    else if ( tag == SOLID_TAG )
    {
        [[maybe_unused]] ObjectID id = in.readUInt32();
        in.readFloat();
    }
    else if ( tag == CLAMP_TAG )
    {
#if NEW_SOLID_CLAMP
        clamp_stiff = in.readFloat();
        in.readFloats(clamp_place, DIM);
#else
        Vector tmp;
        in.readFloat();
        in.readFloats(tmp, DIM);
        std::clog << "Ignored " << reference() << " clamp @ " << tmp << "\n";
#endif
    }
}


void Solid::print(std::ostream& os, bool write_shape) const
{
    std::streamsize p = os.precision(4);
    os << "new solid " << reference() << '\n';
    os << "{\n";
    os << " nb_points = " << nPoints << '\n';
    for ( index_t n = 0; n < nPoints ; ++n )
    {
        os << " point" << n+1 << " = ";
        if ( write_shape )
            os << std::setw(8) << std::fixed << Vector(soShape+DIM*n);
        else
            os << std::setw(8) << std::fixed << Vector(pPos+DIM*n);
        if ( radius(n) > 0 )
            os << ", " << radius(n);
        os << '\n';
    }
    os << "}" << '\n';
    os.precision(p);
}


std::ostream& operator << (std::ostream& os, Solid const& arg)
{
    arg.print(os, false);
    return os;
}

