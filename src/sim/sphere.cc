// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "cymdef.h"
#include "primitives.h"
#include "assert_macro.h"
#include "exceptions.h"
#include "messages.h"
#include "glossary.h"
#include "mecapoint.h"
#include "sphere_prop.h"
#include "object_set.h"
#include "space_prop.h"
#include "space.h"
#include "sphere.h"
#include "meca.h"
#include "modulo.h"
#include "simul.h"


//------------------- construction and destruction ---------------------------

/**
 The Sphere is made with no point
 */
Sphere::Sphere(SphereProp const* p)
: spRadius(0), spDrag(0), spDragRot(0), sDir(nullptr), prop(p)
{
}


/*
 This will initialize with a center point and references
 */
Sphere::Sphere(SphereProp const* p, real rad)
: spRadius(rad), spDrag(0), spDragRot(0), sDir(nullptr), prop(p)
{
    if ( !prop )
        throw InvalidParameter("Sphere:prop should be specified");
    
    if ( rad <= 0 )
        throw InvalidParameter("sphere:radius should be > 0");
    
    // center point
    assert_true( nPoints == 0 );
    addPoint( Vector(0,0,0) );
    
    // reference points to track the orientation of the sphere
    if ( DIM >= 2 )
        addPoint( Vector(rad,0,0) );
    if ( DIM >= 3 ) {
        addPoint( Vector(0,rad,0) );
        addPoint( Vector(0,0,rad) );
    }
    
    setDragCoefficient();
}


Sphere::Sphere(const Sphere & o)
: Mecable(o), sDir(nullptr)
{
    prop     = o.prop;
    spRadius = o.spRadius;
}


Sphere & Sphere::operator = (const Sphere & o)
{
    prop     = o.prop;
    spRadius = o.spRadius;
    return *this;
}


Sphere::~Sphere()
{
    if ( objset() )
        simul().singles.deleteWrists(this);
    release();
    prop = nullptr;
}

//------------------------------------------------------------------------------
#pragma mark -

/*
 here 'cp' is the vector from the center to the point to be added,
 in other words, the position of the point in the local reference frame.
 */
index_t Sphere::addSurfacePoint(Vector const& cp)
{
    return addPoint(posP(0)+cp.normalized(spRadius));
}


/**
 @ingroup NewObject
 
 Specify radius and number of surface points of a Sphere:

     new sphere NAME
     {
        radius = VALUE, DEVIATION, MINIMUM
        point1 = INTEGER, POSITION [, SINGLE]
     }
 
 Variability around the mean radius is added if 'DEVIATION' and 'MINIMUM' are specified.
 All values must be positive. When specifying the points on the surface of the Sphere,
 an `INTEGER` set the number of points created, and `POSITION` can be a `VECTOR`,
 or the string 'surface'.  Multiple points commands can be given: point1, point2, etc.
 
 <h3> Add Singles to a Sphere </h3>
 
 The parameter 'attach' can be used to add Single to the points of a Solid:
 
     new sphere NAME
     {
        radius   = ...
        point1   = ...
        etc.
        attach   = SINGLE [, SINGLE] ...   //this will attach to all points
        attach1  = SINGLE [, SINGLE] ...   //this will attach to point1
        etc.
     }
 
 Where `SINGLE` is string containing at most 3 words: `[INTEGER] NAME [each]`,
 where the `INTEGER` specifies the number of Singles, `NAME` specifies their name,
 and the optional word `each` species that the command applies to every point.
 
 The command `attach` applies to all the points of the Solid, while `attach1`,
 `attach2`, etc. apply to the points specified by `point1`, `point2`, etc.
 With `attach`, the Singles are distributed randomly on all the points,
 and if `each` is specified, the specification is repeated for each point.
 
 For example if `grafted` is the name of a Single, one can use:
 
     new solid NAME
     {
        attach1 = 1 grafted each
        attach2 = 10 grafted
     }
 */
ObjectList Sphere::build(Glossary & opt, Simul& sim)
{
    ObjectList objs(this);
    std::string str;
    index_t inp = 1, inx = 0, nbp = 1;

    if ( opt.has_key("point0") )
        throw InvalidParameter("point indices start at 1 (use `point1`, `point2`, etc.)");

    // interpret each instruction as a command to add points:
    std::string var = "point1";
    while ( opt.has_key(var) )
    {
        inx = 0;
        nbp = 1;
        if ( opt.set_positive_integer(nbp, var) )
            ++inx;
        
        if ( nbp > 0 )
        {
            index_t fip = nPoints;
            str = opt.value(var, inx);
            // add 'nbp' points:
            for ( size_t n = 0; n < nbp; ++n )
            {
                Vector vec(0,0,0);
                if ( str == "surface" )
                    vec = Vector::randU(radius());
                else
                {
                    vec = Cytosim::findPosition(str, nullptr);
                    if ( 8 * vec.norm() < spRadius )
                        throw InvalidParameter(var+" cannot be brought to the Sphere surface");
                }
                addSurfacePoint(vec);
            }
            
            // attach Single to this set of points:
            while ( opt.set(str, var, ++inx) )
                sim.singles.makeWrists(objs, this, fip, nbp, str);
            
            // attach Single to this set of points:
            inx = 0;
            var = "attach" + std::to_string(inp);
            while ( opt.set(str, var, inx++) )
                sim.singles.makeWrists(objs, this, fip, nbp, str);
        }
        
        // set next keyword:
        var = "point" + std::to_string(++inp);
    }
    
    
    // attach Singles distributed over the surface points:
    inx = 0;
    while ( opt.set(str, "attach", inx++) )
        sim.singles.makeWrists(objs, this, nbRefPoints, nbSurfacePoints(), str);

    
    // final verification of the number of points:
    nbp = 0;
    if ( opt.set(nbp, "nb_points")  &&  nbp != nPoints )
    {
        throw InvalidParameter("could not find the number of points specified in solid:nb_points");
    }
    return objs;
}


//------------------------------------------------------------------------------

void Sphere::resize(const real R)
{
    //std::clog << "Sphere::resize " << R << '\n';
    if ( R > 0 )
    {
        spRadius = R;
        reshape();
    }
}


/**
 Expect higher friction due to flow around the sphere in a narrow tube.
 This is only valid if (r -a)/a << 1, where r = radius of the tube, and
 a = radius of the sphere.
 
 The formula are taken from:
     The Motion of a Closely-Fitting Sphere in a Fluid-Filled Tube
     P. Bungay and H. Brenner, Int. J. Multiphase Flow
     Vol 1, pp. 25-56, 1973 (see 3.6, 4.68a and 5.11)
 @latex
     \gamma = \frac{9 \pi^2 \sqrt{2} }{ 4\epsilon^{5/2}} \,\eta \, r
     \gamma^{rot} = 2 \pi^2 \sqrt{\frac{2}{\epsilon}} \,\eta\, r^3
 @end
 */
void Sphere::setDragCoefficientPiston(real thick)
{
    const real rad = radius();
    real eps = ( 0.5 * thick - rad ) / rad;
    
    if ( eps <= 0 )
        throw InvalidParameter("Sphere's piston_effect would yield invalid values");
    if ( eps > 1 )
        throw InvalidParameter("Sphere's piston_effect is inapplicable");

    real vr = M_PI * M_PI * prop->viscosity * rad;
    spDrag    = vr * ( 2.25 * M_SQRT2 ) * std::pow(eps,-2.5);
    spDragRot = vr * square(rad) * std::sqrt(8/eps);
        
    //report the reduced mobility of the sphere:
    Cytosim::log.print("Sphere of radius %.3f has piston mobility %.2e\n", spRadius, spDrag);
}


/**
 Except if the `piston_effect` is enabled,
 the mobility is that of a sphere in an infinite fluid (Stokes law):
 
 Translation:
     dposition/dtime = mu_T * force
     mu_T = 6 * PI * viscosity * radius
 
 Rotation:
     dangle/dtime = mu_R * torque
     mu_R = 8 * PI * viscosity * radius^3
 
 */
void Sphere::setDragCoefficient()
{
    assert_true( spRadius > 0 );
    const real rad = radius();
    
    real vr = prop->viscosity * rad;
    // hydrodynamic not corrected: infinite fluid is assumed
    spDrag    = ( 6 * M_PI ) * vr;
    spDragRot = ( 8 * M_PI ) * ( vr * square(rad) );

    if ( prop->piston_effect )
    {
        if ( prop->confine_space )
            setDragCoefficientPiston(prop->confine_space->thickness());
        else
            Cytosim::warn("Piston effect ignored for lack of confine space\n");
    }
}


//------------------------------------------------------------------------------
#pragma mark - Mecable

void Sphere::allocateMecable(index_t nbp)
{
    real * ptr = Mecable::allocateMemory(nbp, DIM);
    if ( ptr )
        sDir = ptr;
}


void Sphere::release()
{
}


void Sphere::prepareMecable()
{
    setDragCoefficient();
    assert_true( spDrag > 0 );
    assert_true( spDragRot > 0 );
    
    makeProjection();
}


void Sphere::setInteractions(Meca& meca) const
{
    if ( prop->confine != CONFINE_OFF )
    {
        Space const* spc = prop->confine_space;
        
        switch ( prop->confine )
        {
            case CONFINE_INSIDE:
            {
                Vector cen(pPos);
                if ( ! spc->inside(cen) )
                    spc->setConfinement(cen, Mecapoint(this, 0), meca, prop->confine_stiff);
            } break;
            
            case CONFINE_ALL_INSIDE:
            {
                Vector cen(pPos);
                if ( ! spc->allInside(cen, spRadius) )
                    spc->setConfinement(cen, Mecapoint(this, 0), spRadius, meca, prop->confine_stiff);
            } break;
            
            case CONFINE_ON:
                spc->setConfinement(posP(0), Mecapoint(this, 0), meca, prop->confine_stiff);
            
            default:
                throw InvalidParameter("Invalid sphere::confine");
        }
    }
}

//------------------------------------------------------------------------------

real Sphere::addBrownianForces(real const* fce, real alpha, real* rhs) const
{
    real bT = 0;
    if ( ! std::isinf(spDrag) )
        bT = std::sqrt( alpha * spDrag );
    
    // magnitude of diffusion of surface points:
    real bS = 0;
    if ( prop->point_mobility > 0 )
        bS = std::sqrt( alpha / prop->point_mobility );

    Vector F(0, 0, 0);
    Torque T(nullTorque);

    // position of center:
    real cX = pPos[0];
    real cY = pPos[1];
    real cZ = pPos[2];

    /*
     Stage 1: add random forces to the surface points, calculating resulting force
     and torque in F and T. This will be subtracted in stage 3.
     */
    for ( size_t p = nbRefPoints; p < nPoints; ++p )
    {
        real * res = rhs + DIM * p;
        real * pos = pPos + DIM * p;
        Vector fp = bS * Vector(res);
        
        F += fp;
        
        res[0] = fce[DIM*p+0] + fp.XX;
#if   ( DIM == 2 )
        res[1] = fce[DIM*p+1] + fp.YY;
        T += cross(Vector(pos[0]-cX, pos[1]-cY), fp);
#elif ( DIM >= 3 )
        res[1] = fce[DIM*p+1] + fp.YY;
        res[2] = fce[DIM*p+2] + fp.ZZ;
        T += cross(Vector(pos[0]-cX, pos[1]-cY, pos[2]-cZ), fp);
#endif
    }

    /*
     Stage 2: The Torque is distributed to reference points on the surface of the Sphere.
     In 2D, there is one such point, and the coefficient is therefore 1.
     in 3D, there are 3 reference points, but always one is parallel to the axis of the torque,
     and the decomposition over these 3 points gives a factor 2.
     */
    T /= ( DIM - 1 ) * spRadius * spRadius;
    Vector R = cross(Vector(cX,cY,cZ), T);

    for ( size_t p = 1; p < nbRefPoints; ++p )
    {
        real * res = rhs + DIM * p;
        real const* pos = pPos + DIM * p;
        Vector fp = bT * Vector(res);
        assert_true( fce[DIM*p+0] == 0 );
#if   ( DIM == 2 )
        assert_true( fce[DIM*p+1] == 0 );
        res[0] = fp.XX - R.XX + T * pos[1];
        res[1] = fp.YY - R.YY - T * pos[0];
        F += fp + cross(Vector(pos[0]-cX, pos[1]-cY), T);
#elif ( DIM >= 3 )
        assert_true( fce[DIM*p+1] == 0 );
        assert_true( fce[DIM*p+2] == 0 );
        res[0] = fp.XX - R.XX - T.YY * pos[2] + T.ZZ * pos[1];
        res[1] = fp.YY - R.YY - T.ZZ * pos[0] + T.XX * pos[2];
        res[2] = fp.ZZ - R.ZZ - T.XX * pos[1] + T.YY * pos[0];
        F += fp + cross(Vector(pos[0]-cX, pos[1]-cY, pos[2]-cZ), T);
#endif
    }
    
    /*
     Stage 3: add random displacement to the center of the sphere, but subtract
     the force generated in stage 1, such as to not overestimate diffusion
     */
#if   ( DIM == 2 )
    rhs[0] = fce[0] + bT * rhs[0] - F.XX;
    rhs[1] = fce[1] + bT * rhs[1] - F.YY;
#elif ( DIM >= 3 )
    rhs[0] = fce[0] + bT * rhs[0] - F.XX;
    rhs[1] = fce[1] + bT * rhs[1] - F.YY;
    rhs[2] = fce[2] + bT * rhs[2] - F.ZZ;
#endif

    return std::max(bT/spDrag, bS*prop->point_mobility);
}


/**
 Here we start from the i-th Vector and make the other ones orthogonal.
 There must be a better way to do this...
 */
void Sphere::orthogonalize(index_t i)
{
#if ( DIM >= 3 )
    const index_t ix = 1 + i;
    const index_t iy = 1 + (i+1)%3;
    const index_t iz = 1 + (i+2)%3;
    
    Vector cen(pPos);
    assert_true( nPoints >= nbRefPoints );
    
    // reduce to the center of mass an normalize
    Vector tmpX = posP(ix) - cen;
    Vector tmpY = posP(iy) - cen;
    Vector tmpZ = normalize( posP(iz) - cen );
    
    // make tmpY orthogonal to tmpZ, and normalized
    tmpY = normalize(tmpY - dot(tmpZ, tmpY) * tmpZ);
    
    // make tmpX orthogonal to tmpZ and tmpY
    tmpX = normalize(tmpX - dot(tmpZ, tmpX) * tmpZ - dot(tmpY, tmpX) * tmpY);
    
    // store corrected vectors back into the array
    ( cen + spRadius * tmpX ).store(pPos+DIM*ix);
    ( cen + spRadius * tmpY ).store(pPos+DIM*iy);
    ( cen + spRadius * tmpZ ).store(pPos+DIM*iz);
#endif
}


/**
 we get rid of finite-step errors but conserve the shape
 by projecting back onto the sphere,
 without changing the position of point zero (the center)
*/
void Sphere::reshape()
{
    assert_true( nPoints > 0 );
    assert_true( spRadius > 0 );
    Vector cen(pPos);
    
    for ( index_t i = 1; i < nPoints; ++i )
    {
        Vector off = ( posP(i) - cen ).normalized(spRadius);
        setPoint(i, cen + off);
    }
    
#if ( DIM >= 3 )
    orthogonalize(RNG.pint32(3));
#endif
}


/// adjust position, projecting the surface point on the sphere
void Sphere::getPoints(real const* arg)
{
    assert_true( nPoints > 0 );
    assert_true( spRadius > 0 );

    Vector cen(arg);
    copy_real(DIM, arg, pPos);
    
    for ( index_t i = 1; i < nPoints; ++i )
    {
        Vector off = ( Vector(arg+DIM*i) - cen ).normalized(spRadius);
        setPoint(i, cen + off);
    }
    
#if ( DIM >= 3 )
    orthogonalize(RNG.pint32(3));
#endif
}


//------------------------------------------------------------------------------
//------------------- methods for the projection -------------------------------
#pragma mark -

#if ( DIM == 1 )

//this is unsafe, don't use Sphere in 1D!
void Sphere::makeProjection() { ABORT_NOW("Sphere is not implemented in 1D"); }
void Sphere::projectForces(const real* X, real* Y) const {}

#else

/**
 prepare variables for the projection 
 */
void Sphere::makeProjection()
{
    real cen[DIM];
    real curv = 1.0 / spRadius;
    assert_true( nPoints >= nbRefPoints );
    for ( size_t d = 0; d < DIM; ++d )
        cen[d] = curv * pPos[d];
    // calculate radial vectors from center:
    for ( size_t p = nbRefPoints; p < nPoints; ++p )
    {
        real * dir = sDir + DIM * p;
        real * pos = pPos + DIM * p;
        for ( int d = 0; d < DIM; ++d )
            dir[d] = curv * pos[d] - cen[d];
    }
}


// The function should set Y <- sc * mobility * X.
void Sphere::projectForces(const real* X, real* Y) const
{
    // total force:
    Vector F(0,0,0);
    
    // total torque:
#if   ( DIM == 2 )
    real T = 0;
#elif ( DIM >= 3 )
    Vector T(0,0,0);
#endif
    
    for ( size_t p = 0; p < nPoints; ++p )
    {
        real * pos = pPos + DIM * p;
        real const* xxx = X + DIM * p;
        
        F.XX += xxx[0];
        F.YY += xxx[1];
#if   ( DIM == 2 )
        T    += pos[0] * xxx[1] - pos[1] * xxx[0];
#elif ( DIM >= 3 )
        F.ZZ += xxx[2];
        T.XX += pos[1] * xxx[2] - pos[2] * xxx[1];
        T.YY += pos[2] * xxx[0] - pos[0] * xxx[2];
        T.ZZ += pos[0] * xxx[1] - pos[1] * xxx[0];
#endif
    }
    
    Vector cen(pPos);
    
    T -= cross(cen, F); // reduce the torque to the center of mass
    T *= 1.0/spDragRot; // multiply by the mobility
    // scale force and add component due to shift of origin:
    F  = F * (1.0/spDrag) + cross(cen, T);
    
    for ( size_t p = 0; p < nbRefPoints; ++p )
    {
        real * yyy = Y + DIM * p;
        real * pos = pPos + DIM * p;

        for ( size_t d = 0; d < DIM; ++d )
#if   ( DIM == 2 )
        yyy[0] = F.XX - T * pos[1];
        yyy[1] = F.YY + T * pos[0];
#elif ( DIM >= 3 )
        yyy[0] = F.XX + T.YY * pos[2] - T.ZZ * pos[1];
        yyy[1] = F.YY + T.ZZ * pos[0] - T.XX * pos[2];
        yyy[2] = F.ZZ + T.XX * pos[1] - T.YY * pos[0];
#endif
    }
    
    //scale by point mobility:
    const real mob = prop->point_mobility;

    for ( size_t p = nbRefPoints; p < nPoints; ++p )
    {
        real * yyy = Y + DIM * p;
        real * pos = pPos + DIM * p;
        real * rad = sDir + DIM * p;
        real const* xxx = X + DIM * p;
#if   ( DIM == 2 )
        real a = rad[0] * xxx[0] + rad[1] * xxx[1];
        yyy[0] = F.XX - T * pos[1] + mob * ( xxx[0] - a * rad[0] );
        yyy[1] = F.YY + T * pos[0] + mob * ( xxx[1] - a * rad[1] );
#elif ( DIM >= 3 )
        real a = rad[0] * xxx[0] + rad[1] * xxx[1] + rad[2] * xxx[2];
        yyy[0] = F.XX + T.YY * pos[2] - T.ZZ * pos[1] + mob * ( xxx[0] - a * rad[0] );
        yyy[1] = F.YY + T.ZZ * pos[0] - T.XX * pos[2] + mob * ( xxx[1] - a * rad[1] );
        yyy[2] = F.ZZ + T.XX * pos[1] - T.YY * pos[0] + mob * ( xxx[2] - a * rad[2] );
#endif
    }
}

#endif


//------------------------------------------------------------------------------
#pragma mark -

void Sphere::write(Outputter& out) const
{
    writeMarker(out, TAG);
    out.writeFloat(radius());
    Mecable::write(out);
}


void Sphere::read(Inputter& in, Simul& sim, ObjectTag tag)
{
    try
    {
        real rad;
        rad = in.readFloat();
        Mecable::read(in, sim, tag);
        resize(rad);
    }
    catch( Exception & e )
    {
        clearPoints();
        throw;
    }
}


void Sphere::print(std::ostream& os) const

{
    os << "new sphere " << reference() << '\n';
    os << "{\n";
    os << " nb_points = " << nbPoints() << '\n';
    for ( index_t i = 0; i < nbPoints() ; ++i )
        os << " point" << i+1 << " = " << posP(i) << '\n';
    os << "}" << '\n';
}


std::ostream& operator << (std::ostream& os, Sphere const& arg)
{
    arg.print(os);
    return os;
}

