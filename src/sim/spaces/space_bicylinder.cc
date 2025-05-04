// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University.
#include "space_bicylinder.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "mecapoint.h"
#include "glossary.h"
#include "meca.h"
#include "project_ellipse.h"

SpaceBicylinder::SpaceBicylinder(SpaceProp const* p)
: Space(p)
{
    if ( DIM < 3 )
        throw InvalidParameter("bicylinder is only valid in 3D: use rectangle instead");
    radius_ = 0;
}


void SpaceBicylinder::resize(Glossary& opt)
{
    real rad = radius_;
    
    if ( opt.set(rad, "diameter") )
        rad *= 0.5;
    else opt.set(rad, "radius");

    if ( rad < 0 )
        throw InvalidParameter("cylinder:radius must be >= 0");
    
    radius_ = rad;
}


void SpaceBicylinder::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(-radius_,-radius_,-radius_);
    sup.set( radius_, radius_, radius_);
}


bool SpaceBicylinder::inside(Vector const& W) const
{
#if ( DIM > 2 )
    return max_real(square(W.XX), square(W.YY)) <= square(radius_) - square(W.ZZ);
#else
    return false;
#endif
}


bool SpaceBicylinder::allInside(Vector const& W, const real rad) const
{
    assert_true( rad >= 0 );
#if ( DIM > 2 )
    return max_real(square(W.XX), square(W.YY)) <= square(radius_-rad) - square(W.ZZ);
#else
    return false;
#endif
}


Vector SpaceBicylinder::place() const
{
#if ( DIM >= 3 )
    Vector W;
    do {
        W = Vector::randS(radius_);
    } while (max_real(square(W.XX), square(W.YY)) > square(radius_) - square(W.ZZ) );
    return W;
#else
    return Vector(0, 0, 0);
#endif
}

//------------------------------------------------------------------------------
Vector SpaceBicylinder::project(Vector const& W) const
{
#if ( DIM >= 3 )
    real Y = radius_ / W.normXZ();
    real X = radius_ / W.normYZ();

    if ( 1 <= X ) // inside Cylinder of axis X
    {
        // inside both cylinders, we project on the closest one:
        if ( X < Y )
            return Vector(W.XX, X*W.YY, X*W.ZZ);
        // inside cylinder of axis X, we project on the other cylinder
        return Vector(Y*W.XX, W.YY, Y*W.ZZ);
    }
    else if ( 1 <= Y )
    {
        // inside Cylinder of axis Y, we projection on X-Cylinder
        return Vector(W.XX, X*W.YY, X*W.ZZ);
    }
    
    // projection on the cylinder:
    if ( abs_real(Y) > abs_real(X) )
    {
        // projection on cylinder aligned with X:
        Vector pX(W.XX, X*W.YY, X*W.ZZ);
        if ( pX.normXZSqr() < square(radius_) )
            return Vector(W.XX, X*W.YY, X*W.ZZ);
    }
    else
    {
        // projection on cylinder aligned with Y:
        Vector pY(Y*W.XX, W.YY, Y*W.ZZ);
        if ( pY.normYZSqr() < square(radius_) )
            return Vector(Y*W.XX, W.YY, Y*W.ZZ);
    }
    
    // projection on the junction between the two cylinders:
    if ( W.XX * W.YY > 0 )
    {
        // in plane X=Y
        projectEllipse(X, Y, M_SQRT1_2*(W.XX+W.YY), W.ZZ, M_SQRT2*radius_, radius_);
        return Vector(M_SQRT1_2*X, M_SQRT1_2*X, Y);
    }
    else
    {
        // in plane X=-Y
        projectEllipse(X, Y, M_SQRT1_2*(W.XX-W.YY), W.ZZ, M_SQRT2*radius_, radius_);
        return Vector(M_SQRT1_2*X, -M_SQRT1_2*X, Y);
    }

#endif
    return Vector(W);
}

//------------------------------------------------------------------------------

void SpaceBicylinder::setConfinement(Vector const& W, Mecapoint const& mp, Meca& meca,
                                     real stiff, const real rad)
{
#if ( DIM >= 3 )
    if ( min_real(square(W.XX), square(W.YY)) > square(rad) - square(W.ZZ) )
    {
        real Y = rad / W.normXZ();
        real X = rad / W.normYZ();

        // projection on the cylinder:
        if ( abs_real(Y) > abs_real(X) )
        {
            // projection on cylinder aligned with X:
            Vector pX(W.XX, X*W.YY, X*W.ZZ);
            if ( pX.normXZSqr() < square(rad) )
                meca.addCylinderClampX(mp, rad, stiff);
        }
        else
        {
            // projection on cylinder aligned with Y:
            Vector pY(Y*W.XX, W.YY, Y*W.ZZ);
            if ( pY.normYZSqr() < square(rad) )
                meca.addCylinderClampY(mp, rad, stiff);
        }

        Vector E, T;
        if ( W.XX * W.YY > 0 )
        {
            // in plane X=Y
            projectEllipse(X, Y, M_SQRT1_2*(W.XX+W.YY), W.ZZ, M_SQRT2*rad, rad);
            X *= M_SQRT1_2;
            E.set(X, X, Y);
            // X = -M_SQRT2 * Y; Y = M_SQRT1_2 * X
            real s = 1.0 / std::sqrt(2*square(Y) + square(X));
            T.set(s*Y, s*Y, -s*X);
        }
        else
        {
            // in plane X=-Y
            projectEllipse(X, Y, M_SQRT1_2*(W.XX-W.YY), W.ZZ, M_SQRT2*rad, rad);
            X *= M_SQRT1_2;
            E.set(X, -X, Y);
            real s = 1.0 / std::sqrt(2*square(Y) + square(X));
            T.set(-s*Y, s*Y, s*X);
        }
        meca.addLineClamp(mp, E, T, stiff);
        //meca.addPointClamp(mp, E, stiff);
    }
    else if ( abs_real(W.YY) > abs_real(W.XX) )
    {
        meca.addCylinderClampX(mp, rad, stiff);
    }
    else
    {
        meca.addCylinderClampY(mp, rad, stiff);
    }
#endif
}


void SpaceBicylinder::setConfinement(Vector const& pos, Mecapoint const& mp, Meca& meca, real stiff) const
{
    setConfinement(pos, mp, meca, stiff, radius_);
}


void SpaceBicylinder::setConfinement(Vector const& pos, Mecapoint const& mp,
                                     real rad, Meca& meca, real stiff) const
{
    real R = max_real(0, radius_ - rad);
    
    setConfinement(pos, mp, meca, stiff, R);
}

//------------------------------------------------------------------------------

void SpaceBicylinder::write(Outputter& out) const
{
    writeMarker(out, TAG);
    writeShape(out, "RR");
    out.writeUInt16(2);
    out.writeFloat(radius_);
    out.writeFloat(radius_);
}


void SpaceBicylinder::setLengths(const real len[8])
{
    radius_ = len[0];
    radius_ = len[1];
}

void SpaceBicylinder::read(Inputter& in, Simul&, ObjectTag)
{
    real len[8] = { 0 };
    readShape(in, 8, len, "RR");
    setLengths(len);
}

//------------------------------------------------------------------------------
#pragma mark - OpenGL display

#ifdef DISPLAY

#include "gle.h"
#include "gym_flute.h"
#include "gym_draw.h"


void SpaceBicylinder::draw3D() const
{
    const float R(radius_);
    const unsigned inc = 1;
    unsigned i = 0;
    flute6* flu = gym::mapBufferV3N3(2*gle::pi_4half+4);
    for ( unsigned n = gle::pi_1half; n < gle::pi_3half; n += inc )
    {
        float C = gle::cos_(n), S = gle::sin_(n);
        flu[i++] = {R*C, +R*C, R*S, C, 0, S };
        flu[i++] = {R*C, -R*C, R*S, C, 0, S };
    }
    for ( unsigned n = gle::pi_3half; n <= gle::pi_5half; n += inc )
    {
        float C = gle::cos_(n), S = gle::sin_(n);
        flu[i++] = {R*C, -R*C, R*S, C, 0, S };
        flu[i++] = {R*C, +R*C, R*S, C, 0, S };
    }
    gym::unmapBufferV3N3();
    gym::drawTriangleStrip(0, i);
    i = 0;
    flu = gym::mapBufferV3N3(2*gle::pi_4half+4);
    for ( unsigned n = gle::pi_1half; n < gle::pi_3half; n += inc )
    {
        float C = gle::cos_(n), S = gle::sin_(n);
        flu[i++] = {-R*C, R*C, R*S, 0, C, S };
        flu[i++] = {+R*C, R*C, R*S, 0, C, S };
    }
    for ( unsigned n = gle::pi_3half; n <= gle::pi_5half; n += inc )
    {
        float C = gle::cos_(n), S = gle::sin_(n);
        flu[i++] = {+R*C, R*C, R*S, 0, C, S };
        flu[i++] = {-R*C, R*C, R*S, 0, C, S };
    }
    gym::unmapBufferV3N3();
    gym::drawTriangleStrip(0, i);
}

#else

void SpaceBicylinder::draw3D() const {}

#endif

