// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "dim.h"
#include "space_tee.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "glossary.h"
#include "quartic_solver.h"
#include "random.h"


SpaceTee::SpaceTee(SpaceProp const* p)
: Space(p)
{
    if ( DIM == 1 )
        throw InvalidParameter("tee cannot be used in 1D");
}


void SpaceTee::resize(Glossary& opt)
{
    real arm = tArmLength, jun = tJunction;
    real len = tLength, rad = tRadius;
    
    if ( opt.set(rad, "diameter") )
        rad *= 0.5;
    else opt.set(rad, "radius");
    opt.set(len, "length");
    opt.set(jun, "junction");
    opt.set(arm, "arm");

    if ( len <= 0 || rad <= 0 || arm < 0 )
        throw InvalidParameter("Space tee can't have negative length, arm length or radius.");
    if ( abs_real(jun)+rad > len )
        throw InvalidParameter("Space tee: the position of the branch plus the radius must lie within the length of the T.");
    
    tLength = len;
    tRadius = rad;
    tArmLength = arm;
    tJunction = jun;
    update();
}

//------------------------------------------------------------------------------

void SpaceTee::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(-tRadius-tLength,-tRadius, -tRadius );
    sup.set( tRadius+tLength, tArmLength+2*tRadius, tRadius );
}


real SpaceTee::volume() const
{
#if ( DIM == 1 )
    return 0;
#elif ( DIM == 2 )
    real base = 4 * tLength * tRadius + M_PI * tRadiusSq;
    real arm  = 2 * tArmLength * tRadius + M_PI_2 * tRadiusSq;
    return( base + arm );
#else
    //the complete base cylinder
    real base  = 2*M_PI * tLength * tRadiusSq + 4*M_PI/3.0 * tRadius * tRadiusSq;
    //the part of the arm with y > tRadius
    real arm   =  tArmLength * M_PI * tRadiusSq + 2*M_PI/3.0 * tRadius * tRadiusSq;
    //the part of the arm with y < tRadius without the intersection with the base 
    real extra = ( M_PI - 8./3. )*tRadius*tRadiusSq;
    return( base + arm + extra );
#endif
}


//------------------------------------------------------------------------------
bool SpaceTee::inside(Vector const& W) const
{
#if ( DIM > 1 )
    real nrmSq      = 0;
    const real x    = abs_real(W.XX);
    const real xRel = (W.XX - tJunction);
    
    //check if w is inside the base cylinder
    if ( x > tLength )
        nrmSq = square(x - tLength);
#if ( DIM == 2 )
    nrmSq += square(W.YY);
#elif ( DIM > 2 )
    nrmSq += square(W.YY) + square(W.ZZ);
#endif
    if ( nrmSq <= tRadiusSq ) return( true );
    
    //check if w is inside the arm
    if ( W.YY >= 0 )
    {
        nrmSq = 0;
        if ( W.YY > tArmLength+tRadius )
            nrmSq = square(W.YY - (tArmLength+tRadius));
#if ( DIM == 2 )
        nrmSq += square(xRel);
#elif ( DIM > 2 )
        nrmSq += square(xRel) + square(W.ZZ);
#endif
        return( nrmSq <= tRadiusSq );
    }
#endif
    return false;
}


//------------------------------------------------------------------------------
real SpaceTee::projectOnBase(const Vector& W, Vector& P) const
{
real scale, nrm = 0;
#if ( DIM > 1 )
#if ( DIM == 2 )
    nrm = square(W.YY);
#elif ( DIM > 2 )
    nrm = square(W.YY) + square(W.ZZ);
#endif
    if ( W.XX >  tLength )
        nrm += square(W.XX - tLength);
    else if ( W.XX < -tLength )
        nrm += square(W.XX + tLength);
    
    if ( nrm > 0 ) {
        nrm   = sqrt(nrm);
        scale = tRadius/nrm;
    }
    else {
        nrm   = 0;
        scale = 0;
    }
    
    real pX, pY = 0, pZ = 0;
    
    if ( W.XX >  tLength )
        pX =  tLength + scale*(W.XX - tLength);
    else if ( W.XX < -tLength )
        pX = -tLength + scale*(W.XX + tLength);
    else
        pX = W.XX;
    
    if ( scale != 0 )
        pY = scale*W.YY;
    else
        pY = tRadius;
#if ( DIM > 2 )
    pZ = scale*W.ZZ;
#endif
    
    P.set(pX, pY, pZ);
#endif
    return( abs_real(nrm - tRadius) );
}


//------------------------------------------------------------------------------
real SpaceTee::projectOnArm(const Vector& W, Vector& P) const
{
    real  scale, nrm = 0;
#if ( DIM > 1 )
    const real totArmLength = tArmLength+tRadius;
    const real xRel         = (W.XX - tJunction);
    
    //this projection is only valid for W.YY >= 0
    assert_true( W.YY >= 0 );
    
#if ( DIM == 2 )
    nrm = square(xRel);
#elif ( DIM > 2 )
    nrm = square(xRel) + square(W.ZZ);
#endif
    if ( W.YY > totArmLength )
        nrm += square(W.YY-totArmLength);
    if ( nrm > 0 ) {
        nrm   = sqrt(nrm);
        scale = tRadius/nrm;
    }
    else {
        nrm   = 0;
        scale = 0;
    }
    
    real pX, pY = 0, pZ = 0;
    if ( scale != 0 )
        pX = tJunction + scale*xRel;
    else
        pX = tJunction + tRadius;
    
    if ( W.YY > totArmLength )
        pY = totArmLength + scale*(W.YY-totArmLength);
    else
        pY = W.YY;

#if ( DIM > 2 )
    pZ = scale*W.ZZ;
#endif
    P.set(pX, pY, pZ);
#endif
    return( abs_real(nrm - tRadius) );
}


//------------------------------------------------------------------------------
void SpaceTee::projectOnInter(const Vector& W, Vector& P) const
{
    const real xRel = (W.XX - tJunction);
    real pX, pY, pZ;

#if ( DIM == 2 )
    //Points in the intersection area are projected to the corners or to the bottom.
    //The parameterisation of the line of equal distance between a line and a point
    //given by    xl(t) = t       xp(t) = xp         with parameter t
    //            yl(t) = yl      yp(t) = yp
    //is          xi(t) = t
    //            yi(t) = (xp - t)^2 / 2(yp - yl) + (yp^2 - yl^2) / 2(yp - yl)
    //For yl = -tRadius, yp = tRadius and xp = +-tRadius we get
    //            yi(t) = (+-tRadius - t)^2 / 4tWidth
    
    if ( W.XX <= tJunction ) {
        if ( W.YY >= (square(-tRadius - xRel)/(4.*tRadius)) ) {
            //w is projected on the corner
            pX =  tJunction-tRadius;
            pY =  tRadius;
        }
        else {
            //w is projected on the bottom of the base cylinder
            pX =  W.XX;
            pY = -tRadius;
        }
    }
    else {
        if ( W.YY >= (square(tRadius - xRel)/(4.*tRadius)) ) {
            //w is projected on the corner
            pX =  tJunction+tRadius;
            pY =  tRadius;
        }
        else {
            //w is projected on the bottom of the base cylinder
            pX =  W.XX;
            pY = -tRadius;
        }
    }
#endif
    
#if ( DIM >= 3 )
    //w is in the intersection area and projected on the intersection line,
    //which is an ellipse in 3D. The two halfaxis of the ellipse are given
    //by    a = tRadius * sqrt(2)
    //      b = tRadius
    
    //check for pathological cases
    if ( W.XX == 0 ) {
        //The point lies on the short half axis "b" and is
        //always projected to x=0 and z=b or z=-b
        P.set(0,0,std::copysign(tRadius, W.ZZ));
        return;
    }
    M_PI;
    //turn the point, so that the intersection ellipse is in the xz-plane
    real xTurned = ( xRel + std::copysign(W.YY, xRel) ) * M_SQRT1_2;
    real xTurnedSq = square(xTurned);
    
    if ( W.ZZ == 0 ) {
        //The point lies on the long halfaxis "a".
        //In this case the quartic has exactly three solutions, two of which are
        //trivially known: x1=+a or x2=-a, since the half axis are perpendicular
        //to the ellipse. The third solution can be easily found by polynomial
        //division: x3=2*x
        
        //if |x| > a/2, the closest perpendicular projection is on the tips,
        //otherwise the closest projection is solution x3
        if ( abs_real(xTurned)* M_SQRT2 > tRadius ) {
            //we set the final points, already turned back
            pX = std::copysign(tRadius, xTurned) + tJunction;
            pY = tRadius;
            pZ = 0;
        }
        else {
            pX = xTurned*M_SQRT2 + tJunction;
            pY = abs_real(xTurned)*M_SQRT2;
            //we randomly distribute the points to +z or -z
            pZ = RNG.sflip()*sqrt(tRadiusSq - 2*xTurnedSq);
        }
    }
    else {
        real s1, s2, s3, s4;     // solutions of the quartic
        int    nSol;             // number of real solutions
        real   xSol, xSolTurned; // the correct solutions of the quartic and of x
        
        // solve the quartic
        nSol = QuarticSolver::solveQuartic(1, 6, (13-   (2*xTurnedSq +   W.ZZ*W.ZZ) / tRadiusSq),
                                                 (12- 4*(  xTurnedSq +   W.ZZ*W.ZZ) / tRadiusSq),
                                                   4- 2*(  xTurnedSq + 2*W.ZZ*W.ZZ) / tRadiusSq,
                                                  s1, s2, s3, s4);
        
        if ( nSol < 1 )
            ABORT_NOW("Failed to solve quartic for the intersection area.");
        
        // calculate x from t
        xSolTurned = 2.*xTurned/(s1 + 2.);
        
        // turn the point back to it's original position
        xSol = xSolTurned * M_SQRT1_2;
        pX = xSol + tJunction;
        pY = abs_real(xSol);
        pZ = std::copysign(sqrt(tRadiusSq-xSol*xSol), W.ZZ);
    }
    P.set(pX,pY,pZ);
#endif
}


//------------------------------------------------------------------------------
Vector SpaceTee::project(Vector const& W) const
{
    Vector P(W);
#if ( DIM > 1 )
    const real xRel = (W.XX - tJunction); //the x coordinate of w
                                          //relative to tJunction
    if ( inside(W) )
    {
#if ( DIM == 2 )
        if ( W.YY > tRadius ) {
            //w is inside the arm
            projectOnArm(W, P);
        } else if ( (xRel >= -tRadius) && (xRel <= tRadius) && (W.YY >= 0) ) {
            //w is inside the intersection area
            projectOnInter(W, P);
        }
        else {
            // w is inside the base cylinder
            projectOnBase(W, P);
        }
#endif
        
#if ( DIM > 2 )
        if ( (xRel >  tRadius)
            || (xRel < -tRadius)
            || (W.YY < 0)
            || (W.YY*W.YY*(tRadiusSq - xRel*xRel) < square(xRel*W.ZZ)) )
        {
            //w is projected on the base cylinder, if
            //    the point is on the right side of the arm
            //or  the point is on the left side of the arm
            //or  the point is in the lower half of the base cylinder
            //or  the y coordinate of the point is low enough, so that it can
            //    be projected perpendicularly on the base cylinder:
            //    y < z |xRel| / sqrt( tRadius^2 - xRel^2 )
            projectOnBase(W, P);
        } else if ( square(W.YY*xRel) + square(W.YY*W.ZZ) > square(xRel*tRadius) ) {
            //w is projected on the arm, if
            //the y coordinate of the point is greater than the y coordinate of
            //the corresponding point on the intersection ellipse:
            //y > r |xRel| / sqrt( xRel^2 + z^2 )
            projectOnArm(W, P);
        }
        else {
            //w is projected on the intersection ellipse
            projectOnInter(W, P);
        }
#endif
    }
    else 
    {
        //point w is outside the tee
        Vector pArm;                      //projection of w on the arm
        real dBase = projectOnBase(W, P); //distance of w to the base cylinder
        //all points with y<0 are projected on the base cylinder
        if ( W.YY >= 0 )
        {
            //check if w is closer to base or arm
            if ( dBase <= projectOnArm(W, pArm) )
                return P;
            else
                return pArm;
        }
    }
#endif
    return P;
}


//------------------------------------------------------------------------------

void SpaceTee::write(Outputter& out) const
{
    out.put_characters("tee", 16);
    out.writeUInt16(4);
    out.writeFloat(tLength);
    out.writeFloat(tRadius);
    out.writeFloat(tJunction);
    out.writeFloat(tArmLength);
}


void SpaceTee::setLengths(const real len[])
{
    tLength    = len[0];
    tRadius    = len[1];
    tJunction  = len[2];
    tArmLength = len[3];
    update();
}

void SpaceTee::read(Inputter& in, Simul&, ObjectTag)
{
    real len[8] = { 0 };
    read_data(in, len, "tee");
    setLengths(len);
}


//------------------------------------------------------------------------------
//                         OPENGL  DISPLAY
//------------------------------------------------------------------------------

#ifdef DISPLAY

#include "gle.h"
using namespace gle;

bool SpaceTee::draw() const
{
#if ( DIM == 2 )
    constexpr size_t fin = 8 * gle::finesse;
    GLfloat c[2*fin+1], s[2*fin+1];
    gle::circle(2*fin, c, s, tRadius);
    
    glBegin(GL_LINE_LOOP);
    //the upper side from the tJunction to the left
    gleVertex(  tJunction-tRadius, tRadius, 0 );
    gleVertex( -tLength,          tRadius, 0 );
    
    //the left cap
    for ( size_t n = 0; n < fin; ++n )
        gleVertex( -tLength-s[n], c[n], 0 );
    
    //the lower side from left to right
    gleVertex( -tLength, -tRadius, 0 );
    gleVertex(  tLength, -tRadius, 0 );
    
    //the right cap
    for ( size_t n = 0; n < fin; ++n )
        gleVertex( tLength+s[n], -c[n], 0 );
    
    //the upper side from the right to the tJunction
    gleVertex( tLength,          tRadius, 0 );
    gleVertex( tJunction+tRadius, tRadius, 0 );
    
    //the right side of the arm
    gleVertex( tJunction+tRadius, tRadius+tArmLength, 0 );
    
    //the cap of the arm
    for ( size_t n = 0; n < fin; ++n )
        gleVertex( tJunction+c[n], tArmLength+tRadius+s[n], 0 );
    
    //the left side of the arm
    gleVertex( tJunction-tRadius, tRadius+tArmLength, 0 );

    glEnd();
#endif
    
#if ( DIM > 2 )

    const GLenum glp1 = GL_CLIP_PLANE4;
    const GLenum glp2 = GL_CLIP_PLANE5;
    const GLdouble planeZ[] = { 0, 0, 1, 0 };
    const GLdouble plane1[] = { 0, -sqrt(0.5), sqrt(0.5), 0 };
    const GLdouble plane2[] = { 0,  sqrt(0.5), sqrt(0.5), 0 };

    GLdouble L = tLength/tRadius;
    GLdouble J = tJunction/tRadius;
    GLdouble A = tArmLength/tRadius;

    glEnable(glp1);
    glPushMatrix();
    gleScale(tRadius);

    //right side:
    glPushMatrix();
    glRotated(90, 0, 1, 0);
    glTranslated(0, 0, J);
    glClipPlane(glp1, plane1);
    gleTube0(0, L-J, 1);
    glTranslated(0, 0, L-J);
    glClipPlane(glp1, planeZ);
    gleSphere8B();
    glPopMatrix();

    //left side:
    glPushMatrix();
    glRotated(90, 0, -1, 0);
    glTranslated(0, 0, -J);
    glClipPlane(glp1, plane1);
    gleTube0(0, L+J, 1);
    glTranslated(0, 0, L+J);
    glClipPlane(glp1, planeZ);
    gleSphere8B();
    glTranslated(0, 0, -L-J);
    glPopMatrix();

    //the arm:
    glPushMatrix();
    glTranslated(J, 0, 0);
    glRotated(-90, 1, 0, 0);
    glRotated( 90, 0, 0, 1);
    glEnable(glp2);
    glClipPlane(glp1, plane1);
    glClipPlane(glp2, plane2);
    gleTube0(0, A+1, 1);
    glTranslated(0, 0, A+1);
    glDisable(glp2);
    glClipPlane(glp1, planeZ);
    gleSphere8B();
    glPopMatrix();

    glPopMatrix();
    glDisable(glp1);
    glDisable(glp2);
#endif
    
    return true;
}

#else

bool SpaceTee::draw() const
{
    return false;
}


#endif


