// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "dim.h"
#include "space_tee.h"
#include "exceptions.h"
#include "quartic_solver.h"
#include "random.h"



SpaceTee::SpaceTee(const SpaceProp* p)
: Space(p)
{
    if ( DIM == 1 )
        throw InvalidParameter("tee cannot be used in 1D");
}


void SpaceTee::resize()
{
    tLength    = length(0);
    tWidth     = length(1);
    tJunction  = length(2);
    tArmLength = length(3);
    tWidthSq   = tWidth * tWidth;
    if ( tLength <= 0 || tWidth <= 0 || tArmLength < 0 )
        throw InvalidParameter("Space tee can't have negative length, arm length or radius.");
    if ( fabs(tJunction)+tWidth > tLength )
        throw InvalidParameter("Space tee: the position of the branch plus the radius must lie within the length of the T.");
}

//------------------------------------------------------------------------------

void SpaceTee::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(-tWidth-tLength,-tWidth, -tWidth );
    sup.set( tWidth+tLength, tArmLength+2*tWidth, tWidth );
}


real SpaceTee::volume() const
{
#if ( DIM == 1 )
    return 0;
#elif ( DIM == 2 )
    real base = 4 * tLength * tWidth + M_PI * tWidthSq;
    real arm  = 2 * tArmLength * tWidth + M_PI_2 * tWidthSq;
    return( base + arm );
#else
    //the complete base cylinder
    real base  = 2 * tLength * M_PI * tWidthSq + 4./3. * M_PI * tWidth * tWidthSq;
    //the part of the arm with y > tWidth
    real arm   =  tArmLength * M_PI * tWidthSq + 2./3. * M_PI * tWidth * tWidthSq;
    //the part of the arm with y < tWidth without the intersection with the base 
    real inter = ( M_PI - 8./3. )*tWidth*tWidthSq;
    return( base + arm + inter );
#endif
}


//------------------------------------------------------------------------------
bool SpaceTee::inside(Vector const& w) const
{
    real nrmSq      = 0;
    const real x    = fabs( w[0] );
    const real xRel = (w[0] - tJunction);
    
    //check if w is inside the base cylinder
    if ( x > tLength )
        nrmSq = square(x - tLength);
#if ( DIM == 2 )
    nrmSq += w[1]*w[1];
#elif ( DIM >= 3 )
    nrmSq += w[1]*w[1] + w[2]*w[2];
#endif
    if ( nrmSq <= tWidthSq ) return( true );
    
    //check if w is inside the arm
    if ( w[1] >= 0)
    {
        nrmSq = 0;
        if ( w[1] > tArmLength+tWidth )
            nrmSq = square(w[1] - (tArmLength+tWidth));
#if ( DIM == 2 )
        nrmSq += xRel*xRel;
#elif ( DIM >= 3 )
        nrmSq += xRel*xRel + w[2]*w[2];
#endif
        return( nrmSq <= tWidthSq );
    }
    
    return( false );
}


//------------------------------------------------------------------------------
real SpaceTee::projectOnBase(const Vector w, Vector& p) const
{
    real scale, nrm = 0;
#if ( DIM > 2 )
    nrm = w[1]*w[1] + w[2]*w[2];
#elif ( DIM > 1 )
    nrm = w[1]*w[1];
#endif
    if ( w[0] >  tLength )
        nrm += square(w[0] - tLength);
    else if ( w[0] < -tLength )
        nrm += square(w[0] + tLength);
    
    if ( nrm > 0 ) {
        nrm   = sqrt(nrm);
        scale = tWidth/nrm;
    }
    else {
        nrm   = 0;
        scale = 0;
    }
    
    if ( w[0] >  tLength )
        p[0] =  tLength + scale*(w[0] - tLength);
    else if ( w[0] < -tLength )
        p[0] = -tLength + scale*(w[0] + tLength);
    else
        p[0] = w[0];
    
    if ( scale != 0 )
        p[1] = scale*w[1];
    else
        p[1] = tWidth;
    
#if ( DIM > 2 )
    p[2] = scale*w[2];
#endif
    
    return( fabs(nrm - tWidth) );
}


//------------------------------------------------------------------------------
real SpaceTee::projectOnArm(const Vector w, Vector& p) const
{
    real  scale, nrm = 0;
    const real totArmLength = tArmLength+tWidth;
    const real xRel         = (w[0] - tJunction);
    
    //this projection is only valid for w[1] >= 0
    assert_true( w[1] >= 0 );
    
#if ( DIM == 2 )
    nrm = xRel*xRel;
#elif ( DIM > 2 )
    nrm = xRel*xRel + w[2]*w[2];
#endif
    if ( w[1] > totArmLength )
        nrm += square(w[1]-totArmLength);
    
    if ( nrm > 0 ) {
        nrm   = sqrt(nrm);
        scale = tWidth/nrm;
    }
    else {
        nrm   = 0;
        scale = 0;
    }
    
    if ( scale != 0 )
        p[0] = tJunction + scale*xRel;
    else
        p[0] = tJunction + tWidth;
    
    if ( w[1] > totArmLength )
        p[1] = totArmLength + scale*(w[1]-totArmLength);
    else
        p[1] = w[1];
    
#if ( DIM > 2 )
    p[2] = scale*w[2];
#endif
    
    return( fabs(nrm - tWidth) );
}


//------------------------------------------------------------------------------
void SpaceTee::projectOnInter(const Vector w, Vector& p) const
{
    const real xRel = (w[0] - tJunction);
    
#if ( DIM == 2 )
    //Points in the intersection area are projected to the corners or to the bottom.
    //The parameterisation of the line of equal distance between a line and a point
    //given by    xl(t) = t       xp(t) = xp         with parameter t
    //            yl(t) = yl      yp(t) = yp
    //is          xi(t) = t
    //            yi(t) = (xp - t)^2 / 2(yp - yl) + (yp^2 - yl^2) / 2(yp - yl)
    //For yl = -tWidth, yp = tWidth and xp = +-tWidth we get
    //            yi(t) = (+-tWidth - t)^2 / 4tWidth
    
    if ( w[0] <= tJunction ) {
        if ( w[1] >= (square(-tWidth - xRel)/(4.*tWidth)) ) {
            //w is projected on the corner
            p[0] =  tJunction-tWidth;
            p[1] =  tWidth;
        }
        else {
            //w is projected on the bottom of the base cylinder
            p[0] =  w[0];
            p[1] = -tWidth;
        }
    }
    else {
        if ( w[1] >= (square(tWidth - xRel)/(4.*tWidth)) ) {
            //w is projected on the corner
            p[0] =  tJunction+tWidth;
            p[1] =  tWidth;
        }
        else {
            //w is projected on the bottom of the base cylinder
            p[0] =  w[0];
            p[1] = -tWidth;
        }
    }
#endif
    
#if ( DIM >= 3 )
    //w is in the intersection area and projected on the intersection line,
    //which is an ellipse in 3D. The two halfaxis of the ellipse are given
    //by    a = tWidth * sqrt(2)
    //      b = tWidth
    
    //check for pathological cases
    if ( w[0] == 0 ) {
        //The point lies on the short half axis "b" and is
        //always projected to x=0 and z=b or z=-b
        p[0] = 0;
        p[1] = 0;
        p[2] = std::copysign(tWidth, w[2]);
        return;
    }
    
    //turn the point, so that the intersection ellipse is in the xz-plane
    real   xTurned;
    if ( xRel >= 0 )
        xTurned   =  (xRel + w[1]) / sqrt(2);
    else
        xTurned   =  (xRel - w[1]) / sqrt(2);
    
    real xTurnedSq = square(xTurned);
    
    if ( w[2] == 0 ) {
        //The point lies on the long halfaxis "a".
        //In this case the quartic has exactly three solutions, two of which are
        //trivially known: x1=+a or x2=-a, since the half axis are perpendicular
        //to the ellipse. The third solution can be easily found by polynomial
        //division: x3=2*x
        
        //if |x| > a/2, the closest perpendicular projection is on the tips,
        //otherwise the closest projection is solution x3
        if ( fabs(xTurned)*sqrt(2) > tWidth ) {
            //we set the final points, already turned back
            p[0] = std::copysign(tWidth, xTurned) + tJunction;
            p[1] = tWidth;
            p[2] = 0;
        }
        else {
            p[0] = xTurned*sqrt(2) + tJunction;
            p[1] = fabs(xTurned)*sqrt(2);
            //we randomly distribute the points to +z or -z
            p[2] = RNG.sflip()*sqrt(tWidthSq - 2.*xTurnedSq);
        }
        return;
        
    }
    else {
        real s1, s2, s3, s4;     // solutions of the quartic
        int    nSol;             // number of real solutions
        real   xSol, xSolTurned; // the correct solutions of the quartic and of x
        
        // solve the quartic
        nSol = QuarticSolver::solveQuartic(1, 6, (13-   (2*xTurnedSq +   w[2]*w[2]) / tWidthSq),
                                                 (12- 4*(  xTurnedSq +   w[2]*w[2]) / tWidthSq),
                                                   4- 2*(  xTurnedSq + 2*w[2]*w[2]) / tWidthSq,
                                                  s1, s2, s3, s4);
        
        if ( nSol < 1 )
            ABORT_NOW("Failed to solve quartic for the intersection area.");
        
        // calculate x from t
        xSolTurned = 2.*xTurned/(s1 + 2.);
        
        // turn the point back to it's original position
        xSol = xSolTurned / sqrt(2);
        p[0] = xSol + tJunction;
        p[1] = fabs(xSol);
        if ( w[2] > 0 )
            p[2] =  sqrt( tWidthSq - xSol*xSol );
        else
            p[2] = -sqrt( tWidthSq - xSol*xSol );
    }
#endif
}


//------------------------------------------------------------------------------
Vector SpaceTee::project(Vector const& w) const
{
    Vector p;
    const real xRel = (w[0] - tJunction); //the x coordinate of w
                                          //relative to tJunction
    if ( inside(w) )
    {
#if ( DIM == 2 )
        if ( w[1] > tWidth ) {
            //w is inside the arm
            projectOnArm(w, p);
        } else if ( (xRel >= -tWidth) && (xRel <= tWidth) && (w[1] >= 0) ) {
            //w is inside the intersection area
            projectOnInter(w, p);
        }
        else {
            // w is inside the base cylinder
            projectOnBase(w, p);
        }
#endif
        
#if ( DIM > 2 )
        if (  (xRel >  tWidth)
              || (xRel < -tWidth)
              || (w[1] < 0)
              || (w[1]*w[1]*(tWidthSq - xRel*xRel) < xRel*xRel*w[2]*w[2]) )
        {
            //w is projected on the base cylinder, if 
            //    the point is on the right side of the arm
            //or  the point is on the left side of the arm
            //or  the point is in the lower half of the base cylinder
            //or  the y coordinate of the point is low enough, so that it can
            //    be projected perpendicularly on the base cylinder:
            //    y < z |xRel| / sqrt( tWidth^2 - xRel^2 )
            projectOnBase(w, p);
        } else if ( w[1]*w[1]*(xRel*xRel + w[2]*w[2]) > xRel*xRel*tWidthSq ) {
            //w is projected on the arm, if
            //the y coordinate of the point is greater than the y coordinate of
            //the corresponding point on the intersection ellipse:
            //y > r |xRel| / sqrt( xRel^2 + z^2 )
            projectOnArm(w, p);
        }
        else {
            //w is projected on the intersection ellipse
            projectOnInter(w, p);
        }
#endif
        
    }
    else 
    {
        //point w is outside the tee
        
        Vector pArm;                      //projection of w on the arm
        real dBase = projectOnBase(w, p); //distance of w to the base cylinder
        
        //all points with y<0 are projected on the base cylinder
        if ( w[1] >= 0 )
        {
            //check if w is closer to base or arm
            if ( dBase <= projectOnArm(w, pArm) )
                return p;
            else
                return pArm;
        }
    }
    return p;
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
    gle::circle(2*fin, c, s, tWidth);
    
    glBegin(GL_LINE_LOOP);
    //the upper side from the tJunction to the left
    gleVertex(  tJunction-tWidth, tWidth, 0 );
    gleVertex( -tLength,          tWidth, 0 );
    
    //the left cap
    for ( size_t n = 0; n < fin; ++n )
        gleVertex( -tLength-s[n], c[n], 0 );
    
    //the lower side from left to right
    gleVertex( -tLength, -tWidth, 0 );
    gleVertex(  tLength, -tWidth, 0 );
    
    //the right cap
    for ( size_t n = 0; n < fin; ++n )
        gleVertex( tLength+s[n], -c[n], 0 );
    
    //the upper side from the right to the tJunction
    gleVertex( tLength,          tWidth, 0 );
    gleVertex( tJunction+tWidth, tWidth, 0 );
    
    //the right side of the arm
    gleVertex( tJunction+tWidth, tWidth+tArmLength, 0 );
    
    //the cap of the arm
    for ( size_t n = 0; n < fin; ++n )
        gleVertex( tJunction+c[n], tArmLength+tWidth+s[n], 0 );
    
    //the left side of the arm
    gleVertex( tJunction-tWidth, tWidth+tArmLength, 0 );

    glEnd();
#endif
    
#if ( DIM > 2 )

    const GLenum glp1 = GL_CLIP_PLANE4;
    const GLenum glp2 = GL_CLIP_PLANE5;
    const GLdouble planeZ[] = { 0, 0, 1, 0 };
    const GLdouble plane1[] = { 0, -sqrt(0.5), sqrt(0.5), 0 };
    const GLdouble plane2[] = { 0,  sqrt(0.5), sqrt(0.5), 0 };

    GLdouble L = tLength/tWidth;
    GLdouble J = tJunction/tWidth;
    GLdouble A = tArmLength/tWidth;

    glEnable(glp1);
    glPushMatrix();
    gleScale(tWidth);

    //right side:
    glPushMatrix();
    glRotated(90, 0, 1, 0);
    glTranslated(0, 0, J);
    glClipPlane(glp1, plane1);
    gleTubeF(0, L-J, 1);
    glTranslated(0, 0, L-J);
    glClipPlane(glp1, planeZ);
    gleSphere8B();
    glPopMatrix();

    //left side:
    glPushMatrix();
    glRotated(90, 0, -1, 0);
    glTranslated(0, 0, -J);
    glClipPlane(glp1, plane1);
    gleTubeF(0, L+J, 1);
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
    gleTubeF(0, A+1, 1);
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


