// Cytosim was created by Francois Nedelec. Copyright Cambridge University 2020

#include "dim.h"
#include "sim.h"
#include "simul.h"
#include "display2.h"
#include "display_color.h"
#include "modulo.h"

#include "fake.h"
#include "line_disp.h"
#include "point_disp.h"

#include "opengl.h"
#include "gle.h"
#include "gle_color_list.h"
#include "glut.h"

using namespace gle;
extern Modulo const* modulo;

//------------------------------------------------------------------------------

Display2::Display2(DisplayProp const* dp) : Display(dp)
{
}


void Display2::drawSimul(Simul const& sim)
{
    glDepthMask(GL_FALSE);
    glDisable(GL_LIGHTING);
    glDisable(GL_CULL_FACE);
    drawFields(sim.fields);
    
    glEnable(GL_LIGHTING);
#if ( DIM > 2 )
    glEnable(GL_CULL_FACE);
    glDepthMask(GL_TRUE);
#endif
    drawSpaces(sim.spaces);
    
    glDisable(GL_LIGHTING);
    glDisable(GL_CULL_FACE);

    if ( prop->couple_select & 1 )
        drawCouplesF(sim.couples);
    
    if ( prop->couple_select & 2 )
        drawCouplesA(sim.couples);
    
    if ( prop->single_select & 1 )
        drawSinglesF(sim.singles);
    
    drawFibers(sim.fibers);
    
#if ( DIM == 3 )
    glEnable(GL_LIGHTING);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
#else
    glDisable(GL_LIGHTING);
#endif
    
    drawBeads(sim.beads);
    drawSolids(sim.solids);
    drawSpheres(sim.spheres);
    
#if ( DIM == 3 )
    glDisable(GL_LIGHTING);
    glDisable(GL_CULL_FACE);
#endif

    if ( prop->couple_select & 4 )
        drawCouplesB(sim.couples);
    
    if ( prop->single_select & 2 )
        drawSinglesA(sim.singles);
    
#if ( DIM == 3 )
    glEnable(GL_LIGHTING);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
#endif

    drawOrganizers(sim.organizers);
    drawMisc(sim);
}

//------------------------------------------------------------------------------
#pragma mark -

void Display2::drawBead(Bead const& obj)
{
    const PointDisp * disp = obj.prop->disp;

    // display center:
    if ( disp->style & 2  && disp->perceptible )
    {
        bodyColor(obj);
        drawObject(obj.position(), disp->size, gle::sphere1);
    }
    
#if ( DIM == 2 )
    // display outline:
    if ( disp->style & 4 )
    {
        bodyColorF(obj).load();
        lineWidth(disp->width);
        drawFlat(obj.position(), obj.radius(), gle::circle);
    }
#endif
}


/**
 Display a semi-transparent disc / sphere
 */
void Display2::drawBeadT(Bead const& obj)
{
    const PointDisp * disp = obj.prop->disp;
    
    if ( disp->style & 1 )
    {
        drawBallT(obj.position(), obj.radius(), bodyColorF(obj));
    }
}

//------------------------------------------------------------------------------
#pragma mark -

void Display2::drawSolid(Solid const& obj)
{
    const PointDisp * disp = obj.prop->disp;
    
    //display points:
    if ( disp->style & 2  &&  disp->size > 0  && disp->perceptible )
    {
        bodyColor(obj);
        for ( size_t i = 0; i < obj.nbPoints(); ++i )
            drawObject(obj.posP(i), disp->size, gle::hedron(obj.radius(i)>0));
    }
    
    //display outline of spheres
    if ( disp->style & 4 )
    {
#if ( DIM == 2 )
        glDisable(GL_LIGHTING);
        lineWidth(disp->width);
        bodyColorF(obj).load();
        for ( size_t i = 0; i < obj.nbPoints(); ++i )
        {
            if ( obj.radius(i) > 0 )
                drawFlat(obj.posP(i), obj.radius(i), gle::circle);
        }
#elif ( DIM == 3 )
        //special display for ParM simulations (DYCHE 2006; KINETOCHORES 2019)
        if ( obj.mark()  &&  obj.nbPoints() > 1 )
        {
            glEnable(GL_LIGHTING);
            bodyColor(obj);
            //gle::gleObject(obj.posP(0), obj.diffPoints(1, 0), obj.radius(0), gle::circle);
            glPushMatrix();
            Vector A = obj.posP(0), B = obj.posP(1);
            gle::transAlignZ(0.5*(A+B), obj.radius(0), A-B);
            gle::cylinderZ();
            glPopMatrix();
        }
#endif
    }
    
    //print the number for each Solid
    if ( disp->style & 8 )
    {
        char tmp[8];
        bodyColorF(obj).load();
        snprintf(tmp, sizeof(tmp), "%u", obj.identity());
        gle::drawText(obj.posP(0), tmp, GLUT_BITMAP_HELVETICA_10);
    }
    
    //draw polygon around vertices of Solid
    if ( disp->style & 16 )
    {
        lineWidth(disp->width);
        bodyColorF(obj).load();
        glBegin(GL_LINE_LOOP);
        for ( size_t i = 0; i < obj.nbPoints(); ++i )
            gleVertex(obj.posPoint(i));
        glEnd();
    }
}

/**
 Display a semi-transparent disc / sphere
 */
void Display2::drawSolidT(Solid const& obj, size_t inx)
{
    const PointDisp * disp = obj.prop->disp;

    if (( disp->style & 1 ) & ( obj.radius(inx) > 0 ))
    {
        Vector X = obj.posP(inx);
        size_t near[3];
        size_t num = obj.closestSpheres(inx, near[0], near[1], near[2]);
        //printf("nearest balls to %lu / %lu are %lu %lu %lu\n", inx, obj.nbPoints(), near[0], near[1], near[2]);
        // set clipping planes with nearest balls
        for ( size_t i = 0; i < num; ++i )
        {
            size_t J = near[i];
            Vector P = obj.posPoint(J);
            real A = ( square(obj.radius(inx)) - square(obj.radius(J)) ) / distanceSqr(X, P);
            GLenum glp = GL_CLIP_PLANE5 - i;
            glEnable(glp);
            gle::setClipPlane(glp, normalize(X-P), (0.5-0.5*A)*X+(0.5+0.5*A)*P);
        }
        drawBallT(X, obj.radius(inx), bodyColorF(obj));
        glDisable(GL_CLIP_PLANE3);
        glDisable(GL_CLIP_PLANE4);
        glDisable(GL_CLIP_PLANE5);
    }
}

//------------------------------------------------------------------------------
#pragma mark -

void Display2::drawSphere(Sphere const& obj)
{
    const PointDisp * disp = obj.prop->disp;
    
    //display center and surface points
    if ( disp->style & 2  &&  disp->perceptible )
    {
        bodyColor(obj);
        drawObject(obj.posP(0), disp->size, gle::star);
        for ( size_t i = obj.nbRefPoints; i < obj.nbPoints(); ++i )
            drawObject(obj.posP(i), disp->size, gle::sphere1);
    }
    
    //display reference points
    if ( disp->style & 8  &&  disp->perceptible )
    {
        bodyColor(obj);
        for ( size_t i = 1; i < obj.nbRefPoints; ++i )
            drawObject(obj.posP(i), disp->size, gle::tetrahedron);
    }
}


void Display2::drawSphereT(Sphere const& obj)
{
    const PointDisp * disp = obj.prop->disp;

    if ( disp->style & 5 )
    {
        const Vector C = obj.posP(0);
#if ( DIM < 3 )
        bodyColorF(obj).load();
        if ( disp->style & 1 )
            drawFlat(C, obj.radius(), gle::circle);
        if ( disp->style & 2 )
            drawFlat(C, obj.radius(), gle::discUp);
#else
        /* Note: The rotation matrix for the sphere calculated below from the
         reference points, includes scaling by the radius of the sphere.
         We then use a primitive for a sphere of radius 1.
         */
        bodyColorF(obj).load_both();
        Display::drawSphereT(C, obj.posP(1)-C, obj.posP(2)-C, obj.posP(3)-C, disp->style);
#endif
    }
}

//------------------------------------------------------------------------------
void Display2::drawOrganizer(Organizer const& obj) const
{
    PointDisp const* disp = obj.disp();
    
    if ( !disp )
        return;

    if ( disp->style & 2 )
    {
        Vector P, Q;
        glDisable(GL_LIGHTING);
        
        for ( size_t i = 0; obj.getLink(i, P, Q); ++i )
            disp->drawI(P);

        lineWidth(disp->width);
        bodyColorF(disp, obj.signature()).load();
        glBegin(GL_LINES);
        for ( size_t i = 0; obj.getLink(i, P, Q); ++i )
        {
            if ( modulo ) modulo->fold(Q, P);
            gleVertex(P);
            gleVertex(Q);
        }
        glEnd();
    }

    /**
     This display the Solid connecting two Aster as a spindle.
     Used for Cleo Kozlowski simulation of C. elegans (2007)
     */
    if ( disp->style & 1 && obj.tag() == Fake::TAG )
    {
        Solid const* sol = Solid::toSolid(obj.core());
        if ( sol && sol->nbPoints() >= 4 )
        {
#if ( DIM == 3 )
            bodyColor(*sol);
            glEnable(GL_LIGHTING);
            glPushMatrix();
            Vector3 a = 0.5*(sol->posP(0) + sol->posP(2));
            Vector3 b = 0.5*(sol->posP(1) + sol->posP(3));
            stretchAlignZ(a, b, 1);
            gle::dualPass(gle::barrel);
            glPopMatrix();
            glDisable(GL_LIGHTING);
#else
            glDisable(GL_LIGHTING);
            bodyColorF(*sol).load();
            glBegin(GL_LINES);
            for ( size_t i = 0; i < sol->nbPoints(); ++i )
                gleVertex(sol->posPoint(i));
            glEnd();
#endif
        }
    }
}

//------------------------------------------------------------------------------
#pragma mark -

/**
 display the attached position of free singles
 */
void Display2::drawSinglesF(const SingleSet & set) const
{
    for ( Single * obj=set.firstF(); obj ; obj=obj->next() )
        obj->disp()->drawF(obj->posFoot());
}


void Display2::drawSinglesA(const SingleSet & set) const
{
    // display the Hands
    for ( Single * obj=set.firstA(); obj ; obj=obj->next() )
    {
        const PointDisp * disp = obj->disp();
        if ( disp->perceptible  &&  obj->fiber()->disp->visible )
        {
            Vector ph = obj->posHand();
            
            disp->drawA(ph);
            
            if ( obj->hasForce() && disp->width > 0 )
            {
                Vector ps = obj->sidePos();
                Vector pf = obj->posFoot();
                if ( modulo )
                {
                    modulo->fold(pf, ph);
                    modulo->fold(ps, ph);
                }
                
                disp->color.load();
#if ( DIM >= 3 )
                gleTube(pf, ph, disp->width*sFactor, gle::cone);
                //drawCone(pf, ph-pf, disp->width*sFactor);
#else
                gle::drawBand(ph, disp->width*sFactor, ps, disp->width*sFactor);
                gle::drawBand(ps, disp->width*sFactor, disp->color, pf, disp->width*sFactor, disp->color.alpha_scaled(0.5));
#endif
            }
        }
    }
}

//------------------------------------------------------------------------------
#pragma mark -

/**
 Always display Hand1 of the Couple.
 */
void Display2::drawCouplesF1(CoupleSet const& set) const
{
    for ( Couple * cx = set.firstFF() ; cx ; cx=cx->next() )
    {
        if ( cx->active() )
            cx->disp1()->drawF(cx->posFree());
        else
            cx->disp1()->drawI(cx->posFree());
    }
}


PointDisp const* Couple::disp12() const
{
    if ( disp1()->visible )
        return disp1();
    else
        return disp2();
}


PointDisp const* Couple::disp21() const
{
    if ( disp2()->visible )
        return disp2();
    else
        return disp1();
}


/**
 Display either Hand1 or Hand2, exposing both sides with equal chances.
 This gives the impression that Couple flicker randomly between frames,
 as if they were two-sided balls 'rotating' very fast.
 */
void Display2::drawCouplesF2(CoupleSet const& set) const
{
    Couple * nxt;
    Couple * obj = set.firstFF();
    // this loop is unrolled, processing objects 2 by 2:
    if ( set.sizeFF() & 1 )
    {
        nxt = obj->next();
        if ( obj->active() )
            obj->disp12()->drawF(obj->posFree());
        else
            obj->disp12()->drawI(obj->posFree());
        obj = nxt;
    }
    while ( obj )
    {
        nxt = obj->next();
        if ( obj->active() )
            obj->disp21()->drawF(obj->posFree());
        else
            obj->disp21()->drawI(obj->posFree());
        obj = nxt->next();
        if ( nxt->active() )
            nxt->disp12()->drawF(nxt->posFree());
        else
            nxt->disp12()->drawI(nxt->posFree());
    }
}


void Display2::drawCouplesA(CoupleSet const& set) const
{
    // display bound couples
    for ( Couple * cx=set.firstAF(); cx ; cx=cx->next() )
    {
        if ( cx->fiber1()->disp->visible )
        {
            if ( cx->active() )
                cx->disp1()->drawF(cx->posHand1());
            else
                cx->disp1()->drawI(cx->posHand1());
        }
    }
    
    for ( Couple * cx=set.firstFA(); cx ; cx=cx->next() )
    {
        if ( cx->fiber2()->disp->visible )
        {
            if ( cx->active() )
                cx->disp2()->drawF(cx->posHand2());
            else
                cx->disp2()->drawI(cx->posHand2());
            
        }
    }
}


void Display2::drawCoupleB(Couple const* cx) const
{
    const PointDisp * pd1 = cx->disp1();
    const PointDisp * pd2 = cx->disp2();
    
    Vector p1 = cx->posHand1();
    Vector p2 = cx->posHand2();
    
    if ( modulo )
        modulo->fold(p2, p1);
    
    if ( pd1 == pd2 )
    {
        if ( pd1->perceptible )
        {
            pd1->color.load();
#if ( 1 )
            /*
             Can shift positions towards the minus-end by couple's length
             to create an effect to highlight the configuration:
             ///// on antiparallel fibers
             >>>>> on parallel fibers
             */
            real len = 0.5 * cx->prop->length;
            lineWidth(pd1->width);
            glBegin(GL_LINE_STRIP);
            Vector d1 = len * cx->dirFiber1();
            Vector d2 = len * cx->dirFiber2();
            gleVertex(p1);
            gleVertex(0.5*((p1+d1)+(p2+d2)));
            gleVertex(p2);
            glEnd();
#else
            lineWidth(pd1->width);
            glBegin(GL_LINES);
            gleVertex(p1);
            gleVertex(p2);
            glEnd();
#endif
        }
    }
    else
    {
        if ( pd1->perceptible || pd2->perceptible )
        {
#if ( 1 )
            /*
             Can shift positions towards the minus-end by couple's length
             to create an effect to highlight the configuration:
             ///// on antiparallel fibers
             >>>>> on parallel fibers
             */
            real len = 0.5 * cx->prop->length;
            lineWidth(pd1->width);
            glBegin(GL_LINE_STRIP);
            Vector d1 = len * cx->dirFiber1();
            Vector d2 = len * cx->dirFiber2();
            gleVertex(p1);
            gleVertex(0.5*((p1+d1)+(p2+d2)));
            gleVertex(p2);
            glEnd();
#else
            lineWidth(pd1->width);
            glBegin(GL_LINES);
            gleVertex(p1);
            gleVertex(p2);
            glEnd();
#endif
        }
    }
    
    if ( cx->active() )
    {
        pd1->drawA(p1);
        pd2->drawA(p2);
    }
    else
    {
        pd1->drawI(p1);
        pd2->drawI(p2);
    }
}

