// Cytosim was created by Francois Nedelec. Copyright Cambridge University 2020

#include "dim.h"
#include "sim.h"
#include "simul.h"
#include "display1.h"
#include "modulo.h"

#include "opengl.h"
#include "gle_color_list.h"
#include "glut.h"
#include "gle.h"

#include "fake.h"
#include "fiber_disp.h"
#include "line_disp.h"
#include "point_disp.h"
#include "display_color.h"
#include "vector_float.h"

using namespace gle;
extern Modulo const* modulo;


#define ENABLE_EXPLODE_DISPLAY ( DIM < 3 )


Display1::Display1(DisplayProp const* dp) : Display(dp)
{
}


void Display1::drawSimul(Simul const& sim)
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
    glEnableClientState(GL_COLOR_ARRAY);

    if (( prop->couple_select & 1 ) && ( sim.couples.sizeFF() > 0 ))
        drawCouplesF(sim.couples);

    if (( prop->single_select & 1 ) && ( sim.singles.sizeF() > 0 ))
        drawSinglesF(sim.singles);
    
    glDisableClientState(GL_COLOR_ARRAY);

    drawFibers(sim.fibers);

    glEnableClientState(GL_COLOR_ARRAY);
    if (( prop->couple_select & 2 ) && ( sim.couples.sizeA() > 0 ))
        drawCouplesA(sim.couples);
    glDisableClientState(GL_COLOR_ARRAY);

#if ( DIM == 3 )
    
    glEnable(GL_LIGHTING);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);

#endif
    
    drawBeads(sim.beads);
    drawSolids(sim.solids);
    drawSpheres(sim.spheres);
    
#if ( DIM == 3 )
    
    glDisable(GL_LIGHTING);
    glDisable(GL_CULL_FACE);

#endif
    
    glEnableClientState(GL_COLOR_ARRAY);
    if (( prop->couple_select & 4 ) && ( sim.couples.sizeAA() > 0 ))
        drawCouplesB(sim.couples);
    
    if (( prop->single_select & 2 ) && ( sim.singles.sizeA() > 0 ))
        drawSinglesA(sim.singles);
    glDisableClientState(GL_COLOR_ARRAY);

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


void Display1::drawFiber(Fiber const& fib)
{
#if ENABLE_EXPLODE_DISPLAY
    //translate whole display to display the Fiber
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    gle::translate(0, fib.disp->explode_shift, 0);
#endif
    
    Display::drawFiber(fib);

#if ENABLE_EXPLODE_DISPLAY
    glPopMatrix();
#endif
}

//------------------------------------------------------------------------------
#pragma mark -

void Display1::drawBead(Bead const& obj)
{
    const PointDisp * disp = obj.prop->disp;

    // display center:
    if ( disp->style & 2 )
    {
        bodyColor(obj);
        drawObject(obj.position(), disp->size, gle::tetrahedron);
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
void Display1::drawBeadT(Bead const& obj)
{
    const PointDisp * disp = obj.prop->disp;
    
    if ( disp->style & 1 )
    {
        drawBallT(obj.position(), obj.radius(), bodyColorF(obj));
    }
}

//------------------------------------------------------------------------------
#pragma mark -

void Display1::drawSolid(Solid const& obj)
{
    const PointDisp * disp = obj.prop->disp;
    
    //display points:
    if ( disp->style & 2  &&  disp->size > 0 )
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
            gle::transAlignZ(A, obj.radius(0), B-A);
            gle::cylinder1();
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
        drawStrip(obj.nbPoints(), obj.addrPoints(), GL_LINE_LOOP);
    }
}

/**
 Display a semi-transparent disc / sphere
 */
void Display1::drawSolidT(Solid const& obj, size_t inx)
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

void Display1::drawSphere(Sphere const& obj)
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


void Display1::drawSphereT(Sphere const& obj)
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
            drawFlat(C, obj.radius(), gle::disc);
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
void Display1::drawOrganizer(Organizer const& obj) const
{
    PointDisp const* disp = obj.disp();
    size_t i = 0, cnt = obj.nbLinks();

    if ( disp && ( disp->style & 2 ))
    {
        Vector P, Q;
        fluteD* pts = gle::mapVertexBuffer(2*cnt);
        while ( obj.getLink(i, P, Q) & ( i < cnt ) )
        {
            if ( modulo ) modulo->fold(Q, P);
            pts[  2*i] = P;
            pts[1+2*i] = Q;
            ++i;
        }
        gle::unmapVertexBuffer();
        glDisable(GL_LIGHTING);
        bodyColorF(disp, obj.signature()).load();
        lineWidth(disp->width);
        glDrawArrays(GL_LINES, 0, 2*i);

        gle::bindVertexBuffer((DIM>1?DIM:2));
        pointSize(disp->size);
        glDrawArrays(GL_POINTS, 0, i);
    }

    /**
     This display the Solid connecting two Aster as a spindle.
     Used for Cleo Kozlowski simulation of C. elegans (2007)
     */
    if ( disp && ( disp->style & 1 ) && obj.tag() == Fake::TAG )
    {
        Solid const* sol = Solid::toSolid(obj.core());
        if ( sol && sol->nbPoints() >= 4 )
        {
#if ( DIM == 3 )
            glEnable(GL_LIGHTING);
            bodyColor(*sol);
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
            drawStrip(sol->nbPoints(), sol->addrPoints(), GL_LINES);
#endif
        }
    }
}

//------------------------------------------------------------------------------
#pragma mark - Graphics primitives with shift

template < typename OBJ >
inline fluteD setVertex(Vector const& pos, const OBJ& obj)
{
#  if ( DIM == 1 )
    return flute2::cast(pos.XX, obj->signature()*0x1p-28-4);
#  else
    return fluteD{pos};
#  endif
}

#if ENABLE_EXPLODE_DISPLAY
inline fluteD setVertex(Vector const& pos, const Fiber * fib)
{
    GLfloat shift = fib->disp->explode_shift;
#  if ( DIM == 1 )
    return flute2::cast(pos.XX, shift);
#  elif ( DIM == 2 )
    return flute2::cast(pos.XX, pos.YY+shift);
#  else
    return flute3::cast(pos.XX, pos.YY+shift, pos.ZZ);
#  endif
}
#else
inline fluteD setVertex(Vector const& pos, const Fiber*)
{
    return fluteD{pos};
}
#endif

//------------------------------------------------------------------------------
#pragma mark -

void Display1::drawSinglesF(const SingleSet & set) const
{
    if ( prop->point_size > 0 )
    {
        size_t i = 0, cnt = set.sizeF();
        fluteD* pts = gle::mapVertexBuffer(cnt);
        flute4* col = gle::mapColorBuffer(cnt);
        for ( Single * obj=set.firstF(); obj ; obj=obj->next() )
        {
            if ( obj->disp()->perceptible )
            {
                pts[i] = setVertex(obj->posFoot(), obj);
                col[i++] = flute4{obj->disp()->color2};
            }
        }
        gle::unmapVertexBuffer();
        gle::unmapColorBuffer();
        pointSize(prop->point_size);
        glDrawArrays(GL_POINTS, 0, i);
    }
}


void Display1::drawSinglesA(const SingleSet & set) const
{
    size_t i = 0, cnt = 2*set.sizeF();
    if ( prop->point_size > 0 )
    {
        fluteD* pts = gle::mapVertexBuffer(cnt);
        flute4* col = gle::mapColorBuffer(cnt);
        for ( Single * obj=set.firstA(); obj ; obj=obj->next() )
        {
            Fiber const* fib = obj->fiber();
            if ( obj->disp()->perceptible & fib->disp->visible )
            {
                pts[i] = setVertex(obj->posHand(), fib);
                col[i++] = flute4{obj->disp()->color2};
            }
        }
        gle::unmapVertexBuffer();
        gle::unmapColorBuffer();
        pointSize(prop->point_size);
        glDrawArrays(GL_POINTS, 0, i);
    }
    
    // display links to anchor points
    if ( prop->link_width > 0 )
    {
        fluteD* pts = gle::mapVertexBuffer(cnt);
        flute4* col = gle::mapColorBuffer(cnt);
        i = 0;
        for ( Single * obj=set.firstA(); obj ; obj=obj->next() )
        {
            if ( obj->hasForce() )
            {
                Fiber const* fib = obj->fiber();
                if ( obj->disp()->perceptible & fib->disp->visible )
                {
                    Vector P = obj->posHand();
                    Vector Q = obj->posFoot();
                    if ( modulo ) modulo->fold(P, Q);
                    pts[i] = setVertex(P, fib);
                    col[i++] = flute4{obj->disp()->color};
                    pts[i] = setVertex(Q, obj);
                    col[i++] = flute4{obj->disp()->color};
                }
            }
        }
        gle::unmapVertexBuffer();
        gle::unmapColorBuffer();
        lineWidth(prop->link_width);
        glDrawArrays(GL_LINES, 0, i);
    }
}

//------------------------------------------------------------------------------
#pragma mark -
/**
Always display Hand1 of Couple
 */
void Display1::drawCouplesF1(CoupleSet const& set) const
{
    if ( prop->point_size > 0 )
    {
        size_t i = 0, cnt = set.sizeFF();
        fluteD* pts = gle::mapVertexBuffer(cnt);
        flute4* col = gle::mapColorBuffer(cnt);
        for ( Couple * obj = set.firstFF(); obj ; obj=obj->next() )
        {
            if ( obj->active() & obj->disp1()->perceptible )
            {
                pts[i] = setVertex(obj->posFree(), obj);
                col[i++] = flute4{obj->disp1()->color2};
            }
        }
        gle::unmapVertexBuffer();
        gle::unmapColorBuffer();
        pointSize(prop->point_size);
        glDrawArrays(GL_POINTS, 0, i);

#if ( DIM > 1 )
        // display inactive Couples with smaller size:
        i = 0;
        pts = gle::mapVertexBuffer(cnt);
        col = gle::mapColorBuffer(cnt);
        for ( Couple * obj = set.firstFF(); obj ; obj=obj->next() )
        {
            if ( !obj->active() && obj->disp1()->perceptible )
            {
                pts[i] = fluteD{obj->posFree()};
                col[i++] = flute4{obj->disp1()->color2};
            }
        }
        gle::unmapVertexBuffer();
        gle::unmapColorBuffer();
        pointSize(M_SQRT1_2*prop->point_size);
        glDrawArrays(GL_POINTS, 0, i);
#endif
    }
}


/**
 Display either Hand1 or Hand2, exposing both sides with equal chances.
 This gives the impression that Couple flicker randomly between frames,
 as if they were two-sided balls 'rotating' very fast.
 */
void Display1::drawCouplesF2(CoupleSet const& set) const
{
    if ( prop->point_size > 0 )
    {
        size_t i = 0, cnt = set.sizeFF();
        fluteD* pts = gle::mapVertexBuffer(cnt);
        flute4* col = gle::mapColorBuffer(cnt);
        Couple * nxt;
        Couple * obj = set.firstFF();

        if ( set.sizeFF() & 1 )
        {
            nxt = obj->next();
            if ( obj->disp12()->perceptible )
            {
                pts[i] = setVertex(obj->posFree(), obj);
                col[i++] = flute4{obj->disp12()->color2};
            }
            obj = nxt;
        }
        while ( obj )
        {
            nxt = obj->next();
            if ( obj->disp21()->perceptible )
            {
                pts[i] = setVertex(obj->posFree(), obj);
                col[i++] = flute4{obj->disp21()->color2};
            }
            obj = nxt->next();
            if ( nxt->disp12()->perceptible )
            {
                pts[i] = setVertex(nxt->posFree(), nxt);
                col[i++] = flute4{nxt->disp12()->color2};
            }
        }
        gle::unmapVertexBuffer();
        gle::unmapColorBuffer();
        pointSize(prop->point_size);
        glDrawArrays(GL_POINTS, 0, i);
    }
}


void Display1::drawCouplesA(CoupleSet const& set) const
{
    if ( prop->point_size > 0 )
    {
        size_t i = 0, cnt = set.sizeA();
        fluteD* pts = gle::mapVertexBuffer(cnt);
        flute4* col = gle::mapColorBuffer(cnt);
        for ( Couple * obj=set.firstAF(); obj ; obj=obj->next() )
        {
            Fiber const* fib = obj->fiber1();
            if ( obj->disp1()->perceptible & fib->disp->visible )
            {
                pts[i] = setVertex(obj->posHand1(), fib);
                col[i++] = flute4{obj->disp1()->color2};
            }
        }
        for ( Couple * obj=set.firstFA(); obj ; obj=obj->next() )
        {
            Fiber const* fib = obj->fiber2();
            if ( obj->disp2()->perceptible & fib->disp->visible )
            {
                pts[i] = setVertex(obj->posHand2(), fib);
                col[i++] = flute4{obj->disp2()->color2};
            }
        }
        gle::unmapVertexBuffer();
        gle::unmapColorBuffer();
        pointSize(prop->point_size);
        glDrawArrays(GL_POINTS, 0, i);
    }
}


void Display1::drawCouplesB(CoupleSet const& set) const
{
    size_t i = 0, cnt = 2 * set.sizeAA();

    if ( prop->point_size > 0 )
    {
        fluteD* pts = gle::mapVertexBuffer(cnt);
        flute4* col = gle::mapColorBuffer(cnt);
        for ( Couple * obj=set.firstAA(); obj ; obj=obj->next() )
        {
#if ( 0 )
            // only display if bridging two anti-parallel filaments
            if ( prop->couple_select & 8  && obj->cosAngle() > 0 )
                continue;
            // only display if bridging two parallel filaments
            if ( prop->couple_select & 16 && obj->cosAngle() < 0 )
                continue;
#endif
            Fiber const* fib1 = obj->fiber1();
            Fiber const* fib2 = obj->fiber2();
            if ( obj->disp1()->perceptible & fib1->disp->visible )
            {
                pts[i] = setVertex(obj->posHand1(), fib1);
                col[i++] = flute4{obj->disp1()->color};
            }
            if ( obj->disp2()->perceptible & fib2->disp->visible )
            {
                pts[i] = setVertex(obj->posHand2(), fib2);
                col[i++] = flute4{obj->disp2()->color};
            }
        }
        gle::unmapVertexBuffer();
        gle::unmapColorBuffer();
        pointSize(prop->point_size);
        glDrawArrays(GL_POINTS, 0, i);
    }
    
    // display the link for bridging couples
    if ( prop->link_width > 0 )
    {
#if ENABLE_EXPLODE_DISPLAY
        cnt *= 2;
#endif
        fluteD* pts = gle::mapVertexBuffer(cnt);
        flute4* col = gle::mapColorBuffer(cnt);
        i = 0;
        for ( Couple * obj=set.firstAA(); obj ; obj=obj->next() )
        {
#if ( 0 )
            // only display if bridging two anti-parallel filaments
            if ( prop->couple_select & 8  && obj->cosAngle() > 0 )
                continue;
            // only display if bridging two parallel filaments
            if ( prop->couple_select & 16 && obj->cosAngle() < 0 )
                continue;
#endif
            Fiber const* fib1 = obj->fiber1();
            Fiber const* fib2 = obj->fiber2();
            bool vis1 = obj->disp1()->perceptible & fib1->disp->visible;
            bool vis2 = obj->disp2()->perceptible & fib2->disp->visible;
            if ( vis1 | vis2 )
            {
                Vector P = obj->posHand1();
                Vector Q = obj->posHand2();
                if ( modulo ) modulo->fold(Q, P);
#if ENABLE_EXPLODE_DISPLAY
                // display the link twice, attached to each (shifted) fiber
                if ( vis1 )
                {
                    pts[i] = setVertex(P, fib1);
                    col[i++] = flute4{obj->disp1()->color};
                    pts[i] = setVertex(Q, fib1);
                    col[i++] = flute4{obj->disp2()->color};
                }
                if ( vis2 )
                {
                    pts[i] = setVertex(P, fib2);
                    col[i++] = flute4{obj->disp1()->color};
                    pts[i] = setVertex(Q, fib2);
                    col[i++] = flute4{obj->disp2()->color};
                }
#else
                pts[i] = setVertex(P, fib1);
                col[i++] = flute4{obj->disp1()->color};
                pts[i] = setVertex(Q, fib2);
                col[i++] = flute4{obj->disp2()->color};
#endif
            }
        }
        gle::unmapVertexBuffer();
        gle::unmapColorBuffer();
        lineWidth(prop->link_width);
        glDrawArrays(GL_LINES, 0, i);
    }
}

