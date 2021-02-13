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

    if ( prop->couple_select & 1 )
        drawCouplesF(sim.couples);

    if ( prop->single_select & 1 )
        drawSinglesF(sim.singles);
    
    glDisableClientState(GL_COLOR_ARRAY);

    drawFibers(sim.fibers);

    glEnableClientState(GL_COLOR_ARRAY);
    if ( prop->couple_select & 2 )
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
    
    if ( !disp )
        return;

    if ( disp->style & 2 )
    {
        Vector P, Q;
        glDisable(GL_LIGHTING);

        bodyColorF(disp, obj.signature()).load();
        pointSize(0.75f*disp->size);
        glBegin(GL_POINTS);
        for ( size_t i = 0; obj.getLink(i, P, Q); ++i )
            gleVertex(P);
        glEnd();

        lineWidth(disp->width);
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
            drawStrip(sol->nbPoints(), sol->addrPoints(), GL_LINES);
#endif
        }
    }
}

//------------------------------------------------------------------------------
#pragma mark -

template < typename OBJ >
inline floatD setVertex(Vector const& pos, const OBJ& obj)
{
#  if ( DIM == 1 )
    return float2::cast(pos.XX, obj->signature()*0x1p-28-4);
#  else
    return floatD{pos};
#  endif
}

#if ENABLE_EXPLODE_DISPLAY
inline floatD setVertex(Vector const& pos, const Fiber * fib)
{
    GLfloat shift = fib->disp->explode_shift;
#  if ( DIM == 1 )
    return float2::cast(pos.XX, shift);
#  elif ( DIM == 2 )
    return float2::cast(pos.XX, pos.YY+shift);
#  else
    return float3::cast(pos.XX, pos.YY+shift, pos.ZZ);
#  endif
}
#else
inline floatD setVertex(Vector const& pos, const Fiber*)
{
    return floatD{pos};
}
#endif

#if ENABLE_EXPLODE_DISPLAY

inline void shiftedVertex(Vector const& pos, const Fiber * fib)
{
    real shift = fib->disp->explode_shift;
#if ( DIM == 3 )
    glVertex3f(pos.XX, pos.YY+shift, pos.ZZ);
#elif ( DIM == 2 )
    glVertex2f(pos.XX, pos.YY+shift);
#else
    glVertex2f(pos.XX, shift);
#endif
}

inline void drawLine(Vector const& a, Fiber const* fib, PointDisp const* disp, Vector const& b)
{
    if ( disp->visible && fib->disp->visible )
    {
        disp->color.load();
        shiftedVertex(a, fib);
        disp->color2.load();
        shiftedVertex(b, fib);
    }
}


/**
 */
inline void drawLine(Vector const& a, const Fiber * fibA, const PointDisp* dispA,
                     Vector const& b, const Fiber * fibB, const PointDisp* dispB)
{
#if ( 1 )
    //draw two segments if `explode` is enabled
    if ( dispA->visible && fibA->disp->visible )
    {
        dispA->color.load();
        shiftedVertex(a, fibA);
        dispB->color.load();
        shiftedVertex(b, fibA);
    }
    if ( dispB->visible && fibB->disp->visible && fibB->prop->disp->explode_style )
    {
        dispA->color.load();
        shiftedVertex(a, fibB);
        dispB->color.load();
        shiftedVertex(b, fibB);
    }
#else
    if ( fibA->disp->visible || fibB->disp->visible )
    {
        //draw one segment
        dispA->color.load();
        shiftedVertex(a, fibA);
        dispB->color.load();
        shiftedVertex(b, fibB);
    }
#endif
}

#else

inline void drawLine(Vector const& a, const Fiber * fib, const PointDisp* disp, Vector const& b)
{
    if ( disp->perceptible && fib->disp->visible )
    {
        disp->color.load();
        gleVertex(a);
        disp->color2.load();
        gleVertex(b);
    }
}

inline void drawLine(Vector const& a, const Fiber * fibA, const PointDisp* dispA,
                     Vector const& b, const Fiber * fibB, const PointDisp* dispB)
{
    if (   dispA->perceptible && dispB->perceptible
        && ( fibA->disp->visible || fibB->disp->visible ))
    {
        dispA->color.load();
        gleVertex(a);
        dispB->color.load();
        gleVertex(b);
    }
}

#endif

//------------------------------------------------------------------------------
#pragma mark -

void Display1::drawSinglesF(const SingleSet & set) const
{
    if (( prop->point_size > 0 ) & ( set.sizeF() > 0 ))
    {
        size_t i = 0, cnt = set.sizeF();
        floatD* pts = gle::mapVertexBuffer(cnt);
        float4* col = gle::mapColorBuffer(cnt);
        for ( Single * obj=set.firstF(); obj ; obj=obj->next() )
        {
            if ( obj->disp()->perceptible )
            {
                pts[i] = setVertex(obj->posFoot(), obj);
                col[i++] = float4{obj->disp()->color2};
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
    // display positions of Hands
    if (( prop->point_size > 0 ) & ( set.sizeF() > 0 ))
    {
        glEnableClientState(GL_COLOR_ARRAY);
        size_t i = 0, cnt = set.sizeF();
        floatD* pts = gle::mapVertexBuffer(cnt);
        float4* col = gle::mapColorBuffer(cnt);
        for ( Single * obj=set.firstA(); obj ; obj=obj->next() )
        {
            Fiber const* fib = obj->fiber();
            if ( obj->disp()->perceptible & fib->disp->visible )
            {
                pts[i] = setVertex(obj->posHand(), fib);
                col[i++] = float4{obj->disp()->color2};
            }
        }
        gle::unmapVertexBuffer();
        gle::unmapColorBuffer();
        pointSize(prop->point_size);
        glDrawArrays(GL_POINTS, 0, i);
        glDisableClientState(GL_COLOR_ARRAY);
    }
    
    // display links to anchor points
    if ( prop->link_width > 0 )
    {
        lineWidth(prop->link_width);
        glBegin(GL_LINES);
        for ( Single * obj=set.firstA(); obj ; obj=obj->next() )
        {
            if ( obj->hasForce() )
            {
                Vector ph = obj->posHand();
                Vector pf = obj->posFoot();
                if ( modulo ) modulo->fold(pf, ph);
                drawLine(ph, obj->fiber(), obj->disp(), pf);
            }
        }
        glEnd();
    }
}

//------------------------------------------------------------------------------
#pragma mark -
/**
Always display Hand1 of Couple
 */
void Display1::drawCouplesF1(CoupleSet const& set) const
{
    if (( prop->point_size > 0 ) & ( set.sizeFF() > 0 ))
    {
        size_t i = 0, cnt = set.sizeFF();
        floatD* pts = gle::mapVertexBuffer(cnt);
        float4* col = gle::mapColorBuffer(cnt);
        for ( Couple * obj = set.firstFF(); obj ; obj=obj->next() )
        {
            if ( obj->active() & obj->disp1()->perceptible )
            {
                pts[i] = setVertex(obj->posFree(), obj);
                col[i++] = float4{obj->disp1()->color2};
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
                pts[i] = floatD{obj->posFree()};
                col[i++] = float4{obj->disp1()->color2};
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
    if (( prop->point_size > 0 ) & ( set.sizeFF() > 0 ))
    {
        size_t i = 0, cnt = set.sizeFF();
        floatD* pts = gle::mapVertexBuffer(cnt);
        float4* col = gle::mapColorBuffer(cnt);
        Couple * nxt;
        Couple * obj = set.firstFF();

        if ( set.sizeFF() & 1 )
        {
            nxt = obj->next();
            if ( obj->disp12()->perceptible )
            {
                pts[i] = setVertex(obj->posFree(), obj);
                col[i++] = float4{obj->disp12()->color2};
            }
            obj = nxt;
        }
        while ( obj )
        {
            nxt = obj->next();
            if ( obj->disp21()->perceptible )
            {
                pts[i] = setVertex(obj->posFree(), obj);
                col[i++] = float4{obj->disp21()->color2};
            }
            obj = nxt->next();
            if ( nxt->disp12()->perceptible )
            {
                pts[i] = setVertex(nxt->posFree(), nxt);
                col[i++] = float4{nxt->disp12()->color2};
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
    if (( prop->point_size > 0 ) & ( set.sizeAF()+set.sizeFA() > 0 ))
    {
        size_t i = 0, cnt = set.sizeAF()+set.sizeFA();
        floatD* pts = gle::mapVertexBuffer(cnt);
        float4* col = gle::mapColorBuffer(cnt);
        for ( Couple * obj=set.firstAF(); obj ; obj=obj->next() )
        {
            Fiber const* fib = obj->fiber1();
            if ( obj->disp1()->perceptible & fib->disp->visible )
            {
                pts[i] = setVertex(obj->posHand1(), fib);
                col[i++] = float4{obj->disp1()->color2};
            }
        }
        for ( Couple * obj=set.firstFA(); obj ; obj=obj->next() )
        {
            Fiber const* fib = obj->fiber2();
            if ( obj->disp2()->perceptible & fib->disp->visible )
            {
                pts[i] = setVertex(obj->posHand2(), fib);
                col[i++] = float4{obj->disp2()->color2};
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
    if (( prop->point_size > 0 ) & ( set.sizeAA() > 0 ) )
    {
        glEnableClientState(GL_COLOR_ARRAY);
        size_t i = 0, cnt = 2 * set.sizeAA();
        floatD* pts = gle::mapVertexBuffer(cnt);
        float4* col = gle::mapColorBuffer(cnt);
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
                col[i++] = float4{obj->disp1()->color};
            }
            if ( obj->disp2()->perceptible & fib2->disp->visible )
            {
                pts[i] = setVertex(obj->posHand2(), fib2);
                col[i++] = float4{obj->disp2()->color};
            }
        }
        gle::unmapVertexBuffer();
        gle::unmapColorBuffer();
        pointSize(prop->point_size);
        glDrawArrays(GL_POINTS, 0, i);
        glDisableClientState(GL_COLOR_ARRAY);
    }
    
    // display the link for bridging couples
    if ( prop->link_width > 0 )
    {
        lineWidth(prop->link_width);
        glBegin(GL_LINES);
        for ( Couple * cx=set.firstAA(); cx ; cx=cx->next() )
        {
#if ( 0 )
            // only display if bridging two anti-parallel filaments
            if ( prop->couple_select & 8  && cx->cosAngle() > 0 )
                continue;
            // only display if bridging two parallel filaments
            if ( prop->couple_select & 16 && cx->cosAngle() < 0 )
                continue;
#endif
            Vector P = cx->posHand1();
            Vector Q = cx->posHand2();
            if ( modulo ) modulo->fold(Q, P);
            drawLine(P, cx->fiber1(), cx->disp1(), Q, cx->fiber2(), cx->disp2());
        }
        glEnd();
    }
}

