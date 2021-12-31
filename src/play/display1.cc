// Cytosim was created by Francois Nedelec. Copyright Cambridge University 2020

#include "dim.h"
#include "simul.h"
#include "display1.h"
#include "modulo.h"

#include "gle.h"
#include "gle_color_list.h"
#include "glut.h"

#include "fiber_disp.h"
#include "line_disp.h"
#include "point_disp.h"
//#include "display_color.h"
#include "gym_flute.h"

extern Modulo const* modulo;


#define ENABLE_EXPLODED_DISPLAY ( DIM < 2 )


Display1::Display1(DisplayProp const* dp) : Display(dp)
{
}


void Display1::drawObjects(Simul const& sim)
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

#if ( DIM >= 3 )
    
    glEnable(GL_LIGHTING);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);

#endif
    
    drawBeads(sim.beads);
    drawSolids(sim.solids);
    drawSpheres(sim.spheres);
    
#if ( DIM >= 3 )
    
    glDisable(GL_LIGHTING);
    glDisable(GL_CULL_FACE);

#endif
    
    glEnableClientState(GL_COLOR_ARRAY);
    if (( prop->couple_select & 4 ) && ( sim.couples.sizeAA() > 0 ))
        drawCouplesB(sim.couples);
    
    if (( prop->single_select & 2 ) && ( sim.singles.sizeA() > 0 ))
        drawSinglesA(sim.singles);
    glDisableClientState(GL_COLOR_ARRAY);

#if ( DIM >= 3 )
    
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
#if ENABLE_EXPLODED_DISPLAY
    //translate whole display to display the Fiber
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    gle::translate(0, fib.disp->explode_shift, 0);
#endif
    
    Display::drawFiber(fib);

#if ENABLE_EXPLODED_DISPLAY
    glPopMatrix();
#endif
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

#if ENABLE_EXPLODED_DISPLAY
inline fluteD setVertex(Vector const& pos, const Fiber * fib)
{
    real shift = fib->disp->explode_shift;
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
        fluteD4* flu = gym::mapBufferC4VD(cnt);
        for ( Single * obj=set.firstF(); obj ; obj=obj->next() )
        {
            if ( obj->disp()->perceptible )
                flu[i++] = {obj->disp()->color2, setVertex(obj->posFoot(), obj)};
        }
        gym::unmapBufferC4VD();
        pointSize(prop->point_size);
        glDrawArrays(GL_POINTS, 0, i);
    }
}


void Display1::drawSinglesA(const SingleSet & set) const
{
    gle_color air(0,0,0,0);
    size_t i = 0, cnt = 2*set.sizeA();
    fluteD4* flu = gym::mapBufferC4VD(cnt);
    for ( Single * obj=set.firstA(); obj ; obj=obj->next() )
    {
        Fiber const* fib = obj->fiber();
        if ( obj->disp()->perceptible & fib->disp->visible )
        {
            gle_color d, c = obj->disp()->color;
            Vector Q, P = obj->posHand();
            if ( obj->hasLink() ) {
                d = c;
                Q = obj->posFoot();
                if ( modulo ) modulo->fold(Q, P);
            } else {
                Q = P;
                d = air;
            }
            flu[i++] = {c, setVertex(P, fib)};
            flu[i++] = {d, setVertex(Q, obj)};
        }
    }
    gym::unmapBufferC4VD();
    
    if ( prop->link_width > 0 )
    {
        lineWidth(prop->link_width);
        glDrawArrays(GL_LINES, 0, i);
    }
    
    if ( prop->point_size > 0 )
    {
        gym::bindBufferC4VD(2);
        pointSize(prop->point_size);
        glDrawArrays(GL_POINTS, 0, i/2);
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
        fluteD4* flu = gym::mapBufferC4VD(cnt);
        for ( Couple * obj = set.firstFF(); obj ; obj=obj->next() )
        {
            if ( obj->active() && obj->disp1()->perceptible )
                flu[i++] = {obj->disp1()->color2, setVertex(obj->posFree(), obj)};
        }
        gym::unmapBufferC4VD();
        pointSize(prop->point_size);
        glDrawArrays(GL_POINTS, 0, i);

#if ( DIM > 1 )
        // display inactive Couples with square dots:
        i = 0;
        flu = gym::mapBufferC4VD(cnt);
        for ( Couple * obj = set.firstFF(); obj ; obj=obj->next() )
        {
            if ( !obj->active() && obj->disp1()->perceptible )
                flu[i++] = {obj->disp1()->color2, obj->posFree()};
        }
        gym::unmapBufferC4VD();
        glDisable(GL_POINT_SMOOTH);
        pointSize(M_SQRT1_2*prop->point_size);
        glDrawArrays(GL_POINTS, 0, i);
        glEnable(GL_POINT_SMOOTH);
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
        fluteD4* flu = gym::mapBufferC4VD(cnt);
        Couple * nxt;
        Couple * obj = set.firstFF();

        if ( set.sizeFF() & 1 )
        {
            nxt = obj->next();
            if ( obj->disp12()->perceptible )
                flu[i++] = { obj->disp12()->color2, setVertex(obj->posFree(), obj) };
            obj = nxt;
        }
        while ( obj )
        {
            nxt = obj->next();
            if ( obj->disp21()->perceptible )
                flu[i++] = { obj->disp21()->color2, setVertex(obj->posFree(), obj) };
            obj = nxt->next();
            if ( nxt->disp12()->perceptible )
                flu[i++] = { nxt->disp12()->color2, setVertex(nxt->posFree(), nxt) };
        }
        gym::unmapBufferC4VD();
        pointSize(prop->point_size);
        glDrawArrays(GL_POINTS, 0, i);
    }
}


void Display1::drawCouplesA(CoupleSet const& set) const
{
    if ( prop->point_size > 0 )
    {
        size_t i = 0, cnt = set.sizeA();
        fluteD4* flu = gym::mapBufferC4VD(cnt);
        for ( Couple * obj=set.firstAF(); obj ; obj=obj->next() )
        {
            Fiber const* fib = obj->fiber1();
            if ( obj->disp1()->perceptible & fib->disp->visible )
                flu[i++] = { obj->disp1()->color2, setVertex(obj->posHand1(), fib) };
        }
        for ( Couple * obj=set.firstFA(); obj ; obj=obj->next() )
        {
            Fiber const* fib = obj->fiber2();
            if ( obj->disp2()->perceptible & fib->disp->visible )
                flu[i++] = { obj->disp2()->color2, setVertex(obj->posHand2(), fib) };
        }
        gym::unmapBufferC4VD();
        pointSize(prop->point_size);
        glDrawArrays(GL_POINTS, 0, i);
    }
}


void Display1::drawCouplesB(CoupleSet const& set) const
{
    gle_color air(0,0,0,0);
    size_t i = 0, cnt = 2 * set.sizeAA() * (1+ENABLE_EXPLODED_DISPLAY);
    fluteD4* flu = gym::mapBufferC4VD(cnt);
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
        gle_color col1 = vis1 ? obj->disp1()->color : air;
        gle_color col2 = vis2 ? obj->disp2()->color : air;
        Vector P = obj->posHand1();
        Vector Q = obj->posHand2();
        if ( modulo ) modulo->fold(Q, P);
        if ( vis1 | vis2 )
        {
            flu[i++] = { col1, setVertex(P, fib1) };
            flu[i++] = { col2, setVertex(Q, fib1) };
        }
#if ENABLE_EXPLODED_DISPLAY
        if ( vis2 )
        {
            flu[i++] = { col2, setVertex(Q, fib2) };
            flu[i++] = { col1, setVertex(P, fib2) };
        }
#endif
    }
    gym::unmapBufferC4VD();
    
    if ( prop->link_width > 0 )
    {
        lineWidth(prop->link_width);
        glDrawArrays(GL_LINES, 0, i);
    }

    if ( prop->point_size > 0 )
    {
        pointSize(prop->point_size);
#if ENABLE_EXPLODED_DISPLAY
        gym::bindBufferC4VD(2);
        glDrawArrays(GL_POINTS, 0, i/2);
#else
        glDrawArrays(GL_POINTS, 0, i);
#endif
    }
}


