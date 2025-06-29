// Cytosim was created by Francois Nedelec. Copyright Cambridge University 2020

#include "dim.h"
#include "simul.h"
#include "display1.h"
#include "modulo.h"

#include "gle.h"
#include "gym_color_list.h"

#include "fiber_disp.h"
#include "line_disp.h"
#include "point_disp.h"
#include "gym_flute.h"
#include "gym_flute_dim.h"
#include "gym_draw.h"


#define ENABLE_EXPLODED_DISPLAY ( DIM < 2 )


Display1::Display1(DisplayProp const* dp) : Display(dp)
{
    linkWidth = 1;
    pointSize = 1;
    smallSize = 1;
}


void Display1::drawObjects(Simul const& sim)
{
    linkWidth = pixwidth(prop->link_width);
    pointSize = pixwidth(prop->point_size);
    smallSize = pixwidth(prop->point_size);
    
    gym::closeDepthMask();
    gym::disableLighting();
    gym::disableCullFace();
    drawFields(sim.fields);
    
#if ( DIM > 2 )
    gym::openDepthMask();
#endif
    gym::enableLighting();
    drawSpaces(sim.spaces);
    gym::disableLighting();
    gym::disableCullFace();

    if (( prop->couple_select & 1 ) && ( sim.couples.sizeFF() > 0 ))
        drawCouplesF(sim.couples);

    if (( prop->single_select & 1 ) && ( sim.singles.sizeF() > 0 ))
        drawSinglesF(sim.singles);

#if ( DIM >= 3 )
    gym::enableCullFace(GL_BACK);
    gym::enableLighting();
#endif
    
    drawBeads(sim.beads);
    drawSolids(sim.solids);
    drawSpheres(sim.spheres);
    
#if ( DIM >= 3 )
    gym::disableLighting();
    gym::disableCullFace();
#endif
    
    drawFibers(sim.fibers);
#if !ENABLE_EXPLODED_DISPLAY
    drawFiberTexts(sim.fibers);
#endif
    gym::disableLighting();

    if (( prop->couple_select & 2 ) && ( sim.couples.sizeA() > 0 ))
        drawCouplesA(sim.couples);

    if (( prop->couple_select & 4 ) && ( sim.couples.sizeAA() > 0 ))
        drawCouplesB1(sim.couples);
    
    if (( prop->single_select & 2 ) && ( sim.singles.sizeA() > 0 ))
        drawSinglesA(sim.singles);

#if ( DIM >= 3 )
    gym::enableLighting();
    gym::enableCullFace(GL_BACK);
#endif

    drawOrganizers(sim.organizers);
}

//------------------------------------------------------------------------------
#pragma mark -


void Display1::drawFibers(FiberSet const& set)
{
#if ENABLE_EXPLODED_DISPLAY
    //translate whole display to display the Fiber
    GLfloat ref[16];
    gym::get_view(ref);
    for ( Fiber const* fib = set.first(); fib ; fib=fib->next() )
    {
        if ( fib->disp->visible )
        {
            gym::set_view(ref, 0, explodeShift(fib), 0);
            drawFiber(*fib);
        }
    }
    gym::set_view(ref);
#else
    for ( Fiber const* fib = set.first(); fib ; fib=fib->next() )
    {
        if ( fib->disp->visible )
            drawFiber(*fib);
    }
#endif
}

//------------------------------------------------------------------------------
#pragma mark - Methods to shift vertices vertically in exploded display

#if ENABLE_EXPLODED_DISPLAY

template < typename OBJ >
inline float Display1::explodeShift(OBJ const* obj) const
{
    if ( prop->explode_style )
        return prop->explode_range * ( obj->signature() * 0x1p-32 - 0.5 );
    return 0;
}


template < typename FLOATS, typename OBJ >
inline void Display1::shiftVertex(FLOATS * ptr, OBJ const* obj) const
{
    float shift = explodeShift(obj);
#  if ( DIM == 1 )
    ptr->setY(shift);
#  else
    ptr->setY(ptr->xyz[5]+shift);
#  endif
}

template < typename FLOATS >
inline void Display1::shiftVertex(FLOATS * ptr, FLOATS * qrt, const Fiber* fib) const
{
    float shift = explodeShift(fib);
#  if ( DIM == 1 )
    ptr->setY(shift);
    qrt->setY(shift);
#  else
    ptr->setY(ptr->xyz[5]+shift);
    qrt->setY(ptr->xyz[5]+shift);
#  endif
}

#else

template < typename FLOATS, typename OBJ >
inline void Display1::shiftVertex(FLOATS *, OBJ const*) const
{
}

template < typename FLOATS >
inline void Display1::shiftVertex(FLOATS *, FLOATS *, Fiber const*) const
{
}

#endif


//------------------------------------------------------------------------------
#pragma mark -

void Display1::drawSinglesF(const SingleSet & set) const
{
    if ( prop->point_size > 0 )
    {
        size_t cnt = set.sizeF();
        flute4D* flu = gym::mapBufferC4VD(cnt);
        flute4D* ptr = flu;
        for ( Single * obj=set.firstF(); obj ; obj=obj->next() )
        {
            if ( obj->disp()->perceptible )
            {
                ptr[0] = {obj->disp()->color2, obj->posFoot()};
                shiftVertex(ptr, obj);
                ++ptr;
            }
        }
        assert_true( ptr <= flu+cnt );
        gym::unmapBufferC4VD();
        gym::ref_view();
        gym::drawPoints(smallSize, 0, ptr-flu);
        gym::cleanupCV();
    }
}


void Display1::drawSinglesA(const SingleSet & set) const
{
    gym_color air(0,0,0,0);
    size_t cnt = 2 * set.sizeA();
    flute4D* flu = gym::mapBufferC4VD(cnt);
    flute4D* ptr = flu;
    for ( Single * obj=set.firstA(); obj ; obj=obj->next() )
    {
        Fiber const* fib = obj->fiber();
        if ( obj->disp()->perceptible & fib->disp->visible )
        {
            gym_color d, c = obj->disp()->color;
            Vector Q, P = obj->posHand();
            if ( obj->hasLink() ) {
                d = c;
                Q = obj->posFoot();
                if ( modulo ) modulo->fold(Q, P);
            } else {
                Q = P;
                d = air;
            }
            ptr[0] = {c, P};
            ptr[1] = {d, Q};
            shiftVertex(ptr, ptr+1, fib);
            ptr += 2;
        }
    }
    assert_true( ptr <= flu+cnt );
    gym::unmapBufferC4VD();
    gym::ref_view();

    if ( prop->link_width > 0 )
    {
        gym::drawLines(linkWidth, 0, ptr-flu);
    }
    
    if ( prop->point_size > 0 )
    {
        gym::rebindBufferC4VD(2);
        gym::drawPoints(pointSize, 0, (ptr-flu)/2);
    }
    gym::cleanupCV();
}

//------------------------------------------------------------------------------
#pragma mark - Couples


void Display1::drawCouplesF(CoupleSet const& set) const
{
    if ( prop->point_size > 0 )
    {
        if ( prop->couple_flip )
            drawCouplesF2(set);
        else
            drawCouplesF1(set);
    }
}


/**
Always display Hand1 of Couple
 */
void Display1::drawCouplesF1(CoupleSet const& set) const
{
    size_t cnt = set.sizeFF();
    flute4D* flu = gym::mapBufferC4VD(cnt);
    flute4D* ptr = flu;
    flute4D* end = flu + cnt;
    for ( Couple * obj = set.firstFF(); obj ; obj=obj->next() )
    {
        if ( obj->disp1()->perceptible )
        {
            flute4D * f4d = ( obj->active() ? ptr++ : --end );
            *f4d = {obj->disp1()->color2, obj->posFree()};
            shiftVertex(f4d, obj);
        }
    }
    assert_true( ptr <= flu+cnt );
    gym::unmapBufferC4VD();
    gym::ref_view();
    gym::drawPoints(smallSize, 0, ptr-flu);
    // display inactive Couples with square dots:
    gym::drawSquarePoints(0.5f*smallSize, end-flu, cnt-(end-flu));
    gym::cleanupCV();
}


/**
 Display either Hand1 or Hand2, exposing both sides with equal chances.
 This gives the impression that Couple flicker randomly between frames,
 as if they were two-sided balls 'rotating' very fast.
 */
void Display1::drawCouplesF2(CoupleSet const& set) const
{
    size_t cnt = set.sizeFF();
    flute4D* flu = gym::mapBufferC4VD(cnt);
    flute4D* ptr = flu;
    
    Couple * nxt;
    Couple * obj = set.firstFF();
    if ( set.sizeFF() & 1 )
    {
        nxt = obj->next();
        if ( obj->disp12()->perceptible )
        {
            ptr[0] = { obj->disp12()->color2, obj->posFree() };
            shiftVertex(ptr, obj);
            ++ptr;
        }
        obj = nxt;
    }
    while ( obj )
    {
        nxt = obj->next();
        if ( obj->disp21()->perceptible )
        {
            ptr[0] = { obj->disp21()->color2, obj->posFree() };
            shiftVertex(ptr, obj);
            ++ptr;
        }
        obj = nxt->next();
        if ( nxt->disp12()->perceptible )
        {
            ptr[0] = { nxt->disp12()->color2, nxt->posFree() };
            shiftVertex(ptr, nxt);
            ++ptr;
        }
    }
    assert_true( ptr <= flu+cnt );
    gym::unmapBufferC4VD();
    gym::ref_view();
    gym::drawPoints(smallSize, 0, ptr-flu);
    gym::cleanupCV();
}


void Display1::drawCouplesA(CoupleSet const& set) const
{
    if ( prop->point_size > 0 )
    {
        size_t cnt = set.sizeA();
        flute4D* flu = gym::mapBufferC4VD(cnt);
        flute4D* ptr = flu;
        for ( Couple * obj=set.firstAF(); obj ; obj=obj->next() )
        {
            Fiber const* fib = obj->fiber1();
            if ( obj->disp1()->perceptible & fib->disp->visible )
            {
                ptr[0] = { obj->disp1()->color2, obj->posHand1() };
                shiftVertex(ptr, fib);
                ++ptr;
            }
        }
        for ( Couple * obj=set.firstFA(); obj ; obj=obj->next() )
        {
            Fiber const* fib = obj->fiber2();
            if ( obj->disp2()->perceptible & fib->disp->visible )
            {
                ptr[0] = { obj->disp2()->color2, obj->posHand2() };
                shiftVertex(ptr, fib);
                ++ptr;
            }
        }
        assert_true( ptr <= flu+cnt );
        gym::unmapBufferC4VD();
        gym::ref_view();
        gym::drawPoints(pointSize, 0, ptr-flu);
        gym::cleanupCV();
    }
}


void Display1::drawCouplesB1(CoupleSet const& set) const
{
    gym_color air(0,0,0,0);
    size_t cnt = 2 * set.sizeAA() * (1+ENABLE_EXPLODED_DISPLAY);
    flute4D* flu = gym::mapBufferC4VD(cnt);
    flute4D* ptr = flu;
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
        gym_color col1 = vis1 ? obj->disp1()->color : air;
        gym_color col2 = vis2 ? obj->disp2()->color : air;
        Vector P = obj->posHand1();
        Vector Q = obj->posHand2();
        if ( modulo ) modulo->fold(Q, P);
        if ( vis1 | vis2 )
        {
            ptr[0] = { col1, P };
            ptr[1] = { col2, Q };
            shiftVertex(ptr, ptr+1, fib1);
            ptr += 2;
        }
#if ENABLE_EXPLODED_DISPLAY
        if ( vis2 )
        {
            ptr[0] = { col2, Q };
            ptr[1] = { col1, P };
            shiftVertex(ptr, ptr+1, fib2);
            ptr += 2;
        }
#endif
    }
    assert_true( ptr <= flu+cnt );
    gym::unmapBufferC4VD();
    gym::ref_view();

    if ( prop->link_width > 0 )
    {
        gym::drawLines(linkWidth, 0, ptr-flu);
    }

    if ( prop->point_size > 0 )
    {
#if ENABLE_EXPLODED_DISPLAY
        gym::rebindBufferC4VD(2);
        gym::drawPoints(pointSize, 0, (ptr-flu)/2);
#else
        gym::drawPoints(pointSize, 0, ptr-flu);
#endif
    }
    gym::cleanupCV();
}


void Display1::drawCouplesB0(CoupleSet const& set) const
{
    size_t cnt = 2 * set.sizeAA();
    flute4D* flu = gym::mapBufferC4VD(cnt);
    flute4D* ptr = flu;
    for ( Couple * obj=set.firstAA(); obj ; obj=obj->next() )
    {
        Vector P = obj->posHand1();
        Vector Q = obj->posHand2();
        if ( modulo ) modulo->fold(Q, P);
        ptr[0] = { obj->disp1()->color, P };
        ptr[1] = { obj->disp2()->color, Q };
        ptr += 2;
    }
    assert_true( ptr <= flu+cnt );
    gym::unmapBufferC4VD();
    gym::ref_view();
    if ( prop->link_width > 0 )
        gym::drawLines(linkWidth, 0, ptr-flu);
    if ( prop->point_size > 0 )
        gym::drawPoints(pointSize, 0, ptr-flu);
    gym::cleanupCV();
}

