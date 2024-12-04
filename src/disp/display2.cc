// Cytosim was created by Francois Nedelec. Copyright Cambridge University 2020

#include "dim.h"
#include "simul.h"
#include "display2.h"
#include "modulo.h"

#include "gle.h"
#include "gym_color_list.h"

#include "line_disp.h"
#include "point_disp.h"
#include "fiber_disp.h"
#include "gym_flute.h"
#include "gym_flute_dim.h"
#include "gym_draw.h"
#include "gym_check.h"

//------------------------------------------------------------------------------

Display2::Display2(DisplayProp const* dp) : Display(dp)
{
}


void Display2::drawObjects(Simul const& sim)
{
    gym::closeDepthMask();
    gym::disableLighting();
    gym::disableCullFace();
    drawFields(sim.fields);
    
    gym::enableLighting();
#if ( DIM > 2 )
    gym::openDepthMask();
#endif
    drawSpaces(sim.spaces);
    gym::disableLighting();
    gym::disableCullFace();

    if (( prop->couple_select & 1 ) && ( sim.couples.sizeFF() > 0 ))
        drawCouplesF(sim.couples);
    
    if (( prop->couple_select & 2 ) && ( sim.couples.sizeA() > 0 ))
        drawCouplesA(sim.couples);
    
    if (( prop->single_select & 1 ) && ( sim.singles.sizeF() > 0 ))
        drawSinglesF(sim.singles);
    
    glBindBuffer(GL_PIXEL_UNPACK_BUFFER, 0);
    
    drawFibers(sim.fibers);
    drawFiberTexts(sim.fibers);

#if ( DIM >= 3 )
    gym::enableLighting();
    gym::enableCullFace(GL_BACK);
#else
    gym::disableLighting();
#endif
    
    drawBeads(sim.beads);
    drawSolids(sim.solids);
    drawSpheres(sim.spheres);
    
#if ( DIM >= 3 )
    gym::disableCullFace();
#endif
    gym::disableLighting();

    if (( prop->couple_select & 4 ) && ( sim.couples.sizeAA() > 0 ))
        drawCouplesB(sim.couples);
    
    if (( prop->single_select & 2 ) && ( sim.singles.sizeA() > 0 ))
        drawSinglesA(sim.singles);
    glBindBuffer(GL_PIXEL_UNPACK_BUFFER, 0);

#if ( DIM >= 3 )
    gym::enableLighting();
    gym::enableCullFace(GL_BACK);
#endif

    drawOrganizers(sim.organizers);
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
    CHECK_GL_ERROR("in Display::drawSinglesF()");
}


void Display2::drawSinglesA(const SingleSet & set) const
{
    // display the Hands
    for ( Single * obj=set.firstA(); obj ; obj=obj->next() )
    {
        PointDisp const* dis = obj->disp();
        if ( dis->perceptible  &&  obj->fiber()->disp->visible )
        {
            Vector ph = obj->posHand();
            
            dis->drawA(ph);
            
            if ( obj->hasLink() && dis->width > 0 )
            {
                Vector ps = obj->sidePos();
                Vector pf = obj->posFoot();
                if ( modulo )
                {
                    modulo->fold(pf, ph);
                    modulo->fold(ps, ph);
                }
                
                gym::color(dis->color);
#if ( DIM >= 3 )
                gym::stretchAlignZ(pf, ph, pixscale(dis->width));
                gle::cutCone();
                //drawCone(pf, ph-pf, pixscale(dis->width));
#else
                gle::drawBand(ph, pixscale(dis->width), ps, pixscale(dis->width));
                gle::drawBand(ps, pixscale(dis->width), dis->color, pf, pixscale(dis->width), dis->color.alpha_scaled(0.5f));
#endif
            }
        }
    }
}

//------------------------------------------------------------------------------
#pragma mark -

void Display2::drawCouplesF(CoupleSet const& set) const
{
    if ( prop->couple_flip )
        drawCouplesF2(set);
    else
        drawCouplesF1(set);
}

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
    gym_color air(0,0,0,0);
    PointDisp const* pd1 = cx->disp1();
    PointDisp const* pd2 = cx->disp2();
    
    Vector p1 = cx->posHand1();
    Vector p2 = cx->posHand2();
    
    if ( modulo )
        modulo->fold(p2, p1);
    
    if ( pd1->perceptible || pd2->perceptible )
    {
        flute4D* flu = gym::mapBufferC4VD(4);
#if DIM < 2
        flu[0] = { pd1->color, p1 };
        flu[1] = { pd2->color, p2 };
#else
        /*
         Can shift positions towards the minus-end by couple's length
         to create an effect to highlight the configuration:
         ///// on antiparallel fibers
         >>>>> on parallel fibers
         */
        Vector d1 = cx->dirFiber1();
        Vector d2 = cx->dirFiber2();
        Vector pp = 0.5*(p1+p2) + (0.25*cx->prop->length)*(d1+d2);
        gym_color col1 = pd1->visible ? pd1->color : air;
        gym_color col2 = pd2->visible ? pd2->color : air;
        flu[0] = { col1, p1 };
        flu[1] = { col1, pp };
        flu[2] = { col2, pp };
        flu[3] = { col2, p2 };
#endif
        gym::unmapBufferC4VD();
        gym::ref_view();
        gym::drawLines(pd1->widthX, 0, 4);
        gym::cleanupCV();
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

