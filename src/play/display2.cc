// Cytosim was created by Francois Nedelec. Copyright Cambridge University 2020

#include "dim.h"
#include "sim.h"
#include "simul.h"
#include "display2.h"
#include "modulo.h"

#include "opengl.h"
#include "gle.h"
#include "gle_color_list.h"
#include "glut.h"

#include "line_disp.h"
#include "point_disp.h"
#include "fiber_disp.h"
//#include "display_color.h"
#include "gle_flute.h"

using namespace gle;
extern Modulo const* modulo;

//------------------------------------------------------------------------------

Display2::Display2(DisplayProp const* dp) : Display(dp)
{
}


void Display2::drawObjects(Simul const& sim)
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

    if (( prop->couple_select & 1 ) && ( sim.couples.sizeFF() > 0 ))
        drawCouplesF(sim.couples);
    
    if (( prop->couple_select & 2 ) && ( sim.couples.sizeA() > 0 ))
        drawCouplesA(sim.couples);
    
    if (( prop->single_select & 1 ) && ( sim.singles.sizeF() > 0 ))
        drawSinglesF(sim.singles);
    
    glBindBuffer(GL_PIXEL_UNPACK_BUFFER, 0);
    
    drawFibers(sim.fibers);

#if ( DIM >= 3 )
    glEnable(GL_LIGHTING);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
#else
    glDisable(GL_LIGHTING);
#endif
    
    drawBeads(sim.beads);
    drawSolids(sim.solids);
    drawSpheres(sim.spheres);
    
#if ( DIM >= 3 )
    glDisable(GL_LIGHTING);
    glDisable(GL_CULL_FACE);
#endif

    if (( prop->couple_select & 4 ) && ( sim.couples.sizeAA() > 0 ))
        drawCouplesB(sim.couples);
    
    if (( prop->single_select & 2 ) && ( sim.singles.sizeA() > 0 ))
        drawSinglesA(sim.singles);
    glBindBuffer(GL_PIXEL_UNPACK_BUFFER, 0);

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
        const PointDisp * disp = obj->disp();
        if ( disp->perceptible  &&  obj->fiber()->disp->visible )
        {
            Vector ph = obj->posHand();
            
            disp->drawA(ph);
            
            if ( obj->hasLink() && disp->width > 0 )
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
                stretchTube(pf, ph, disp->width*sizeScale, gle::truncatedCone);
                //drawCone(pf, ph-pf, disp->width*sizeScale);
#else
                gle::drawBand(ph, disp->width*sizeScale, ps, disp->width*sizeScale);
                gle::drawBand(ps, disp->width*sizeScale, disp->color, pf, disp->width*sizeScale, disp->color.alpha_scaled(0.5f));
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
    gle_color air(0,0,0,0);
    const PointDisp * pd1 = cx->disp1();
    const PointDisp * pd2 = cx->disp2();
    
    Vector p1 = cx->posHand1();
    Vector p2 = cx->posHand2();
    
    if ( modulo )
        modulo->fold(p2, p1);
    
    if ( pd1->perceptible || pd2->perceptible )
    {
        glEnableClientState(GL_COLOR_ARRAY);
        fluteVC* flu = gle::mapVertexColorBuffer(4);
#if 0
        flu[0] = { p1, pd1->color };
        flu[1] = { p2, pd2->color };
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
        gle_color col1 = pd1->visible ? pd1->color : air;
        gle_color col2 = pd2->visible ? pd2->color : air;
        flu[0] = { p1, col1 };
        flu[1] = { pp, col1 };
        flu[2] = { pp, col2 };
        flu[3] = { p2, col2 };
#endif
        gle::unmapVertexColorBuffer();
        lineWidth(pd1->width);
        glDrawArrays(GL_LINES, 0, 4);
        glDisableClientState(GL_COLOR_ARRAY);
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

