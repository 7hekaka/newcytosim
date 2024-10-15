// Cytosim was created by Francois Nedelec. Copyright Cambridge University 2021

#include "dim.h"
#include "simul.h"
#include "display3.h"
#include "modulo.h"
#include "random_pcg.h"
using namespace PCG32;

#include "gle.h"
#include "gym_color.h"
#include "gym_color_list.h"
#include "gym_view.h"
#include "gym_draw.h"
#include "gym_check.h"

#include "organizer.h"
#include "fiber_disp.h"
#include "line_disp.h"
#include "point_disp.h"
#include "display_color.cc"

/// Cliping planes permit nice junctions to be made between tubes, but are slow
#define USE_CLIP_PLANES 1

//------------------------------------------------------------------------------
/*
 Should use instanced rendering to draw multiple Sphere and Tube, etc.
 Possibly also for the tubes if we transfer 'tube length' to the GPU

 https://learnopengl.com/Advanced-OpenGL/Instancing
 https://www.khronos.org/opengl/wiki/Vertex_Rendering#Instancing
 https://metalbyexample.com/instanced-rendering/
 */

Display3::Display3(DisplayProp const* dp) : Display(dp)
{
}


void Display3::drawObjects(Simul const& sim)
{
    gym::closeDepthMask();
    gym::disableLighting();
    gym::disableCullFace();
    drawFields(sim.fields);
    
    gym::enableLighting();
    gym::openDepthMask();
    drawSpaces(sim.spaces);
    gym::disableCullFace();

    if ( stencil_ )
    {
        /*
         We use here the stencil test to make sure that nothing else is drawn
         where the inner side of the fibers is visible. This improves the
         display with clipping planes, as fibers appear as cut solid objects
         */
        glClearStencil(0);
        glClear(GL_STENCIL_BUFFER_BIT);
        glEnable(GL_STENCIL_TEST);
        glStencilFunc(GL_ALWAYS, 0, ~0);
        
        gym::enableCullFace(GL_FRONT);
        //drawFibers(sim.fibers);
        FiberSet const& set = sim.fibers;
        // set Stencil to 1 for inner surfaces of fibers:
        glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);
        GLint val = 0;
        // display the Fiber always in the same order:
        for( Fiber const* fib = set.firstID(); fib; fib=set.nextID(fib) )
        {
            if ( fib->disp->visible )
            {
                FiberDisp const* dis = fib->prop->disp;
                glStencilFunc(GL_ALWAYS, ++val, ~0);
                drawFiberLines(*fib, dis->line_style, dis->line_width);
            }
        }
        // set Stencil to 0 for outer surfaces:
        gym::switchCullFace(GL_BACK);
        glStencilOp(GL_KEEP, GL_KEEP, GL_ZERO);
        val = 0;
        for( Fiber const* fib = set.firstID(); fib; fib=set.nextID(fib) )
        {
            if ( fib->disp->visible )
            {
                glStencilFunc(GL_ALWAYS, ++val, ~0);
                drawFiber(*fib);
            }
        }

        glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);
        glStencilFunc(GL_EQUAL, 0, ~0);
    }
    else
    {
        /**
         If the display is 'cut', we might see the inner sides,
         but rendering would be faster with Culling enabled
        */
        // gym::enableCullFace(GL_BACK);
        drawFibers(sim.fibers);
    }
    
    drawFiberTexts(sim.fibers);
    gym::enableLighting();
    gym::enableCullFace(GL_BACK);

    drawBeads(sim.beads);
    drawSolids(sim.solids);
    drawSpheres(sim.spheres);
    
    gym::enableLighting();

    if (( prop->single_select & 1 ) && ( sim.singles.sizeF() > 0 ))
        drawSinglesF(sim.singles);
    
    if (( prop->couple_select & 1 ) && ( sim.couples.sizeFF() > 0 ))
        drawCouplesF(sim.couples);

    if (( prop->couple_select & 2 ) && ( sim.couples.sizeA() > 0 ))
        drawCouplesA(sim.couples);

    if (( prop->couple_select & 4 ) && ( sim.couples.sizeAA() > 0 ))
        drawCouplesB(sim.couples);

    if (( prop->single_select & 2 ) && ( sim.singles.sizeA() > 0 ))
        drawSinglesA(sim.singles);

    if ( stencil_ )
    {
        glClearStencil(0);
        glDisable(GL_STENCIL_TEST);
    }

    drawOrganizers(sim.organizers);
    gym::disableCullFace();
}


//------------------------------------------------------------------------------
#pragma mark - Drawing primitives


inline void Display3::drawPoint(Vector const& pos, PointDisp const* dis) const
{
    if ( dis->perceptible )
    {
        assertLighting();
        drawObject(pos, pixscale(dis->size), gle::sphere1);
    }
}


inline void Display3::drawObject3(Vector const& pos, float rad, void(*obj)()) const
{
    drawObject(pos, pixscale(rad), obj);
}


/// draw a point with a small sphere
inline void Display3::drawHand(Vector const& pos, PointDisp const* dis) const
{
    if ( dis->perceptible )
    {
        assertLighting();
        gym::color_both(dis->color);
        drawObject(pos, pixscale(dis->size), gle::blob);
    }
}

/// draw a point with a small sphere
inline void Display3::drawHandF(Vector const& pos, PointDisp const* dis) const
{
    if ( dis->perceptible )
    {
        assertLighting();
        gym::color_both(dis->color2);
        drawObject(pos, pixscale(dis->size), gle::blob);
    }
}


//------------------------------------------------------------------------------
#pragma mark - Nicely Joined Fiber Rendering using clipping planes

/**
This draws the model-segments, using function `select_color` to set display colors
*/
void Display3::drawFiberSegmentsClip(Fiber const& fib, float rad,
                                     gym_color (*select_color)(Fiber const&, index_t)) const
{
    const index_t last = fib.lastSegment();
    Vector old = fib.posPoint(0);
    Vector pos = fib.posPoint(1);
    Vector nxt, dir;
    gym::color_front(select_color(fib, 0));

    if ( last > 0 )
    {
        nxt = fib.posPoint(2);
        dir = normalize(nxt-old);
        gym::enableClipPlane(4);
        gym::setClipPlane(4, -dir, pos);
        gym::transAlignZ(old, rad, pos-old);
        gle::capedTube();
        gym::setClipPlane(4, dir, pos);
        
        // draw inner segments
        gym::enableClipPlane(5);
        for ( index_t i = 1; i < last; ++i )
        {
            old = pos;
            pos = nxt;
            nxt = fib.posPoint(i+2);
            dir = normalize(nxt-old);
            gym::color_front(select_color(fib, i));
            gym::setClipPlane(5, -dir, pos);
            gym::transAlignZ(old, rad, pos-old);
            gle::innerTube();
            gym::setClipPlane(4,  dir, pos);
        }
        gym::disableClipPlane(5);
    }
    else
    {
        dir = (pos-old) / fib.segmentation();
        nxt = pos;
        pos = old;
        old = 0.5 * ( nxt + pos );
        gym::enableClipPlane(4);
        gym::setClipPlane(4, -dir, old);
        gym::transAlignZ(pos, rad, nxt-pos);
        gle::capedTube();
        gym::setClipPlane(4, dir, old);
    }
    // draw last segment:
    gym::color_front(select_color(fib, last));
    gym::transAlignZ(nxt, rad, pos-nxt);
    gle::endedTube();
    gym::disableClipPlane(4);
}


/**
This draws segments of length 'len' which may not correspond to the vertices
used to model the Fiber. All abscissa is relative to the minus end.
The function `select_color` is called to set the color of the segments.
*/
void Display3::drawFiberSectionsClip(Fiber const& fib, float rad,
                                     int inx, const int last,
                                     real abs, const real inc,
                                     gym_color (*select_color)(Fiber const&, int, real),
                                     real fac, real facM, real facP) const
{
    Vector old = fib.displayPosM(abs);
    Vector pos = fib.displayPosM(abs+inc);
    Vector nxt = fib.displayPosM(abs+inc*2);
    Vector dir = normalize(nxt-pos);
    
    gym_color col = select_color(fib, inx++, facM);
    gym::enableClipPlane(4);
    if ( col.visible() )
    {
        gym::color_front(col);
        gym::setClipPlane(4, -dir, pos);
        gym::transAlignZ(old, rad, pos-old);
        if ( abs <= 0 )
            gle::capedTube();
        else
            gle::halfTube2();
    }
    gym::setClipPlane(4, dir, pos);
    
    // keep abs to match to the end of the section to be drawn
    abs += 2 * inc;

    gym::enableClipPlane(5);
    // draw segments
    while ( inx < last )
    {
        abs += inc;
        old = pos;
        pos = nxt;
        nxt = fib.displayPosM(abs);
        dir = normalize(nxt-old);
        gym::setClipPlane(5, -dir, pos);
        col = select_color(fib, inx++, fac);
        if ( col.visible() )
        {
            gym::color_front(col);
            gym::transAlignZ(old, rad, pos-old);
            gle::innerTube();
        }
        gym::setClipPlane(4, dir, pos);
    }
    gym::disableClipPlane(5);
    
    // draw last segment, which may be truncated:
    col = select_color(fib, last, facP);
    if ( col.visible() )
    {
        gym::color_front(col);
        if ( abs + REAL_EPSILON > fib.length() )
        {
            gym::stretchAlignZ1(nxt, -rad, fib.dirEndP(), -rad);
            gle::endedTube();
        }
        else
        {
            gym::transAlignZ(nxt, rad, pos-nxt);
            gle::halfTube2();
        }
    }
    gym::disableClipPlane(4);
}


//------------------------------------------------------------------------------
#pragma mark - draw Fibers using longer tubes to fill the gaps at the junctions

/**
This draws the model-segments, using function `select_color` to set display colors
*/
void Display3::drawFiberSegmentsJoin(Fiber const& fib, float rad,
                                     gym_color (*select_color)(Fiber const&, index_t)) const
{
    const index_t last = fib.lastSegment();
    Vector pos = fib.posPoint(0);
    Vector nxt = fib.posPoint(1);
    
    gym::color_front(select_color(fib, 0));
    gym::transAlignZ(pos, rad, nxt-pos);
    gle::dome();
    gym::scale(1, 1, fib.segmentation()/rad);
    if ( last == 0 )
    {
        gle::tube4();
        gle::discTop2();
        return;
    }
    gle::tubeS();

    for ( index_t i = 1; i < last; ++i )
    {
        pos = nxt;
        nxt = fib.posPoint(i+1);
        gym::color_front(select_color(fib, i));
        gym::stretchAlignZ(pos, nxt, rad);
        gle::tubeM();
    }
    pos = nxt;
    nxt = fib.posPoint(last+1);
    gym::color_front(select_color(fib, last));
    gym::stretchAlignZ(pos, nxt, rad);
    gle::tubeE();
    gle::discTop2();
}


/**
This draws segments of length 'len' which may not correspond to the vertices
used to model the Fiber. All abscissa is relative to the minus end.
The function `select_color` is called to set the color of the segments.
*/
void Display3::drawFiberSectionsJoin(Fiber const& fib, float rad,
                                     int inx, const int last,
                                     real abs, const real inc,
                                     gym_color (*select_color)(Fiber const&, int, real),
                                     real fac, real facM, real facP) const
{
    Vector pos = fib.displayPosM(abs);
    Vector nxt = fib.displayPosM(abs+inc);
    
    gym_color col = select_color(fib, inx++, facM);
    if ( col.visible() )
    {
        gym::color_front(col);
        if ( abs <= 0 )
        {
            real len = (nxt-pos).norm();
            gym::stretchAlignZ1(pos, rad, (nxt-pos)/len, rad);
            gle::dome();
            gym::scale(1, 1, len/rad);
        }
        else
        {
            gym::stretchAlignZ(pos, nxt, rad);
        }
        if ( last == inx )
        {
            gle::tube4();
            gle::discTop2();
            return;
        }
        gle::tubeS();
    }
    // keep abs to match to the end of the section already drawn
    abs += inc;
    
    while ( inx < last )
    {
        abs += inc;
        pos = nxt;
        nxt = fib.displayPosM(abs);
        col = select_color(fib, inx++, fac);
        if ( col.visible() )
        {
            gym::color_front(col);
            gym::stretchAlignZ(pos, nxt, rad);
            gle::tubeM();
        }
    }
    // draw last segment, which may be truncated:
    col = select_color(fib, last, facP);
    if ( col.visible() )
    {
        gym::color_front(col);
        if ( abs+inc >= fib.length() )
        {
            gym::stretchAlignZ1(nxt, rad, fib.dirEndP(), fib.length()-abs);
            gle::discTop2();
        }
        else
        {
            gym::stretchAlignZ(nxt, fib.displayPosM(abs+inc), rad);
        }
        gle::tubeE();
    }
}

//------------------------------------------------------------------------------
#pragma mark -

#if USE_CLIP_PLANES
#  define drawFiberSegments drawFiberSegmentsClip
#  define drawFiberSections drawFiberSectionsClip
#else
#  define drawFiberSegments drawFiberSegmentsJoin
#  define drawFiberSections drawFiberSectionsJoin
#endif

void Display3::drawFiberLines(Fiber const& fib, int style, float width) const
{
    FiberDisp const*const dis = fib.prop->disp;
    const float rad = pixscale(width);

    // set back color:
    if ( dis->coloring )
        gym::color_back(fib.disp->color);
    else
        gym::color_back(dis->back_color);
    gym::enableLighting();
    
    switch ( style )
    {
        case 1:
            drawFiberSegments(fib, rad, color_fiber);
            break;
        case 2:
            drawFiberSegments(fib, rad, color_by_tension);
            break;
        case 3:
            drawFiberSegments(fib, rad, color_by_tension_jet);
            break;
        case 4:
            drawFiberSegments(fib, rad, color_by_direction);
            break;
        case 5:
            drawFiberSegments(fib, rad, color_seg_curvature);
            break;
        case 6: {
            /** This is using transparency with segments that are not depth sorted
             but this code is only used in 2D normally, so it's okay */
            gym::disableCullFace();
            drawFiberSegments(fib, rad, color_by_abscissaM);
            gym::restoreCullFace();
        } break;
        case 8: if ( fib.endStateM() == STATE_GREEN )
        case 7: {
            /** This is using transparency with segments that are not depth sorted
             but this code is only used in 2D normally, so it's okay */
            gym::enableCullFace(GL_BACK);
            drawFiberSegments(fib, rad, color_by_abscissaP);
            gym::restoreCullFace();

        } break;
        case 9:
            drawFiberSegments(fib, rad, color_by_height);
            break;
        case 10:
            drawFiberSegments(fib, rad, color_by_grid);
            break;
    }
}


// displays segment 'inx' with transparency
void Display3::drawFiberSegmentT(Fiber const& fib, unsigned inx) const
{
    FiberDisp const*const dis = fib.prop->disp;
    const int style = dis->line_style;
    const real iseg = fib.segmentationInv();
    real rad = pixscale(dis->line_width);
    Vector A = fib.posP(inx);
    Vector B = fib.posP(inx+1);

    if ( style == 6 )
        gym::color_front(color_by_abscissaM(fib, inx));
    else if ( style == 7 || style == 8 )
        gym::color_front(color_by_abscissaP(fib, inx));
    else if ( style == 2 )
        gym::color_front(color_by_tension(fib, inx));
    else if ( style == 3 )
        gym::color_front(color_by_tension_jet(fib, inx));
    else
        gym::color_front(fib.disp->color);

    // truncate terminal segment according to length_scale
    if ( inx == 0 && style == 6 )
    {
        real x = 3 * dis->length_scale * iseg;
        if ( x < 1.0 )
        {
            B = A + x * ( B - A );
            color_by_abscissaM(fib, inx);
        }
    }
    
    // truncate terminal segment according to length_scale
    if ( inx == fib.lastSegment() && ( style == 7 || style == 8 ) )
    {
        real x = 3 * dis->length_scale * iseg;
        if ( x < 1.0 )
        {
            A = B + x * ( A - B );
            color_by_abscissaP(fib, inx);
        }
    }
    
    gym::enableLighting();
    /* Either CULL_FACE should be enable to hide the back side,
     or every primitive should be renderred with a double pass*/
    gym::enableCullFace(GL_BACK);

#if USE_CLIP_PLANES
    if ( inx == 0 )
    {
        gym::enableClipPlane(5);
        if ( inx == fib.lastSegment() )
        {
            gym::setClipPlane(5, (A-B)*iseg, (A+B)*0.5);
            gym::transAlignZ(A, rad, B-A);
            gle::capedTube();
            gym::setClipPlane(5, (B-A)*iseg, (A+B)*0.5);
            gym::transAlignZ(B, rad, A-B);
            gle::endedTube();
        }
        else
        {
            gym::setClipPlane(5, (A-B)*iseg, B);
            gym::transAlignZ(A, rad, B-A);
            gle::capedTube();
        }
        gym::disableClipPlane(5);
        return;
    }
    else
    {
        gym::setClipPlane(5, normalize(B-fib.posP(inx-1)), A);
    }
    
    if ( inx == fib.lastSegment() )
    {
        gym::enableClipPlane(4);
        gym::setClipPlane(4, (B-A)*iseg, A);
        gym::transAlignZ(B, rad, A-B);
        gle::endedTube();
        gym::disableClipPlane(4);
        return;
    }
    else
    {
        gym::setClipPlane(4, normalize(A-fib.posP(inx+2)), B);
    }
    
    gym::enableClipPlane(5);
    gym::enableClipPlane(4);
    gym::transAlignZ(A, rad, B-A);
    gle::innerTube();
    gym::disableClipPlane(4);
    gym::disableClipPlane(5);
#else
    gym::transAlignZ(A, rad, B-A);
    if ( inx == 0 )
        gle::dome();
    gym::scale(1, 1, fib.segmentation()/rad);
    gle::tube4();
    if ( inx == fib.lastSegment() )
        gle::discTop2();
#endif
    gym::restoreCullFace();
}


//------------------------------------------------------------------------------
#pragma mark - Display Lattice


void Display3::drawFiberLattice(Fiber const& fib, VisibleLattice const& lat, float rad,
                                gym_color (*select_color)(Fiber const&, int, real)) const
{
    FiberDisp const*const dis = fib.prop->disp;

    gym::enableLighting();
    GLfloat blk[] = { 0.f, 0.f, 0.f, 1.f };
    GLfloat bak[] = { 0.f, 0.f, 0.f, 1.f };
    dis->back_color.store(bak);
    glMaterialfv(GL_BACK, GL_AMBIENT,  blk);
    glMaterialfv(GL_BACK, GL_DIFFUSE,  bak);
    glMaterialfv(GL_FRONT, GL_AMBIENT,  blk);
    glMaterialfv(GL_FRONT, GL_DIFFUSE,  blk);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, blk);
    glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, blk);
    glMateriali (GL_FRONT_AND_BACK, GL_SHININESS, 32);

    const real fac = 1 / dis->lattice_scale;
    const real uni = lat.unit();
    const auto inf = lat.indexM();
    const auto sup = lat.indexP();

    real abs = uni * inf - fib.abscissaM();       // should be non-positive!
    assert_true( abs <= 0 && 0 <= abs+uni );

    real facM = fac;
    real facP = fac;
    
    if ( dis->lattice_rescale )
    {
        real lenM = abs + uni;                     // should be positive!
        real lenP = fib.abscissaP() - uni * sup;   // should be positive!
        facM = ( lenM > 0.001*uni ? fac*uni/lenM : fac );
        facP = ( lenP > 0.001*uni ? fac*uni/lenP : fac );
    }

    drawFiberSections(fib, rad, inf, sup, abs, uni, select_color, fac, facM, facP);
}


void Display3::drawFiberLattice1(Fiber const& fib, VisibleLattice const& lat, float rad) const
{
    drawFiberLattice(fib, lat, pixscale(rad), color_by_lattice);
}

void Display3::drawFiberLattice2(Fiber const& fib, VisibleLattice const& lat, float rad) const
{
    drawFiberLattice(fib, lat, pixscale(rad), color_by_lattice_jet);
}

void Display3::drawFiberLattice3(Fiber const& fib, VisibleLattice const& lat, float rad) const
{
    drawFiberLattice(fib, lat, pixscale(rad), color_by_lattice_white);
}

void Display3::drawFiberLatticeEdges(Fiber const& fib, VisibleLattice const& lat, float rad) const
{
#if 0
    gym_color col = fib.dis->color;
    gym_color lor = col.darken(0.75);
    const real uni = lat.unit();
    drawFiberStriped(fib, pixscale(rad), uni, col, uni, lor);
#endif
    drawFiberLattice(fib, lat, pixscale(rad), color_alternate);
}

//------------------------------------------------------------------------------
#pragma mark -


/**
 Display the minus end of a Fiber, according to `style`:
 - 1: sphere
 - 2: cone
 - 3: flat cylinder
 - 4: arrow-head
 - 5: arrow-head in reverse direction
 - 6: cube
 - 7: cylinder placed backward so as to overlap with the fiber
 - 8: hemisphere at the tip
 .
 with 3D objects
 */
void Display3::drawFiberEndMinus(Fiber const& fib, int style, float size) const
{
    const float rad = pixscale(size);
    if ( rad > pixelSize ) switch(style)
    {
        default: break;
        case 1: drawObject(fib.posEndM(), rad, gle::sphere2); break;
        case 2: drawObject(fib.posEndM(),-fib.dirEndM(), rad, gle::longCone); break;
        case 3: drawObject(fib.posEndM(),-fib.dirEndM(), rad, gle::cylinder1); break;
        case 4: drawObject(fib.posEndM(),-fib.dirEndM(), rad, gle::arrowTail); break;
        case 5: drawObject(fib.posEndM(), fib.dirEndM(), rad, gle::arrowTail); break;
        case 6: drawObject(fib.posEndM(),-fib.dirEndM(), rad, gle::cube); break;
        case 7: drawObject(fib.posEndM(), fib.dirEndM(), rad, gle::cylinder1); break;
        case 8: drawObject(fib.posEndM(), fib.dirEndM(), rad, gle::hemisphere2); break;
    }
}


/**
 Display the plus end of a Fiber, according to `style`:
 - 1: sphere
 - 2: cone
 - 3: flat cylinder
 - 4: arrow-head
 - 5: arrow-head in reverse direction
 - 6: cube
 - 7: cylinder placed backward so as to overlap with the fiber
 - 8: hemisphere at the tip
 .
 with 3D objects
 */
void Display3::drawFiberEndPlus(Fiber const& fib, int style, float size) const
{
    const float rad = pixscale(size);
    if ( rad > pixelSize ) switch(style)
    {
        default: break;
        case 1: drawObject(fib.posEndP(), rad, gle::sphere2); break;
        case 2: drawObject(fib.posEndP(), fib.dirEndP(), rad, gle::longCone); break;
        case 3: drawObject(fib.posEndP(), fib.dirEndP(), rad, gle::cylinder1); break;
        case 4: drawObject(fib.posEndP(), fib.dirEndP(), rad, gle::arrowTail); break;
        case 5: drawObject(fib.posEndP(),-fib.dirEndP(), rad, gle::arrowTail); break;
        case 6: drawObject(fib.posEndP(), fib.dirEndP(), rad, gle::cube); break;
        case 7: drawObject(fib.posEndP(),-fib.dirEndP(), rad, gle::cylinder1); break;
        case 8: drawObject(fib.posEndP(),-fib.dirEndP(), rad, gle::hemisphere2); break;
    }
}


void Display3::drawFiberSpeckles(Fiber const& fib) const
{
    FiberDisp const*const dis = fib.prop->disp;
    const float rad = pixscale(dis->speckle_size);
    if ( rad < pixelSize )
        return;
    gym::color_front(fib.disp->color);
    gym::color_back(dis->back_color);
    assertLighting();

    // display random speckles:
    if ( dis->speckle_style == 1 )
    {
        /*
         A simple random number generator seeded by fib.signature()
         is used to distribute points always at the same position
         with respect to the lattice of each fiber.
         */
        
        const real gap = dis->speckle_gap;
        constexpr real TINY = 0x1p-32;
        // draw speckles below the origin of abscissa:
        if ( fib.abscissaM() < 0 )
        {
            uint64_t Z = pcg32_init(fib.signature());
            real a = gap * std::log(pcg32(Z)*TINY);
            while ( a > fib.abscissaP() )
                a += gap * std::log(pcg32(Z)*TINY);
            while ( a >= fib.abscissaM() )
            {
                gym::transScale(fib.pos(a), rad);
                gle::icosahedron();
                a += gap * std::log(pcg32(Z)*TINY);
            }
        }
        // draw speckles above the origin of abscissa:
        if ( fib.abscissaP() > 0 )
        {
            uint64_t Z = pcg32_init(~fib.signature());
            real a = -gap * std::log(pcg32(Z)*TINY);
            while ( a < fib.abscissaM() )
                a -= gap * std::log(pcg32(Z)*TINY);
            while ( a <= fib.abscissaP() )
            {
                gym::transScale(fib.pos(a), rad);
                gle::icosahedron();
                a -= gap * std::log(pcg32(Z)*TINY);
            }
        }
    }
    else if ( dis->speckle_style == 2 )
    {
        //we distribute points regularly along the center line
        const real gap = dis->speckle_gap;
        real a = gap * std::ceil( fib.abscissaM() / gap );
        while ( a <= fib.abscissaP() )
        {
            gym::transScale(fib.pos(a), rad);
            gle::icosahedron();
            a += gap;
        }
    }
}


void Display3::drawFiberPoints(Fiber const& fib) const
{
    FiberDisp const*const dis = fib.prop->disp;
    // diameter of lines and points in space units:
    const float rad = pixscale(dis->point_size);
    int style = dis->point_style;
    gym::color_front(fib.disp->color);
    gym::color_back(dis->back_color);

    if ( rad > pixelSize )
    {
        if ( style == 1 )
        {
            gym::enableLighting();
            // display vertices:
            for ( unsigned i = 0; i < fib.nbPoints(); ++i )
                drawObject(fib.posP(i), rad, gle::cube);
        }
        else if ( style == 2 )
        {
            gym::enableLighting();
            // display arrowheads along the fiber:
            const real gap = dis->point_gap;
            real ab = std::ceil(fib.abscissaM()/gap) * gap;
            for ( ; ab <= fib.abscissaP(); ab += gap )
                gle::drawCone(fib.pos(ab), fib.dir(ab), rad);
        }
        else if ( style == 3 )
        {
            gym::enableLighting();
            // display cones regularly along the fiber:
            const real gap = dis->point_gap;
            real ab = std::ceil(fib.abscissaM()/gap) * gap;
            for ( ; ab <= fib.abscissaP(); ab += gap )
                gle::drawCone(fib.pos(ab), -fib.dir(ab), rad);
        }
        else if ( style == 4 )
        {
            // display middle of fiber:
            drawObject(fib.posMiddle(), 2*rad, gle::sphere2);
        }
    }
}

//------------------------------------------------------------------------------

void Display3::drawOrganizer(Organizer const& obj) const
{
    Solid const* sol = obj.solid();
    Sphere const* sph = obj.sphere();
    if ( !sol && !sph ) return;
    PointDisp const* dis = ( sol ? sol->prop->disp : sph->prop->disp );
    if ( !dis ) return;
    gym_color col;
    if ( sol ) col = bodyColorF(*sol).match_a(dis->color);
    if ( sph ) col = bodyColorF(*sph).match_a(dis->color);

    if ( dis->style & 2 )
    {
        // draw links between Solid/Sphere and Mecables
        Vector P, Q;
        gym::color_front(col);
        const float wid = pixscale(dis->width);

        for ( index_t i = 0; obj.getLink(i, P, Q); ++i )
        {
            drawPoint(P, dis);
            if ( modulo ) modulo->fold(Q, P);
            gym::stretchAlignZ(P, Q, wid);
            gle::tube1();
        }
    }
    /**
     This displays the Solid connecting two Aster as a spindle.
     Used for Cleo Kozlowski simulation of C. elegans (2007)
     */
    if (( dis->style & 1 ) && obj.tag() == Organizer::FAKE_TAG )
    {
        if ( sol && sol->nbPoints() >= 4 )
        {
            gym::color_front(col);
            gym::color_back(col.darken(0.5));
#if ( DIM >= 3 )
            Vector3 a = 0.5 * (sol->posP(0) + sol->posP(2));
            Vector3 b = 0.5 * (sol->posP(1) + sol->posP(3));
            gym::stretchAlignZ(a, b, 1);
            gle::dualPassBarrel();
#else
            const float wid = pixscale(dis->width);
            for ( index_t i = 0; i < sol->nbPoints(); i+=2 )
            {
                gym::stretchAlignZ(sol->posPoint(i), sol->posPoint(i+1), wid);
                gle::hexTube();
            }
#endif
        }
    }
}

//------------------------------------------------------------------------------
#pragma mark -

void Display3::drawSinglesF(SingleSet const& set) const
{
    assertLighting();
    for ( Single * obj=set.firstF(); obj ; obj=obj->next() )
    {
        PointDisp const* dis = obj->disp();
        if ( dis->perceptible )
        {
            //drawHandF(obj->posFoot(), obj->disp());
            if ( obj->base() )
                gym::color_both(dis->color2);
            else
                gym::color_both(dis->color2.tweak(obj->signature()));
#if ( 0 )
            // draw anchored single using a flatten dome, instead of a sphere
            Space const* spc = obj->confineSpace();
            if ( spc )
            {
                const float rad = pixscale(dis->size);
                /// draw a flat dome tangent to the Space:
                Vector pos = obj->posFoot();
                Vector dir = spc->normalToEdge(pos);
                drawObject(pos, dir, rad, gle::dome);
            }
            else
#endif
            drawObject(obj->posFoot(), pixscale(dis->size), gle::blob);
        }
    }
}


void Display3::drawSingleA(Single const* obj) const
{
    PointDisp const* dis = obj->disp();
    Vector ph = obj->posHand();
    gym::color_both(dis->color);
    drawHand(ph, dis);
}


void Display3::drawSingleB(Single const* obj) const
{
    PointDisp const* dis = obj->disp();

    if ( dis->perceptible )
    {
        Vector ph = obj->posHand();
        Vector pf = obj->posFoot();
        if ( modulo ) modulo->fold(pf, ph);
        const float wid = pixscale(dis->width);
        const float rad = pixscale(dis->size);

        gym::color_both(dis->color);
#if ( 0 )
        // draw a dome tangent to the Space:
        Space const* spc = obj->confineSpace();
        if ( spc )
            drawObject(pf, spc->normalToEdge(pf), 1.4142f*rad, gle::dome);
        else
#endif
        {
            gym::transScale(pf, wid);
            gle::blob(); // the foot
        }
        //gym::color_both(dis->color);
#if ( DIM > 2 )
        Vector dif = pf - ph;
        float L = norm(dif);
        
        if ( L > 1e-3 )
        {
            // displace hands away from the fiber's centerline:
            const real R = pixscale(obj->fiber()->prop->disp->line_width + 0.25 * dis->size);
            Vector T = obj->dirFiber();
            ph += ( dif - dot(dif,T) * T ) * min_real(0.45, R/L);
            dif = pf - ph;
            L = norm(dif);
        }
        gym::transAlignZ(ph, rad, dif/L);
        gle::blob();
        gym::scale(wid/rad, wid/rad, L/rad);
        gle::hexTube();
#elif ( 0 )
        Vector dir = normalize( pf - ph );
        gym::enableClipPlane(5);
        gym::setClipPlane(5, -dir, pf);
        gym::transAlignZ(ph, rad, dir);
        gle::needle();
        gym::disableClipPlane(5);
#else
        if ( obj->base() )
            drawObject(ph, rad, gle::octahedron);
        else
            drawHand(ph, dis);
        gle::drawBand(ph, wid, dis->color, pf, wid, dis->color.alpha_scaled(0.5f));
#endif
    }
}

void Display3::drawSinglesA(SingleSet const& set) const
{
    for ( Single * obj=set.firstA(); obj ; obj=obj->next() )
    {
        if ( obj->disp()->perceptible && obj->fiber()->disp->visible )
        {
            if ( obj->hasLink() )
                drawSingleB(obj);
            else
                drawSingleA(obj);
        }
    }
}

//------------------------------------------------------------------------------
#pragma mark - Couple

void Display3::drawCouplesF(CoupleSet const& set) const
{
    if ( prop->couple_flip )
        drawCouplesF2(set);
    else
        drawCouplesF1(set);
}


/**
 Display always Hand1 of Couple
 */
void Display3::drawCouplesF1(CoupleSet const& set) const
{
    assertLighting();
    for ( Couple * obj = set.firstFF(); obj ; obj=obj->next() )
    {
        //drawHandF(cx->posFree(), cx->disp1());
        PointDisp const* dis = obj->disp1();
        if ( dis->perceptible )
        {
            gym::color_both(dis->color2.tweak(obj->signature()));
            drawObject(obj->posFree(), pixscale(dis->size), gle::blob);
        }
    }
}


/**
 Display either Hand1 or Hand2, exposing both sides with equal chances.
 This gives the impression that Couple flicker randomly between frames,
 as if they were two-sided balls 'rotating' very fast.
 */
void Display3::drawCouplesF2(CoupleSet const& set) const
{
    Couple * nxt;
    Couple * obj = set.firstFF();
    
    if ( set.sizeFF() & 1 )
    {
        nxt = obj->next();
        drawHandF(obj->posFree(), obj->disp12());
        obj = nxt;
    }
    while ( obj )
    {
        nxt = obj->next();
        drawHandF(obj->posFree(), obj->disp21());
        obj = nxt->next();
        drawHandF(nxt->posFree(), nxt->disp12());
    }
}


void Display3::drawCouplesA(CoupleSet const& set) const
{
    //const float R = pixscale(0.5f);
    for ( Couple * cx=set.firstAF(); cx ; cx=cx->next() )
    {
        Hand const* h = cx->hand1();
        PointDisp const* dis = h->disp();

        if ( h->fiber()->disp->visible && dis->perceptible )
        {
#if ( 0 )  // ENDOCYTOSIS 2015
            if ( cx->fiber1()->disp->color.transparent() )
            {
                gym::color_both(dis->color, cx->fiber1()->disp->color.transparency());
                drawPoint(h->pos(), dis);
                continue;
            }
#endif
            drawHand(h->pos(), dis);
            //drawObject(h->outerPos(), h->pos()-h->outerPos(), R*dis->size, gle::tetrahedron);
        }
    }
    
    for ( Couple * cx=set.firstFA(); cx ; cx=cx->next() )
    {
        Hand const* h = cx->hand2();
        PointDisp const* dis = h->disp();

        if ( cx->fiber2()->disp->visible && cx->disp2()->perceptible )
        {
#if ( 0 )  // ENDOCYTOSIS 2015
            if ( cx->fiber2()->disp->color.transparent() )
            {
                gym::color_both(dis->color, cx->fiber2()->disp->color.transparency());
                drawPoint(h->pos(), dis);
                continue;
            }
#endif
            drawHand(h->pos(), dis);
            //drawObject(h->outerPos(), h->pos()-h->outerPos(), R*dis->size, gle::tetrahedron);
        }
    }
}

/** two blobs joined by a tube */
void Display3::drawCoupleBcrude(Couple const* cx) const
{
    PointDisp const* pd1 = cx->disp1();
    PointDisp const* pd2 = cx->disp2();
    
    Vector p1 = cx->posHand1();
    Vector p2 = cx->posHand2();

    gym::color_both(pd1->color);
    gym::stretchAlignZ(p1, p2, pixscale(pd1->width));
    gle::hexTube();
    
    if ( pd1->visible ) drawHand(p1, pd1);
    if ( pd2->visible ) drawHand(p2, pd2);
}


/** each hand is represented by two feet and the position of these feet
 are set as a function of hand's abscissa, to immitate a walk:
 one foot is on the surface of the fiber, and the other one away from it.
 In ICE back from Bad Honnef, 4.10.2024
 */
void Display3::drawCoupleBwalk(Couple const* cx) const
{
    const float DK = 0.625; // factor by which one hand is darkened
    PointDisp const* pd1 = cx->disp1();
    PointDisp const* pd2 = cx->disp2();
    
    Vector P1 = cx->posHand1();
    Vector P2 = cx->posHand2();
    
    Vector dif = P2 - P1;
    real dns = dif.normSqr();
    
    if ( dns > 1e-6 )
    {
        // displace hands away from the fiber's centerline:
        dns = 1.0 / std::sqrt(dns);
        // position the heads at the surface of the filaments:
        const real R1 = pixscale(cx->fiber1()->prop->disp->line_width + 0.25 * pd1->size);
        const real R2 = pixscale(cx->fiber2()->prop->disp->line_width + 0.25 * pd2->size);
        // move points orthogonal to the fiber's axis
        Vector T1 = cx->dirFiber1();
        Vector T2 = cx->dirFiber2();
        Vector h1 = ( dif - dot(dif,T1) * T1 ) * min_real(0.45, R1*dns);
        Vector h2 = ( dot(dif,T2) * T2 - dif ) * min_real(0.45, R2*dns);
        
        gym::color_both(pd2->color.mix(pd1->color));
        gym::stretchAlignZ(P1+h1, P2+h2, pixscale(pd1->width));
        gle::thinTube();

#if ( DIM > 2 )
        // the axes around which uplifted hands revolve are tilted w/r dif
        T1 = dif + cross(dif, T1);
        T2 = cross(T2, dif) - dif;
#endif
        // grounded foot is at position 'b' and uplifted one at (x, y):
        const real R = 0.008;
        Fiber const* fib = cx->fiber1();
        real a = cx->hand1()->abscissa();
        real i = std::round(a/R);
        real b = i * R;
        real x = 2 * a - b;
        real y = std::sqrt((R*R - square(x-b))) * ( dns * 0.4 );
        
        float S = pixscale(0.5*pd1->size);
        gym_color col = pd1->color;
        gym_color lor = pd1->color;
        if ( int(i) & 1 ) col = col.darken(DK); else lor = lor.darken(DK);
        gym::color_both(col);
        gym::transScale(h1+fib->pos(b), S);
        gle::blob();
        gym::color_both(lor);
        gym::transScale(h1+fib->pos(x)+y*T1, S);
        gle::blob();

        // grounded foot is at position 'b' and uplifted one at (x, y):
        fib = cx->fiber2();
        a = cx->hand2()->abscissa();
        i = std::round(a/R);
        b = i * R;
        x = 2 * a - b;
        y = std::sqrt((R*R - square(x-b))) * ( dns * 0.4 );
        
        S = pixscale(0.5*pd2->size);
        col = pd2->color;
        lor = pd2->color;
        if ( int(i) & 1 ) col = col.darken(DK); else lor = lor.darken(DK);
        gym::color_both(col);
        gym::transScale(h2+fib->pos(b), S);
        gle::blob();
        gym::color_both(lor);
        gym::transScale(h2+fib->pos(x)+y*T2, S);
        gle::blob();
    }
}


void Display3::drawCoupleBhomo(Couple const* cx, PointDisp const* dis) const
{
    Vector p1 = cx->posHand1();
    Vector p2 = cx->posHand2();

#if !FIBER_HAS_FAMILY
    Vector dif = p2 - p1;
    real dns = dif.normSqr();
    if ( dns > 1e-6 )
    {
        // displace hands away from the fiber's centerline:
        dns = 1.0 / std::sqrt(dns);
        const real R1 = pixscale(cx->fiber1()->prop->disp->line_width + 0.4 * dis->size);
        const real R2 = pixscale(cx->fiber2()->prop->disp->line_width + 0.4 * dis->size);
        // move points orthogonal to the fiber's axis
        Vector t1 = cx->dirFiber1();
        Vector t2 = cx->dirFiber2();
        p1 += ( dif - dot(dif,t1) * t1 ) * min_real(0.45, R1*dns);
        p2 -= ( dif - dot(dif,t2) * t2 ) * min_real(0.45, R2*dns);
    }
#endif
    float R = pixscale(dis->size);
    gym::color_both(dis->color);
    //gym::stretchAlignZ(p1, p2, R);
    gym::transAlignZ(p1, R, p2-p1);
    gle::droplet();
    //gym::stretchAlignZ(p2, p1, R);
    gym::transAlignZ(p2, R, p1-p2);
    gle::droplet();
}


void Display3::drawCoupleBside(Couple const* cx) const
{
    PointDisp const* pd1 = cx->disp1();
    PointDisp const* pd2 = cx->disp2();
    
    Vector p1 = cx->posHand1();
    Vector pS = cx->sidePos1();
    Vector p2 = cx->posHand2();
    if ( modulo ) modulo->fold(p2, p1);
    
#if FIBER_HAS_FAMILY
    if ( pd1 == pd2 )
    {
        // semi-accurate rendering of Couple's side-side link
        gym::color_both(pd1->color);
        Vector mid = 0.5 * ( cx->sidePos1() + cx->sidePos2() );
        drawPoint(mid, pd1);
        gym::stretchAlignZ(p2, mid, pixscale(pd2->width));
        gle::hexTube();
        gym::stretchAlignZ(p1, mid, pixscale(pd1->width));
        gle::hexTube();
        drawPoint(p1, pd1);
        drawPoint(p2, pd2);
        return;
    }
#endif
    
    float rad = pixscale(pd1->size);
    float Lr = cx->prop->length / rad;
    float iLr = ( pd1->width / pd1->size );
    
    //if ( cx->cosAngle() > 0 ) gym::color_both(1, 0.5, 0.25, 1)); else gym::color_both(0, 1, 0, 1);
    if ( pd1->visible )
    {
        float Z = cx->fiber1()->prop->disp->line_width / pd1->size + 0.4;
        gym::color_both(pd1->color);
        gym::transAlignZ(p1, rad, pS-p1);
        gym::translate(0, 0, Z);
        gle::blob();
        gym::translate(0, 0,Lr-Z);
        gle::cuboid();
        gym::shift(0, 0,-Lr);
        gym::scale(iLr, iLr, Lr);
        gle::hexTube();
    }

    // draw a link between pS and p2
    gym::transAlignZ(p2, rad, pS-p2);
    if ( pd2->visible )
    {
        gym::color_both(pd2->color);
        gle::blob();
    }
    gym::scale(iLr, iLr, norm(pS-p2) / rad);
    gle::hexTube();
}



void Display3::drawCoupleBori(Couple const* cx) const
{
    PointDisp const* pd1 = cx->disp1();
    PointDisp const* pd2 = cx->disp2();
    
    Vector p1 = cx->posHand1();
    Vector p2 = cx->posHand2();
    if ( modulo ) modulo->fold(p2, p1);
    Vector dif = p2 - p1;
    real dns = dif.normSqr();
    
    if ( dns > 1e-6 )
    {
#if !FIBER_HAS_FAMILY
        // displace hands away from the fiber's centerline:
        dns = 1.0 / std::sqrt(dns);
        const real R1 = pixscale(cx->fiber1()->prop->disp->line_width + 0.4 * pd1->size);
        const real R2 = pixscale(cx->fiber2()->prop->disp->line_width + 0.4 * pd2->size);
        // move points orthogonal to the fiber's axis
        Vector dir1 = cx->dirFiber1();
        Vector dir2 = cx->dirFiber2();
        p1 += ( dif - dot(dif,dir1) * dir1 ) * min_real(0.45, R1*dns);
        p2 -= ( dif - dot(dif,dir2) * dir2 ) * min_real(0.45, R2*dns);
        dif = p2 - p1;
#endif
        if ( pd1->visible )
        {
            float R1 = pixscale(pd1->size);
            gym::color_both(pd1->color);
            //gym::stretchAlignZ(p1, p2, R1);
            gym::transAlignZ(p1, R1, dif);
            gle::droplet();
        }
        if ( pd2->visible )
        {
            float R2 = pixscale(pd2->size);
            gym::color_both(pd2->color);
            //gym::stretchAlignZ(p2, p1, R2);
            gym::transAlignZ(p2, -R2, dif);
            gle::droplet();
        }
    }
    else
    {
        if ( pd1->visible ) drawHand(p1, pd1);
        if ( pd2->visible ) drawHand(p2, pd2);
    }
}


void Display3::drawCoupleBalt(Couple const* cx) const
{
    PointDisp const* pd1 = cx->disp1();
    PointDisp const* pd2 = cx->disp2();
    
    Vector p1 = cx->posHand1();
    Vector p2 = cx->posHand2();
    if ( modulo ) modulo->fold(p2, p1);
    
#if !FIBER_HAS_FAMILY
    Vector dif = p2 - p1;
    real dns = dif.normSqr();
    if ( dns > 1e-6 )
    {
        // displace hands away from the fiber's centerline:
        dns = 1.0 / std::sqrt(dns);
        const real R1 = pixscale(cx->fiber1()->prop->disp->line_width + 0.4 * pd1->size);
        const real R2 = pixscale(cx->fiber2()->prop->disp->line_width + 0.4 * pd2->size);
        // move points orthogonal to the fiber's axis
        Vector dir1 = cx->dirFiber1();
        Vector dir2 = cx->dirFiber2();
        p1 += ( dif - dot(dif,dir1) * dir1 ) * min_real(0.45, R1*dns);
        p2 -= ( dif - dot(dif,dir2) * dir2 ) * min_real(0.45, R2*dns);
    }
#endif
    
    float wid = pixscale(pd1->width);
    float R = pixscale(pd1->size) / wid;
    float Lr = wid / norm( p2 - p1 );
    
    gym::stretchAlignZ(p1, p2, wid);
    gym::color_both(pd1->color);
    gle::hexTube();
    if ( pd1->visible )
    {
        gym::scale(R, R, R*Lr);
        gle::blob();
    }
    if ( pd2->visible )
    {
        gym::translate(0, 0, 1/(R*Lr));
        gym::color_both(pd2->color);
        gle::blob();
    }
}


void Display3::drawCoupleB(Couple const* cx) const
{
    PointDisp const* pd1 = cx->disp1();
    PointDisp const* pd2 = cx->disp2();
    
    if ( pd1 == pd2 )
    {
        if ( pd1->visible )
            drawCoupleBwalk(cx);
    }
    else
    {
        Vector p1 = cx->posHand1();
        Vector p2 = cx->posHand2();
        if ( modulo ) modulo->fold(p2, p1);
        
        Vector dif = p2 - p1;
        
        if ( pd1->visible )
        {
            float R1 = pixscale(pd1->size);
            gym::color_both(pd1->color);
            //gym::stretchAlignZ(p1, p2, R1);
            gym::transAlignZ(p1, R1, dif);
            gle::droplet();
        }
        if ( pd2->visible )
        {
            float R2 = pixscale(pd2->size);
            gym::color_both(pd2->color);
            //gym::stretchAlignZ(p2, p1, R2);
            gym::transAlignZ(p2, -R2, dif);
            gle::droplet();
        }
    }
}
