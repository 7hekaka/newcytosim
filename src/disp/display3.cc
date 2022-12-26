// Cytosim was created by Francois Nedelec. Copyright Cambridge University 2021

#include "dim.h"
#include "simul.h"
#include "display3.h"
#include "modulo.h"
#include "random_pcg.h"

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

extern Modulo const* modulo;

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
        glStencilFunc(GL_EQUAL, 0, ~0);
        
        // set Stencil to 1 for inner surfaces of fibers:
        gym::enableCullFace(GL_FRONT);
        //drawFibers(sim.fibers);
        FiberSet const& set = sim.fibers;
        // display the Fiber always in the same order:
        // set Stencil to 0 for outer surfaces:
        glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);
        GLint val = 0;
        for( Fiber const* fib = set.firstID(); fib; fib=set.nextID(fib) )
        {
            if ( fib->disp->visible )
            {
                glStencilFunc(GL_ALWAYS, ++val, ~0);
                drawFiber(*fib);
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
                glStencilFunc(GL_EQUAL, ++val, ~0);
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
    drawMisc(sim);
}


//------------------------------------------------------------------------------
#pragma mark - Drawing primitives


inline void Display3::drawPoint(Vector const& pos, float rad) const
{
    assert_enabled(GL_LIGHTING);
    drawObject(pos, rad, gle::sphere1);
}


inline void Display3::drawPoint(Vector const& pos, PointDisp const* dis) const
{
    if ( dis->perceptible )
    {
        assert_enabled(GL_LIGHTING);
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
        assert_enabled(GL_LIGHTING);
        gym::color_both(dis->color);
        drawObject(pos, pixscale(dis->size), gle::blob);
    }
}

/// draw a point with a small sphere
inline void Display3::drawHandF(Vector const& pos, PointDisp const* dis) const
{
    if ( dis->perceptible )
    {
        assert_enabled(GL_LIGHTING);
        gym::color_both(dis->color2);
        drawObject(pos, pixscale(dis->size), gle::blob);
    }
}


//------------------------------------------------------------------------------
#pragma mark - Nicely Joined Fiber Rendering using clipping planes

/**
This draws the model-segments, using function `select_color` to set display colors
*/
void Display3::drawFiberSegmentsClip(Fiber const& fib, real rad,
                                     gym_color (*select_color)(Fiber const&, size_t)) const
{
    const size_t last = fib.lastSegment();
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
        gle::drawTube(old, rad, pos, gle::capedTube);
        gym::setClipPlane(4, dir, pos);
        
        // draw inner segments
        gym::enableClipPlane(5);
        for ( size_t i = 1; i < last; ++i )
        {
            old = pos;
            pos = nxt;
            nxt = fib.posPoint(i+2);
            dir = normalize(nxt-old);
            gym::color_front(select_color(fib, i));
            gym::setClipPlane(5, -dir, pos);
            gle::drawTube(old, rad, pos, gle::centralTube);
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
        gle::drawTube(pos, rad, nxt, gle::capedTube);
        gym::setClipPlane(4, dir, old);
    }
    // draw last segment:
    gym::color_front(select_color(fib, last));
    gle::drawTube(nxt, rad, pos, gle::endedTube);
    gym::disableClipPlane(4);
}


/**
This draws segments of length 'len' which may not correspond to the vertices
used to model the Fiber. All abscissa is relative to the MINUS_END.
The function `set_color` is called to set the color of the segments.
*/
void Display3::drawFiberSectionsClip(Fiber const& fib, real rad,
                                     long inx, const long last,
                                     real abs, const real inc,
                                     gym_color (*select_color)(Fiber const&, long, real),
                                     real fac, real facM, real facP) const
{
    Vector old = fib.displayPosM(abs);
    Vector pos = fib.displayPosM(abs+inc);
    Vector nxt = fib.displayPosM(abs+inc*2);
    Vector dir = normalize(nxt-pos);
    
    gym::color_front(select_color(fib, inx++, facM));
    gym::enableClipPlane(4);
    gym::setClipPlane(4, -dir, pos);
    if ( abs <= 0 )
        gle::drawTube(old, rad, pos, gle::capedTube);
    else
        gle::drawTube(old, rad, pos, gle::halfTube2);
    gym::setClipPlane(4, dir, pos);
    
    // keep abs to match to the end of the section already drawn
    abs += inc;

    gym::enableClipPlane(5);
    // draw segments
    while ( inx < last )
    {
        abs += inc;
        old = pos;
        pos = nxt;
        nxt = fib.displayPosM(abs+inc);
        dir = normalize(nxt-old);
        gym::color_front(select_color(fib, inx++, fac));
        gym::setClipPlane(5, -dir, pos);
        // could add a disc to close the tube: gle::endedTube
        gle::drawTube(old, rad, pos, gle::centralTube);
        gym::setClipPlane(4, dir, pos);
    }
    gym::disableClipPlane(5);

    // draw last segment, which may be truncated:
    gym::color_front(select_color(fib, last, facP));
    if ( abs+inc >= fib.length() )
    {
        gym::stretchAlignZ1(nxt, rad, -fib.dirEndP(), fib.length()-abs);
        gle::endedTube();
    }
    else
    {
        gle::drawTube(nxt, rad, pos, gle::halfTube2);
    }
    gym::disableClipPlane(4);
}


//------------------------------------------------------------------------------
#pragma mark - draw Fibers using longer tubes to fill the gaps at the junctions

/**
This draws the model-segments, using function `select_color` to set display colors
*/
void Display3::drawFiberSegmentsJoin(Fiber const& fib, real rad,
                                     gym_color (*select_color)(Fiber const&, size_t)) const
{
    const size_t last = fib.lastSegment();
    Vector pos = fib.posPoint(0);
    Vector nxt = fib.posPoint(1);
    
    gym::color_front(select_color(fib, 0));
    gym::transAlignZ(pos, rad, nxt-pos);
    gle::hemisphere4();
    gym::scale(1, 1, fib.segmentation()/rad);
    if ( last == 0 )
    {
        gle::tube4();
        gle::discTop2();
        return;
    }
    gle::tubeS();

    for ( size_t i = 1; i < last; ++i )
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
used to model the Fiber. All abscissa is relative to the MINUS_END.
The function `set_color` is called to set the color of the segments.
*/
void Display3::drawFiberSectionsJoin(Fiber const& fib, real rad,
                                     long inx, const long last,
                                     real abs, const real inc,
                                     gym_color (*select_color)(Fiber const&, long, real),
                                     real fac, real facM, real facP) const
{
    Vector pos = fib.displayPosM(abs);
    Vector nxt = fib.displayPosM(abs+inc);
    
    gym::color_front(select_color(fib, inx++, facM));
    if ( abs <= 0 )
    {
        real len = (nxt-pos).norm();
        gym::stretchAlignZ1(pos, rad, (nxt-pos)/len, rad);
        gle::hemisphere4();
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

    // keep abs to match to the end of the section already drawn
    abs += inc;
    
    while ( inx < last )
    {
        abs += inc;
        pos = nxt;
        nxt = fib.displayPosM(abs);
        gym::color_front(select_color(fib, inx++, fac));
        gym::stretchAlignZ(pos, nxt, rad);
        gle::tubeM();
    }
    // draw last segment, which may be truncated:
    gym::color_front(select_color(fib, last, facP));
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

//------------------------------------------------------------------------------
#pragma mark -

#if USE_CLIP_PLANES
#  define drawFiberSegments drawFiberSegmentsClip
#  define drawFiberSections drawFiberSectionsClip
#else
#  define drawFiberSegments drawFiberSegmentsJoin
#  define drawFiberSections drawFiberSectionsJoin
#endif

void Display3::drawFiberLines(Fiber const& fib, int style) const
{
    FiberDisp const*const disp = fib.prop->disp;
    const real rad = pixscale(disp->line_width);

    // set back color:
    if ( disp->coloring )
        gym::color_back(fib.disp->color);
    else
        gym::color_back(disp->back_color);
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
            drawFiberSegments(fib, rad, color_seg_curvature);
            break;
        case 5:
            drawFiberSegments(fib, rad, color_by_direction);
            break;
        case 6: {
            /** This is using transparency with segments that are not depth sorted
             but this code is only used in 2D normally, so it's okay */
            gym::disableCullFace();
            drawFiberSegments(fib, rad, color_by_abscissaM);
            gym::restoreCullFace();
        } break;
        case 7: {
            /** This is using transparency with segments that are not depth sorted
             but this code is only used in 2D normally, so it's okay */
            gym::enableCullFace(GL_BACK);
            drawFiberSegments(fib, rad, color_by_abscissaP);
            gym::restoreCullFace();

        } break;
        case 8:
            drawFiberSegments(fib, rad, color_by_height);
            break;
        case 9:
            drawFiberSegments(fib, rad, color_by_grid);
            break;
    }
}


// displays segment 'inx' with transparency
void Display3::drawFiberSegmentT(Fiber const& fib, size_t inx) const
{
    FiberDisp const*const disp = fib.prop->disp;
    const real rad = pixscale(disp->line_width);
    real iseg = fib.segmentationInv();

    Vector A = fib.posP(inx);
    Vector B = fib.posP(inx+1);

    gym::enableLighting();
    /* Either CULL_FACE should be enable to hide the back side,
     or every primitive should be renderred with a double pass*/
    gym::enableCullFace(GL_BACK);

    if ( disp->line_style == 6 )
        gym::color_front(color_by_abscissaM(fib, inx));
    else if ( disp->line_style == 7 )
        gym::color_front(color_by_abscissaP(fib, inx));
    else if ( disp->line_style == 2 )
        gym::color_front(color_by_tension(fib, inx));
    else if ( disp->line_style == 3 )
        gym::color_front(color_by_tension_jet(fib, inx));
    else
        gym::color_front(fib.disp->color);
    
    // truncate terminal segment according to length_scale
    if ( inx == 0 && disp->line_style == 6 )
    {
        real x = 3 * disp->length_scale * iseg;
        if ( x < 1.0 )
        {
            B = A + x * ( B - A );
            color_by_abscissaM(fib, inx);
        }
    }
    
    // truncate terminal segment according to length_scale
    if ( inx == fib.lastSegment() && disp->line_style == 7 )
    {
        real x = 3 * disp->length_scale * iseg;
        if ( x < 1.0 )
        {
            A = B + x * ( A - B );
            color_by_abscissaP(fib, inx);
        }
    }

#if USE_CLIP_PLANES
    if ( inx == 0 )
    {
        gym::enableClipPlane(5);
        if ( inx == fib.lastSegment() )
        {
            gym::setClipPlane(5, (A-B)*iseg, (A+B)*0.5);
            gle::drawTube(A, rad, B, gle::capedTube);
            gym::setClipPlane(5, (B-A)*iseg, (A+B)*0.5);
            gle::drawTube(B, rad, A, gle::endedTube);
        }
        else
        {
            gym::setClipPlane(5, (A-B)*iseg, B);
            gle::drawTube(A, rad, B, gle::capedTube);
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
        gle::drawTube(B, rad, A, gle::endedTube);
        gym::disableClipPlane(4);
        return;
    }
    else
    {
        gym::setClipPlane(4, normalize(A-fib.posP(inx+2)), B);
    }
    
    gym::enableClipPlane(5);
    gym::enableClipPlane(4);
    gle::drawTube(A, rad, B, gle::centralTube);
    gym::disableClipPlane(4);
    gym::disableClipPlane(5);
#else
    gym::stretchAlignZ(A, B, rad);
    gle::tube4();
    if ( inx == 0 )
        gle::hemisphere();
    if ( inx == fib.lastSegment() )
        gle::discTop2();
#endif
    gym::restoreCullFace();
}


//------------------------------------------------------------------------------
#pragma mark - Display Lattice


void Display3::drawFiberLattice(Fiber const& fib, VisibleLattice const& lat, real rad,
                                gym_color (*select_color)(Fiber const&, long, real)) const
{
    FiberDisp const*const disp = fib.prop->disp;

    gym::enableLighting();
    GLfloat blk[] = { 0.f, 0.f, 0.f, 1.f };
    GLfloat bak[] = { 0.f, 0.f, 0.f, 1.f };
    disp->back_color.store(bak);
    glMaterialfv(GL_BACK, GL_AMBIENT,  blk);
    glMaterialfv(GL_BACK, GL_DIFFUSE,  bak);
    glMaterialfv(GL_FRONT, GL_AMBIENT,  blk);
    glMaterialfv(GL_FRONT, GL_DIFFUSE,  blk);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, blk);
    glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, blk);
    glMateriali (GL_FRONT_AND_BACK, GL_SHININESS, 32);

    const real fac = 1 / disp->lattice_scale;
    const real uni = lat.unit();
    const auto inf = lat.indexM();
    const auto sup = lat.indexP();

    real abs = uni * inf - fib.abscissaM();       // should be non-positive!
    assert_true( abs <= 0 && 0 <= abs+uni );

    real facM = fac;
    real facP = fac;
    
    if ( disp->lattice_rescale )
    {
        real lenM = abs + uni;                     // should be positive!
        real lenP = fib.abscissaP() - uni * sup;   // should be positive!
        facM = ( lenM > 0.001*uni ? fac*uni/lenM : fac );
        facP = ( lenP > 0.001*uni ? fac*uni/lenP : fac );
    }

    drawFiberSections(fib, rad, inf, sup, abs, uni, select_color, fac, facM, facP);
}


void Display3::drawFiberLattice1(Fiber const& fib, VisibleLattice const& lat, real rad) const
{
    drawFiberLattice(fib, lat, rad, color_by_lattice);
}

void Display3::drawFiberLattice2(Fiber const& fib, VisibleLattice const& lat, real rad) const
{
    drawFiberLattice(fib, lat, rad, color_by_lattice_jet);
}

void Display3::drawFiberLattice3(Fiber const& fib, VisibleLattice const& lat, real rad) const
{
    drawFiberLattice(fib, lat, rad, color_by_lattice_white);
}

void Display3::drawFiberLatticeEdges(Fiber const& fib, VisibleLattice const& lat, real rad) const
{
    drawFiberLattice(fib, lat, rad, color_alternate);
}

//------------------------------------------------------------------------------
#pragma mark -


/**
 Display the MINUS_END of a Fiber, according to `style`:
 - 1: draw a sphere
 - 2: draw a cone
 - 3: draw a flat cylinder
 - 4: draw an arrow-head
 - 5: arrow-head in reverse direction
 .
 with 3D objects
 */
void Display3::drawFiberMinusEnd(Fiber const& fib, int style, float size) const
{
    const float rad = pixscale(size);
    if ( rad > 0 ) switch(style)
    {
        default: break;
        case 1: drawObject(fib.posEndM(), rad, gle::sphere2); break;
        case 2: drawObject(fib.posEndM(),-fib.dirEndM(), rad, gle::longCone); break;
        case 3: drawObject(fib.posEndM(),-fib.dirEndM(), rad, gle::cylinder2); break;
        case 4: drawObject(fib.posEndM(),-fib.dirEndM(), rad, gle::arrowTail); break;
        case 5: drawObject(fib.posEndM(), fib.dirEndM(), rad, gle::arrowTail); break;
        case 6: drawObject(fib.posEndM(),-fib.dirEndM(), rad, gle::cube); break;
        case 7: drawObject(fib.posEndM(), fib.dirEndM(), rad, gle::cylinder2); break;
    }
}


/**
 Display the PLUS_END of a Fiber, according to `style`:
 - 1: draw a sphere
 - 2: draw a cone
 - 3: draw a flat cylinder
 - 4: draw an arrow-head
 - 5: arrow-head in reverse direction
 .
 with 3D objects
 */
void Display3::drawFiberPlusEnd(Fiber const& fib, int style, float size) const
{
    const float rad = pixscale(size);
    if ( rad > 0 ) switch(style)
    {
        default: break;
        case 1: drawObject(fib.posEndP(), rad, gle::sphere2); break;
        case 2: drawObject(fib.posEndP(), fib.dirEndP(), rad, gle::longCone); break;
        case 3: drawObject(fib.posEndP(), fib.dirEndP(), rad, gle::cylinder2); break;
        case 4: drawObject(fib.posEndP(), fib.dirEndP(), rad, gle::arrowTail); break;
        case 5: drawObject(fib.posEndP(),-fib.dirEndP(), rad, gle::arrowTail); break;
        case 6: drawObject(fib.posEndP(), fib.dirEndP(), rad, gle::cube); break;
        case 7: drawObject(fib.posEndP(),-fib.dirEndP(), rad, gle::cylinder2); break;
    }
}


void Display3::drawFiberSpeckles(Fiber const& fib) const
{
    FiberDisp const*const disp = fib.prop->disp;
    const float rad = pixscale(disp->speckle_size);
    gym::color_front(fib.disp->color);
    gym::color_back(disp->back_color);

    // display random speckles:
    if ( disp->speckle_style == 1 )
    {
        /*
         A simple random number generator seeded by fib.signature()
         is used to distribute points always at the same position
         with respect to the lattice of each fiber.
         */
        
        const real gap = disp->speckle_gap;
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
                drawPoint(fib.pos(a), rad);
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
                drawPoint(fib.pos(a), rad);
                a -= gap * std::log(pcg32(Z)*TINY);
            }
        }
    }
    else if ( disp->speckle_style == 2 )
    {
        //we distribute points regularly along the center line
        const real gap = disp->speckle_gap;
        real a = gap * std::ceil( fib.abscissaM() / gap );
        while ( a <= fib.abscissaP() )
        {
            drawPoint(fib.pos(a), rad);
            a += gap;
        }
    }
}


void Display3::drawFiberPoints(Fiber const& fib) const
{
    FiberDisp const*const disp = fib.prop->disp;
    // diameter of lines and points in space units:
    const float rad = pixscale(disp->point_size);
    int style = disp->point_style;
    gym::color_front(fib.disp->color);
    gym::color_back(disp->back_color);

    if ( style == 1 )
    {
        gym::enableLighting();
        // display vertices:
        for ( size_t i = 0; i < fib.nbPoints(); ++i )
            drawObject(fib.posP(i), rad, gle::cube);
    }
    else if ( style == 2 )
    {
        gym::enableLighting();
        // display arrowheads along the fiber:
        const real gap = disp->point_gap;
        real ab = std::ceil(fib.abscissaM()/gap) * gap;
        for ( ; ab <= fib.abscissaP(); ab += gap )
            gle::drawCone(fib.pos(ab), fib.dir(ab), rad);
    }
    else if ( style == 3 )
    {
        gym::enableLighting();
        // display chevrons along the fiber:
        const real gap = disp->point_gap;
        real ab = std::ceil(fib.abscissaM()/gap) * gap;
        for ( ; ab <= fib.abscissaP(); ab += gap )
            gle::drawCone(fib.pos(ab), -fib.dir(ab), rad);
    }
    else if ( style == 4 )
    {
        // display middle of fiber:
        drawPoint(fib.posMiddle(), 2*rad);
    }
}

//------------------------------------------------------------------------------

void Display3::drawOrganizer(Organizer const& obj) const
{
    PointDisp const* disp = obj.disp();

    if ( disp && ( disp->style & 2 ))
    {
        Vector P, Q;
        bodyColor(disp, obj.signature());
        const float wid = pixscale(disp->width);

        for ( size_t i = 0; obj.getLink(i, P, Q); ++i )
        {
            drawPoint(P, disp);
            if ( modulo ) modulo->fold(Q, P);
            gle::stretchTube(P, wid, Q, gle::tube1);
        }
    }
    /**
     This displays the Solid connecting two Aster as a spindle.
     Used for Cleo Kozlowski simulation of C. elegans (2007)
     */
    if ( disp && ( disp->style & 1 ) && obj.tag() == Organizer::TAG_FAKE )
    {
        Solid const* sol = Solid::toSolid(obj.core());
        if ( sol && sol->nbPoints() >= 4 )
        {
            bodyColor(*sol);
#if ( DIM >= 3 )
            Vector3 a = 0.5 * (sol->posP(0) + sol->posP(2));
            Vector3 b = 0.5 * (sol->posP(1) + sol->posP(3));
            gym::stretchAlignZ(a, b, 1);
            gym::color(0.6f,0.6f,0.6f);
            gle::dualPassBarrel();
#else
            const float wid = pixscale(disp->width);
            for ( size_t i = 0; i < sol->nbPoints(); i+=2 )
                gle::stretchTube(sol->posPoint(i), wid, sol->posPoint(i+1), gle::hexTube);
#endif
        }
    }
}

//------------------------------------------------------------------------------
#pragma mark -

void Display3::drawSinglesF(SingleSet const& set) const
{
    for ( Single * obj=set.firstF(); obj ; obj=obj->next() )
    {
        if ( obj->disp()->perceptible )
        {
            gym::color_both(obj->disp()->color2);
#if ( 0 )
            if ( obj->disp()->style == 2 )
            {
                Space const* spc = obj->confineSpace();
                if ( spc )
                {
                    const float rad = pixscale(obj->disp()->size);
                    /// draw a disc tangent to the Space:
                    Vector pos = obj->posFoot();
                    Vector dir = spc->normalToEdge(pos);
                    drawObject(pos, dir, rad, gle::disc);
                    continue;
                }
            }
#endif
            drawHandF(obj->posFoot(), obj->disp());
        }
    }
}


void Display3::drawSingleA(Single const* obj) const
{
    const PointDisp * disp = obj->disp();
    Vector ph = obj->posHand();
    gym::color_both(disp->color);
    drawHand(ph, disp);
}


void Display3::drawSingleB(Single const* obj) const
{
    const PointDisp * disp = obj->disp();
    Vector ph = obj->posHand();
    Vector pf = obj->posFoot();
    if ( modulo ) modulo->fold(pf, ph);
    const float wid = pixscale(disp->width);
    const float rad = pixscale(disp->size);

    gym::color_both(disp->color2);
#if ( 0 )
    if ( obj->disp()->style == 2 && obj->confineSpace() )
    {
        // draw a disc tangent to the Space:
        drawObject(pf, obj->confineSpace()->normalToEdge(pf), rad, gle::disc);
    }
    else
#endif
    {
        gym::transScale(pf, wid);
        gle::blob(); // the foot
    }
    gym::color_both(disp->color);
#if ( DIM > 2 )
    Vector diff = pf - ph;
    float L = norm(diff);
    gym::transAlignZ(ph, rad, diff/L);
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
        drawHand(ph, disp);
    gle::drawBand(ph, wid, disp->color, pf, wid, disp->color.alpha_scaled(0.5f));
#endif
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
    for ( Couple * cx = set.firstFF(); cx ; cx=cx->next() )
        drawHandF(cx->posFree(), cx->disp1());
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
        PointDisp const* disp = h->prop->disp;

        if ( h->fiber()->disp->visible && disp->visible )
        {
#if ( 0 )  // ENDOCYTOSIS 2015
            if ( cx->fiber1()->disp->color.transparent() )
            {
                gym::color_both(disp->color, cx->fiber1()->disp->color.transparency());
                drawPoint(h->pos(), disp);
                continue;
            }
#endif
            drawHand(h->pos(), disp);
            //drawObject(h->outerPos(), h->pos()-h->outerPos(), R*disp->size, gle::tetrahedron);
        }
    }
    
    for ( Couple * cx=set.firstFA(); cx ; cx=cx->next() )
    {
        Hand const* h = cx->hand2();
        PointDisp const* disp = h->prop->disp;

        if ( cx->fiber2()->disp->visible && cx->disp2()->visible )
        {
#if ( 0 )  // ENDOCYTOSIS 2015
            if ( cx->fiber2()->disp->color.transparent() )
            {
                gym::color_both(disp->color, cx->fiber2()->disp->color.transparency());
                drawPoint(h->pos(), disp);
                continue;
            }
#endif
            drawHand(h->pos(), disp);
            //drawObject(h->outerPos(), h->pos()-h->outerPos(), R*disp->size, gle::tetrahedron);
        }
    }
}

void Display3::drawCoupleBplain(Couple const* cx) const
{
    const PointDisp * pd1 = cx->disp1();
    const PointDisp * pd2 = cx->disp2();
    
    Vector p1 = cx->posHand1();
    Vector p2 = cx->posHand2();

    gym::color_both(pd1->color);
    gle::stretchTube(p1, pixscale(pd1->width), p2, gle::hexTube);
    if ( pd1->visible ) drawHand(p1, pd1);
    if ( pd2->visible ) drawHand(p2, pd2);
}


void Display3::drawCoupleBside(Couple const* cx) const
{
    const PointDisp * pd1 = cx->disp1();
    const PointDisp * pd2 = cx->disp2();
    
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
        drawPoint(mid, pixscale(pd1->width));
        gle::stretchTube(p2, pixscale(pd2->width), mid, gle::hexTube);
        gle::stretchTube(p1, pixscale(pd1->width), mid, gle::hexTube);
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
        gym::translate(0, 0,-Lr);
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



void Display3::drawCoupleB(Couple const* cx) const
{
    const PointDisp * pd1 = cx->disp1();
    const PointDisp * pd2 = cx->disp2();
    
    Vector p1 = cx->posHand1();
    Vector p2 = cx->posHand2();
    if ( modulo ) modulo->fold(p2, p1);
    
    Vector dif = p2 - p1;
    real dns = dif.normSqr();
    
    if ( dns > 1e-6 )
    {
#if !FIBER_HAS_FAMILY
        // moving the 'hands' to the surface of the fiber:
        dns = sizeScale / std::sqrt(dns);
        // position the heads at the surface of the filaments:
        const real rad1 = cx->fiber1()->prop->disp->line_width + 0.4 * pd1->size;
        const real rad2 = cx->fiber2()->prop->disp->line_width + 0.4 * pd2->size;
        // move points along the link
        //p1 += dif * min_real(0.45, rad1*dns);
        //p2 -= dif * min_real(0.45, rad2*dns);
        // move points orthogonal to the fiber's axis
        Vector dir1 = cx->dirFiber1();
        Vector dir2 = cx->dirFiber2();
        p1 += ( dif - dot(dif,dir1) * dir1 ) * min_real(0.45, rad1*dns);
        p2 -= ( dif - dot(dif,dir2) * dir2 ) * min_real(0.45, rad2*dns);
#endif
        if ( pd1->visible )
        {
            real R1 = pixscale(pd1->size);
            gym::color_both(pd1->color);
            //gym::stretchAlignZ(p1, p2, R1);
            gym::transAlignZ(p1, R1, p2-p1);
            gle::droplet();
        }
        if ( pd2->visible )
        {
            real R2 = pixscale(pd2->size);
            gym::color_both(pd2->color);
            //gym::stretchAlignZ(p2, p1, R2);
            gym::transAlignZ(p2, R2, p1-p2);
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
    const PointDisp * pd1 = cx->disp1();
    const PointDisp * pd2 = cx->disp2();
    
    Vector p1 = cx->posHand1();
    Vector p2 = cx->posHand2();
    if ( modulo ) modulo->fold(p2, p1);
    
#if !FIBER_HAS_FAMILY
    Vector dif = p2 - p1;
    real dns = dif.normSqr();
    if ( dns > 1e-6 )
    {
        // moving the 'hands' to the surface of the fiber:
        dns = sizeScale / std::sqrt(dns);
        // position the heads at the surface of the filaments:
        const real rad1 = cx->fiber1()->prop->disp->line_width + 0.4 * pd1->size;
        const real rad2 = cx->fiber2()->prop->disp->line_width + 0.4 * pd2->size;
        // move points along the link
        //p1 += dif * min_real(0.45, rad1*dns);
        //p2 -= dif * min_real(0.45, rad2*dns);
        // move points orthogonal to the fiber's axis
        Vector dir1 = cx->dirFiber1();
        Vector dir2 = cx->dirFiber2();
        p1 += ( dif - dot(dif,dir1) * dir1 ) * min_real(0.45, rad1*dns);
        p2 -= ( dif - dot(dif,dir2) * dir2 ) * min_real(0.45, rad2*dns);
    }
#endif

#if ( 0 ) // ENDOCYTOSIS 2015
    if ( cx->fiber1()->disp->color.transparent() )
    {
        gym::color_both(pd1->color, cx->fiber1()->disp->color.transparency());
        gym::closeDepthMask();
        gle::stretchTube(p1, pixscale(pd2->width), p2, gle::hexTube);
        gym::openDepthMask();
        return;
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

