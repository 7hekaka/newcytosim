// Cytosim was created by Francois Nedelec. Copyright Cambridge University 2021

#include "dim.h"
#include "simul.h"
#include "display3.h"
#include "modulo.h"

#include "gle.h"
#include "gle_color_list.h"
#include "glut.h"

#include "fake.h"
#include "fiber_disp.h"
#include "line_disp.h"
#include "point_disp.h"
#include "display_color.h"

/// Cliping planes permit nice junctions to be made between tubes, but are slow
#define USE_CLIP_PLANES 0

using namespace gle;
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
    glDepthMask(GL_FALSE);
    glDisable(GL_LIGHTING);
    glDisable(GL_CULL_FACE);
    drawFields(sim.fields);
    
    glEnable(GL_LIGHTING);
    glEnable(GL_CULL_FACE);
    glDepthMask(GL_TRUE);
    drawSpaces(sim.spaces);
    
    glDisable(GL_CULL_FACE);

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
        glEnable(GL_CULL_FACE);
        //drawFibers(sim.fibers);
        FiberSet const& set = sim.fibers;
        // display the Fiber always in the same order:
        // set Stencil to 0 for outer surfaces:
        glCullFace(GL_FRONT);
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
        glCullFace(GL_BACK);
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
        //glEnable(GL_CULL_FACE);
        //glCullFace(GL_BACK);
        drawFibers(sim.fibers);
    }
    
    glEnable(GL_LIGHTING);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    //gym::bindBuffer();

    drawBeads(sim.beads);
    drawSolids(sim.solids);
    drawSpheres(sim.spheres);
    
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

    //gle::unbindBuffer();

    if ( stencil_ )
    {
        glClearStencil(0);
        glDisable(GL_STENCIL_TEST);
    }

    drawOrganizers(sim.organizers);
    glDisable(GL_CULL_FACE);
    drawMisc(sim);
}


//------------------------------------------------------------------------------
#pragma mark - Drawing primitives


inline void Display3::drawPoint(Vector const& pos, float rad) const
{
    assert_true(glIsEnabled(GL_LIGHTING));
    drawObject(pos, rad, gle::sphere1);
}


inline void Display3::drawPoint(Vector const& pos, PointDisp const* dis) const
{
    if ( dis->perceptible )
    {
        assert_true(glIsEnabled(GL_LIGHTING));
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
        assert_true(glIsEnabled(GL_LIGHTING));
        dis->color.load_both();
        drawObject(pos, pixscale(dis->size), gle::blob);
    }
}

/// draw a point with a small sphere
inline void Display3::drawHandF(Vector const& pos, PointDisp const* dis) const
{
    if ( dis->perceptible )
    {
        assert_true(glIsEnabled(GL_LIGHTING));
        dis->color2.load_both();
        drawObject(pos, pixscale(dis->size), gle::blob);
    }
}


//------------------------------------------------------------------------------
#pragma mark - Nicely Joined Fiber Rendering using clipping planes

#if USE_CLIP_PLANES
/**
This draws the model-segments, using function `select_color` to set display colors
*/
void Display3::drawFiberSegments(Fiber const& fib, real rad,
                                 gle_color (*select_color)(Fiber const&, size_t)) const
{
    const size_t last = fib.lastSegment();
    Vector old = fib.posPoint(0);
    Vector pos = fib.posPoint(1);
    Vector nxt;
    Vector dir = (pos-old) / fib.segmentation();
    select_color(fib, 0).load_front();

    if ( last > 0 )
    {
        nxt = fib.posPoint(2);
        glEnable(GL_CLIP_PLANE4);
        setClipPlane(GL_CLIP_PLANE4, -dir, nxt);
        drawTube(old, rad, pos, gle::capedTube2);
        setClipPlane(GL_CLIP_PLANE4, dir, pos);
        
        // draw inner segments
        glEnable(GL_CLIP_PLANE5);
        for ( size_t i = 1; i < last; ++i )
        {
            old = pos;
            pos = nxt;
            nxt = fib.posPoint(i+2);
            dir = normalize(nxt-old);
            select_color(fib, i).load_front();
            setClipPlane(GL_CLIP_PLANE5, -dir, pos);
            drawTube(old, rad, pos, gle::longTube2);
            setClipPlane(GL_CLIP_PLANE4,  dir, pos);
        }
        glDisable(GL_CLIP_PLANE5);
    }
    else
    {
        nxt = pos;
        pos = old;
        old = 0.5 * ( nxt + pos );
        glEnable(GL_CLIP_PLANE4);
        setClipPlane(GL_CLIP_PLANE4, -dir, old);
        drawTube(pos, rad, nxt, gle::capedTube2);
        setClipPlane(GL_CLIP_PLANE4, dir, old);
    }
    // draw last segment:
    select_color(fib, last).load_front();
    drawTube(nxt, rad, pos, gle::endedTube2);
    glDisable(GL_CLIP_PLANE4);
}


/**
This draws segments of length 'len' which may not correspond to the vertices
used to model the Fiber. All abscissa is relative to the MINUS_END.
The function `set_color` is called to set the color of the segments.
*/
void Display3::drawFiberSubSegments(Fiber const& fib, real rad,
                                    long inx, const long last,
                                    real abs, const real inc,
                                    gle_color (*select_color)(Fiber const&, long, real),
                                    real fac, real facM, real facP) const
{
    Vector old = fib.displayPosM(abs);
    Vector pos = fib.displayPosM(abs+inc);
    Vector nxt = fib.displayPosM(abs+inc*2);
    Vector dir = normalize(nxt-pos);
    
    select_color(fib, inx++, facM).load_front();
    glEnable(GL_CLIP_PLANE4);
    setClipPlane(GL_CLIP_PLANE4, -dir, pos);
    if ( abs <= 0 )
        drawTube(old, rad, pos, gle::capedTube2);
    else
        drawTube(old, rad, pos, gle::halfTube2);
    setClipPlane(GL_CLIP_PLANE4, dir, pos);
    
    // keep abs to match to the end of the section already drawn
    abs += inc;

    glEnable(GL_CLIP_PLANE5);
    // draw segments
    while ( inx < last )
    {
        abs += inc;
        old = pos;
        pos = nxt;
        nxt = fib.displayPosM(abs+inc);
        dir = normalize(nxt-old);
        select_color(fib, inx++, fac).load_front();
        setClipPlane(GL_CLIP_PLANE5, -dir, pos);
        // could add a disc to close the tube: gle::endedTube2
        drawTube(old, rad, pos, gle::longTube2);
        setClipPlane(GL_CLIP_PLANE4, dir, pos);
    }
    glDisable(GL_CLIP_PLANE5);

    // draw last segment, which may be truncated:
    select_color(fib, last, facP).load_front();
    if ( abs+inc >= fib.length() )
    {
        glPushMatrix();
        gle::stretchAlignZ1(nxt, rad, -fib.dirEndP(), fib.length()-abs);
        gle::endedTube2();
        glPopMatrix();
    }
    else
    {
        drawTube(nxt, rad, pos, gle::halfTube2);
    }
    glDisable(GL_CLIP_PLANE4);
}

#else

//------------------------------------------------------------------------------
#pragma mark - draw Fibers using longer tubes to fill the gaps at the junctions

/**
This draws the model-segments, using function `select_color` to set display colors
*/
void Display3::drawFiberSegments(Fiber const& fib, real rad,
                                 gle_color (*select_color)(Fiber const&, size_t)) const
{
    const size_t last = fib.lastSegment();
    Vector pos = fib.posPoint(0);
    Vector nxt = fib.posPoint(1);
    
    select_color(fib, 0).load_front();
    glPushMatrix();
    gle::transAlignZ(pos, rad, nxt-pos);
    gle::hemisphere4();
    glScalef(1, 1, fib.segmentation()/rad);
    if ( last == 0 )
    {
        gle::tube4();
        gle::discTop2();
        glPopMatrix();
        return;
    }
    gle::tubeS();
    glPopMatrix();

    for ( size_t i = 0; i < last; ++i )
    {
        pos = nxt;
        nxt = fib.posPoint(i+1);
        select_color(fib, i).load_front();
        glPushMatrix();
        gle::stretchAlignZ(pos, nxt, rad);
        gle::tubeM();
        glPopMatrix();
    }
    glPushMatrix();
    pos = nxt;
    nxt = fib.posPoint(last+1);
    select_color(fib, last).load_front();
    gle::stretchAlignZ(pos, nxt, rad);
    gle::tubeE();
    gle::discTop2();
    glPopMatrix();
}


/**
This draws segments of length 'len' which may not correspond to the vertices
used to model the Fiber. All abscissa is relative to the MINUS_END.
The function `set_color` is called to set the color of the segments.
*/
void Display3::drawFiberSubSegments(Fiber const& fib, real rad,
                                    long inx, const long last,
                                    real abs, const real inc,
                                    gle_color (*select_color)(Fiber const&, long, real),
                                    real fac, real facM, real facP) const
{
    Vector pos = fib.displayPosM(abs);
    Vector nxt = fib.displayPosM(abs+inc);
    
    select_color(fib, inx++, facM).load_front();
    glPushMatrix();
    if ( abs <= 0 )
    {
        real len = (nxt-pos).norm();
        gle::stretchAlignZ1(pos, rad, (nxt-pos)/len, rad);
        gle::hemisphere4();
        glScalef(1, 1, len/rad);
    }
    else
    {
        gle::stretchAlignZ(pos, nxt, rad);
    }
    if ( last == inx )
    {
        gle::tube4();
        gle::discTop2();
        glPopMatrix();
        return;
    }
    gle::tubeS();
    glPopMatrix();
    
    // keep abs to match to the end of the section already drawn
    abs += inc;
    
    while ( inx < last )
    {
        abs += inc;
        pos = nxt;
        nxt = fib.displayPosM(abs);
        select_color(fib, inx++, fac).load_front();
        glPushMatrix();
        gle::stretchAlignZ(pos, nxt, rad);
        gle::tubeM();
        glPopMatrix();
    }
    // draw last segment, which may be truncated:
    select_color(fib, last, facP).load_front();
    glPushMatrix();
    if ( abs+inc >= fib.length() )
    {
        gle::stretchAlignZ1(nxt, rad, fib.dirEndP(), fib.length()-abs);
        gle::discTop2();
    }
    else
    {
        gle::stretchAlignZ(nxt, fib.displayPosM(abs+inc), rad);
    }
    gle::tubeE();
    glPopMatrix();
}

#endif

//------------------------------------------------------------------------------
#pragma mark -

void Display3::drawFiberLines(Fiber const& fib, int style) const
{
    FiberDisp const*const disp = fib.prop->disp;
    const real rad = pixscale(disp->line_width);

    // set back color:
    if ( disp->coloring )
        fib.disp->color.load_back();
    else
        disp->back_color.load_back();
    glEnable(GL_LIGHTING);
    
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
            GLboolean cull = glIsEnabled(GL_CULL_FACE);
            glEnable(GL_CULL_FACE);
            drawFiberSegments(fib, rad, color_by_abscissaM);
            if ( !cull ) glDisable(GL_CULL_FACE);
        } break;
        case 7: {
            /** This is using transparency with segments that are not depth sorted
             but this code is only used in 2D normally, so it's okay */
            GLboolean cull = glIsEnabled(GL_CULL_FACE);
            glEnable(GL_CULL_FACE);
            glCullFace(GL_BACK);
            drawFiberSegments(fib, rad, color_by_abscissaP);
            if ( !cull ) glDisable(GL_CULL_FACE);
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

    glEnable(GL_LIGHTING);
    /* Either CULL_FACE should be enable to hide the back side,
     or every primitive should be renderred with a double pass*/
    GLboolean cull = glIsEnabled(GL_CULL_FACE);
    glEnable(GL_CULL_FACE);
    
    if ( disp->line_style == 6 )
        color_by_abscissaM(fib, inx).load_front();
    else if ( disp->line_style == 7 )
        color_by_abscissaP(fib, inx).load_front();
    else if ( disp->line_style == 2 )
        color_by_tension(fib, inx).load_front();
    else if ( disp->line_style == 3 )
        color_by_tension_jet(fib, inx).load_front();
    else
        fib.disp->color.load_front();
    
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
        glEnable(GL_CLIP_PLANE5);
        if ( inx == fib.lastSegment() )
        {
            setClipPlane(GL_CLIP_PLANE5, (A-B)*iseg, (A+B)*0.5);
            drawTube(A, rad, B, gle::capedTube2);
            setClipPlane(GL_CLIP_PLANE5, (B-A)*iseg, (A+B)*0.5);
            drawTube(B, rad, A, gle::endedTube2);
        }
        else
        {
            setClipPlane(GL_CLIP_PLANE5, (A-B)*iseg, B);
            drawTube(A, rad, B, gle::capedTube2);
        }
        glDisable(GL_CLIP_PLANE5);
        return;
    }
    else
    {
        setClipPlane(GL_CLIP_PLANE5, normalize(B-fib.posP(inx-1)), A);
    }
    
    if ( inx == fib.lastSegment() )
    {
        glEnable(GL_CLIP_PLANE4);
        setClipPlane(GL_CLIP_PLANE4, (B-A)*iseg, A);
        drawTube(B, rad, A, gle::endedTube2);
        glDisable(GL_CLIP_PLANE4);
        return;
    }
    else
    {
        setClipPlane(GL_CLIP_PLANE4, normalize(A-fib.posP(inx+2)), B);
    }
    
    glEnable(GL_CLIP_PLANE5);
    glEnable(GL_CLIP_PLANE4);
    drawTube(A, rad, B, gle::longTube2);
    glDisable(GL_CLIP_PLANE4);
    glDisable(GL_CLIP_PLANE5);
#else
    glPushMatrix();
    gle::stretchAlignZ(A, B, rad);
    gle::tube4();
    if ( inx == 0 )
        gle::hemisphere4();
    if ( inx == fib.lastSegment() )
        gle::discTop2();
    glPopMatrix();
#endif
    if ( !cull ) glDisable(GL_CULL_FACE);
}


//------------------------------------------------------------------------------
#pragma mark - Display Lattice


void Display3::drawFiberLattice(Fiber const& fib, VisibleLattice const& lat, real width,
                                gle_color (*select_color)(Fiber const&, long, real)) const
{
    FiberDisp const*const disp = fib.prop->disp;

    glEnable(GL_LIGHTING);
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
    const real rad = pixscale(width);
    
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

    drawFiberSubSegments(fib, rad, inf, sup, abs, uni, select_color, fac, facM, facP);
}


void Display3::drawFiberLattice1(Fiber const& fib, VisibleLattice const& lat, real width) const
{
    drawFiberLattice(fib, lat, width, color_by_lattice);
}

void Display3::drawFiberLattice2(Fiber const& fib, VisibleLattice const& lat, real width) const
{
    drawFiberLattice(fib, lat, width, color_by_lattice_jet);
}

void Display3::drawFiberLattice3(Fiber const& fib, VisibleLattice const& lat, real width) const
{
    drawFiberLattice(fib, lat, width, color_by_lattice_striped);
}

void Display3::drawFiberLatticeEdges(Fiber const& fib, VisibleLattice const& lat, real width) const
{
    drawFiberLattice(fib, lat, width, color_alternate);
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

    // display random speckles:
    if ( disp->speckle_style == 1 )
    {
        /*
         A simple random number generator seeded by fib.signature()
         is used to distribute points always at the same position
         with respect to the lattice of each fiber.
         */
        
        const real spread = disp->speckle_gap;
        constexpr real TINY = 0x1p-32;
        // draw speckles below the origin of abscissa:
        if ( fib.abscissaM() < 0 )
        {
            uint32_t z = fib.signature();
            real a = spread * std::log(z*TINY);
            while ( a > fib.abscissaP() )
            {
                z = lcrng2(z);
                a += spread * std::log(z*TINY);
            }
            while ( a >= fib.abscissaM() )
            {
                drawPoint(fib.pos(a), rad);
                z = lcrng2(z);
                a += spread * std::log(z*TINY);
            }
        }
        // draw speckles above the origin of abscissa:
        if ( fib.abscissaP() > 0 )
        {
            uint32_t z = ~fib.signature();
            real a = -spread * std::log(z*TINY);
            while ( a < fib.abscissaM() )
            {
                z = lcrng1(z);
                a -= spread * std::log(z*TINY);
            }
            while ( a <= fib.abscissaP() )
            {
                drawPoint(fib.pos(a), rad);
                z = lcrng1(z);
                a -= spread * std::log(z*TINY);
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
    int style = disp->point_style & 3;

    if ( style == 1 )
    {
        // display vertices:
        for ( size_t i = 0; i < fib.nbPoints(); ++i )
            drawObject(fib.posP(i), rad, gle::cube);
    }
    else if ( style == 2 )
    {
        glEnable(GL_LIGHTING);
        // display arrowheads along the fiber:
        const real gap = disp->point_gap;
        real ab = std::ceil(fib.abscissaM()/gap) * gap;
        for ( ; ab <= fib.abscissaP(); ab += gap )
            gle::drawCone(fib.pos(ab), fib.dir(ab), rad);
    }
    else if ( style == 3 )
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
            stretchTube(P, Q, wid, gle::tube1);
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
            glPushMatrix();
            Vector3 a = 0.5 * (sol->posP(0) + sol->posP(2));
            Vector3 b = 0.5 * (sol->posP(1) + sol->posP(3));
            stretchAlignZ(a, b, 1);
            glColor3f(0.6f,0.6f,0.6f);
            gle::dualPass(gle::barrel);
            glPopMatrix();
#else
            const float wid = pixscale(disp->width);
            for ( size_t i = 0; i < sol->nbPoints(); i+=2 )
                stretchTube(sol->posPoint(i), sol->posPoint(i+1), wid, gle::hexTube);
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
            obj->disp()->color2.load_both();
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
    disp->color.load_both();
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

    disp->color2.load_both();
#if ( 0 )
    if ( obj->disp()->style == 2 && obj->confineSpace() )
    {
        // draw a disc tangent to the Space:
        drawObject(pf, obj->confineSpace()->normalToEdge(pf), rad, gle::disc);
    }
    else
#endif
    {
        glPushMatrix();
        gle::transScale(pf, wid);
        gle::blob(); // the foot
        glPopMatrix();
    }
    disp->color.load_both();
#if ( DIM > 2 )
    Vector diff = pf - ph;
    float L = norm(diff);
    glPushMatrix();
    transAlignZ(ph, rad, diff/L);
    gle::blob();
    glScalef(wid/rad, wid/rad, L/rad);
    gle::hexTube();
    glPopMatrix();
#elif ( 0 )
    Vector dir = normalize( pf - ph );
    glEnable(GL_CLIP_PLANE5);
    setClipPlane(GL_CLIP_PLANE5, -dir, pf);
    glPushMatrix();
    transAlignZ(ph, rad, dir);
    gle::needle();
    glPopMatrix();
    glDisable(GL_CLIP_PLANE5);
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
#pragma mark -
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
                disp->color.load_both(cx->fiber1()->disp->color.transparency());
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
                disp->color.load_both(cx->fiber2()->disp->color.transparency());
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

    pd1->color.load_both();
    stretchTube(p1, p2, pixscale(pd1->width), gle::hexTube);
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
        pd1->color.load_both();
        Vector mid = 0.5 * ( cx->sidePos1() + cx->sidePos2() );
        drawPoint(mid, pixscale(pd1->width));
        stretchTube(p2, mid, pixscale(pd2->width), gle::hexTube);
        stretchTube(p1, mid, pixscale(pd1->width), gle::hexTube);
        drawPoint(p1, pd1);
        drawPoint(p2, pd2);
        return;
    }
#endif
    
    float rad = pixscale(pd1->size);
    float Lr = cx->prop->length / rad;
    float iLr = ( pd1->width / pd1->size );
    
    //if ( cx->cosAngle() > 0 ) gle_color(1, 0.5, 0.25).load_both(); else gle_color(0, 1, 0).load_both();
    if ( pd1->visible )
    {
        float Z = cx->fiber1()->prop->disp->line_width / pd1->size + 0.4;
        pd1->color.load_both();
        glPushMatrix();
        gle::transAlignZ(p1, rad, pS-p1);
        glTranslatef(0, 0, Z);
        gle::blob();
        glTranslatef(0, 0,Lr-Z);
        gle::cuboid();
        glTranslatef(0, 0,-Lr);
        glScalef(iLr, iLr, Lr);
        gle::hexTube();
        glPopMatrix();
    }

    // draw a link between pS and p2
    glPushMatrix();
    gle::transAlignZ(p2, rad, pS-p2);
    if ( pd2->visible )
    {
        pd2->color.load_both();
        gle::blob();
    }
    glScalef(iLr, iLr, norm(pS-p2) / rad);
    gle::hexTube();
    glPopMatrix();
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
        dns = pixscale(gle::invsqrt(dns));
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
        Vector mid = 0.5 * ( p1 + p2 );
        Vector dir = normalize( p2 - p1 );
        glEnable(GL_CLIP_PLANE5);
        if ( pd1->visible )
        {
            setClipPlane(GL_CLIP_PLANE5, -dir, mid);
            pd1->color.load_both();
            glPushMatrix();
            transAlignZ(p1, pixscale(pd1->size), dir);
            //gle::needle();
            hemisphere2(); gle::cone3();
            glPopMatrix();
        }
        if ( pd2->visible )
        {
            setClipPlane(GL_CLIP_PLANE5,  dir, mid);
            pd2->color.load_both();
            glPushMatrix();
            transAlignZ(p2, pixscale(pd2->size), -dir);
            //gle::needle();
            hemisphere2(); gle::cone3();
            glPopMatrix();
        }
        glDisable(GL_CLIP_PLANE5);
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
        dns = pixscale(gle::invsqrt(dns));
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
        pd1->color.load_both(cx->fiber1()->disp->color.transparency());
        glDepthMask(GL_FALSE);
        stretchTube(p1, p2, pixscale(pd2->width), gle::hexTube);
        glDepthMask(GL_TRUE);
        return;
    }
#endif
    
    float wid = pixscale(pd1->width);
    float R = pixscale(pd1->size) / wid;
    float Lr = wid / norm( p2 - p1 );
    
    glPushMatrix();
    gle::stretchAlignZ(p1, p2, wid);
    pd1->color.load_both();
    gle::hexTube();
    if ( pd1->visible )
    {
        glScalef(R, R, R*Lr);
        gle::blob();
    }
    if ( pd2->visible )
    {
        glTranslatef(0, 0, 1/(R*Lr));
        pd2->color.load_both();
        gle::blob();
    }
    glPopMatrix();
}

