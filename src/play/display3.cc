// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "dim.h"
#include "simul.h"
#include "display3.h"
#include "display_color.h"
#include "modulo.h"

#include "fake.h"
#include "fiber_disp.h"
#include "line_disp.h"
#include "point_disp.h"

#include "opengl.h"
#include "gle_color_list.h"
#include "glut.h"
#include "gle.h"

using namespace gle;
extern Modulo const* modulo;


Display3::Display3(DisplayProp const* dp) : Display(dp)
{
}

//------------------------------------------------------------------------------
#pragma mark -

void Display3::drawSimul(Simul const& sim)
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

    drawBeads(sim.beads);
    drawSolids(sim.solids);
    drawSpheres(sim.spheres);
    
    if ( prop->single_select & 1 )
        drawSinglesF(sim.singles);
    
    if ( prop->couple_select & 1 )
        drawCouplesF(sim.couples);

    if ( prop->couple_select & 2 )
        drawCouplesA(sim.couples);

    if ( prop->couple_select & 4 )
        drawCouplesB(sim.couples);

    if ( prop->single_select & 2 )
        drawSinglesA(sim.singles);
    
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
#pragma mark -

inline void Display3::drawBallT(Vector const& pos, float radius, gle_color const& col) const
{
    glEnable(GL_LIGHTING);
    col.load_both();
    glPushMatrix();
    gle::translate(pos);
    gle::scale(radius);
    gle::dualPassSphere4();
    glPopMatrix();
}


inline void Display3::drawPoint(Vector const& pos, float size) const
{
    glEnable(GL_LIGHTING);
    glPushMatrix();
    gle::translate(pos);
    gle::scale(size*sFactor);
    gle::sphere2();
    glPopMatrix();
}


inline void Display3::drawPoint(Vector const& pos, PointDisp const* disp) const
{
    if ( disp->perceptible )
    {
        glEnable(GL_LIGHTING);
        glPushMatrix();
        gle::translate(pos);
        gle::scale(disp->size*sFactor);
        gle::sphere1();
#if ( 0 )
        if ( disp->symbol )
        {
            glDisable(GL_LIGHTING);
            disp->symbol_color.load();
            glRasterPos2f(0,0);
            glBitmap(0,0,0,0,-5,-4,0);
            glutBitmapCharacter(GLUT_BITMAP_9_BY_15, disp->symbol);
            glEnable(GL_LIGHTING);
        }
#endif
        glPopMatrix();
    }
}


//------------------------------------------------------------------------------
#pragma mark - Nicely Joined Fiber Rendering using clipping planes

#if ( 1 )
/**
This draws the model-segments, using function `set_color` to set display colors
*/
void Display3::drawFiberSegments(Fiber const& fib, real rad,
                                 void (*set_color)(Fiber const&, size_t, real), real beta) const
{
    Vector pos = fib.posPoint(0), old;
    Vector nxt = fib.posPoint(1);
    Vector dir = normalize(nxt-pos);  // could use _mm_rsqrt_ss
    
    glEnable(GL_CLIP_PLANE4);
    setClipPlane(GL_CLIP_PLANE4, -dir, nxt);
    set_color(fib, 0, beta);
    drawTube(pos, rad, nxt, gle::capedTube2);
    setClipPlane(GL_CLIP_PLANE4, dir, nxt);
    glEnable(GL_CLIP_PLANE5);

    // draw inner segments
    const size_t last = fib.lastSegment();
    for ( size_t inx = 1; inx <= last; ++inx )
    {
        old = pos;
        pos = nxt;
        nxt = fib.posPoint(inx+1);
        dir = normalize(nxt-old);
        set_color(fib, inx, beta);
        setClipPlane(GL_CLIP_PLANE5, -dir, pos);
        drawTube(old, rad, pos, gle::longTube2);
        setClipPlane(GL_CLIP_PLANE4,  dir, pos);
    }

    glDisable(GL_CLIP_PLANE5);
    // draw last segment:
    set_color(fib, last, beta);
    drawTube(nxt, rad, pos, gle::endedTube2);
    glDisable(GL_CLIP_PLANE4);
}


/**
This draws segments of length 'len' which may not correspond to the vertices
used to model the Fiber. All abscissa is relative to the MINUS_END.
The function `set_color` is called to set the color of the segments.
*/
void Display3::drawFiberSubSegments(Fiber const& fib, real rad,
                                    VisibleLattice::lati_t inx, const VisibleLattice::lati_t last,
                                    real abs, const real inc,
                                    void (*set_color)(Fiber const&, long, real),
                                    real fac, real facM, real facP) const
{
    Vector old = fib.displayPosM(abs);
    Vector pos = fib.displayPosM(abs+inc);
    Vector nxt = fib.displayPosM(abs+inc*2);
    Vector dir = normalize(nxt-pos);
    
    set_color(fib, inx++, facM);
    glEnable(GL_CLIP_PLANE4);
    setClipPlane(GL_CLIP_PLANE4, -dir, pos);
    if ( abs <= 0 )
        drawTube(old, rad, pos, gle::capedTube2);
    else
        drawTube(old, rad, pos, gle::halfTube2);
    setClipPlane(GL_CLIP_PLANE4, dir, pos);

    glEnable(GL_CLIP_PLANE5);
    abs += inc*2;
    // draw segments
    while ( inx < last )
    {
        abs += inc;
        old = pos;
        pos = nxt;
        nxt = fib.displayPosM(abs);
        dir = normalize(nxt-old);
        set_color(fib, inx++, fac);
        setClipPlane(GL_CLIP_PLANE5, -dir, pos);
        // can draw a disc to close the tube:
        drawTube(old, rad, pos, gle::longTube2);  //gle::obturatedTube2
        setClipPlane(GL_CLIP_PLANE4, dir, pos);
    }
    glDisable(GL_CLIP_PLANE5);

    // draw last segment, which may be truncated:
    set_color(fib, inx, facP);
    if ( abs >= fib.length() )
        drawTube(nxt, rad, pos, gle::endedTube2);
    else
        drawTube(nxt, rad, pos, gle::halfTube2);
    glDisable(GL_CLIP_PLANE4);
}

#else

//------------------------------------------------------------------------------
#pragma mark - Fiber Rendering using Spherocylinders

/**
This draws the model-segments, using function `set_color` to set display colors
*/
void Display3::drawFiberSegments(Fiber const& fib, real rad,
                                 void (*set_color)(Fiber const&, size_t, real), real beta) const
{
    Vector nxt = fib.posPoint(0);

    // rendering tubes as spherocylinders that join properly at any angle!
    const size_t last = fib.lastSegment();
    for ( size_t inx = 0; inx <= last; ++inx )
    {
        Vector pos = nxt;
        nxt = fib.posPoint(inx+1);
        set_color(fib, inx, beta);
        glPushMatrix();
        gle::transAlignZ(pos, 1.0, nxt-pos);
        gle::capsuleZ(norm(nxt-pos), rad);
        glPopMatrix();
    }
}


/**
This draws segments of length 'len' which may not correspond to the vertices
used to model the Fiber. All abscissa is relative to the MINUS_END.
The function `set_color` is called to set the color of the segments.
*/
void Display3::drawFiberSubSegments(Fiber const& fib, real rad,
                                    VisibleLattice::lati_t inx, const VisibleLattice::lati_t last,
                                    real abs, const real inc,
                                    void (*set_color)(Fiber const&, long, real),
                                    real fac, real facM, real facP) const
{
    Vector nxt = fib.displayPosM(abs);
    
    // rendering tubes as spherocylinders that join properly at any angle!
    for ( ; inx <= last; ++inx )
    {
        abs += inc;
        Vector pos = nxt;
        nxt = fib.displayPosM(abs);
        set_color(fib, inx, fac);
        glPushMatrix();
        gle::transAlignZ(pos, 1.0, nxt-pos);
        gle::capsuleZ(norm(nxt-pos), rad);
        glPopMatrix();
    }
}

#endif

//------------------------------------------------------------------------------
#pragma mark -

void color_not(Fiber const&, size_t, real)
{
}

// defined in display.cc
extern gle_color color_by_tension(Fiber const& fib, size_t seg);
extern gle_color color_by_tension_jet(Fiber const& fib, size_t seg);

void color_seg_tension(Fiber const& fib, size_t seg, real beta)
{
    real x = beta * fib.tension(seg);
    fib.disp->color.load_front(x);
}

void color_seg_tension_jet(Fiber const& fib, size_t seg, real beta)
{
    real x = beta * fib.tension(seg);
    // use jet coloring, where Lagrange multipliers are negative under compression
    gle_color::jet_color_alpha(x).load_front();
}

void color_seg_curvature(Fiber const& fib, size_t seg, real beta)
{
    if ( fib.nbPoints() > 2 )
    {
        real c = fib.curvature(std::max(seg, 1LU));
        real d = fib.curvature(std::min(seg+1, fib.lastSegment()));
        gle_color::jet_color(beta*(c+d)).load_front();
    }
    else
        gle_color::jet_color(0).load_front();
}

void color_seg_direction(Fiber const& fib, size_t seg, real)
{
    gle::radial_color(fib.dirSegment(seg)).load_front();
}

/// using distance from the minus end to the start of segment `seg`
void color_seg_distanceM(Fiber const& fib, size_t seg, real beta)
{
    real x = std::min(seg*beta, (real)32.0);
    fib.disp->color.load_front(std::exp(-x));
}

/// using distance from the plus end to the end of segment `seg`
void color_seg_distanceP(Fiber const& fib, size_t seg, real beta)
{
    real x = std::min((fib.lastSegment()-seg)*beta, (real)32.0);
    fib.disp->color.load_front(std::exp(-x));
}

/// color set according to distance to the confining Space
void color_seg_height(Fiber const& fib, size_t seg, real beta)
{
    real Z = 0;
    Space const* spc = fib.prop->confine_space_ptr;
    if ( spc )
        Z = -spc->signedDistanceToEdge(fib.posPoint(seg, 0.5));
#if ( DIM > 2 )
    else
        Z = fib.posPoint(seg,0.5).ZZ;
#endif
    gle_color::jet_color(Z*beta).load_front();
}


/// color set according to steric grid
void color_seg_grid(Fiber const& fib, size_t seg, real beta)
{
    Map<DIM> const& map = fib.simul().pointGridF.map();
    Vector w = fib.posPoint(seg, 0.5);
    size_t i = map.index(w);
    gle::alt_color(i).load_front();
}


void Display3::drawFiberLines(Fiber const& fib) const
{
    FiberDisp const*const disp = fib.prop->disp;
    const real rad = disp->line_width * sFactor;

    // set back color:
    if ( fib.prop->disp->coloring )
        fib.disp->color.load_back();
    else
        fib.prop->disp->back_color.load_back();
    glEnable(GL_LIGHTING);
    
    switch ( disp->line_style )
    {
        case 1:
            fib.disp->color.load_front();
            drawFiberSegments(fib, rad, color_not, 1.0);
            break;
        case 2:
        {
            const real beta = 1.0 / disp->tension_scale;
            drawFiberSegments(fib, rad, color_seg_tension, beta);
        } break;
        case 3:
        {
            const real beta = 1.0 / disp->tension_scale;
            drawFiberSegments(fib, rad, color_seg_tension_jet, beta);
        } break;
        case 4:
        {
            const real beta = 0.5 * fib.prop->disp->length_scale;
            drawFiberSegments(fib, rad, color_seg_curvature, beta);
        } break;
        case 5:
            drawFiberSegments(fib, rad, color_seg_direction, 1.0);
            break;
        case 6:
        {
            /** This is using transparency with segments that are not depth sorted
             but this code is only used in 2D normally, so it's okay */
            GLboolean cull = glIsEnabled(GL_CULL_FACE);
            glEnable(GL_CULL_FACE);
            const real beta = fib.segmentation() / disp->length_scale;
            drawFiberSegments(fib, rad, color_seg_distanceM, beta);
            if ( !cull ) glDisable(GL_CULL_FACE);
        } break;
        case 7:
        {
            /** This is using transparency with segments that are not depth sorted
             but this code is only used in 2D normally, so it's okay */
            GLboolean cull = glIsEnabled(GL_CULL_FACE);
            glEnable(GL_CULL_FACE);
            glCullFace(GL_BACK);
            const real beta = fib.segmentation() / disp->length_scale;
            drawFiberSegments(fib, rad, color_seg_distanceP, beta);
            if ( !cull ) glDisable(GL_CULL_FACE);
        } break;
        case 8:
        {
            const real beta = 1.0 / disp->length_scale;
            drawFiberSegments(fib, rad, color_seg_height, beta);
        } break;
        case 9:
            drawFiberSegments(fib, rad, color_seg_grid, 1.0);
            break;
    }
}


// displays segment 'inx' with transparency
void Display3::drawFiberSegmentT(Fiber const& fib, size_t inx) const
{
    FiberDisp const*const disp = fib.prop->disp;
    const real rad = disp->line_width * sFactor;
    real iseg = fib.segmentationInv();

    Vector A = fib.posP(inx);
    Vector B = fib.posP(inx+1);
    
    glEnable(GL_LIGHTING);
    /* Either CULL_FACE should be enable to hide the back side,
     or every primitive should be renderred with a double pass*/
    GLboolean cull = glIsEnabled(GL_CULL_FACE);
    glEnable(GL_CULL_FACE);
    
    if ( disp->line_style == 6 )
        color_seg_distanceM(fib, inx, fib.segmentation()/disp->length_scale);
    else if ( disp->line_style == 7 )
        color_seg_distanceP(fib, inx, fib.segmentation()/disp->length_scale);
    else if ( disp->line_style == 2 )
        color_by_tension(fib, inx).load_front();
    else if ( disp->line_style == 3 )
        color_by_tension_jet(fib, inx).load_front();
    else
        fib.disp->color.load_front();

#if 1
    if ( inx == 0 )
    {
        glEnable(GL_CLIP_PLANE5);
        if ( disp->line_style == 6 )
        {
            // cut the terminal segment according to length_scale
            real x = 2 * disp->length_scale * iseg;
            if ( x < 1.0 )
            {
                B = A + x * ( B - A );
                color_seg_distanceM(fib, 0.0, 1.0);
            }
        }
        setClipPlane(GL_CLIP_PLANE5, (A-B)*iseg, B);
        drawTube(A, rad, B, gle::capedTube2);
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
        if ( disp->line_style == 7 )
        {
            // cut the terminal segment according to length_scale
            real x = 2 * disp->length_scale * iseg;
            if ( x < 1.0 )
            {
                A = B + x * ( A - B );
                color_seg_distanceP(fib, inx, 1.0);
            }
        }
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
    gle::transAlignZ(0.5*(A+B), 1.0, B-A);
    gle::capsuleZ(fib.segmentation(), rad);
    glPopMatrix();
#endif
    if ( !cull ) glDisable(GL_CULL_FACE);
}


//------------------------------------------------------------------------------
#pragma mark -


void color_alternate(Fiber const& fib, long ix, real)
{
    if ( ix & 1 )
        fib.disp->color.load_front();
    else
        fib.disp->color.darken(0.75).load_front();
}


void color_by_lattice(Fiber const& fib, long ix, real scale)
{
    gle_color col = fib.disp->color.darken(scale*fib.visibleLattice()->data(ix));
    if ( col.transparent() )
        col.load_both();
    else
        col.load_front();
}


void color_by_lattice_jet(Fiber const& fib, long ix, real scale)
{
    gle_color::jet_color(scale*fib.visibleLattice()->data(ix)).load_front();
}


void Display3::drawFiberLattice(Fiber const& fib, VisibleLattice const& lat, real width,
                                void (*set_color)(Fiber const&, long, real)) const
{
    FiberDisp const*const disp = fib.prop->disp;

    glEnable(GL_LIGHTING);
    GLfloat blk[] = { 0.0, 0.0, 0.0, 1.0 };
    GLfloat bak[] = { 0.0, 0.0, 0.0, 1.0 };
    disp->back_color.store(bak);
    glMaterialfv(GL_BACK, GL_AMBIENT,  blk);
    glMaterialfv(GL_BACK, GL_DIFFUSE,  bak);
    glMaterialfv(GL_FRONT, GL_AMBIENT,  blk);
    glMaterialfv(GL_FRONT, GL_DIFFUSE,  blk);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, blk);
    glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, blk);
    glMateriali (GL_FRONT_AND_BACK, GL_SHININESS, 32);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);

    const real fac = 1.0 / disp->lattice_scale;
    const real uni = lat.unit();
    const real rad = width * sFactor;
    
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

    drawFiberSubSegments(fib, rad, inf, sup, abs, uni, set_color, fac, facM, facP);
}


void Display3::drawFiberLattice1(Fiber const& fib, VisibleLattice const& lat, real width) const
{
    drawFiberLattice(fib, lat, width, color_by_lattice);
}

void Display3::drawFiberLattice2(Fiber const& fib, VisibleLattice const& lat, real width) const
{
    drawFiberLattice(fib, lat, width, color_by_lattice_jet);
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
void Display3::drawFiberMinusEnd(Fiber const& fib, int style, real rad) const
{
    if ( rad > 0 )
    {
        switch(style)
        {
            default: break;
            case 1: drawObject(fib.posEndM(), rad, gle::sphere2); break;
            case 2: drawObject(fib.posEndM(), -fib.dirEndM(), rad, gle::longCone); break;
            case 3: drawObject(fib.posEndM(), -fib.dirEndM(), rad, gle::cylinder); break;
            case 4: drawObject(fib.posEndM(), -fib.dirEndM(), rad, gle::arrowTail); break;
            case 5: drawObject(fib.posEndM(),  fib.dirEndM(), rad, gle::arrowTail); break;
            case 6: drawObject(fib.posEndM(), -fib.dirEndM(), rad, gle::cube); break;
        }
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
void Display3::drawFiberPlusEnd(Fiber const& fib, int style, real rad) const
{
    if ( rad > 0 )
    {
        switch(style)
        {
            default: break;
            case 1: drawObject(fib.posEndP(), rad, gle::sphere2); break;
            case 2: drawObject(fib.posEndP(),  fib.dirEndP(), rad, gle::longCone); break;
            case 3: drawObject(fib.posEndP(),  fib.dirEndP(), rad, gle::cylinder); break;
            case 4: drawObject(fib.posEndP(),  fib.dirEndP(), rad, gle::arrowTail); break;
            case 5: drawObject(fib.posEndP(), -fib.dirEndP(), rad, gle::arrowTail); break;
            case 6: drawObject(fib.posEndP(),  fib.dirEndP(), rad, gle::cube); break;
        }
    }
}


void Display3::drawFiberSpeckles(Fiber const& fib) const
{
    FiberDisp const*const disp = fib.prop->disp;
    real rad = disp->speckle_size;

    if ( rad * uFactor < 2 )
        return;

    // display random speckles:
    if ( disp->speckle_style == 1 )
    {
        /*
         A simple random number generator seeded by fib.signature()
         is used to distribute points always at the same position
         with respect to the lattice of each fiber.
         */
        
        const real spread = disp->speckle_interval;
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
        const real sep = disp->speckle_interval;
        real a = sep * std::ceil( fib.abscissaM() / sep );
        while ( a <= fib.abscissaP() )
        {
            drawPoint(fib.pos(a), rad);
            a += sep;
        }
    }
}


void Display3::drawFiberPoints(Fiber const& fib) const
{
    FiberDisp const*const disp = fib.prop->disp;
    // diameter of lines and points in space units:
    real rad = disp->point_size;
    int style = disp->point_style & 3;

    if ( rad * uFactor < 2 )
        return;

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
        const float rad = disp->point_size*sFactor;
        const real sep = disp->point_interval;
        real ab = std::ceil(fib.abscissaM()/sep) * sep;
        for ( ; ab <= fib.abscissaP(); ab += sep )
            gle::drawCone(fib.pos(ab), fib.dir(ab), rad);
    }
    else if ( style == 3 )
    {
        // display middle of fiber:
        drawPoint(fib.posMiddle(), 2*rad);
    }
}


//------------------------------------------------------------------------------
#pragma mark -

void Display3::drawBead(Bead const& obj)
{
    const PointDisp * disp = obj.prop->disp;

    // display center:
    if ( disp->style & 2 )
    {
        bodyColor(obj);
        drawPoint(obj.position(), disp);
    }
    
#if ( DIM == 2 )
    // display outline:
    if ( disp->style & 4 )
    {
        glDisable(GL_LIGHTING);
        bodyColorF(obj).load();
        lineWidth(disp->width);
        drawFlat(obj.position(), obj.radius(), gle::circle);
    }
#endif
}

/**
 Display a bead as a sphere
 */
void Display3::drawBeadT(Bead const& obj)
{
    const PointDisp * disp = obj.prop->disp;
    
    if ( disp->style & 1 )
    {
        drawBallT(obj.position(), obj.radius(), bodyColorF(obj));
    }
}

//------------------------------------------------------------------------------
#pragma mark -

void Display3::drawSolid(Solid const& obj)
{
    const PointDisp * disp = obj.prop->disp;
    
    //display points:
    if ( disp->style & 2  &&  disp->size > 0 )
    {
        bodyColor(obj);
        for ( size_t p = 0; p < obj.nbPoints(); ++p )
        {
            if ( obj.radius(p) > 0 )
                drawPoint(obj.posP(p), disp);
            else
                drawObject(obj.posP(p), disp->size, gle::cube);
        }
    }
    
#if ( DIM == 3 )
    //special display for ParM simulations (DYCHE 2006; KINETOCHORES 2019)
    if ( obj.mark()  &&  disp->style & 4  &&  obj.nbPoints() > 1 )
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

    //display a signature for each Solid
    if ( disp->style & 8 )
    {
        char tmp[8];
        bodyColor(obj);
        snprintf(tmp, sizeof(tmp), "%u", obj.identity());
        gle::drawText(obj.posP(0), tmp, GLUT_BITMAP_HELVETICA_10);
    }
    
    //draw polygon around vertices of Solid
    if ( disp->style & 16 )
    {
        glDisable(GL_LIGHTING);
        lineWidth(1.0);
        bodyColorF(obj).load();
        glBegin(GL_LINE_LOOP);
        for ( size_t ii = 0; ii < obj.nbPoints(); ++ii )
            gleVertex(obj.posPoint(ii));
        glEnd();
    }
}


/**
 Display a semi-transparent disc / sphere
 */
void Display3::drawSolidT(Solid const& obj, size_t inx)
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

void Display3::drawSphere(Sphere const& obj)
{
    const PointDisp * disp = obj.prop->disp;
    
    //display center and surface points
    if ( disp->size > 0  &&  disp->style & 2 )
    {
        bodyColor(obj);
        drawObject(obj.posP(0), disp->size, gle::octahedron);
        for ( size_t i = obj.nbRefPoints; i < obj.nbPoints(); ++i )
            drawObject(obj.posP(i), disp->size, gle::cube);
    }

    //display reference points
    if ( disp->size > 0  &&  disp->style & 8 )
    {
        bodyColor(obj);
        for ( size_t i = 1; i < obj.nbRefPoints; i++ )
            drawPoint(obj.posP(i), disp);
    }
}

void Display3::drawSphereT(Sphere const& obj)
{
    const PointDisp * disp = obj.prop->disp;

    if ( disp->style & 5 )
    {
        const Vector C = obj.posP(0);
        bodyColorF(obj).load_front();
#if ( DIM < 3 )
        if ( disp->style & 1 )
            drawFlat(C, obj.radius(), gle::circle);
        if ( disp->style & 2 )
            drawBallT(C, obj.radius(), bodyColorF(obj));
#else
        Display::drawSphereT(C, obj.posP(1)-C, obj.posP(2)-C, obj.posP(3)-C, disp->style);
#endif
    }
}

//------------------------------------------------------------------------------

void Display3::drawOrganizer(Organizer const& obj) const
{
    PointDisp const* disp = obj.disp();
    
    if ( !disp )
        return;

    const real w = disp->width*sFactor;

    if ( disp->style & 2 )
    {
        Vector P, Q;
        bodyColor(disp, obj.signature());
        
        for ( size_t i = 0; obj.getLink(i, P, Q); ++i )
        {
            drawPoint(P, disp);
            if ( modulo ) modulo->fold(Q, P);
            gleTube(P, Q, w, gle::tube1);
        }
    }
    /**
     This displays the Solid connecting two Aster as a spindle.
     Used for Cleo Kozlowski simulation of C. elegans (2007)
     */
    if ( disp->style & 1 && obj.tag() == Fake::TAG )
    {
        Solid const* sol = Solid::toSolid(obj.core());
        if ( sol && sol->nbPoints() >= 4 )
        {
            bodyColor(*sol);
#if ( DIM == 3 )
            glPushMatrix();
            Vector3 a = 0.5 * (sol->posP(0) + sol->posP(2));
            Vector3 b = 0.5 * (sol->posP(1) + sol->posP(3));
            stretchAlignZ(a, b, 1);
            glColor3f(0.6f,0.6f,0.6f);
            gle::dualPass(gle::barrel);
            glPopMatrix();
#else
            for ( size_t ii = 0; ii < sol->nbPoints(); ii+=2 )
                gleTube(sol->posPoint(ii), sol->posPoint(ii+1), w, gle::hexTube);
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
#if ( 1 )
            if ( obj->disp()->style == 2 )
            {
                Space const* spc = obj->confineSpace();
                if ( spc )
                {
                    /// draw a disc tangent to the Space:
                    Vector pos = obj->posFoot();
                    Vector dir = spc->normalToEdge(pos);
                    drawObject(pos, dir, obj->disp()->size, gle::discUp);
                    continue;
                }
            }
#endif
            if ( obj->base() )
                drawObject(obj->posFoot(), obj->disp()->size, gle::octahedron);
            else
                drawPoint(obj->posFoot(), obj->disp());
        }
    }
}


void Display3::drawSingleA(Single const* obj) const
{
    const PointDisp * disp = obj->disp();
    Vector ph = obj->posHand();
    disp->color.load_both();
    
    if ( obj->hasForce() && disp->width > 0 )
    {
        Vector pf = obj->posFoot();
        if ( modulo ) modulo->fold(pf, ph);
#if ( 0 )
        if ( obj->disp()->style == 2 )
        {
            Space const* spc = obj->confineSpace();
            if ( spc )
            {
                // draw a disc tangent to the Space:
                disp->color2.load_both();
                drawObject(pf, spc->normalToEdge(pf), disp->size, gle::discUp);
                disp->color.load_both();
            }
        }
#endif
#if ( DIM > 2 )
        gleTube(pf, ph, disp->width*sFactor, gle::cone);
#else
        gle::drawBand(ph, disp->width*sFactor, disp->color, pf, disp->width*sFactor, disp->color.alpha_scaled(0.5));
#endif
    }
    if ( obj->base() )
        drawObject(ph, disp->size, gle::octahedron);
    else
        drawPoint(ph, disp);
}

void Display3::drawSinglesA(SingleSet const& set) const
{
    for ( Single * obj=set.firstA(); obj ; obj=obj->next() )
    {
        if ( obj->disp()->perceptible && obj->fiber()->disp->visible )
            drawSingleA(obj);
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
        drawHand2(cx->posFree(), cx->disp1());
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
        drawHand2(obj->posFree(), obj->disp12());
        obj = nxt;
    }
    while ( obj )
    {
        nxt = obj->next();
        drawHand2(obj->posFree(), obj->disp21());
        obj = nxt->next();
        drawHand2(nxt->posFree(), nxt->disp12());
    }
}


void Display3::drawCouplesA(CoupleSet const& set) const
{
    for ( Couple * cx=set.firstAF(); cx ; cx=cx->next() )
    {
        if ( cx->fiber1()->disp->visible && cx->disp1()->visible )
        {
#if ( 0 )
            // ENDOCYTOSIS 2015
            if ( cx->fiber1()->disp->color.transparent() )
                cx->disp1()->color.load_both(cx->fiber1()->disp->color.transparency());
            else
#endif
            drawHand(cx->posHand1(), cx->disp1());
        }
    }
    
    for ( Couple * cx=set.firstFA(); cx ; cx=cx->next() )
    {
        if ( cx->fiber2()->disp->visible && cx->disp2()->visible )
        {
#if ( 0 )
            // ENDOCYTOSIS 2015
            if ( cx->fiber2()->disp->color.transparent() )
                cx->disp1()->color.load_both(cx->fiber2()->disp->color.transparency());
            else
#endif
            drawHand(cx->posHand2(), cx->disp2());
        }
    }
}

void Display3::drawCoupleBfast(Couple const* cx) const
{
    const PointDisp * pd1 = cx->disp1();
    const PointDisp * pd2 = cx->disp2();
    
    Vector p1 = cx->posHand1();
    Vector p2 = cx->posHand2();
    
    if ( pd1 == pd2 )
    {
        if ( pd1->visible )
        {
            pd1->color.load_both();
            gleTube(p1, p2, pd2->width*sFactor, gle::hexTube);
            drawHand(p1, pd1);
            drawHand(p2, pd2);
        }
    }
    else if ( pd1->visible || pd2->visible )
    {
        pd1->color.load_both();
        gleTube(p1, p2, pd1->width*sFactor, gle::tube1);
        if ( pd1->visible ) drawHand(p1, pd1);
        if ( pd2->visible ) drawHand(p2, pd2);
    }
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
    
#if !FIBER_HAS_FAMILY
    if ( dns > 1e-6 )
    {
        // moving the 'hands' to the surface of the fiber:
        dns = sFactor / std::sqrt(dns);
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
    
    if ( pd1 == pd2 )
    {
#if ( 0 )
        // ENDOCYTOSIS 2015
        if ( cx->fiber1()->disp->color.transparent() )
        {
            pd1->color.load_both(cx->fiber1()->disp->color.transparency());
            glDepthMask(GL_FALSE);
            gleTube(p1, p2, pd2->width*sFactor, gle::hexTube);
            glDepthMask(GL_TRUE);
            continue;
        }
#endif
        if ( pd1->visible )
        {
            pd1->color.load_both();
            drawPoint(p1, pd1);
            drawPoint(p2, pd2);
#if FIBER_HAS_FAMILY
            // accurate rendering of Couple's composite link
            Vector mid = 0.5 * ( cx->sidePos1() + cx->sidePos2() );
            drawPoint(mid, pd1->width);
            gleTube(p2, mid, pd2->width*sFactor, gle::tube2); //gle::hexTube);
            gleTube(p1, mid, pd1->width*sFactor, gle::tube2); //gle::hexTube);
#elif ( 0 )
            drawPoint(cx->sidePos1(), pd1->width);
            drawPoint(cx->sidePos2(), pd1->width);
            gleTube(p2, cx->sidePos2(), pd2->width*sFactor, gle::hexTube);
            gleTube(p1, cx->sidePos1(), pd1->width*sFactor, gle::hexTube);
            gleTube(cx->sidePos1(), cx->sidePos2(), pd2->width*sFactor, gle::hexTube);
#else
            gleTube(p1, p2, pd2->width*sFactor, gle::hexTube);
#endif
        }
    }
    else if ( dns > 1e-6 )
    {
        Vector mid = 0.5 * ( p1 + p2 );
        
        glEnable(GL_CLIP_PLANE5);
        if ( pd1->visible )
        {
            setClipPlane(GL_CLIP_PLANE5, -dif, mid);
            pd1->color.load_front();
            gleTube(p1, p2, pd1->width*sFactor, gle::tube1);
            drawHand(p1, pd1);
        }
        
        if ( pd2->visible )
        {
            setClipPlane(GL_CLIP_PLANE5,  dif, mid);
            pd2->color.load_front();
            gleTube(p2, p1, pd2->width*sFactor, gle::tube1);
            drawHand(p2, pd2);
        }
        glDisable(GL_CLIP_PLANE5);
    }
    else
    {
        if ( pd1->visible ) drawHand(p1, pd1);
        if ( pd2->visible ) drawHand(p2, pd2);
    }
}


