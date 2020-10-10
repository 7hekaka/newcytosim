// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "display.h"
#include "organizer.h"
#include "property_list.h"
#include "hand_prop.h"
#include "sphere_prop.h"
#include "fiber_prop.h"
#include "point_disp.h"
#include "fiber_disp.h"
#include "line_disp.h"
#include "opengl.h"
#include "modulo.h"
#include "simul.h"
#include "field.h"
#include "gle.h"
#include "gle_color.h"
#include "gle_color_list.h"
#include "glapp.h"
#include "glut.h"

extern Modulo const* modulo;


//------------------------------------------------------------------------------
#pragma mark -


Display::Display(DisplayProp const* dp)
: pixelSize(1), uFactor(1), sFactor(1), prop(dp)
{
    assert_true(dp);
    
    prep_time = -1;
}

void Display::setPixelFactors(GLfloat ps, GLfloat u)
{
    pixelSize = ps;
    uFactor   = u;
    /*
     the 0.5 below comes from the fact that glPointSize uses diameter
     while most gle::primitives use radius as arguments
     */
    sFactor = 0.5f * u * ps;
}


void Display::drawSimul(Simul const& sim)
{
    gle::gleDrawText(Vector(0,0,0), "Empty Display::display", GLUT_BITMAP_8_BY_13);
}


void Display::display(Simul const& sim)
{
    // clear list of transparent objects
    zObjects.clear();
    
    /*
     Draw opaque objects:
     - depth buffer is writable
     - glColor specifies the Front material color
     */

#if ( DIM >= 3 )
    glEnable(GL_LIGHTING);
    glDepthMask(GL_TRUE);
#else
    glDisable(GL_LIGHTING);
#endif
    
    drawSimul(sim);
    
    /*
     Draw translucent objects:
     - make depth buffer readable only
     - objects are depth-sorted, from far to near
     - Dual pass is used to display back before front
     */

#if ( DIM >= 3 )
    
    glDisable(GL_CULL_FACE);
    glDepthMask(GL_FALSE);
    
    if ( zObjects.size() )
        drawTransparentObjects(zObjects);
    
    glEnable(GL_CULL_FACE);
    drawTransparentSpaces(sim.spaces);

#endif
    
    glDepthMask(GL_TRUE);
    
#ifndef NDEBUG
    gle::gleReportErrors(stderr, "in Display::display()");
#endif
}


/**
 To get correct display, it would be necessary to display all opaque objects first,
 and then all transparent objects for all tiles. Here, we calls Display::display()
 a number of times, and objects are sorted within each tile. The result is not perfect.
 */
void Display::displayTiled(Simul const& sim, int arg)
{
    assert_true(modulo);
    
    int l[3] = { 0 };
    int u[3] = { 0 };
    
    for ( size_t d = 0; d < DIM; ++d )
    {
        if ( modulo->isPeriodic(d) )
        {
            l[d] = (arg & (1<<d)) ? -1 : 0;
            u[d] = +1;
        }
    }
    
    glMatrixMode(GL_MODELVIEW);
    
    const Vector px = modulo->period(0);
    const Vector py = modulo->period(1);
    const Vector pz = modulo->period(2);

    for ( int dx = l[0]; dx <= u[0]; ++dx )
    for ( int dy = l[1]; dy <= u[1]; ++dy )
    for ( int dz = l[2]; dz <= u[2]; ++dz )
    {
        Vector T = dx * px + dy * py + dz * pz;
        gle::gleTranslate( T);
        display(sim);
        gle::gleTranslate(-T);
    }
}

//------------------------------------------------------------------------------
#pragma mark -


/**
 Create a FiberDisp for this Property if necessary
 */
void Display::prepareFiberDisp(FiberProp* fp, PropertyList& alldisp, gle_color col)
{
    FiberDisp *& disp = fp->disp;
    
    // recover existing property:
    if ( !disp )
        disp = static_cast<FiberDisp*>(alldisp.find("fiber:display", fp->name()));

    // create new property with default values:
    if ( disp == nullptr )
    {
        disp = new FiberDisp(fp->name());
        alldisp.push_back(disp);
        // set default:
        disp->color       = col;
        disp->back_color  = col.darken(0.5);
        disp->point_size  = prop->point_size;
        disp->line_width  = prop->line_width;
    }
    
    // parse user-provided values:
    if ( fp->display_fresh )
    {
        disp->read_string(fp->display, " in "+fp->name()+":display");
        fp->display_fresh = false;
    }
    
    if ( disp->coloring == FiberDisp::COLORING_CLUSTER )
        fiber_prep |= 1;
    
    if ( disp->line_style == 2 || disp->line_style == 3 )
        fiber_prep |= 2;

    if ( disp->coloring == FiberDisp::COLORING_AGE )
        fiber_prep |= 4;
}


/**
 set LineDisp for given Fiber
 */
void Display::prepareLineDisp(const Fiber * fib)
{
    assert_true(fib->prop);
    FiberDisp const*const disp = fib->prop->disp;
    
    if ( !fib->disp )
        fib->disp = new LineDisp();
    
    gle_color col = disp->color;
    
    // change body color depending on coloring mode:
    switch ( disp->coloring )
    {
        default:
        case FiberDisp::COLORING_OFF:
            col = disp->color;
            break;
        case FiberDisp::COLORING_RANDOM:
            col = gle::bright_color(fib->signature()).match_a(disp->color);
            break;
        case FiberDisp::COLORING_DIRECTION:
            col = gle::radial_color(fib->avgDirection());
            break;
        case FiberDisp::COLORING_MARK:
            col = gle::nice_color(fib->mark());
            break;
        case FiberDisp::COLORING_FLAG:
            col = gle::std_color(fib->flag());
            break;
#if FIBER_HAS_FAMILY
        case FiberDisp::COLORING_FAMILY:
            if ( fib->family_ )
                col = gle::nice_color(fib->family_->signature());
            else
                col = disp->color;
            break;
#endif
        case FiberDisp::COLORING_CLUSTER:
            col = gle::alt_color(fib->flag());
            break;
        case FiberDisp::COLORING_AGE:
            col = gle_color::jet_color((fib->age()-age_start)*age_scale, 1.0);
            break;
        case FiberDisp::COLORING_PSTATE:
            if ( fib->endStateP() > 0 )
                col = disp->end_color[std::min(fib->endStateP(),5U)];
            break;
    }
    
    LineDisp * self = fib->disp;
    self->color = col;
    
#if ( 0 )
    // colors of ends set to match body color:
    self->end_color[0] = col;
    self->end_color[1] = col;
#else
    // colors of ends for non-dynamic filaments:
    self->end_color[0] = disp->end_color[0];
    self->end_color[1] = disp->end_color[0];
#endif
    
    // For dynamic Fibers, change colors of tips according to state:
    if ( fib->endStateP() > 0 )
        self->end_color[0] = disp->end_color[std::min(fib->endStateP(),5U)];
    
    if ( fib->endStateM() > 0 )
        self->end_color[1] = disp->end_color[std::min(fib->endStateM(),5U)];

    bool hide = false;
    // hide right or left-pointing fibers:
    if (( disp->hide & 1 )  &&  dot(fib->diffPoints(0), disp->hide_axis) < 0 )
        hide = true;
    if (( disp->hide & 2 )  &&  dot(fib->diffPoints(0), disp->hide_axis) > 0 )
        hide = true;
    
#if ( DIM == 2 )
    // hide clockwise or counter-clockwise orientated fibers:
    if (( disp->hide & 4 )  &&  cross(fib->posP(0), fib->diffPoints(0)) < 0 )
        hide = true;
    if (( disp->hide & 8 )  &&  cross(fib->posP(0), fib->diffPoints(0)) > 0 )
        hide = true;
#elif ( DIM == 3 )
    // hide clockwise or counter-clockwise orientated fibers in the XY plane
    if (( disp->hide & 4 )  &&  cross(fib->posP(0), fib->diffPoints(0)).ZZ < 0 )
        hide = true;
    if (( disp->hide & 8 )  &&  cross(fib->posP(0), fib->diffPoints(0)).ZZ > 0 )
        hide = true;
#endif
    
#if ( 1 )
    // hide fibers depending on mask
    if ( fib->signature() & disp->mask_bitfield )
        hide = true;
#else
    if ( fib->mark() & disp->mask_bitfield )
        hide = true;
#endif

    // hide fibers in a specified state
    if ( fib->endStateP() == disp->hide_state )
        hide = true;
    
    // change color of 'hidden' filament:
    if ( hide )
        self->color = disp->hide_color;
    
    // default visibility set from class:
    if ( disp->visible )
    {
        // change visibility flag according to body color:
        if ( !self->color.visible() )
            self->visible = 0;
        else if ( self->color.transparent() )
            self->visible = -1;
        else
            self->visible = 1;
    }
    else
        self->visible = 0;
    
    // set parameters for exploded display
    if ( disp->explode )
        self->explode_shift = ( lcrng3(fib->signature()) * 0x1p-32 - 0.5 ) * disp->explode_range;
    else
        self->explode_shift = 0;
}


/**
 Create a PointDisp for this Property if necessary
 */
template < typename T >
void Display::preparePointDisp(T * p, PropertyList& alldisp, gle_color col)
{
    assert_true(p);
        
    PointDisp *& disp = p->disp;
    
    // search for matching property:
    if ( !disp )
        disp = static_cast<PointDisp*>(alldisp.find(p->category()+":display", p->name()));
    
    // create new property:
    if ( !disp )
    {
        //std::clog <<" new " << p->category() << ":display " << p->name() << "\n";
        disp = new PointDisp(p->category()+":display", p->name());
        disp->clear();
        alldisp.push_back(disp);
        // set default:
        disp->color  = col;
        disp->color2 = col.alpha_scaled(0.25);
        disp->size   = prop->point_size;
        if ( p->category() == "hand" )
            disp->width = prop->link_width;
        else
            disp->width = prop->line_width;
    }
    
    // parse display string once:
    if ( p->display_fresh )
    {
        disp->read_string(p->display, " in "+p->name()+":display");
        p->display_fresh = false;
    }
    
    disp->prepare(uFactor, sFactor, prop->style==1);
}

/**
 Perform the operations that are necessary to display the simulation:
 - create FiberDisp, HandDisp, SphereDisp, etc. (one per Property)
 - create LineDisp (one per Fiber)
 - set default values,
 - parse display strings
 .
*/
void Display::prepareForDisplay(Simul const& sim, PropertyList& alldisp)
{
    if ( prop->fold )
        sim.foldPosition();
    
    // counter to give different colors to the objects
    size_t idx = 0;
    
    PropertyList plist = sim.properties.find_all("fiber");
    
    fiber_prep = 0;
    // create a FiberDisp for each FiberProp:
    for ( Property* p : plist )
        prepareFiberDisp(static_cast<FiberProp*>(p), alldisp, gle::nice_color(idx++));

    if ( prep_time != sim.time() )
    {
        // the cluster analysis only needs to be done once per state:
        //prep_time = sim.time();
        if ( fiber_prep & 1 )
            sim.flagClusters(0, 0, 1);
        
        // if fiber tensions are used for display, recompute them now:
        if ( fiber_prep & 2 )
            sim.computeForces();
        
        // calculate Fiber::age() range and set color scaling factor:
        if ( fiber_prep & 4 )
        {
            size_t cnt;
            real avg, dev, mn, mx;
            FiberSet::infoBirthtime(sim.fibers.collect(), cnt, avg, dev, mn, mx);
            if ( mx > mn )
            {
                //std::clog << "=Fiber:age range [" << mn << " " << mx << " ]\n";
                age_start = sim.time() - mx;
                age_scale = 5.0 / ( mx - mn );
            }
            else
            {
                age_start = 0;
                age_scale = 1;
            }
        }
    }
    
    // attribute LineDisp, and set individual display values for all fibers
    for ( Fiber const* fib = sim.fibers.first(); fib; fib = fib->next() )
        prepareLineDisp(fib);
    
    //create a PointDisp for each HandProp:
    for ( Property * i : sim.properties.find_all("hand") )
        preparePointDisp(static_cast<HandProp*>(i), alldisp, gle::nice_color(idx++));
    
    //create a PointDisp for each SphereProp:
    for ( Property * i : sim.properties.find_all("sphere") )
        preparePointDisp(static_cast<SphereProp*>(i), alldisp, gle::bright_color(idx++));
    
    //create a PointDisp for each SolidProp:
    for ( Property * i : sim.properties.find_all("solid", "bead") )
        preparePointDisp(static_cast<SolidProp*>(i), alldisp, gle::bright_color(idx++));
    
    //create a PointDisp for each SpaceProp:
    gle_color col(DIM==3?0x00000044:0xAAAAAAFF);
    for ( Property * i : sim.properties.find_all("space") )
        preparePointDisp(static_cast<SpaceProp*>(i), alldisp, col);
}


//------------------------------------------------------------------------------
#pragma mark -


/**
 Draw a transparent Spaces (3D only)
 */
void Display::drawSpace(Space const* obj, bool opaque)
{
    const PointDisp * disp = obj->prop->disp;
    
    bool back = ( disp->visible & 2 ) & ( disp->color2.opaque() == opaque );
    bool front = ( disp->visible & 1 ) & ( disp->color.opaque() == opaque );
    
    if ( back | front )
    {
        lineWidth(disp->width);
        glEnable(GL_LIGHTING);
        GLboolean cull = glIsEnabled(GL_CULL_FACE);
        glEnable(GL_CULL_FACE);
        if ( back )
        {
            disp->color2.load_back();
            glCullFace(GL_FRONT);
            obj->draw();
        }
        if ( front )
        {
            disp->color.load_front();
            glCullFace(GL_BACK);
            obj->draw();
        }
        if ( !cull ) glDisable(GL_CULL_FACE);
    }
}


void Display::drawSpaces(SpaceSet const& set)
{
#if ( DIM == 3 )
    
    // draw non-transparent Spaces first:
    for ( Space * obj = set.first(); obj; obj=obj->next() )
    {
        if ( obj->prop->disp->visible )
        {
            obj->prop->disp->color.load_load();
            drawSpace(obj, true);
        }
    }

#else
    
    for ( Space * obj = set.first(); obj; obj=obj->next() )
    {
        const PointDisp * disp = obj->prop->disp;
        if ( disp->visible )
        {
            glDisable(GL_LIGHTING);
            lineWidth(disp->width);
            disp->color.load_load();
            obj->draw();
        }
    }
    
#endif
}


/**
 Draw transparent Spaces
 */
void Display::drawTransparentSpaces(SpaceSet const& set)
{
    for ( Space * obj = set.first(); obj; obj=obj->next() )
    {
        if ( obj->prop->disp->visible )
            drawSpace(obj, false);
    }
}


/**
 This displays only one Field, specified by DisplayProp:field_number
 
 GL_CULL_FACE and GL_LIGHTING should be disabled
 */
void Display::drawFields(FieldSet const& set)
{
#if ( DIM == 3 )
    // get current modelview transformation:
    GLfloat mat[16];
    glGetFloatv(GL_MODELVIEW_MATRIX, mat);
    
    // extract axis corresponding to vertical direction:
    Vector3 dir = normalize(Vector3(mat[2], mat[6], mat[10]));
#else
    Vector3 dir(0,0,1);
#endif
    
#if ( 1 )
    Field * obj = set.first();
#else
    for ( Field * obj = set.first(); obj; obj=obj->next() )
#endif
    
    if ( obj && obj->hasField() )
    {
        if ( obj->prop->visible == 1 )
            obj->draw(obj->prop->field_space_ptr, dir, 0);
        else if ( obj->prop->visible == 2 )
            obj->draw();
    }
}


//------------------------------------------------------------------------------
#pragma mark -

void Display::drawAverageFiber(ObjectList const& objs)
{
    Vector G, D, M, P;
    real S = FiberSet::infoPosition(objs, M, G, P);
    
    if ( S > REAL_EPSILON )
    {
        Vector MP = normalize( P - M );
        gle::gleCylinder(M, MP, 10*pixelSize);
        gle::gleCone(P, MP, 10*pixelSize);
        gle::gleObject(G, 10*pixelSize, gle::gleSphere2B);
    }
}


bool selectR(Object const* obj, void const* arg)
{
    Fiber const* fib = static_cast<Fiber const*>(obj);
    return fib->prop == arg  &&  dot(fib->diffPoints(0), fib->prop->disp->hide_axis) > 0;
}

bool selectL(Object const* obj, void const* arg)
{
    Fiber const* fib = static_cast<Fiber const*>(obj);
    return fib->prop == arg  &&  dot(fib->diffPoints(0), fib->prop->disp->hide_axis) < 0;
}

void Display::drawAverageFiber1(FiberSet const& fibers, void const* arg)
{
    ObjectList objs = fibers.collect(match_property, arg);

#if ( 1 )
    // highlight with a black outline
    glLineWidth(3);
    glColor3f(0,0,0);
    glDepthMask(GL_FALSE);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    drawAverageFiber(objs);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glDepthMask(GL_TRUE);
#endif
    
    glColor3f(1,1,1);
    drawAverageFiber(objs);
}


void Display::drawAverageFiber2(FiberSet const& fibers, void const* arg)
{
    ObjectList objsR = fibers.collect(selectR, arg);
    ObjectList objsL = fibers.collect(selectL, arg);

#if ( 1 )
    // highlight with a black outline
    glLineWidth(3);
    glColor3f(0,0,0);
    glDepthMask(GL_FALSE);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    drawAverageFiber(objsR);
    drawAverageFiber(objsL);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glDepthMask(GL_TRUE);
#endif
    
    // display right-pointing fibers in Red
    glColor3f(1,0,0);
    drawAverageFiber(objsR);
    
    // display left-pointing fibers in Green
    glColor3f(0,1,0);
    drawAverageFiber(objsL);
}    


void Display::drawMisc(Simul const& sim)
{
#if ( 0 )
    // display Steric Grid for visual debugging:
    glLineWidth(0.5);
    glColor3f(0, 0, 1);
    sim.pointGrid.draw();
#endif
    
    for ( Property const* i : sim.properties.find_all("fiber") )
    {
        FiberProp const* fp = static_cast<FiberProp const*>(i);
        if ( fp->disp->draw_average == 1 )
            drawAverageFiber1(sim.fibers, fp);
        else if ( fp->disp->draw_average == 2 )
            drawAverageFiber2(sim.fibers, fp);
    }
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
 */
void Display::drawFiberMinusEnd(Fiber const& fib, int style, real size) const
{
    real width = size * sFactor;
    if ( width > 0 )
    {
        switch(style)
        {
            default:
                break;
            case 1:
                gle::gleObject(fib.posEndM(), width, gle::gleSphere2B);
                break;
            case 2:
                gle::gleCone(fib.posEndM(), -fib.dirEndM(), width);
                break;
            case 3:
                gle::gleCylinder(fib.posEndM(), -fib.dirEndM(), width);
                break;
            case 4:
                gle::gleArrowTail(fib.posEndM(), fib.dirEndM(), width);
                break;
            case 5:
                gle::gleArrowTail(fib.posEndM(), -fib.dirEndM(), width);
                break;
            case 6:
                gle::gleObject(fib.posEndM(), width, gle::cube);
                break;
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
 */
void Display::drawFiberPlusEnd(Fiber const& fib, int style, real size) const
{
    real width = size * sFactor;
    if ( width > 0 )
    {
        switch(style)
        {
            default:
                break;
            case 1:
                gle::gleObject(fib.posEndP(), width, gle::gleSphere2B);
                break;
            case 2:
                gle::gleCone(fib.posEndP(), fib.dirEndP(), width);
                break;
            case 3:
                gle::gleCylinder(fib.posEndP(), fib.dirEndP(), width);
                break;
            case 4:
                gle::gleArrowTail(fib.posEndP(), fib.dirEndP(), width);
                break;
            case 5:
                gle::gleArrowTail(fib.posEndP(), -fib.dirEndP(), width);
                break;
            case 6:
                gle::gleObject(fib.posEndP(), width, gle::cube);
                break;
        }
    }
}


inline gle_color color_by_tension(Fiber const& fib, size_t seg, real beta)
{
    real x = beta * fib.tension(seg);
    if ( x > 0 )  // use normal color for extension
        return fib.disp->color.inverted().alpha(x);
    else          // invert color for compression
        return fib.disp->color.alpha(-x);
}

inline gle_color color_by_tension_jet(Fiber const& fib, size_t seg, real beta)
{
    real x = fib.tension(seg) * beta;
    // use rainbow coloring, where Lagrange multipliers are negative under compression
    return gle_color::jet_color(1-x);
}

inline gle_color color_by_curvature(real val)
{
    return gle_color::jet_color(val);
}

inline gle_color color_by_direction(Fiber const& fib, size_t seg)
{
    return gle::radial_color(fib.dirSegment(seg));
}

/// using the distance from the minus end to vertex `pti`
inline gle_color color_by_distanceM(Fiber const& fib, real pti, real beta)
{
    real x = std::min(pti*beta, (real)32.0);
    return fib.disp->color.alpha(std::exp(-x));
}

/// using the distance from the plus end to vertex `pti`
inline gle_color color_by_distanceP(Fiber const& fib, real pti, real beta)
{
    // using the distance at the vertex
    real x = std::min((fib.lastPoint()-pti)*beta, (real)32.0);
    return fib.disp->color.alpha(std::exp(-x));
}

/// color set according to distance to the confining Space
inline gle_color color_by_height(Fiber const& fib, size_t pti, real beta)
{
    real Z = 0;
    Space const* spc = fib.prop->confine_space_ptr;
    if ( spc )
        Z = -spc->signedDistanceToEdge(fib.posPoint(pti));
#if ( DIM > 2 )
    else
        Z = fib.posPoint(pti).ZZ;
#endif
    return gle_color::jet_color(Z*beta);
}



void Display::drawFiberLines(Fiber const& fib) const
{
    FiberDisp const*const disp = fib.prop->disp;
    GLfloat alpha = disp->color.transparency();

    switch ( disp->line_style )
    {
        case 1:
        {
            fib.disp->color.load();
            // display plain lines:
            lineWidth(disp->line_width);
#if ( DIM > 1 ) && REAL_IS_DOUBLE
            glEnableClientState(GL_VERTEX_ARRAY);
            glVertexPointer(DIM, GL_DOUBLE, 0, fib.addrPoints());
            glDrawArrays(GL_LINE_STRIP, 0, fib.nbPoints());
            glDisableClientState(GL_VERTEX_ARRAY);
#else
            glBegin(GL_LINE_STRIP);
            for ( size_t n = 0; n < fib.nbPoints(); ++n )
                gle::gleVertex(fib.posP(n));
            glEnd();
#endif
        } break;
        case 2:
        {
            // display segments with color indicating internal tension
            const real beta = 1.0 / disp->tension_scale;
            lineWidth(disp->line_width);
            glBegin(GL_LINES);
            for ( size_t n = 0; n < fib.lastPoint(); ++n )
            {
                color_by_tension(fib, n, beta).load();
                gle::gleVertex(fib.posP(n));
                gle::gleVertex(fib.posP(n+1));
            }
            glEnd();
        } break;
        case 3:
        {
            // display segments with color indicating internal tension
            const real beta = 1.0 / disp->tension_scale;
            lineWidth(disp->line_width);
            glBegin(GL_LINE_STRIP);
            for ( size_t n = 0; n < fib.lastPoint(); ++n )
            {
                color_by_tension_jet(fib, n, beta).load(alpha);
                gle::gleVertex(fib.posP(n));
                gle::gleVertex(fib.posP(n+1));
            }
            glEnd();
        } break;
        case 4:
        {
            // display segments with color indicating the curvature
            const real beta = 0.5 * disp->length_scale;
            lineWidth(disp->line_width);
            glBegin(GL_LINE_STRIP);
            real C = ( fib.nbPoints() > 2 ? fib.curvature(1) : 0.0 );
            color_by_curvature(beta*C).load(alpha);
            gle::gleVertex(fib.posEndM());
            for ( size_t n = 1; n < fib.lastPoint(); ++n )
            {
                color_by_curvature(beta*fib.curvature(n)).load(alpha);
                gle::gleVertex(fib.posP(n));
            }
            gle::gleVertex(fib.posEndP());
            glEnd();
        } break;
        case 5:
        {
            // color according to the angle with respect to the XY-plane:
            lineWidth(disp->line_width);
            glBegin(GL_LINES);
            for ( size_t n = 0; n < fib.lastPoint(); ++n )
            {
                color_by_direction(fib, n).load(alpha);
                gle::gleVertex(fib.posP(n));
                gle::gleVertex(fib.posP(n+1));
            }
            glEnd();
        } break;
        case 6:
        {
            // color according to the distance from the minus end
            const real beta = fib.segmentation() / disp->length_scale;
            lineWidth(disp->line_width);
            glBegin(GL_LINE_STRIP);
            color_by_distanceM(fib, 0, beta).load();
            gle::gleVertex(fib.posEndM());
            color_by_distanceM(fib, 0.5, beta).load();
            gle::gleVertex(fib.posPoint(0, 0.5));
            for ( size_t n = 1; n < fib.nbPoints(); ++n )
            {
                color_by_distanceM(fib, n, beta).load();
                gle::gleVertex(fib.posP(n));
            }
            glEnd();
        } break;
        case 7:
        {
            // color according to the distance from the plus end
            const real beta = fib.segmentation() / disp->length_scale;
            lineWidth(disp->line_width);
            glBegin(GL_LINE_STRIP);
            const size_t last = fib.lastSegment();
            for ( size_t n = 0; n < last; ++n )
            {
                color_by_distanceP(fib, n, beta).load();
                gle::gleVertex(fib.posP(n));
            }
            color_by_distanceP(fib, last+0.5, beta).load();
            gle::gleVertex(fib.posPoint(last, 0.5));
            color_by_distanceP(fib, last+0.75, beta).load();
            gle::gleVertex(fib.posPoint(last, 0.75));
            color_by_distanceP(fib, last+1.0, beta).load();
            gle::gleVertex(fib.posEndP());
            glEnd();
        } break;
        case 8:
        {
            // color according to distance to the confining Space
            const real beta = 1.0 / disp->length_scale;
            lineWidth(disp->line_width);
            glBegin(GL_LINE_STRIP);
            for ( size_t n = 0; n < fib.nbPoints(); ++n )
            {
                color_by_height(fib, n, beta).load(alpha);
                gle::gleVertex(fib.posP(n));
            }
            glEnd();
        } break;
    }
}


void Display::drawFiberSegmentT(Fiber const& fib, size_t inx) const
{
    FiberDisp const*const disp = fib.prop->disp;
    lineWidth(disp->line_width);
    
    glDisable(GL_LIGHTING);
    if ( disp->line_style == 6 )
    {
        const real beta = fib.segmentation() / disp->length_scale;
        glBegin(GL_LINE_STRIP);
        color_by_distanceM(fib, inx, beta).load();
        gle::gleVertex(fib.posP(inx));
        if ( inx == 0 )
        {
            color_by_distanceM(fib, 0.25, beta).load();
            gle::gleVertex(fib.posPoint(0, 0.25));
            color_by_distanceM(fib, 0.5, beta).load();
            gle::gleVertex(fib.posPoint(0, 0.5));
        }
        color_by_distanceM(fib, inx+1, beta).load();
        gle::gleVertex(fib.posP(inx+1));
        glEnd();
    }
    else if ( disp->line_style == 7 )
    {
        const real beta = fib.segmentation() / disp->length_scale;
        glBegin(GL_LINE_STRIP);
        color_by_distanceP(fib, inx, beta).load();
        gle::gleVertex(fib.posP(inx));
        if ( inx == fib.lastSegment() )
        {
            color_by_distanceP(fib, inx+0.5, beta).load();
            gle::gleVertex(fib.posPoint(inx, 0.5));
            color_by_distanceP(fib, inx+0.75, beta).load();
            gle::gleVertex(fib.posPoint(inx, 0.75));
        }
        color_by_distanceP(fib, inx+1, beta).load();
        gle::gleVertex(fib.posP(inx+1));
        glEnd();
    }
    else
    {
        fib.disp->color.load();
        glBegin(GL_LINES);
        gle::gleVertex(fib.posP(inx));
        gle::gleVertex(fib.posP(inx+1));
        glEnd();
    }
}


void Display::drawFiberSpeckles(Fiber const& fib) const
{
    FiberDisp const*const disp = fib.prop->disp;
    
    // display random speckles:
    if ( disp->speckle_style == 1 )
    {
        /*
         A simple random number generator seeded by fib.signature()
         is used to distribute points always at the same position
         with respect to the lattice of each fiber.
         */
        pointSize(disp->speckle_size);
        glBegin(GL_POINTS);
        
        const real spread = disp->speckle_interval;
        constexpr real TINY = 0x1p-32;
        // draw speckles below the origin of abscissa
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
                gle::gleVertex(fib.pos(a));
                z = lcrng2(z);
                a += spread * std::log(z*TINY);
            }
        }
        // draw speckles above the origin of abscissa
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
                gle::gleVertex(fib.pos(a));
                z = lcrng1(z);
                a -= spread * std::log(z*TINY);
            }
        }
        glEnd();
    }
    else if ( disp->speckle_style == 2 )
    {
        // display regular speckles
        pointSize(disp->speckle_size);
        glBegin(GL_POINTS);
        //we distribute points regularly along the center line
        const real grad = disp->speckle_interval;
        real ab = grad * std::ceil( fib.abscissaM() / grad );
        while ( ab <= fib.abscissaP() ) {
            gle::gleVertex(fib.pos(ab));
            ab += grad;
        }
        glEnd();
    }
}


void Display::drawFiberPoints(Fiber const& fib) const
{
    FiberDisp const*const disp = fib.prop->disp;

    if ( disp->point_style == 1 )
    {
        // display vertices:
        pointSize(disp->point_size);
        glBegin(GL_POINTS);
        for ( size_t ii = 0; ii < fib.nbPoints(); ++ii )
            gle::gleVertex(fib.posP(ii));
        glEnd();
    }
    else if ( disp->point_style == 2 )
    {
        // display arrowheads along the fiber:
        const real siz = disp->point_size*sFactor;
        const real sep = disp->point_interval;
        real ab = std::ceil(fib.abscissaM()/sep) * sep;
        for ( ; ab <= fib.abscissaP(); ab += sep )
            gle::gleCone(fib.pos(ab), fib.dir(ab), siz);
    }
    else if ( disp->point_style == 3 )
    {
        // display only middle of fiber:
        gle::gleObject(fib.posMiddle(), 2*disp->point_size, gle::gleSphere2B);
    }
}


//------------------------------------------------------------------------------
#pragma mark - Lattice

void set_lattice_color(Fiber const& fib, real val)
{
    fib.prop->disp->color.darken(val).load();
    /*
       FiberDisp const*const disp = fib.prop->disp;
       const gle_color col = disp->color;
       if ( disp->lattice_rescale )
           // use this if the lattice cells hold a quantity:
           col.darken( val * len / ( uni * disp->lattice_scale )).load();
       else // use this if the lattice cells hold a concentration:
           col.darken( val / disp->lattice_scale ).load();
    */
}

/**
 This style uses one vertex for each site, positionned at the center of the range
 OpenGL will interpolate the colors, and each site will be covered by a gradient.
 */
void Display::drawFiberLattice1(Fiber const& fib, VisibleLattice const& lat, real width) const
{
    const real uni = lat.unit();
    const auto inf = lat.indexM();
    const auto sup = lat.indexP();
    assert_true( inf <= sup );
    
    FiberDisp const*const disp = fib.prop->disp;
    const real fac = 1.0 / disp->lattice_scale;

    lineWidth(width);
    glBegin(GL_LINE_STRIP);
    if ( inf == sup )
    {
        //the Fiber is entirely covered by one site!
        real len = fib.abscissaP() - fib.abscissaM();
        set_lattice_color(fib, (fac*lat.data(inf))*(uni/len));
        gle::gleVertex(fib.posEndM());
        gle::gleVertex(fib.posEndP());
    }
    else
    {
        const real lenM = uni * (inf+1) - fib.abscissaM();  // should be positive!
        const real lenP = fib.abscissaP() - uni * sup;      // should be positive!

        real facM = fac;
        real facP = fac;
        if ( disp->lattice_rescale )
        {
            facM = ( lenM > 0.01*uni ? fac*uni/lenM : fac );
            facP = ( lenP > 0.01*uni ? fac*uni/lenP : fac );
        }

        // the terminal site may be truncated
        set_lattice_color(fib, facM*lat.data(inf));
        gle::gleVertex(fib.posEndM());
        if ( uni*(inf+0.5) > fib.abscissaM() )
            gle::gleVertex(fib.pos(uni*(inf+0.5)));
        
        for ( auto h = inf+1; h < sup; ++h )
        {
            set_lattice_color(fib, fac*lat.data(h));
            gle::gleVertex(fib.pos(uni*(h+0.5)));
        }
        
        // the terminal site may be truncated
        set_lattice_color(fib, facP*lat.data(sup));
        if ( uni*(sup+0.5) < fib.abscissaP() )
            gle::gleVertex(fib.pos(uni*(sup+0.5)));
        gle::gleVertex(fib.posEndP());
    }
    glEnd();
}


/**
 This style, uses two vertices for each site, positionned at the extremity of the range,
 and each site is entirely covered by the color corresponding to the value.
 */
void Display::drawFiberLattice2(Fiber const& fib, VisibleLattice const& lat, real width) const
{
    const real uni = lat.unit();
    const auto inf = lat.indexM();
    const auto sup = lat.indexP();
    assert_true( inf <= sup );
    
    FiberDisp const*const disp = fib.prop->disp;
    const real fac = 1.0 / disp->lattice_scale;

    lineWidth(width);
    glBegin(GL_LINE_STRIP);
    if ( inf == sup )
    {
        //the Fiber is entirely covered by one site!
        real len = fib.abscissaP() - fib.abscissaM();
        set_lattice_color(fib, (fac*lat.data(inf))*(uni/len));
        gle::gleVertex(fib.posEndM());
        gle::gleVertex(fib.posEndP());
    }
    else
    {
        const real lenM = uni * (inf+1) - fib.abscissaM();  // should be positive!
        const real lenP = fib.abscissaP() - uni * sup;      // should be positive!

        real facM = fac;
        real facP = fac;
        if ( disp->lattice_rescale )
        {
            facM = ( lenM > 0.01*uni ? fac*uni/lenM : fac );
            facP = ( lenP > 0.01*uni ? fac*uni/lenP : fac );
        }

        // the terminal site may be truncated
        set_lattice_color(fib, facM*lat.data(inf));
        gle::gleVertex(fib.posEndM());
        
        for ( auto h = inf+1; h < sup; ++h )
        {
            gle::gleVertex(fib.pos(uni*h));
            set_lattice_color(fib, fac*lat.data(h));
            gle::gleVertex(fib.pos(uni*h));
        }
        
        // the terminal site may be truncated
        set_lattice_color(fib, facP*lat.data(sup));
        gle::gleVertex(fib.pos(uni*sup));
        gle::gleVertex(fib.posEndP());
    }
    glEnd();
}


/**
 Indicate the edges between sites with small dots
 */
void Display::drawFiberLatticeEdges(Fiber const& fib, VisibleLattice const& lat, real) const
{
    const real uni = lat.unit();
    const auto inf = lat.indexM();
    const auto sup = lat.indexP();

    fib.disp->color.load();
    pointSize(fib.prop->disp->point_size);
    glBegin(GL_POINTS);
    for ( auto h = inf+1; h <= sup; ++h )
        gle::gleVertex(fib.posM(uni*h-fib.abscissaM()));
    glEnd();
}


void Display::drawFiberLabels(Fiber const& fib, void* font) const
{
    FiberDisp const*const disp = fib.prop->disp;
    char str[32];

    glDisable(GL_LIGHTING);
    if ( disp->label_style & 1 )
    {
        // draw fiber name and vertex indices
        int n = snprintf(str, sizeof(str), " %u ", fib.identity());
        for ( size_t ii = 0; ii < fib.nbPoints(); ++ii )
        {
            snprintf(str+n, sizeof(str)-n, "%lu", ii);
            gle::gleDrawText(fib.posP(ii), str, font);
        }
    }
    if ( disp->label_style & 2 )
    {
        // draw fiber name and abscissa value at vertices
        int n = snprintf(str, sizeof(str), " %u ", fib.identity());
        for ( size_t ii = 0; ii < fib.nbPoints(); ++ii )
        {
            snprintf(str+n, sizeof(str)-n, "%.3f", fib.abscissaPoint(ii));
            gle::gleDrawText(fib.posP(ii), str, font);
        }
    }
    if ( disp->label_style & 4 )
    {
        // display integral abscissa along the fiber
        snprintf(str, sizeof(str), "%.3f", fib.abscissaM());
        gle::gleDrawText(fib.posEndM(), str, font);
        
        int s = (int)std::ceil(fib.abscissaM());
        int e = (int)std::floor(fib.abscissaP());
        for ( int a = s; a <= e; ++a )
        {
            snprintf(str, sizeof(str), "%i", a);
            gle::gleDrawText(fib.pos(a), str, font);
        }
        
        snprintf(str, sizeof(str), "%.3f", fib.abscissaP());
        gle::gleDrawText(fib.posEndP(), str, font);
    }
}


/// display forces acting on the vertices, with lines
void Display::drawFiberForces(Fiber const& fib, real scale) const
{
    glDisable(GL_LIGHTING);
    glBegin(GL_LINES);
    for ( size_t ii = 0; ii < fib.nbPoints(); ++ii )
    {
        Vector p = fib.posP(ii);
        Vector q = p + scale * fib.netForce(ii);
        gle::gleVertex(p);
        gle::gleVertex(q);
    }
    glEnd();
}

//------------------------------------------------------------------------------
#pragma mark - Specific styles

/// used for drawFilament
inline void drawMonomer(Vector3 const& pos, real rad)
{
    gle::gleObject(pos, rad, gle::gleSphere2B);
}


/**
 This renders a protofilament by drawing spheres of alternating colors,
 along the backbone of the `Fiber` at distance 4nm from each other.
 */
void Display::drawFilament(Fiber const& fib,
                           gle_color const& color1,
                           gle_color const& color2,
                           gle_color const& colorE) const
{
    // axial translation between two sucessive monomers:
    const real dab = 0.004;
    // enlarge radius of monomers to make them overlap
    const real rad = 0.65 * dab;
    
    real ab = 0;
    
    glEnable(GL_CLIP_PLANE4);
    
    int cnt = 0;
    // increment until we reach the MINUS_END
    while ( ab <= fib.abscissaM() )
    {
        ++cnt;
        ab += dab;
    }
    Vector3 p(fib.pos(ab)), q;
    // draw the monomers until the PLUS_END:
    while ( ab < fib.abscissaP() )
    {
        q = p;
        ab += dab;
        p = Vector3(fib.pos(ab));

        // use different tones to individualize the two strands:
        if ( ++cnt & 1 )
            color1.load_load();
        else
            color2.load_load();
        
        // change color for the last monomer:
        if ( ab + dab > fib.abscissaP() )
        {
            colorE.load_load();
            glDisable(GL_CLIP_PLANE4);
        }
        
        // set clipping plane with the next monomer
        gle::setClipPlane(GL_CLIP_PLANE4, normalize(q-p), (p+q)*0.5);
        
        drawMonomer(q, rad);
        
        // set cliping plane with the previous:
        gle::setClipPlane(GL_CLIP_PLANE5, normalize(p-q), (p+q)*0.5);
        
        glEnable(GL_CLIP_PLANE5);
    }
    glDisable(GL_CLIP_PLANE4);
    glDisable(GL_CLIP_PLANE5);
}


/**
 This renders 26 spheres positionned on a right-handed helix,
 making one turn every 74nm, with a max width of ~ 9nm.
 This is roughly Ken Holmes' model of F-actin:
 Nature 347, 44 - 49 (06 September 1990); doi:10.1038/347044a0
 which shows half a turn in 37nm containing 13 monomers.
 */
void Display::drawActin(Fiber const& fib,
                        gle_color const& color1,
                        gle_color const& color2,
                        gle_color const& colorE) const
{    
    // axial translation between two sucessive monomers:
    const real dab = 0.00275;
    // enlarge radius of monomers to make them overlap
    const real rad = 1.3 * dab;
    // distance from central axis to center of monomers
    real off = 0.0045 - dab;
    
    /*
     The filamentous actin structure can be considered to be a single stranded
     levorotatory helix with a rotation of 166° around the helical axis
     and an axial translation of 27.5 Å
    */
    // rotation angle between consecutive monomers
    const real dan = -166 * M_PI / 180;
    const real cs = std::cos(dan);
    const real sn = std::sin(dan);

    real ab = 0;
    Vector3 d(fib.dirEndM());  // unit tangent to centerline
    Vector3 n = fib.adjustedNormal(d);
    //std::clog << fib.reference() << " " << n << "    " << n.normSqr() << '\n';
    
    Vector3 p, q;
    
    glEnable(GL_CLIP_PLANE4);
    
    int cnt = 0;
    // rotate until we reach the MINUS_END
    while ( ab <= fib.abscissaM() )
    {
        ++cnt;
        ab += dab;
        n = d.rotateOrtho(n, cs, sn);
    }
    p = Vector3(fib.pos(ab)) + off * n;
    // draw the monomers until the PLUS_END:
    while ( ab < fib.abscissaP() )
    {
        q = p;
        ab += dab;
        d = Vector3(fib.dir(ab));
        n = d.rotateOrtho(n, cs, sn);
        p = Vector3(fib.pos(ab)) + off * n;
        
        // use different tones to individualize the two strands:
        if ( ++cnt & 1 )
            color1.load_load();
        else
            color2.load_load();

        // change color for the last monomer:
        if ( ab + dab > fib.abscissaP() )
        {
            colorE.load_load();
            glDisable(GL_CLIP_PLANE4);
        }
        
        // set clipping plane with the next monomer
        gle::setClipPlane(GL_CLIP_PLANE4, normalize(q-p), (p+q)*0.5);
        
        drawMonomer(q, rad);
        
        // set cliping plane with the previous:
        gle::setClipPlane(GL_CLIP_PLANE5, normalize(p-q), (p+q)*0.5);
        
        glEnable(GL_CLIP_PLANE5);
    }
    glDisable(GL_CLIP_PLANE4);
    glDisable(GL_CLIP_PLANE5);
}


/**
 This renders a Microtubule using spheres of alternating colors
 colorA for alpha-tubulin
 colorB for beta-tubulin
 */
void Display::drawMicrotubule(Fiber const& fib,
                              gle_color const& colorA,
                              gle_color const& colorB,
                              gle_color const& colorE) const
{
    // precalculated 3-start helical trajectory, for 13 monomers:
    //real dx[] = {0,0.000923,0.001846,0.002769,0.003692,0.004615,0.005538,0.006461,0.007384,0.008308,0.009231,0.010154,0.011077};
    // some protofilaments are shifted by 8 nm backward:
    real dx[] = {0,0.000923,0.001846,0.002769,0.003692,0.004615,0.005538,0.006461,0.007384,0.000308,0.001231,0.002154,0.003077};
    real dy[] = {0.8855,0.5681,0.1205,-0.3546,-0.7485,-0.9709,-0.9709,-0.7485,-0.3546,0.1205,0.5681,0.8855,1.0000};
    real dz[] = {-0.4647,-0.8230,-0.9927,-0.9350,-0.6631,-0.2393,0.2393,0.6631,0.9350,0.9927,0.8230,0.4647,0};

    // axial translation (one monomer)
    const real sa = 0.004;
    // axial translation (one heterodimer)
    const real dab = 0.008;
    // enlarged radius of monomers makes them overlap slighlty
    const real rad = 0.003;
    // distance from central axis to center of monomers, such that diameter is 25nm
    real off = 0.025 / 2 - rad;

    const real abmax = fib.abscissaP();
    real ab = dab * std::ceil( fib.abscissaM() / dab );
    Vector3 d(fib.dir(ab));   // unit tangent vector
    Vector3 n = fib.adjustedNormal(d);
    
    while ( ab+6*sa < abmax )
    {
        d = Vector3(fib.dir(ab));
        Vector3 p(fib.pos(ab));
        // adjust 'n' to keep it orthogonal to 'd':
        n = d.orthogonal(n, 1.0);
        // set two vectors orthogonal to 'd' of norm 'off':
        Vector3 e = n * off;
        Vector3 f = cross(d, e);

        colorA.load_load();
        for ( int i = 0; i < 13; ++i )
            drawMonomer(p+dx[i]*d+dy[i]*e+dz[i]*f, rad);

        colorB.load_load();
        for ( int i = 0; i < 13; ++i )
            drawMonomer(p+(sa+dx[i])*d+dy[i]*e+dz[i]*f, rad);

        ab += dab;
    }
    // at the plus-end, only draw monomers below the end
    while ( ab+sa < abmax )
    {
        d = Vector3(fib.dir(ab));
        Vector3 p(fib.pos(ab));
        // adjust 'n' to keep it orthogonal to 'd':
        n = d.orthogonal(n, 1.0);
        // set two vectors orthogonal to 'd' of norm 'off':
        Vector3 e = n * off;
        Vector3 f = cross(d, e);

        for ( int i = 0; i < 13; ++i )
        {
            if ( ab+sa+dx[i] < abmax )
            {
                colorA.load_load();
                drawMonomer(p+dx[i]*d+dy[i]*e+dz[i]*f, rad);
                if ( ab+5.2*sa+dx[i] < abmax )
                    colorB.load_load();
                else
                    colorE.load_load();
                drawMonomer(p+(sa+dx[i])*d+dy[i]*e+dz[i]*f, rad);
            }
        }
        ab += dab;
    }
}


void Display::drawFiber(Fiber const& fib)
{
    FiberDisp const*const disp = fib.prop->disp;
    int line_style = disp->line_style;

    VisibleLattice const* lat = fib.visibleLattice();
    if ( lat && lat->ready() )
    {
        // if the Lattice is displayed, do not draw backbone:
        switch ( disp->lattice_style )
        {
            case 1:
                drawFiberLattice1(fib, *lat, disp->line_width);
                line_style = 0;
                break;
            case 2:
                drawFiberLattice2(fib, *lat, disp->line_width);
                line_style = 0;
                break;
            case 3:
                drawFiberLatticeEdges(fib, *lat, disp->line_width);
                line_style = 0;
                break;
        }
    }
    
#if ( DIM == 3 )
    // transparent styles in 3D
    if (( line_style==1 && fib.disp->color.transparent()))
    {
        for ( size_t i = 0; i < fib.lastPoint(); ++i )
            zObjects.push_back(zObject(&fib, i));
        line_style = 0;
    }
    else if ( line_style == 6 )
    {
        // color according to the distance from the minus end
        const real beta = fib.segmentation() / disp->length_scale;
        for ( size_t i = 0; i < fib.lastPoint(); ++i )
            if ( color_by_distanceM(fib, i, beta).visible() )
                zObjects.push_back(zObject(&fib, i));
        line_style = 0;
    }
    else if ( line_style == 7 )
    {
        // color according to the distance from the plus end
        const real beta = fib.segmentation() / disp->length_scale;
        for ( size_t i = 0; i < fib.lastPoint(); ++i )
            if ( color_by_distanceP(fib, i+1, beta).visible() )
                zObjects.push_back(zObject(&fib, i));
        line_style = 0;
    }
#endif
    
    if ( line_style )
    {
        gle_color col1 = fib.disp->color;
        gle_color col2 = fib.disp->color.darken(0.625);
        gle_color colE = fib.disp->end_color[0];
        
        // load backface color:
        if ( fib.prop->disp->coloring )
            col1.load_back();
        else
            fib.prop->disp->back_color.load_back();

        if ( disp->line_style != 1 || disp->style == 0 )
            drawFiberLines(fib);
        else if ( disp->style == 1 )
            drawFilament(fib, col1, col2, colE);
        else if ( disp->style == 2 )
            drawActin(fib, col1, col2, colE);
        else if ( disp->style == 3 )
            drawMicrotubule(fib, col1, col2, colE);
    }

    if ( disp->point_style > 0 )
    {
        fib.disp->color.load_load();
        disp->back_color.load_back();
        drawFiberPoints(fib);
    }
    
    if ( disp->speckle_style > 0 )
    {
        fib.disp->color.load_load();
        disp->back_color.load_back();
        drawFiberSpeckles(fib);
    }

    // draw other fiber elements only if fiber is fully visible:
    //if ( fib.disp->visible > 0 )
    {
        if ( disp->label_style )
        {
            fib.disp->color.load(0.5);
            drawFiberLabels(fib, GLUT_BITMAP_HELVETICA_10);
        }
        
        if ( disp->end_style[1] )
        {
            fib.disp->end_color[1].load_load();
            //fib.disp->color.load_load();
            disp->back_color.load_back();
            drawFiberMinusEnd(fib, disp->end_style[1], disp->end_size[1]);
        }
        
        if ( disp->end_style[0] )
        {
            fib.disp->end_color[0].load_load();
            //fib.disp->color.load_load();
            disp->back_color.load_back();
            drawFiberPlusEnd(fib, disp->end_style[0], disp->end_size[0]);
        }
        
        if ( disp->force_scale > 0 )
        {
            disp->force_color.load();
            lineWidth(disp->point_size);
            drawFiberForces(fib, disp->force_scale);
        }
    }
}


void Display::drawFibers(FiberSet const& set)
{
#if ( 1 )
    // display Fibers in a random (ever changing) order:
    for ( Fiber const* fib = set.first(); fib ; fib=fib->next() )
#else
    // display the Fiber always in the same order:
    for( Fiber const* fib = set.firstID(); fib; fib=set.nextID(fib) )
#endif
    {
        if ( fib->disp->visible )
            drawFiber(*fib);
    }
}


//------------------------------------------------------------------------------
#pragma mark -

void Display::drawCouplesF(CoupleSet const& set) const
{
    if ( prop->couple_flip )
        drawCouplesF2(set);
    else
        drawCouplesF1(set);
}


void Display::drawCouplesB(CoupleSet const& set) const
{
    for ( Couple * cx=set.firstAA(); cx ; cx=cx->next() )
    {
        // do not display Couple if the associated Fibers are both hidden
        if ( !cx->fiber1()->disp->visible && !cx->fiber2()->disp->visible )
            continue;
        
        // only display if bridging two anti-parallel filaments
        if ( prop->couple_select & 8  && cx->cosAngle() > 0 )
            continue;
        
        // only display if bridging two parallel filaments
        if ( prop->couple_select & 16 && cx->cosAngle() < 0 )
            continue;
        
        drawCoupleB(cx);
    }
}

//------------------------------------------------------------------------------
#pragma mark -

void Display::drawSolids(SolidSet const& set)
{
    for ( Solid * obj = set.first(); obj; obj=obj->next() )
    {
        if ( obj->prop->disp->visible )
        {
            drawSolid(*obj);
#if ( DIM == 3 )
            if ( obj->prop->disp->color.transparent() )
            {
                for ( size_t ii = 0; ii < obj->nbPoints(); ++ii )
                    zObjects.push_back(zObject(obj, ii));
            }
            else
#endif
            {
                for ( size_t ii = 0; ii < obj->nbPoints(); ++ii )
                    drawSolidT(*obj, ii);
            }
        }
    }
}


void Display::drawBeads(BeadSet const& set)
{
    for ( Bead * obj = set.first(); obj; obj=obj->next() )
    {
        if ( obj->prop->disp->visible )
        {
            drawBead(*obj);
#if ( DIM == 3 )
            if ( obj->prop->disp->color.transparent() )
                zObjects.push_back(zObject(obj));
            else
#endif
                drawBeadT(*obj);
        }
    }
}


void Display::drawSpheres(SphereSet const& set)
{
    for ( Sphere * obj=set.first(); obj ; obj=obj->next() )
    {
        if ( obj->prop->disp->visible )
        {
            drawSphere(*obj);
#if ( DIM == 3 )
            if ( obj->prop->disp->color.transparent() )
                zObjects.push_back(zObject(obj));
            else
#endif
                drawSphereT(*obj);
        }
    }
}


void Display::drawOrganizers(OrganizerSet const& set)
{
    for ( Organizer * obj=set.first(); obj ; obj=obj->next() )
        drawOrganizer(*obj);
}


//------------------------------------------------------------------------------
#pragma mark -


/// display sub-part `inx` of object `obj`
void Display::zObject::draw(Display * disp) const
{
    Mecable const * mec = point_.mecable();
    switch( mec->tag() )
    {
        case Fiber::TAG:
            disp->drawFiberSegmentT(*static_cast<const Fiber*>(mec), point_.point());
            break;
            
        case Solid::TAG:
            disp->drawSolidT(*static_cast<const Solid*>(mec), point_.point());
            break;
            
        case Bead::TAG:
            disp->drawBeadT(*static_cast<const Bead*>(mec));
            break;
            
        case Sphere::TAG:
            disp->drawSphereT(*static_cast<const Sphere*>(mec));
            break;
            
        default:
            std::cerr << "Internal error: unknown zObject pointer" << std::endl;
    }
}


#if ( DIM >= 3 )

/**
 This display objects in `zObjects` from back to front

 Depth-sorting is used in 3D to display transparent objects
 from the furthest to the nearest.
*/
void Display::drawTransparentObjects(Array<zObject>& list)
{
    // get current modelview transformation:
    GLfloat mat[16];
    glGetFloatv(GL_MODELVIEW_MATRIX, mat);
    
    // extract axis corresponding to vertical direction:
    Vector3 vertical(mat[2], mat[6], mat[10]);
    
    for ( zObject & i : list )
        i.depth(dot(i.position(), vertical));
    
    // depth-sort objects:
    list.sort(compareZObject);

    /*
     Enable polygon offset to avoid artifacts with objects of same size,
     particularly the ends of filaments with their tubular shaft.
     */
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1.0, 1.0);

    for ( zObject const& i : list )
        i.draw(this);
    
    glDisable(GL_POLYGON_OFFSET_FILL);
}

#endif


