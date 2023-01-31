// Cytosim was created by Francois Nedelec. Copyright Cambridge University 2021

#include "display.h"
#include "organizer.h"
#include "property_list.h"
#include "hand_prop.h"
#include "sphere_prop.h"
#include "fiber_prop.h"
#include "random_pcg.h"

#define DISPLAY 1
#include "gym_color.h"
#include "gym_color_list.h"

#include "modulo.h"
#include "simul.h"
#include "field.h"

#include "gle.h"
#include "gym_draw.h"
#include "gym_view.h"
#include "gym_flute.h"
#include "gym_flute_dim.h"
#include "gym_check.h"
#include "fg_stroke.h"

#include "point_disp.h"
#include "fiber_disp.h"
#include "line_disp.h"
#include "display_color.cc"

extern Modulo const* modulo;


Display::Display(DisplayProp const* dp)
: pixelSize(1), unitValue(1), sizeScale(1), depthAxis(0,0,1), prop(dp)
{
    assert_true(dp);
    age_start = 0;
    age_scale = 1.0;
    allLineDisp = nullptr;
    numLineDisp = 0;
}

void Display::setPixelFactors(float ps, float uv)
{
    pixelSize = ps;
    unitValue = uv;
    /*
     the 0.5 below comes from the fact that glPointSize uses diameter
     while most gle::primitives have a radius of 1
     */
    sizeScale = 0.5f * uv * ps;
    //printf(" pixelSize %6.3f unitValue %6.3f : %6.3f\n", ps, uv, sizeScale);
}

Display::~Display()
{
    delete[] allLineDisp;
}


void Display::drawText(Vector const& vec, const char str[]) const
{
#if ( DIM == 3 )
    gym::face_view(vec.XX, vec.YY, vec.ZZ);
#elif ( DIM == 2 )
    gym::face_view(vec.XX, vec.YY, 0.f);
#else
    gym::face_view(vec.XX, 0.f, 0.f);
#endif
    fgStrokeString(0, 0, pixelSize, 1, str, 1);
}

//------------------------------------------------------------------------------
#pragma mark - Couple::disp()

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


//------------------------------------------------------------------------------
#pragma mark - drawObject

void Display::drawObject(Vector const& pos, float rad, void(*obj)()) const
{
    gym::transScale(pos, rad);
    obj();
}


void Display::drawObject(Vector const& pos, Vector const& dir, float rad, void(*obj)()) const
{
    gym::transAlignZ(pos, rad, dir);
    obj();
}


void Display::drawBallT(Vector const& pos, real rad, gym_color const& col, unsigned mark) const
{
    gym::transScale(pos, rad);
    //gym::enableLighting();
    gym::color_both(col);
    gle::dualPassSphere2();
    if ( mark )
    {
        gym::color_front(0,0,0);
        gle::footballPentagons();
    }
}


void Display::drawDiscT(Vector const& pos, real rad, gym_color const& col) const
{
    gym::transScale(pos, rad);
    gym::disableLighting();
    gym::color(col);
    gle::disc();
}


/// used for drawFilament
inline void drawMonomer(Vector3 const& pos, float rad)
{
    gym::transScale(pos, rad);
    gle::sphere2();
}

//------------------------------------------------------------------------------
#pragma mark -

/** This is only one version of the display function, see display1.cc, etc. */
void Display::drawObjects(Simul const& sim)
{
    gym::closeDepthMask();
    gym::disableLighting();
    gym::disableCullFace();
    drawFields(sim.fields);
    
    gym::enableLighting();
    gym::enableCullFace(GL_BACK);
    gym::openDepthMask();
    drawSpaces(sim.spaces);
    
    gym::disableCullFace();

    /**
     If the display is 'cut', we might see the inner sides,
     but rendering would be faster with Culling enabled
     */
    //gym::enableCullFace(GL_BACK);
    drawFibers(sim.fibers);
    
    gym::enableLighting();
    gym::enableCullFace(GL_BACK);

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
    
    drawOrganizers(sim.organizers);
    gym::disableCullFace();
    drawMisc(sim);
}


void Display::drawSimul(Simul const& sim)
{
    // clear list of transparent objects
    zObjects.clear();
    CHECK_GL_ERROR("at drawSimul()");

#if ( DIM >= 3 )
    // Draw opaque objects with writable depth buffer
    drawObjects(sim);
    // Draw transparent objects with read-only depth buffer
    gym::closeDepthMask();
    drawTransparentObjects(sim);
    gym::openDepthMask();
#else
    drawObjects(sim);
#endif
    
    CHECK_GL_ERROR("in Display::drawSimul()");
}


/// qsort function comparing the 4th component of two vectors
static int compareVector4(const void * a, const void * b)
{
    real az = ((Vector4 const*)(a))->TT;
    real bz = ((Vector4 const*)(b))->TT;
    return ( az > bz ) - ( bz > az );
}

/**
 To get correct display, it would be necessary to display all opaque objects first,
 and then all transparent objects for all tiles.
 However, we call Display::drawSimul() a number of times,
 and objects are only sorted within each tile.
 So we depth-sort the views, but the result is still imperfect.
 */
void Display::drawTiled(Simul const& sim, int tile)
{
    assert_true(modulo);
    
    int l[3] = { 0 };
    int u[3] = { 0 };
    
    for ( size_t d = 0; d < DIM; ++d )
    {
        if ( modulo->isPeriodic(d) )
        {
            l[d] = -((tile>>d)&1);
            u[d] = +1;
        }
    }
    
    const Vector px = modulo->period(0);
    const Vector py = modulo->period(1);
    const Vector pz = modulo->period(2);

    Vector4 pos[32];
    int cnt = 0;
    
    for ( int dx = l[0]; dx <= u[0]; ++dx )
    for ( int dy = l[1]; dy <= u[1]; ++dy )
    for ( int dz = l[2]; dz <= u[2]; ++dz )
    {
        Vector P = dx * px + dy * py + dz * pz;
        pos[cnt] = Vector4(P);
        pos[cnt++].TT = dot(depthAxis, Vector3(P));
    }
    
    // depth-sort positions:
    qsort(pos, cnt, sizeof(Vector4), &compareVector4);

    float ref[16], mat[16];
    gym::get_view(ref);
    for ( int i = 0; i < cnt; ++i )
    {
        gym::mat_translate(mat, ref, pos[i].XX, pos[i].YY, pos[i].ZZ);
        gym::set_view(mat);
        drawSimul(sim);
    }
}

//------------------------------------------------------------------------------
#pragma mark -


/**
 Create a FiberDisp for this Property if necessary
 */
void Display::prepareFiberDisp(FiberProp* fp, PropertyList& alldisp, gym_color col)
{
    FiberDisp *& disp = fp->disp;
    
    // recover existing property:
    if ( !disp )
        disp = static_cast<FiberDisp*>(alldisp.find("fiber:display", fp->name()));

    // create new property with default values:
    if ( !disp )
    {
        disp = new FiberDisp(fp->name());
        alldisp.push_back(disp);
        // set default:
        disp->color      = col;
        disp->back_color = col.darken(0.5);
        disp->point_size = prop->point_size;
        disp->line_width = prop->line_width;
    }
    
    // parse user-provided values:
    if ( fp->display_fresh )
    {
        disp->read_string(fp->display, fp->name()+":display");
        fp->display_fresh = false;
    }
    
    if ( disp->coloring == FiberDisp::COLORING_CLUSTER )
        prep_flag |= 1;
    
    if ( disp->line_style == 2 || disp->line_style == 3 )
        prep_flag |= 2;

    if ( disp->coloring == FiberDisp::COLORING_AGE )
        prep_flag |= 4;
    
    disp->setPixels(pixelSize, unitValue);
}


/**
 set LineDisp for given Fiber
 */
void Display::prepareLineDisp(const Fiber * fib, FiberDisp const* disp, LineDisp * self)
{
    bool hide = false;
    assert_true(fib->prop);
    gym_color col = disp->color;
    
    // change body color depending on coloring mode:
    switch ( disp->coloring )
    {
        default:
        case FiberDisp::COLORING_OFF:
            col = disp->color;
            break;
        case FiberDisp::COLORING_RANDOM:
            col = gym::bright_color(fib->signature()).match_a(disp->color);
            break;
        case FiberDisp::COLORING_DIRECTION:
            col = radial_color(fib->direction(), disp->color.alpha());
            break;
        case FiberDisp::COLORING_MARK:
            col = gym::nice_color(fib->mark());
            break;
        case FiberDisp::COLORING_FLAG:
            col = gym::std_color(fib->flag());
            break;
#if FIBER_HAS_FAMILY
        case FiberDisp::COLORING_FAMILY:
            if ( fib->family_ )
                col = gym::nice_color(fib->family_->signature());
            else
                col = disp->color;
            break;
#endif
        case FiberDisp::COLORING_CLUSTER:
            col = gym::alt_color(fib->flag());
            break;
        case FiberDisp::COLORING_AGE:
            col = gym_color::jet_color((fib->age()-age_start)*age_scale, 1.0);
            break;
        case FiberDisp::COLORING_PSTATE:
            if ( fib->endStateP() > 0 )
                col = disp->end_colors[std::min(fib->endStateP(),5U)];
            break;
    }
    
    self->color = col;
    self->color_scale = color_scale(fib, disp->line_style);
    //std::cerr << fib->reference() << ":color_scale " << self->color_scale << "\n";
    
#if ( 1 )
    // colors of ends set to match body color:
    self->end_color[0] = col;
    self->end_color[1] = col;
#else
    // colors of ends for non-dynamic filaments:
    self->end_color[0] = disp->end_colors[0];
    self->end_color[1] = disp->end_colors[0];
#endif
    
    // For dynamic Fibers, change colors of tips according to state:
    if ( fib->endStateP() > 0 )
        self->end_color[0] = disp->end_colors[std::min(fib->endStateP(),5U)];
    
    if ( fib->endStateM() > 0 )
        self->end_color[1] = disp->end_colors[std::min(fib->endStateM(),5U)];

    // hide right or left-pointing fibers:
    if (( disp->hide & 1 )  &&  dot(fib->diffPoints(0), Vector(disp->hide_axis)) < 0 )
        hide = true;
    if (( disp->hide & 2 )  &&  dot(fib->diffPoints(0), Vector(disp->hide_axis)) > 0 )
        hide = true;
    
#if ( DIM == 2 )
    // hide clockwise or counter-clockwise orientated fibers:
    if (( disp->hide & 4 )  &&  cross(fib->posP(0), fib->diffPoints(0)) < 0 )
        hide = true;
    if (( disp->hide & 8 )  &&  cross(fib->posP(0), fib->diffPoints(0)) > 0 )
        hide = true;
#elif ( DIM >= 3 )
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
    if ( disp->explode_style )
        self->explode_shift = ( fib->signature() * 0x1p-32 - 0.5 ) * disp->explode_range;
    else
        self->explode_shift = 0;
}


/**
 Create a PointDisp for this Property if necessary
 */
template < typename T >
void Display::preparePointDisp(T * p, PropertyList& alldisp, gym_color col)
{
    PointDisp *& disp = p->disp;
    
    // create new property:
    if ( !disp )
    {
        // search for matching property:
        disp = static_cast<PointDisp*>(alldisp.find(p->category()+":display", p->name()));
        if ( !disp )
        {
            //std::clog <<" new " << p->category() << ":display " << p->name() << "\n";
            disp = new PointDisp(p->category()+":display", p->name());
            disp->clear();
            alldisp.push_back(disp);
            // set default:
            disp->color  = col;
            disp->color2 = col.alpha_scaled(0.25f);
            disp->size   = prop->point_size;
            if ( p->category() == "hand" )
                disp->width = prop->link_width;
            else
                disp->width = prop->line_width;
        }
    }
    
    // parse display string once:
    if ( p->display_fresh )
    {
        disp->read_string(p->display, p->name()+":display");
        p->display_fresh = false;
    }
    
    disp->setPixels(pixelSize, unitValue, prop->style==2);
}


/// attribute LineDisp to all fibers, and set individual display values
void Display::attributeLineDisp(FiberSet const& fibers)
{
    if ( numLineDisp < fibers.size() )
    {
        constexpr size_t chunk = 32;
        numLineDisp = ( fibers.size() + chunk ) & ~ ( chunk - 1 );
        delete[] allLineDisp;
        allLineDisp = new LineDisp[numLineDisp];
#if 0
        printf("sizeof gym_color %lu bytes\n", sizeof(gym_color));
        printf("sizeof LineDisp  %lu bytes\n", sizeof(LineDisp));
        std::clog << " new allLineDisp(" << numLineDisp << ")\n";
#endif
    }
    size_t i = 0;
    for ( Fiber * fib = fibers.first(); fib; fib = fib->next() )
    {
        fib->disp = &allLineDisp[i++];
        prepareLineDisp(fib, fib->prop->disp, fib->disp);
    }
    assert_true( i <= numLineDisp );
}


/**
 Perform the operations that are necessary to display the simulation:
 - create FiberDisp, HandDisp, SphereDisp, etc. (one per Property)
 - create LineDisp (one per Fiber)
 - set default values,
 - parse display strings
 .
*/
void Display::prepareDrawing(Simul const& sim, PropertyList& alldisp, Vector3 const& axis)
{
    depthAxis = axis;
    
    // counter to give different colors to the objects
    size_t idx = 0;

    PropertyList plist = sim.properties.find_all("fiber");
    
    prep_flag = 0;
    // create a FiberDisp for each FiberProp:
    for ( Property* p : plist )
        prepareFiberDisp(static_cast<FiberProp*>(p), alldisp, gym::nice_color(idx++));

    // create a LineDisp for each Fiber:
    attributeLineDisp(sim.fibers);

    if ( prep_flag )
    {
        // compute clusters:
        if ( prep_flag & 1 )
            sim.flagClusters(1, 1, 0);
        
        // if fiber tensions are used for display, recompute them now:
        if ( prep_flag & 2 )
            sim.computeForces();
        
        // calculate Fiber::age() range and set color scaling factor:
        if ( prep_flag & 4 )
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
    
    //create a PointDisp for each HandProp:
    for ( Property * i : sim.properties.find_all("hand") )
        preparePointDisp(static_cast<HandProp*>(i), alldisp, gym::nice_color(idx++));
    
    //create a PointDisp for each SphereProp:
    for ( Property * i : sim.properties.find_all("sphere") )
        preparePointDisp(static_cast<SphereProp*>(i), alldisp, gym::bright_color(idx++));
    
    //create a PointDisp for each SolidProp:
    for ( Property * i : sim.properties.find_all("solid", "bead") )
        preparePointDisp(static_cast<SolidProp*>(i), alldisp, gym::bright_color(idx++));
    
    //create a PointDisp for each SpaceProp:
    gym_color col(DIM==3?0x00000044:0xAAAAAAFF);
    for ( Property * i : sim.properties.find_all("space") )
        preparePointDisp(static_cast<SpaceProp*>(i), alldisp, col);
}


//------------------------------------------------------------------------------
#pragma mark -


/**
 Draw the back and front sides of Spaces in 3D
 This function is called twice: at the start and at the end of drawSimul()
 */
void Display::drawSpace3D(Space const* obj, bool back)
{
    const PointDisp * disp = obj->prop->disp;
    bool front = back ^ ( disp->color.transparent() );
    
    back = ( disp->visible & 2 ) && back;
    front = ( disp->visible & 1 ) && front;
    
    if ( back | front )
    {
        gym::ref_view();
        gym::enableLighting();
        gym::enableCullFace(GL_FRONT);
        if ( back )
        {
            gym::color_back(disp->color2);
            obj->draw3D();
        }
        if ( front )
        {
            gym::switchCullFace(GL_BACK);
            gym::color_load(disp->color);
            obj->draw3D();
        }
        gym::restoreCullFace();
    }
}


void Display::drawSpaces(SpaceSet const& set)
{
#if ( DIM >= 3 )
    
    // draw non-transparent Spaces first:
    //for ( Space const* obj = set.first(); obj; obj=obj->next() )
    for ( Space const* obj = set.firstID(); obj; obj = set.nextID(obj) )
    {
        if ( obj->prop->disp->visible )
            drawSpace3D(obj, true);
    }

#else
    
    for ( Space * obj = set.first(); obj; obj=obj->next() )
    {
        const PointDisp * disp = obj->prop->disp;
        if ( disp->visible )
        {
            gym::disableLighting();
            gym::color(disp->color);
            gym::ref_view();
            obj->draw2D(disp->widthX);
        }
    }
    
#endif
    CHECK_GL_ERROR("in Display::drawSpaces()");
}


/**
 Draw transparent Spaces
 */
void Display::drawTransparentSpaces(SpaceSet const& set)
{
    for ( Space * obj = set.first(); obj; obj=obj->next() )
    {
        if ( obj->prop->disp->visible )
            drawSpace3D(obj, false);
    }
}


/**
 This displays only one Field, specified by DisplayProp:field_number
 
 CULL_FACE and LIGHTING should be disabled
 */
void Display::drawFields(FieldSet const& set)
{
#if ( DIM >= 3 )
    Vector3 dir = depthAxis;
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
    CHECK_GL_ERROR("at Display::drawFields()");
}


//------------------------------------------------------------------------------
#pragma mark -

void Display::drawAverageFiber(ObjectList const& objs, gym_color col) const
{
    Vector G, D, M, P;
    real S = FiberSet::infoPosition(objs, M, G, P);
    
    if ( S > REAL_EPSILON )
    {
        Vector MP = normalize( P - M );
        const float rad = pixscale(10);
        if ( 1 )
        {
            // a black outline
            float blk[4] = { 0, 0, 0, 1 };
            float RAD = rad * 1.1;
            gym::color_front(blk);
            gym::closeDepthMask();
            gle::drawCylinder(M, MP, RAD);
            gle::drawCone(P, MP, RAD);
            drawObject(G, RAD, gle::sphere2);
            gym::openDepthMask();
        }
        gym::color_front(col);
        gle::drawCylinder(M, MP, rad);
        gle::drawCone(P, MP, rad);
        drawObject(G, rad, gle::sphere2);
    }
}


bool selectR(Object const* obj, void const* arg)
{
    Fiber const* fib = static_cast<Fiber const*>(obj);
    return fib->prop==arg  &&  dot(fib->diffPoints(0), Vector(fib->prop->disp->hide_axis)) > 0;
}

bool selectL(Object const* obj, void const* arg)
{
    Fiber const* fib = static_cast<Fiber const*>(obj);
    return fib->prop==arg  &&  dot(fib->diffPoints(0), Vector(fib->prop->disp->hide_axis)) < 0;
}

void Display::drawAverageFiber1(FiberSet const& fibers, void const* arg) const
{
    ObjectList objs = fibers.collect(match_property, arg);
    drawAverageFiber(objs, gym_color(1,1,1));
}


void Display::drawAverageFiber2(FiberSet const& fibers, void const* arg) const
{
    ObjectList objsR = fibers.collect(selectR, arg);
    ObjectList objsL = fibers.collect(selectL, arg);
    
    // display right-pointing fibers in Red
    drawAverageFiber(objsR, gym_color(1,0,0));
    
    // display left-pointing fibers in Green
    drawAverageFiber(objsL, gym_color(0,1,0));
}    


void Display::drawMisc(Simul const& sim)
{
#if ( 0 )
    // display Grids for visual inspection:
    sim.fiberGrid.drawGrid();
    sim.sMeca.locusGrid.drawGrid();
    sim.sMeca.pointGrid.drawGrid();
#endif
#if ( 0 )
    for ( Property const* i : sim.properties.find_all("fiber") )
    {
        FiberProp const* fp = static_cast<FiberProp const*>(i);
        if ( fp->disp->draw_average == 1 )
            drawAverageFiber1(sim.fibers, fp);
        else if ( fp->disp->draw_average == 2 )
            drawAverageFiber2(sim.fibers, fp);
    }
#endif
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
void Display::drawFiberMinusEnd(Fiber const& fib, int style, float size) const
{
    const float rad = pixscale(size);
    if ( rad > pixelSize ) switch(style)
    {
        default: break;
        case 1: drawObject(fib.posEndM(), rad, gle::sphere2); break;
        case 2: gle::drawCone(fib.posEndM(), -fib.dirEndM(), rad); break;
        case 3: gle::drawCylinder(fib.posEndM(), fib.dirEndM(), rad); break;
        case 4: gle::drawArrowTail(fib.posEndM(), -fib.dirEndM(), rad); break;
        case 5: gle::drawArrowTail(fib.posEndM(), fib.dirEndM(), rad); break;
        case 6: drawObject(fib.posEndM(), -fib.dirEndM(), rad, gle::cube); break;
        case 7: gle::drawCylinder(fib.posEndM(), -fib.dirEndM(), rad); break;
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
void Display::drawFiberPlusEnd(Fiber const& fib, int style, float size) const
{
    gym::ref_view();
    const float rad = pixscale(size);
    if ( rad > pixelSize ) switch(style)
    {
        default: break;
        case 1: drawObject(fib.posEndP(), rad, gle::sphere2); break;
        case 2: gle::drawCone(fib.posEndP(), fib.dirEndP(), rad); break;
        case 3: gle::drawCylinder(fib.posEndP(), -fib.dirEndP(), rad); break;
        case 4: gle::drawArrowTail(fib.posEndP(), fib.dirEndP(), rad); break;
        case 5: gle::drawArrowTail(fib.posEndP(), -fib.dirEndP(), rad); break;
        case 6: drawObject(fib.posEndP(), fib.dirEndP(), rad, gle::cube); break;
        case 7: gle::drawCylinder(fib.posEndP(), fib.dirEndP(), rad); break;
    }
}


void Display::drawFiberBackbone(Fiber const& fib) const
{
    gym::ref_view();
    gym::disableLighting();
    gym::color(fib.disp->color);
    gym::loadPoints(fib.nbPoints(), fib.addrPoints());
    gym::drawLineStrip(fib.prop->disp->line_width, 0, fib.nbPoints());
    gym::cleanup();
}


void Display::drawFiberLines(Fiber const& fib, int style) const
{
    size_t cnt = 2 * fib.nbSegments();
    flute4D* flu = gym::mapBufferC4VD(cnt+4);
    flute4D* ptr = flu;
    bool strip = 1;
    
    switch ( style )
    {
        case 1: { // display plain lines:
            gym_color c = fib.disp->color;
            for ( size_t i = 0; i < fib.nbPoints(); ++i )
                flu[i] = {c, fib.posP(i)};
            ptr = flu + fib.nbPoints();
        } break;
        case 2: // display segments with color indicating internal tension
            strip = 0;
            for ( size_t n = 0; n < fib.nbSegments(); ++n )
            {
                gym_color c = color_by_tension(fib, n);
                ptr[0] = {c, fib.posP(n)};
                ptr[1] = {c, fib.posP(n+1)};
                ptr += 2;
            }
            break;
        case 3: // display segments with color indicating internal tension
            strip = 0;
            for ( size_t n = 0; n < fib.nbSegments(); ++n )
            {
                gym_color c = color_by_tension_jet(fib, n);
                ptr[0] = {c, fib.posP(n)};
                ptr[1] = {c, fib.posP(n+1)};
                ptr += 2;
            }
            break;
        case 4: // display segments with color indicating the curvature
            for ( size_t i = 0; i < fib.nbPoints(); ++i )
                flu[i] = {color_by_curvature(fib, i), fib.posP(i)};
            ptr = flu + fib.nbPoints();
            break;
        case 5: // color according to the angle with respect to the XY-plane:
            strip = 0;
            for ( size_t n = 0; n < fib.nbSegments(); ++n )
            {
                gym_color c = color_by_direction(fib, n);
                ptr[0] = {c, fib.posP(n)};
                ptr[1] = {c, fib.posP(n+1)};
                ptr += 2;
            }
            break;
        case 6: // color according to the distance from the minus end
            *ptr++ = {color_by_distanceM(fib, 0), fib.posP(0)};
            for ( real a = 0.0625; a < 0.6; a *= 2 )
                *ptr++ = {color_by_distanceM(fib, a), fib.midPoint(0, a)};
            for ( size_t n = 1; n < fib.nbPoints(); ++n )
                *ptr++ = {color_by_distanceM(fib, n), fib.posP(n)};
            break;
        case 7: { // color according to the distance from the plus end
            const size_t last = fib.lastPoint();
            for ( size_t n = 0; n < last; ++n )
                *ptr++ = {color_by_distanceP(fib, n), fib.posP(n)};
            for ( real a = 0.5; a > 0.06; a /= 2 )
                *ptr++ = {color_by_distanceP(fib, last-a), fib.midPoint(last-1, 1-a)};
            *ptr++ = {color_by_distanceP(fib, last), fib.posP(last)};
        } break;
        case 8: // color according to distance to the confining Space
            for ( size_t i = 0; i < fib.nbPoints(); ++i )
                flu[i] = {color_by_height(fib, i), fib.posP(i)};
            ptr = flu + fib.nbPoints();
            break;
    }
    gym::ref_view();
    gym::disableLighting();
    gym::unmapBufferC4VD();
    if ( strip )
        gym::drawLineStrip(fib.prop->disp->line_widthX, 0, ptr-flu);
    else
        gym::drawLines(fib.prop->disp->line_widthX, 0, ptr-flu);
    gym::cleanup();
}


void Display::drawFiberSegmentT(Fiber const& fib, size_t inx) const
{
    gym::ref_view();
    FiberDisp const*const disp = fib.prop->disp;
    size_t cnt = 8;
    flute4D* flu = gym::mapBufferC4VD(cnt);
    flute4D* ptr = flu;

    if ( disp->line_style == 6 )
    {
        // color by distance to Minus end
        *ptr++ = {color_by_distanceM(fib, inx), fib.posP(inx)};
        if ( inx == 0 )
        {
            for ( real dx = 0.125; dx < 0.6; dx *= 2 )
                *ptr++ = {color_by_distanceM(fib, dx), fib.midPoint(0, dx)};
        }
        *ptr++ = {color_by_distanceM(fib, inx+1), fib.posP(inx+1)};
    }
    else if ( disp->line_style == 7 )
    {
        // color by distance to Plus end
        *ptr++ = {color_by_distanceP(fib, inx), fib.posP(inx)};
        if ( inx == fib.lastSegment() )
        {
            for ( real dx = 0.5; dx > 0.12; dx /= 2 )
                *ptr++ = {color_by_distanceP(fib, inx+1-dx), fib.midPoint(inx, 1-dx)};
        }
        *ptr++ = {color_by_distanceP(fib, inx+1), fib.posP(inx+1)};
    }
    else
    {
        gym_color c;
        if ( disp->line_style == 2 )
            c = color_by_tension(fib, inx);
        else if ( disp->line_style == 3 )
            c = color_by_tension_jet(fib, inx);
        else
            c = fib.disp->color;
        // the whole segment is painted with the same color:
        ptr[0] = {c, fib.posP(inx)};
        ptr[1] = {c, fib.posP(inx+1)};
        ptr += 2;
    }
    gym::disableLighting();
    gym::unmapBufferC4VD();
    gym::drawLineStrip(disp->line_widthX, 0, ptr-flu);
    gym::restoreLighting();
    gym::cleanup();
}


void Display::drawFiberSpeckles(Fiber const& fib) const
{
    FiberDisp const*const disp = fib.prop->disp;
    const real gap = disp->speckle_gap;

    size_t i = 0, cnt = 8 + 4 * std::ceil(fib.length()/gap);
    fluteD* pts = gym::mapBufferVD(cnt);

    // display random speckles:
    if ( disp->speckle_style == 1 )
    {
        /*
         A simple random number generator seeded by fib.signature()
         is used to distribute points always at the same position
         with respect to the lattice of each fiber.
         */
        constexpr real TINY = 0x1p-32;
        // draw speckles below the origin of abscissa
        if ( fib.abscissaM() < 0 )
        {
            uint64_t Z = pcg32_init(fib.signature());
            real a = gap * std::log(pcg32(Z)*TINY);
            while ( a > fib.abscissaP() )
                a += gap * std::log(pcg32(Z)*TINY);
            while ( a >= fib.abscissaM() )
            {
                if ( i < cnt ) pts[i++] = fib.pos(a);
                a += gap * std::log(pcg32(Z)*TINY);
            }
        }
        // draw speckles above the origin of abscissa
        if ( fib.abscissaP() > 0 )
        {
            uint64_t Z = pcg32_init(~fib.signature());
            real a = -gap * std::log(pcg32(Z)*TINY);
            while ( a < fib.abscissaM() )
                a -= gap * std::log(pcg32(Z)*TINY);
            while ( a <= fib.abscissaP() )
            {
                if ( i < cnt ) pts[i++] = fib.pos(a);
                a -= gap * std::log(pcg32(Z)*TINY);
            }
        }
    }
    else if ( disp->speckle_style == 2 )
    {
        // display regular speckles
        real a = gap * std::ceil( fib.abscissaM() / gap );
        while ( a <= fib.abscissaP() )
        {
            if ( i < cnt ) pts[i++] = fib.pos(a);
            a += gap;
        }
    }
    gym::unmapBufferVD();
    gym::ref_view();
    gym::disableLighting();
    gym::color_load(fib.disp->color);
    gym::drawPoints(disp->speckle_size, 0, i);
}


void Display::drawFiberPoints(Fiber const& fib) const
{
    FiberDisp const*const disp = fib.prop->disp;
    int style = disp->point_style;

    if ( style == 1 )
    {
        // display vertices:
        gym::ref_view();
        gym::disableLighting();
        gym::color(fib.disp->color);
        gym::loadPoints(fib.nbPoints(), fib.addrPoints());
        gym::drawSquarePoints(disp->point_sizeX, 0, fib.nbPoints());
    }
    else if ( style == 2 )
    {
        gym::color_load(fib.disp->color);
        gym::color_back(disp->back_color);
        // display arrowheads along the fiber:
        const float rad = pixscale(disp->point_size);
        const real gap = disp->point_gap;
        real ab = std::ceil(fib.abscissaM()/gap) * gap;
        for ( ; ab <= fib.abscissaP(); ab += gap )
            gle::drawCone(fib.pos(ab), fib.dir(ab), rad);
    }
    else if ( style == 3 )
    {
        // display chevrons along the fiber:
        const real gap = disp->point_gap;
        const real rad = pixscale(disp->point_size);
        real beta = fib.segmentationInv() * rad;
        size_t cnt = 4 * fib.length() / gap + 8;
        fluteD* flu = gym::mapBufferVD(cnt);
        fluteD* ptr = flu;
        real a = std::ceil(fib.abscissaM()/gap) * gap - fib.abscissaM();
        Vector axis(depthAxis);
        for ( ; a <= fib.length(); a += gap )
        {
            Interpolation i = fib.interpolateM(a);
            Vector pos = i.pos();
            Vector dir = i.diff() * beta;
#if ( DIM == 3 )
            Vector nor = cross(axis, dir).normalized(rad);
#else
            Vector nor = dir.orthogonal(rad);
#endif
            ptr[0] = pos + dir - nor;
            ptr[1] = pos - dir;
            ptr[2] = pos - dir;
            ptr[3] = pos + dir + nor;
            ptr += 4;
        }
        assert_true(ptr-flu < cnt);
        gym::unmapBufferVD();
        gym::ref_view();
        gym::disableLighting();
        gym::color(fib.disp->color);
        gym::drawLines(disp->line_widthX, 0, ptr-flu);
    }
    else if ( style == 4 )
    {
        gym::color_load(fib.disp->color);
        gym::color_back(disp->back_color);
        // display only middle of fiber:
        drawObject(fib.posMiddle(), pixscale(2*disp->point_size), gle::sphere2);
    }
}


//------------------------------------------------------------------------------
#pragma mark - Lattice

/**
 This style uses one vertex for each site, positionned at the center of the range
 OpenGL will interpolate the colors, and each site will be covered by a gradient.
 */
void Display::drawFiberLattice1(Fiber const& fib, VisibleLattice const& lat, float rad) const
{
    rad = pixwidth(rad);
    const real uni = lat.unit();
    const auto inf = lat.indexM();
    const auto sup = lat.indexP();
    assert_true( inf <= sup );
    
    gym::ref_view();
    FiberDisp const*const disp = fib.prop->disp;
    gym_color c, col = disp->color;
    const real fac = 1 / disp->lattice_scale;
    size_t cnt = 2 * ( sup - inf );
    flute4D* flu = gym::mapBufferC4VD(cnt+4);
    flute4D* ptr = flu;
    
    if ( inf == sup )
    {
        //the Fiber is entirely covered by one site!
        real len = fib.abscissaP() - fib.abscissaM();
        c = lattice_color(col, (fac*lat.data(inf))*(uni/len));
        ptr[0] = {c, fib.posEndM()};
        ptr[1] = {c, fib.posEndP()};
        ptr += 2;
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
        c = lattice_color(col, facM*lat.data(inf));
        *ptr++ = {c, fib.posEndM()};
        if ( uni*(inf+0.5) > fib.abscissaM() )
            *ptr++ = {c, fib.pos(uni*(inf+0.5))};
        
        for ( auto h = inf+1; h < sup; ++h )
        {
            c = lattice_color(col, fac*lat.data(h));
            *ptr++ = {c, fib.pos(uni*(h+0.5))};
        }
        
        // the terminal site may be truncated
        c = lattice_color(col, facP*lat.data(sup));
        if ( uni*(sup+0.5) < fib.abscissaP() )
            *ptr++ = {c, fib.pos(uni*(sup+0.5))};
        *ptr++ = {c, fib.posEndP()};
    }
    gym::unmapBufferC4VD();
    gym::disableLighting();
    gym::drawLineStrip(rad, 0, ptr-flu);
    gym::cleanup();
}


/**
 This style, uses two vertices for each site, positionned at the extremity of the range,
 and each site is entirely covered by the color corresponding to the value.
 */
void Display::drawFiberLattice2(Fiber const& fib, VisibleLattice const& lat, float rad) const
{
    rad = pixwidth(rad);
    const real uni = lat.unit();
    const auto inf = lat.indexM();
    const auto sup = lat.indexP();
    assert_true( inf <= sup );
    
    gym::ref_view();
    FiberDisp const*const disp = fib.prop->disp;
    gym_color c, col = disp->color;
    const real fac = 1 / disp->lattice_scale;
    size_t cnt = 2 * ( sup - inf );
    flute4D* flu = gym::mapBufferC4VD(cnt+4);
    flute4D* ptr = flu;
    
    if ( inf == sup )
    {
        //the Fiber is entirely covered by one site!
        real len = fib.abscissaP() - fib.abscissaM();
        c = lattice_color(col, (fac*lat.data(inf))*(uni/len));
        ptr[0] = {c, fib.posEndM()};
        ptr[1] = {c, fib.posEndP()};
        ptr += 2;
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
        c = lattice_color(col, facM*lat.data(inf));
        *ptr++ = {c, fib.posEndM()};

        for ( auto h = inf+1; h < sup; ++h )
        {
            Vector P = fib.pos(uni*h);
            ptr[0] = {c, P};
            c = lattice_color(col, fac*lat.data(h));
            ptr[1] = {c, P};
            ptr += 2;
        }
        
        // the terminal site may be truncated
        c = lattice_color(col, facP*lat.data(sup));
        ptr[0] = {c, fib.pos(uni*sup)};
        ptr[1] = {c, fib.posEndP()};
        ptr += 2;
    }
    gym::unmapBufferC4VD();
    gym::disableLighting();
    gym::drawLines(rad, 0, ptr-flu);
    gym::cleanup();
}


void Display::drawFiberLattice3(Fiber const& fib, VisibleLattice const& lat, float rad) const
{
    drawFiberLattice2(fib, lat, rad);
    drawFiberLatticeEdges(fib, lat, rad*0.5f);
}


/**
 Indicate the edges between sites with small dots
 */
void Display::drawFiberLatticeEdges(Fiber const& fib, VisibleLattice const& lat, float rad) const
{
    rad = pixwidth(rad);
    const real uni = lat.unit();
    const auto inf = lat.indexM();
    const auto sup = lat.indexP();
    
    size_t cnt = sup - inf + 4;
    fluteD* flu = gym::mapBufferVD(cnt);
    fluteD* ptr = flu;
    real abs = (inf+1) * uni - fib.abscissaM();
    for ( auto h = inf+1; h <= sup; ++h, abs += uni )
        *ptr++ = fib.posM(abs);
    gym::unmapBufferVD();
    gym::ref_view();
    gym::disableLighting();
    gym::color(fib.disp->color);
    gym::drawPoints(rad, 0, ptr-flu);
    gym::cleanup();
}


void Display::drawFiberLabels(Fiber const& fib, int style, gym_color const& col) const
{
    char str[32];
    
    gym::ref_view();
    gym::color(col);
    gym::disableLighting();
    gym::disableAlphaTest();
    if ( style & 1 )
    {
        // draw fiber identity and vertex indices
        int C = snprintf(str, sizeof(str), " %u ", fib.identity());
        for ( size_t i = 0; i < fib.nbPoints(); ++i )
        {
            snprintf(str+C, sizeof(str)-C, "%lu", i);
            drawText(fib.posP(i), str);
        }
    } 
    else if ( style & 2 )
    {
        // draw fiber identity and abscissa value at vertices
        int C = snprintf(str, sizeof(str), " %u ", fib.identity());
        for ( size_t i = 0; i < fib.nbPoints(); ++i )
        {
            snprintf(str+C, sizeof(str)-C, "%.3f", fib.abscissaPoint(i));
            drawText(fib.posP(i), str);
        }
    }
    if ( style & 4 )
    {
        // display integral abscissa along the fiber
        snprintf(str, sizeof(str), "%.3f", fib.abscissaM());
        drawText(fib.posEndM(), str);
        
        int s = (int)std::ceil(fib.abscissaM());
        int e = (int)std::floor(fib.abscissaP());
        for ( int a = s; a <= e; ++a )
        {
            snprintf(str, sizeof(str), "%i", a);
            drawText(fib.pos(a), str);
        }
        
        snprintf(str, sizeof(str), "%.3f", fib.abscissaP());
        drawText(fib.posEndP(), str);
    }
    if ( style & 8 )
    {
        // indicate tension values in the segments
        Vector a = fib.posEndM();
        for ( size_t i = 1; i < fib.nbPoints(); ++i )
        {
            Vector b = fib.posP(i);
            snprintf(str, sizeof(str), "%+4.1f", fib.tension(i-1));
            drawText(0.5*(a+b), str);
            a = b;
        }
    }
    gym::restoreAlphaTest();
    CHECK_GL_ERROR("in Display::drawFiberLabels()");
}


/// display forces acting on the fiber's vertices, using lines scaled by 'mag'
void Display::drawFiberForces(Fiber const& fib, real mag, float size) const
{
    gym::ref_view();
    gym_color col = fib.prop->disp->force_color;
    gym_color lor = col.alpha_scaled(0.5f);
    size_t cnt = 2 * fib.nbPoints();
    flute4D* flu = gym::mapBufferC4VD(cnt);
    for ( size_t i = 0; i < fib.nbPoints(); ++i )
    {
        Vector P = fib.posP(i);
        Vector F = mag * fib.netForce(i);
        flu[  2*i] = { col, P };
        flu[1+2*i] = { lor, P+F };
    }
    gym::unmapBufferC4VD();
    gym::disableLighting();
    gym::drawLines(size, 0, cnt);
    gym::cleanup();
}

//------------------------------------------------------------------------------
#pragma mark - Specific styles


/**
 This renders a protofilament by drawing spheres of alternating colors,
 along the backbone of the `Fiber` at distance 4nm from each other.
 */
void Display::drawFilament(Fiber const& fib,
                           gym_color const& color1,
                           gym_color const& color2,
                           gym_color const& colorE) const
{
    // axial translation between two sucessive monomers:
    const real dab = 0.004;
    // enlarge radius of monomers to make them overlap
    const float rad = 0.65 * dab;
    
    real ab = 0;
    
    gym::enableClipPlane(4);
    
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
            gym::color_load(color1);
        else
            gym::color_load(color2);
        
        // change color for the last monomer:
        if ( ab + dab > fib.abscissaP() )
        {
            gym::color_load(colorE);
            gym::disableClipPlane(4);
        }
        
        // set clipping plane with the next monomer
        gym::setClipPlane(4, normalize(q-p), (p+q)*0.5);
        
        drawMonomer(q, rad);
        
        // set clipping plane with the previous:
        gym::setClipPlane(5, normalize(p-q), (p+q)*0.5);
        gym::enableClipPlane(5);
    }
    gym::disableClipPlane(4);
    gym::disableClipPlane(5);
}


/**
 This renders 26 spheres positionned on a right-handed helix,
 making one turn every 74nm, with a max width of ~ 9nm.
 This is roughly Ken Holmes' model of F-actin:
 Nature 347, 44 - 49 (06 September 1990); doi:10.1038/347044a0
 which shows half a turn in 37nm containing 13 monomers.
 */
void Display::drawActin(Fiber const& fib,
                        gym_color const& color1,
                        gym_color const& color2,
                        gym_color const& colorE) const
{    
    // axial translation between two sucessive monomers:
    const real dab = 0.00275;
    // enlarge radius of monomers to make them overlap
    const float rad = 1.3 * dab;
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
    
    gym::enableClipPlane(4);
    
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
            gym::color_load(color1);
        else
            gym::color_load(color2);

        // change color for the last monomer:
        if ( ab + dab > fib.abscissaP() )
        {
            gym::color_load(colorE);
            gym::disableClipPlane(4);
        }
        
        // set clipping plane with the next monomer
        gym::setClipPlane(4, normalize(q-p), (p+q)*0.5);
        
        drawMonomer(q, rad);
        
        // set clipping plane with the previous:
        gym::setClipPlane(5, normalize(p-q), (p+q)*0.5);
        
        gym::enableClipPlane(5);
    }
    gym::disableClipPlane(4);
    gym::disableClipPlane(5);
}


/**
 This renders a Microtubule using spheres of alternating colors
 colorA for alpha-tubulin
 colorB for beta-tubulin
 */
void Display::drawMicrotubule(Fiber const& fib,
                              gym_color const& colorA,
                              gym_color const& colorB,
                              gym_color const& colorE) const
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
    const float rad = 0.003;
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

        gym::color_load(colorA);
        for ( int i = 0; i < 13; ++i )
            drawMonomer(p+dx[i]*d+dy[i]*e+dz[i]*f, rad);

        gym::color_load(colorB);
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
                gym::color_load(colorA);
                drawMonomer(p+dx[i]*d+dy[i]*e+dz[i]*f, rad);
                if ( ab+5.2*sa+dx[i] < abmax )
                    gym::color_load(colorB);
                else
                    gym::color_load(colorE);
                drawMonomer(p+(sa+dx[i])*d+dy[i]*e+dz[i]*f, rad);
            }
        }
        ab += dab;
    }
}


int Display::drawFiberLattice(Fiber const& fib, int style, float rad) const
{
    VisibleLattice const* lat = fib.visibleLattice();
    if ( lat ) switch ( style )
    {
        case 1:
            drawFiberLattice1(fib, *lat, rad);
            return 1;
        case 2:
            drawFiberLattice2(fib, *lat, rad);
            return 1;
        case 3:
            drawFiberLattice3(fib, *lat, rad);
            return 1;
        case 4:
            drawFiberLatticeEdges(fib, *lat, rad);
            return 0;
    }
    return 0;
}


void Display::drawFiber(Fiber const& fib)
{
    FiberDisp const*const disp = fib.prop->disp;
    int style = disp->line_style;
    
    if ( disp->style )
    {
        if ( disp->style & 4 )
            return drawFiberBackbone(fib);

        gym_color col1 = fib.disp->color;
        gym_color col2 = fib.disp->color.darken(0.75);
        gym_color colE = fib.disp->end_color[0];
        
        // load backface color:
        if ( fib.prop->disp->coloring )
            gym::color_back(col1);
        else
            gym::color_back(fib.prop->disp->back_color);

        switch( disp->style )
        {
            case 1: drawFilament(fib, col1, col2, colE); break;
            case 2: drawActin(fib, col1, col2, colE); break;
            case 3: drawMicrotubule(fib, col1, col2, colE); break;
        }
        style = 0;
    }

#if FIBER_HAS_LATTICE || FIBER_HAS_MESH
    if ( disp->lattice_style )
    {
        // if the Lattice is displayed, do not draw backbone:
        if ( drawFiberLattice(fib, disp->lattice_style, disp->line_width) )
            style = 0;
    }
#endif
#if ( DIM >= 3 )
    /*
     Handle styles in 3D that are using transparency to draw fiber's segments
     */
    if (( style==1 && fib.disp->color.transparent())
        || ( style==2 || style==3 ))
    {
        for ( size_t i = 0; i < fib.lastPoint(); ++i )
            zObjects.push_back(zObject(&fib, i));
        style = 0;
    }
    else if ( style == 6 )
    {
        // color according to the distance from the minus end
        for ( size_t i = 0; i < fib.lastPoint(); ++i )
            if ( color_by_distanceM(fib, i).visible() )
                zObjects.push_back(zObject(&fib, i));
        style = 0;
    }
    else if ( style == 7 )
    {
        // color according to the distance from the plus end
        for ( size_t i = 0; i < fib.lastPoint(); ++i )
            if ( color_by_distanceP(fib, i+1).visible() )
                zObjects.push_back(zObject(&fib, i));
        style = 0;
    }
#endif

    if ( style )
        drawFiberLines(fib, style);

    if ( disp->point_style > 0 )
        drawFiberPoints(fib);
    
    if ( disp->speckle_style > 0 )
        drawFiberSpeckles(fib);

    // draw other fiber elements only if fiber is fully visible:
    //if ( fib.disp->visible > 0 )
    {
        if ( disp->label_style )
        {
            drawFiberLabels(fib, disp->label_style, fib.disp->color);
        }
        
        if ( disp->end_style[1] )
        {
            gym::color_load(fib.disp->end_color[1]);
            //gym::color_load(fib.disp->color);
            gym::color_back(disp->back_color);
            drawFiberMinusEnd(fib, disp->end_style[1], disp->end_size[1]);
        }
        
        if ( disp->end_style[0] )
        {
            gym::color_load(fib.disp->end_color[0]);
            //gym::color_load(fib.disp->color);
            gym::color_back(disp->back_color);
            drawFiberPlusEnd(fib, disp->end_style[0], disp->end_size[0]);
        }
        
        if ( disp->force_style )
        {
            drawFiberForces(fib, disp->force_scale, disp->point_sizeX);
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
    CHECK_GL_ERROR("in Display::drawFibers()");
}


//------------------------------------------------------------------------------
#pragma mark - Couples


void Display::drawCouplesB(CoupleSet const& set) const
{
    for ( Couple * cx=set.firstAA(); cx ; cx=cx->next() )
    {
        // only display if bridging two anti-parallel filaments
        if (( prop->couple_select & 8 ) && cx->cosAngle() > 0 )
            continue;
        
        // only display if bridging two parallel filaments
        if (( prop->couple_select & 16 ) && cx->cosAngle() < 0 )
            continue;
        
        // do not display Couple if the associated Fibers are both hidden
        bool vis1 = cx->fiber1()->disp->visible | cx->fiber2()->disp->visible;
        
        // do not display Couple if both hands are hidden
        bool vis2 = cx->disp1()->visible | cx->disp2()->visible;
        
        if ( vis1 & vis2 )
            drawCoupleB(cx);
    }
}

//------------------------------------------------------------------------------
#pragma mark - Solid

void Display::drawSolid(Solid const& obj)
{
    const PointDisp * disp = obj.prop->disp;
    
    //display points:
    if (( disp->style & 4 ) && disp->perceptible )
    {
        real rad = pixscale(disp->size);
        bodyColor(obj);
        gym::enableLighting();
        for ( size_t i = 0; i < obj.nbPoints(); ++i )
            drawObject(obj.posP(i), rad, gle::hedron(obj.radius(i)>0));
        gym::cleanup(1);

#if NEW_SOLID_HAS_TWIN
        // draw links between 'obj' and its Twin
        rad /= 2;
        const size_t inx = 0;
        Solid const* twi = obj.twin();
        if ( twi && obj.radius(inx) > 0 && obj.nbPoints() > inx + 3 )
        {
            gym::color(bodyColorF(obj));
            for ( size_t i = inx+1; i <= inx+DIM; ++i )
            {
                gym::stretchAlignZ(obj.posPoint(i), twi->posPoint(i), rad);
                gle::tube2();
            }
        }
#endif
    }
    
    //display outline of spheres in 2D
    if ( disp->style & 8 )
    {
#if ( DIM == 2 )
        gym::disableLighting();
        gym::color(bodyColorF(obj));
        for ( size_t i = 0; i < obj.nbPoints(); ++i )
        {
            if ( obj.radius(i) > pixelSize )
            {
                gym::transScale(obj.posP(i), obj.radius(i));
                gle::circle1(disp->width);
            }
        }
        gym::enableLighting();
#elif ( DIM >= 3 )
        //special display for ParM simulations (DYCHE 2006; KINETOCHORES 2019)
        if ( obj.mark()  &&  obj.nbPoints() > 1 )
        {
            gym::enableLighting();
            bodyColor(obj);
            //drawObject(obj.posP(0), obj.diffPoints(1, 0), obj.radius(0), gle::circle);
            gym::stretchAlignZ(obj.posP(0), obj.posP(1), obj.radius(0));
            gle::cylinder1();
        }
#endif
    }
    
    //print the number for each Solid
    if ( disp->style & 16 )
    {
        char tmp[8];
        gym::color(bodyColorF(obj));
        snprintf(tmp, sizeof(tmp), "0:%u", obj.identity());
        drawText(obj.posP(0), tmp);
        for ( size_t i = 1; i < obj.nbPoints(); ++i )
        {
            snprintf(tmp, sizeof(tmp), "%lu", i);
            drawText(obj.posP(i), tmp);
        }
    }
    
    //draw polygon line joining vertices of Solid
    if ( disp->style & 32 )
    {
        gym::disableLighting();
        gym::color(bodyColorF(obj));
        gym::ref_view();
        gym::loadPoints(obj.nbPoints(), obj.addrPoints());
        gym::drawLineStrip(disp->widthX, 0, obj.nbPoints());
        gym::enableLighting();
    }
}

/**
 Display a semi-transparent disc / sphere
 */
void Display::drawSolidT(Solid const& obj, size_t inx) const
{
    Vector X = obj.posP(inx);
    // using clipping planes to cleanup overlapping Spheres
    size_t near[3];
    size_t num = obj.closestSpheres(inx, near[0], near[1], near[2]);
    //printf("nearest Spheres to %lu / %lu are %lu %lu %lu\n", inx, obj.nbPoints(), near[0], near[1], near[2]);
    // set clipping planes with nearest Spheres
    for ( size_t i = 0; i < num; ++i )
    {
        size_t J = near[i];
        Vector P = obj.posP(J);
        real A = ( square(obj.radius(inx)) - square(obj.radius(J)) ) / distanceSqr(X, P);
        gym::enableClipPlane(5-i);
        gym::setClipPlane(5-i, normalize(X-P), (0.5-0.5*A)*X+(0.5+0.5*A)*P);
    }
#if ( DIM > 2 )
    drawBallT(X, obj.radius(inx), bodyColorF(obj), obj.mark()*(inx==0));
#else
    drawDiscT(X, obj.radius(inx), bodyColorF(obj));
#endif
    for ( size_t i = 0; i < num; ++i )
        gym::disableClipPlane(5-i);
}


static void drawFootball(Solid const& obj, size_t inx, gym_color col, gym_color bak, bool flip)
{
    Vector X = obj.posP(inx);
#if ( DIM >= 3 )
    Vector A = obj.posP(inx+1) - X;
    Vector B = obj.posP(inx+2) - X;
    Vector C = obj.posP(inx+3) - X;
    gym::transRotate(X, A, B, C);
    //flip = ( dot(cross(A,B), C) < 0 );
#else
    gym::transScale(X, obj.radius(inx));
#endif
    if ( flip ) glFrontFace(GL_CW);
    gym::color_front(col, 1.0);
    gle::sphere1();
    gym::color_front(bak);
    gle::footballPentagons();
    if ( flip ) glFrontFace(GL_CCW);
}


void Display::drawSolids(SolidSet const& set)
{
    for ( Solid * obj = set.first(); obj; obj=obj->next() )
    {
        const PointDisp * disp = obj->prop->disp;
        if ( disp->visible )
        {
            drawSolid(*obj);
            if ( disp->style & 1 )
            {
#if ( DIM >= 3 )
#if NEW_SOLID_HAS_TWIN
                const size_t inx = 0;
                Solid const* twi = obj->twin();
                if ( twi && obj->radius(inx) > 0 && obj->nbPoints() > inx + 3 )
                {
                    gym::enableLighting();
                    //real len = std::sqrt(obj->twinTensionSqr());
                    //gym_color col = gym_color::dark_jet_color(disp->scale * len);
                    gym_color col = bodyColorF(*obj);
                    gym_color black(0,0,0,1);
                    drawFootball(*obj, inx, col, black, false);
                    drawFootball(*twi, inx, col, black, true);
                }
#endif
                if ( obj->prop->disp->color.transparent() )
                {
                    for ( size_t i = 0; i < obj->nbPoints(); ++i )
                        if ( obj->radius(i) > 0 )
                            zObjects.push_back(zObject(obj, i));
                }
                else
#endif
                {
                    for ( size_t i = 0; i < obj->nbPoints(); ++i )
                        if ( obj->radius(i) > 0 )
                            drawSolidT(*obj, i);
                }
            }
            gym::cleanup(1);
        }
    }
}

//------------------------------------------------------------------------------
#pragma mark - Beads

void Display::drawBead(Bead const& obj)
{
    const PointDisp * disp = obj.prop->disp;

    // display center:
    if ( disp->style & 2 )
    {
        bodyColor(obj);
        drawObject(obj.position(), pixscale(disp->size), gle::tetrahedron);
    }
    
#if ( DIM == 2 )
    // display outline:
    if (( disp->style & 4 ) && obj.radius() > pixelSize )
    {
        gym::disableLighting();
        gym::color(bodyColorF(obj));
        gym::transScale(obj.position(), obj.radius());
        gle::circle1(disp->width);
    }
#endif
}


/**
 Display a semi-transparent disc / sphere
 */
void Display::drawBeadT(Bead const& obj) const
{
    const PointDisp * disp = obj.prop->disp;
    
    if ( disp->style & 1 )
    {
        gym_color col = bodyColorF(obj);
#if ( DIM > 2 )
        drawBallT(obj.position(), obj.radius(), col, obj.mark());
#else
        drawDiscT(obj.position(), obj.radius(), col);
#endif
    }
}


void Display::drawBeads(BeadSet const& set)
{
    gym::enableLighting();
    for ( Bead * obj = set.first(); obj; obj=obj->next() )
    {
        if ( obj->prop->disp->visible )
        {
            drawBead(*obj);
#if ( DIM >= 3 )
            if ( obj->prop->disp->color.transparent() )
                zObjects.push_back(zObject(obj));
            else
#endif
                drawBeadT(*obj);
        }
    }
}


//------------------------------------------------------------------------------
#pragma mark - Sphere

void Display::drawSphere(Sphere const& obj)
{
    const PointDisp * disp = obj.prop->disp;
    const float rad = pixscale(disp->size);
    
    if ( rad > pixelSize )
    {
        // display surface points
        if ( disp->style & 2 )
        {
            bodyColor(obj);
            for ( size_t i = obj.nbRefPoints; i < obj.nbPoints(); ++i )
                drawObject(obj.posP(i), rad, gle::sphere1);
        }
        
        // display center and reference points
        if ( disp->style & 8 )
        {
            bodyColor(obj);
            drawObject(obj.posP(0), rad, gle::star);
            for ( size_t i = 0; i < obj.nbRefPoints; ++i )
                drawObject(obj.posP(i), rad, gle::cube);
        }
    }
}


void Display::drawSphereT(Sphere const& obj) const
{
    const PointDisp * disp = obj.prop->disp;

    if ( disp->style & 7 )
    {
        const Vector C = obj.posP(0);
#if ( DIM < 3 )
        if ( obj.radius() > pixelSize )
        {
            gym::color(bodyColorF(obj));
            gym::transScale(C, obj.radius());
            if ( disp->style & 1 )
                gle::circle1(disp->widthX);
            if ( disp->style & 2 )
                gle::disc();
            if ( disp->style & 4 )
                drawDiscT(C, obj.radius(), bodyColorF(obj));
        }
#else
        /* Note: The rotation matrix for the sphere calculated below from the
         reference points, includes scaling by the radius of the sphere.
         We then use a primitive for a sphere of radius 1.
         */
        gym::color_both(bodyColorF(obj));
        gym::enableLighting();
        gym::transRotate(C, obj.posP(1)-C, obj.posP(2)-C, obj.posP(3)-C);
        if ( disp->style & 1 )
            gle::dualPassSphere2();
        if ( disp->style & 4 )
            gle::threeArrowStrip(0.5, 1);
#endif
    }
}


void Display::drawSpheres(SphereSet const& set)
{
    for ( Sphere * obj=set.first(); obj ; obj=obj->next() )
    {
        if ( obj->prop->disp->visible )
        {
            if ( obj->prop->disp->perceptible )
                drawSphere(*obj);
#if ( DIM >= 3 )
            if ( obj->prop->disp->color.transparent() )
                zObjects.push_back(zObject(obj));
            else
#endif
                drawSphereT(*obj);
        }
    }
}


//------------------------------------------------------------------------------
#pragma mark - Organizers

void Display::drawOrganizer(Organizer const& obj) const
{
    PointDisp const* disp = obj.disp();
    size_t i = 0, cnt = obj.nbLinks();

    if ( disp && ( disp->style & 2 ))
    {
        Vector P, Q;
        fluteD* flu = gym::mapBufferVD(2*cnt);
        while ( obj.getLink(i, P, Q) )
        {
            if ( modulo ) modulo->fold(Q, P);
            flu[  2*i] = Q;
            flu[1+2*i] = P;
            if ( ++i >= cnt ) break;
        }
        gym::unmapBufferVD();
        gym::ref_view();
        gym::disableLighting();
        gym::color(bodyColorF(disp, obj.signature()));
        gym::drawLines(disp->widthX, 0, 2*i);
        gym::rebindBufferVD(2);
        gym::drawPoints(disp->sizeX, 0, i);
    }

    /**
     This display the Solid connecting two Aster as a spindle.
     Used for Cleo Kozlowski simulation of C. elegans (2007)
     */
    if ( disp && ( disp->style & 1 ) && obj.tag() == Organizer::TAG_FAKE )
    {
        Solid const* sol = Solid::toSolid(obj.core());
        if ( sol && sol->nbPoints() >= 4 )
        {
#if ( DIM >= 3 )
            gym::enableLighting();
            bodyColor(*sol);
            Vector3 a = 0.5*(sol->posP(0) + sol->posP(2));
            Vector3 b = 0.5*(sol->posP(1) + sol->posP(3));
            gym::stretchAlignZ(a, b, 1);
            gle::dualPassBarrel();
#else
            gym::disableLighting();
            gym::color(bodyColorF(*sol));
            gym::ref_view();
            gym::loadPoints(sol->nbPoints(), sol->addrPoints());
            gym::drawLines(disp->widthX, 0, sol->nbPoints());
#endif
        }
    }
}


void Display::drawOrganizers(OrganizerSet const& set)
{
    for ( Organizer * obj=set.first(); obj ; obj=obj->next() )
        drawOrganizer(*obj);
}


//------------------------------------------------------------------------------
#pragma mark - Display of transparent objects sorted by decreasing depth

#if ( DIM >= 3 )

/// display sub-part `inx` of object `obj`
void zObject::draw(Display const* disp) const
{
    assert_false( point_.bad() );
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
            std::cerr << "Internal error: unknown zObject\n";
    }
}


/// qsort function comparing the zObjects::depth()
static int compareZObject(const void * A, const void * B)
{
    real a = static_cast<const zObject*>(A)->depth();
    real b = static_cast<const zObject*>(B)->depth();
    return ( a > b ) - ( a < b );
}

/**
 This display objects in `zObjects` from back to front

 Depth-sorting is used in 3D to display transparent objects
 from the furthest to the nearest.
*/
void Display::drawTransparentObjects(Array<zObject>& list)
{
    for ( zObject & i : list )
        i.depth(dot(i.position(), depthAxis));
    
    // depth-sort objects:
    list.sort(compareZObject);

    gym::enableLighting();
    for ( zObject const& i : list )
        i.draw(this);
}


/**
 Draw translucent objects:
 - make depth buffer readable only
 - objects are depth-sorted, from far to near
 - Dual pass is used to display back before front
 */
void Display::drawTransparentObjects(Simul const& sim)
{
    if ( zObjects.size() )
    {
        /*
         Enable polygon offset to avoid artifacts with objects of same size,
         particularly the ends of filaments with their tubular shaft.
         */
        gym::enableCullFace(GL_BACK);
        glEnable(GL_POLYGON_OFFSET_FILL);
        glPolygonOffset(1.0, 1.0);
        drawTransparentObjects(zObjects);
        glDisable(GL_POLYGON_OFFSET_FILL);
        gym::restoreCullFace();
    }
    
    drawTransparentSpaces(sim.spaces);
}

#endif


