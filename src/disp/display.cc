// Cytosim was created by Francois Nedelec. Copyright Cambridge University 2021

#include "smath.h"
#include "display.h"
#include "organizer.h"
#include "property_list.h"
#include "hand_prop.h"
#include "sphere_prop.h"
#include "fiber_prop.h"
#include "random_pcg.h"

extern uint64_t pcg32_state;

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
#include "gym_flat.h"
#include "fg_font.h"

#include "point_disp.h"
#include "fiber_disp.h"
#include "line_disp.h"
#include "display_color.cc"


Display::Display(DisplayProp const* dp)
: pixelSize(1), unitValue(1), sizeScale(1), depthAxis(0,0,1), prop(dp)
{
    assert_true(dp);
    age_start = 0;
    age_scale = 1.0;
    allLineDisp = nullptr;
    numLineDisp = 0;
}

void Display::setParameters(float ps, float uv, Vector3 const& ax)
{
    pixelSize = ps;
    unitValue = uv;
    /*
     the 0.5 below comes from the fact that glPointSize uses diameter
     while most gle::primitives have a radius of 1
     */
    sizeScale = 0.5f * uv * ps;
    //printf(" pixelSize %6.4f unitValue %6.4f : %6.4f\n", ps, uv, sizeScale);
    depthAxis = ax;
}

Display::~Display()
{
    delete[] allLineDisp;
}


static void strokeText(Vector const& vec, const char str[], float scale)
{
#if ( DIM == 3 )
    gym::translate_ref(vec.XX, vec.YY, vec.ZZ);
#elif ( DIM == 2 )
    gym::translate_ref(vec.XX, vec.YY, 0.f);
#else
    gym::translate_ref(vec.XX, 0.f, 0.f);
#endif
    fgBitmapToken(0, 0, scale, 2, str);
    //fgStrokeString(0, 0, scale, 0, str, 1);
    CHECK_GL_ERROR("strokeText");
}


//------------------------------------------------------------------------------
#pragma mark - drawObject

void Display::drawObject(Vector const& pos, float rad, void(*obj)())
{
    gym::transScale(pos, rad);
    obj();
}


void Display::drawObject(Vector const& pos, Vector const& dir, float rad, void(*obj)())
{
    gym::transAlignZ(pos, rad, dir);
    obj();
}


void drawBallT(Vector const& pos, real rad, ObjectMark mark)
{
    gym::transScale(pos, rad);
    if ( mark & 7 )
        gle::football2(gym_color(0,0,0));
    else
        gle::dualPassSphere2();
}

// using sphere4() for presumably smaller objects
void drawBeadS(Vector const& pos, real rad, ObjectMark mark)
{
    gym::transScale(pos, rad);
    if ( mark & 7 )
        gle::football4(gym_color(0,0,0));
    else
        gle::dualPassSphere4();
}


void drawDiscT(Vector const& pos, real rad)
{
    gym::transScale(pos, rad);
    gym::disableLighting();
    gle::disc();
}


/// used for drawFilament
inline void drawMonomer(Vector3 const& pos, float rad)
{
    gym::transScale(pos, rad);
    gle::sphere2();
}

static void drawFootball(Solid const& obj, size_t inx, gym_color col, gym_color bak, bool flip)
{
    assert_true(inx+DIM < obj.nbPoints());
    Vector X = obj.posP(inx);
#if ( DIM >= 3 )
    Vector A = obj.posP(inx+1) - X;
    Vector B = obj.posP(inx+2) - X;
    Vector C = obj.posP(inx+3) - X;
    gym::transRotate(X, A, B, C);
    //bool flip = ( dot(cross(A,B), C) < 0 );
#else
    gym::transScale(X, obj.radius(inx));
#endif
    if ( flip ) glFrontFace(GL_CW);
    gym::color_front(col, 1.0);
    gle::football1(bak);
    //gle::baseball();
    if ( flip ) glFrontFace(GL_CCW);
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
    drawFiberTexts(sim.fibers);

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
    
    for ( int d = 0; d < DIM; ++d )
    {
        if ( modulo->isPeriodic(d) )
        {
            l[d] = -((tile>>d)&1);
            u[d] = +1;
        }
    }
    
    const Vector3 px = modulo->period(0);
    const Vector3 py = modulo->period(1);
    const Vector3 pz = modulo->period(2);
    
    Vector4 pos[32];
    int cnt = 0;
    
    for ( int dx = l[0]; dx <= u[0]; ++dx )
        for ( int dy = l[1]; dy <= u[1]; ++dy )
            for ( int dz = l[2]; dz <= u[2]; ++dz )
            {
                Vector3 P = dx * px + dy * py + dz * pz;
                pos[cnt] = Vector4(P);
                pos[cnt++].TT = dot(depthAxis, Vector3(P));
            }
    
    // depth-sort positions:
    qsort(pos, cnt, sizeof(Vector4), &compareVector4);
    
    float ref[16];
    gym::get_view(ref);
    for ( int i = 0; i < cnt; ++i )
    {
        gym::set_view(ref, pos[i].XX, pos[i].YY, pos[i].ZZ);
        drawSimul(sim);
    }
    gym::set_view(ref);
}

//------------------------------------------------------------------------------
#pragma mark -


/**
 Create a FiberDisp for this Property if necessary
 */
void Display::initFiberDisp(FiberProp* fp, PropertyList& depot, gym_color col)
{
    FiberDisp *& dis = fp->disp;
    
    // recover existing property:
    if ( !dis )
        dis = static_cast<FiberDisp*>(depot.find("fiber:display", fp->name()));

    // create new property with default values:
    if ( !dis )
    {
        dis = new FiberDisp(fp->name());
        depot.push_back(dis);
        // set default:
        dis->color      = col;
        dis->back_color = col.darken(0.5);
        dis->point_size = prop->point_size;
        dis->line_width = prop->line_width;
    }
    
    // parse user-provided values:
    if ( fp->display_fresh )
    {
        dis->read_string(fp->display, fp->name()+":display");
        if ( dis->mask && !dis->mask_bitfield )
            dis->mask_bitfield = distribute_bits(dis->mask, pcg32_state);
        fp->display_fresh = false;
    }
    
    if ( dis->coloring == FiberDisp::COLORING_CLUSTER )
        prep_flag |= 1;
    
    if ( dis->line_style == 2 || dis->line_style == 3 )
        prep_flag |= 2;

    if ( dis->coloring == FiberDisp::COLORING_AGE )
        prep_flag |= 4;
}


/**
 set LineDisp for given Fiber
 */
void Display::initLineDisp(const Fiber * fib, FiberDisp const* dis, LineDisp * self)
{
    bool hide = false;
    assert_true(fib->prop);
    gym_color col = dis->color;
    
    // change body color depending on coloring mode:
    switch ( dis->coloring )
    {
        default:
        case FiberDisp::COLORING_OFF:
            col = dis->color;
            break;
        case FiberDisp::COLORING_RANDOM:
            col = gym::bright_color(fib->signature()).match_a(dis->color);
            break;
        case FiberDisp::COLORING_DIRECTION:
            col = radial_color(fib->direction(), dis->color.alpha());
            break;
        case FiberDisp::COLORING_MARK:
            col = gym::get_color(fib->mark());
            break;
        case FiberDisp::COLORING_FLAG:
            col = gym::std_color(fib->flag());
            break;
#if FIBER_HAS_FAMILY
        case FiberDisp::COLORING_FAMILY:
            if ( fib->family_ )
                col = gym::get_color(fib->family_->signature());
            else
                col = dis->color;
            break;
#endif
        case FiberDisp::COLORING_CLUSTER:
            col = gym::std_color(fib->flag());
            break;
        case FiberDisp::COLORING_AGE:
            col = gym_color::jet_color((fib->age()-age_start)*age_scale, 1.0);
            break;
        case FiberDisp::COLORING_PSTATE:
            if ( fib->endStateP() > 0 )
                col = dis->end_colors[std::min(fib->endStateP(),5U)];
            break;
    }
    
    self->color = col;
    self->color_scale = color_scale(fib, dis->line_style);
    //std::cerr << fib->reference() << ":color_scale " << self->color_scale << "\n";
    
#if ( 1 )
    // colors of ends set to match body color:
    self->end_color[0] = col;
    self->end_color[1] = col;
#else
    // use colors of non-dynamic ends:
    self->end_color[0] = dis->end_colors[0];
    self->end_color[1] = dis->end_colors[0];
#endif
    
#if ( 1 )
    // For dynamic Fibers, change colors of tips according to state:
    if ( fib->endStateP() > 0 )
        self->end_color[0] = dis->end_colors[std::min(fib->endStateP(),5U)];
    if ( fib->endStateM() > 0 )
        self->end_color[1] = dis->end_colors[std::min(fib->endStateM(),5U)];
#else
    // For dynamic Fibers, change colors of tips according to state:
    if ( fib->freshAssemblyP() > 0 )
        self->end_color[0] = dis->end_colors[std::min(fib->endStateP(),5U)];
    if ( fib->freshAssemblyM() > 0 )
        self->end_color[1] = dis->end_colors[std::min(fib->endStateM(),5U)];
#endif
    
    // hide right or left-pointing fibers:
    if (( dis->hide & 1 )  &&  dot(fib->diffPoints(0), Vector(dis->hide_axis)) < 0 )
        hide = true;
    if (( dis->hide & 2 )  &&  dot(fib->diffPoints(0), Vector(dis->hide_axis)) > 0 )
        hide = true;
    
#if ( DIM == 2 )
    // hide clockwise or counter-clockwise orientated fibers:
    if (( dis->hide & 4 )  &&  cross(fib->posP(0), fib->diffPoints(0)) < 0 )
        hide = true;
    if (( dis->hide & 8 )  &&  cross(fib->posP(0), fib->diffPoints(0)) > 0 )
        hide = true;
#elif ( DIM >= 3 )
    // hide clockwise or counter-clockwise orientated fibers in the XY plane
    if (( dis->hide & 4 )  &&  cross(fib->posP(0), fib->diffPoints(0)).ZZ < 0 )
        hide = true;
    if (( dis->hide & 8 )  &&  cross(fib->posP(0), fib->diffPoints(0)).ZZ > 0 )
        hide = true;
#endif
    
#if ( 1 )
    // hide fibers depending on mask
    if ( fib->signature() & dis->mask_bitfield )
        hide = true;
#else
    if ( fib->mark() & dis->mask_bitfield )
        hide = true;
#endif

    // hide fibers in a specified state
    if ( fib->endStateP() == dis->hide_state )
        hide = true;
    
    // hide fibers except those with a certain mark
    if ( dis->show_marked != ~0U && fib->mark() != dis->show_marked )
        hide = true;

    // change color of 'hidden' filament:
    if ( hide )
        self->color = dis->hide_color;
    
    // default visibility set from class:
    if ( dis->visible )
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
    if ( dis->explode_style )
        self->explode_shift = ( fib->signature() * 0x1p-32 - 0.5 ) * dis->explode_range;
    else
        self->explode_shift = 0;
}


/**
 Create a PointDisp for this Property if necessary
 */
template < typename T >
void Display::initPointDisp(T * p, PropertyList& depot, gym_color col)
{
    PointDisp *& dis = p->disp;
    
    // create new property:
    if ( !dis )
    {
        // search for matching property:
        dis = static_cast<PointDisp*>(depot.find(p->category()+":display", p->name()));
        if ( !dis )
        {
            //std::clog <<" new " << p->category() << ":display " << p->name() << "\n";
            dis = new PointDisp(p->category()+":display", p->name());
            depot.push_back(dis);
            // set default:
            dis->clear();
            dis->color  = col;
            dis->color2 = col.alpha_scaled(0.1875f);
            dis->size   = prop->point_size;
            if ( p->category() == "hand" )
                dis->width = prop->link_width;
            else
                dis->width = prop->line_width;
        }
    }
    
    // parse display string once:
    if ( p->display_fresh )
    {
        dis->read_string(p->display, p->name()+":display");
        //std::clog << dis->color << "  " << p->name() << ":display (" << p->category() << ")\n";
        p->display_fresh = false;
    }
    
    dis->setPixels(pixelSize, unitValue, prop->style==2);
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
        initLineDisp(fib, fib->prop->disp, fib->disp);
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
void Display::prepareDrawing(Simul const& sim, PropertyList& fiberDisp, PropertyList& alldisp)
{
    // counter to give different colors to the objects
    size_t idx = 0;

    PropertyList plist = sim.properties.find_all("fiber");
    
    prep_flag = 0;
    // create a FiberDisp for each FiberProp:
    for ( Property* p : plist )
        initFiberDisp(static_cast<FiberProp*>(p), fiberDisp, gym::get_color(idx++));

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
        initPointDisp(static_cast<HandProp*>(i), alldisp, gym::get_color(idx++));
    
    //create a PointDisp for each SphereProp:
    for ( Property * i : sim.properties.find_all("sphere") )
        initPointDisp(static_cast<SphereProp*>(i), alldisp, gym::bright_color(idx++));
    
    //create a PointDisp for each SolidProp:
    for ( Property * i : sim.properties.find_all("solid", "bead") )
        initPointDisp(static_cast<SolidProp*>(i), alldisp, gym::bright_color(idx++));
    
    //create a PointDisp for each SpaceProp:
    gym_color col(DIM==3?0x00000044:0xAAAAAAFF);
    for ( Property * i : sim.properties.find_all("space") )
        initPointDisp(static_cast<SpaceProp*>(i), alldisp, col);
}


//------------------------------------------------------------------------------
#pragma mark -


/**
 Draw the back and front sides of Spaces in 3D
 This function is called twice: at the start and at the end of drawSimul()
 */
void Display::drawSpace3D(Space const* obj, bool back)
{
    PointDisp const* dis = obj->prop->disp;
    bool front = back ^ ( dis->color.transparent() );
    
    back = ( dis->visible & 2 ) && back;
    front = ( dis->visible & 1 ) && front;
    
    if ( back | front )
    {
        gym::ref_view();
        gym::enableLighting();
        gym::enableCullFace(GL_FRONT);
        if ( back )
        {
            gym::color_back(dis->color2);
            obj->draw3D();
        }
        if ( front )
        {
            gym::switchCullFace(GL_BACK);
            gym::color_load(dis->color);
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
        PointDisp const* dis = obj->prop->disp;
        if ( dis->visible )
        {
            gym::disableLighting();
            gym::color(dis->color);
            gym::ref_view();
            obj->draw2D(dis->widthX);
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
#if ( 1 )
    Field * obj = set.first();
#else
    for ( Field * obj = set.first(); obj; obj=obj->next() )
#endif
    
    if ( obj && obj->hasField() )
    {
        if ( obj->prop->visible == 1 )
            obj->draw(obj->prop->field_space_ptr, depthAxis, 0);
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

void Display::drawAverageFiber1(FiberSet const& fibers, Property const* arg) const
{
    ObjectList objs = fibers.collect(arg);
    drawAverageFiber(objs, gym_color(1,1,1));
}


void Display::drawAverageFiber2(FiberSet const& fibers, Property const* arg) const
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
        FiberProp const* P = static_cast<FiberProp const*>(i);
        if ( P->disp->draw_average == 1 )
            drawAverageFiber1(sim.fibers, P);
        else if ( P->disp->draw_average == 2 )
            drawAverageFiber2(sim.fibers, P);
    }
#endif
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
 .
 */
void Display::drawFiberEndMinus(Fiber const& fib, int style, float size) const
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
 .
 */
void Display::drawFiberEndPlus(Fiber const& fib, int style, float size) const
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
        case 8: drawObject(fib.posEndP(),-fib.dirEndP(), rad, gle::hemisphere2); break;
    }
}


/// draw fresh assembly near the plus ends, using white stripes
void Display::drawFiberGrowth(Fiber const& fib, float size) const
{
    if ( fib.freshAssemblyM() > 0 )
    {
        gym::ref_view();
        const float rad = 1.03125 * pixscale(size);
        gym::color_both(1, 1, 1, 1);
        gym::stretchAlignZ1(fib.posEndM(), rad, fib.dirEndM(), fib.freshAssemblyM());
        gle::stripes(1.f);
    }
    if ( fib.freshAssemblyP() > 0 )
    {
        gym::ref_view();
        const float rad = 1.03125 * pixscale(size);
        gym::color_both(1, 1, 1, 1);
        gym::stretchAlignZ1(fib.posEndP(), rad, -fib.dirEndP(), fib.freshAssemblyP());
        gle::stripes(1.f);
    }
}


void Display::drawFiberBackbone(Fiber const& fib, gym_color col, float width) const
{
    gym::ref_view();
    gym::color(col);
    gym::disableLighting();
    gym::loadPoints(fib.nbPoints(), fib.addrPoints());
    gym::drawLineStrip(width, 0, fib.nbPoints());
    gym::cleanupV();
}


void Display::drawFiberLines(Fiber const& fib, const int style, float width) const
{
    if ( style == 8 && fib.endStateP() != STATE_GREEN )
        return;
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
        case 4: // color according to the angle with respect to the XY-plane:
            strip = 0;
            for ( size_t n = 0; n < fib.nbSegments(); ++n )
            {
                gym_color c = color_by_direction(fib, n);
                ptr[0] = {c, fib.posP(n)};
                ptr[1] = {c, fib.posP(n+1)};
                ptr += 2;
            }
            break;
        case 5: // display segments with color indicating the curvature
            for ( size_t i = 0; i < fib.nbPoints(); ++i )
                flu[i] = {color_by_curvature(fib, i), fib.posP(i)};
            ptr = flu + fib.nbPoints();
            break;
        case 6: // color according to the distance from the minus end
            *ptr++ = {color_by_distanceM(fib, 0), fib.posP(0)};
            for ( real a = 0.0625; a < 0.6; a *= 2 )
                *ptr++ = {color_by_distanceM(fib, a), fib.midPoint(0, a)};
            for ( size_t n = 1; n < fib.nbPoints(); ++n )
                *ptr++ = {color_by_distanceM(fib, n), fib.posP(n)};
            break;
        case 7: case 8: { // color according to the distance from the plus end
            const size_t last = fib.lastPoint();
            for ( size_t n = 0; n < last; ++n )
                *ptr++ = {color_by_distanceP(fib, n), fib.posP(n)};
            for ( real a = 0.5; a > 0.06; a /= 2 )
                *ptr++ = {color_by_distanceP(fib, last-a), fib.midPoint(last-1, 1-a)};
            *ptr++ = {color_by_distanceP(fib, last), fib.posP(last)};
        } break;
        case 9: // color according to distance to the confining Space
            for ( size_t i = 0; i < fib.nbPoints(); ++i )
                flu[i] = {color_by_height(fib, i), fib.posP(i)};
            ptr = flu + fib.nbPoints();
            break;
    }
    gym::ref_view();
    gym::disableLighting();
    gym::unmapBufferC4VD();
    float w = pixwidth(width);
    if ( strip )
        gym::drawLineStrip(w, 0, ptr-flu);
    else
        gym::drawLines(w, 0, ptr-flu);
    gym::cleanupCV();
}


void Display::drawFiberSegmentT(Fiber const& fib, unsigned inx) const
{
    FiberDisp const*const dis = fib.prop->disp;
    const int style = dis->line_style;
    if ( style == 8 && fib.endStateP() != STATE_GREEN )
        return;
    size_t cnt = 8;
    flute4D* flu = gym::mapBufferC4VD(cnt);
    flute4D* ptr = flu;
    
    if ( style == 6 )
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
    else if ( style == 7 || style == 8 )
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
        if ( style == 2 )
            c = color_by_tension(fib, inx);
        else if ( style == 3 )
            c = color_by_tension_jet(fib, inx);
        else
            c = fib.disp->color;
        // the whole segment is painted with the same color:
        ptr[0] = {c, fib.posP(inx)};
        ptr[1] = {c, fib.posP(inx+1)};
        ptr += 2;
    }
    gym::ref_view();
    gym::disableLighting();
    gym::unmapBufferC4VD();
    gym::drawLineStrip(pixwidth(dis->line_width), 0, ptr-flu);
    gym::restoreLighting();
    gym::cleanupCV();
}


void Display::drawFiberSpeckles(Fiber const& fib) const
{
    FiberDisp const*const dis = fib.prop->disp;
    const real gap = dis->speckle_gap;

    size_t i = 0, cnt = 8 + 4 * std::ceil(fib.length()/gap);
    fluteD* pts = gym::mapBufferVD(cnt);

    // display random speckles:
    if ( dis->speckle_style == 1 )
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
    else if ( dis->speckle_style == 2 )
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
    gym::color(fib.disp->color);
    gym::drawPoints(dis->speckle_size, 0, i);
}


void Display::drawFiberPoints(Fiber const& fib) const
{
    FiberDisp const*const dis = fib.prop->disp;
    int style = dis->point_style;

    if ( style == 1 )
    {
        // display vertices:
        gym::ref_view();
        gym::disableLighting();
        gym::color(fib.disp->color);
        gym::loadPoints(fib.nbPoints(), fib.addrPoints());
        gym::drawSquarePoints(pixwidth(dis->point_size), 0, fib.nbPoints());
    }
    else if ( style == 2 )
    {
        gym::color_load(fib.disp->color);
        gym::color_back(dis->back_color);
        // display arrowheads along the fiber:
        const float rad = pixscale(dis->point_size);
        const real gap = dis->point_gap;
        real ab = std::ceil(fib.abscissaM()/gap) * gap;
        for ( ; ab <= fib.abscissaP(); ab += gap )
            gle::drawCone(fib.pos(ab), fib.dir(ab), rad);
    }
    else if ( style == 3 )
    {
        // display chevrons along the fiber:
        const real gap = dis->point_gap;
        const real rad = pixscale(dis->point_size);
        real beta = fib.segmentationInv() * rad;
        size_t cnt = 4 * fib.length() / gap + 8;
        fluteD* flu = gym::mapBufferVD(cnt);
        fluteD* ptr = flu;
        real a = std::ceil(fib.abscissaM()/gap) * gap - fib.abscissaM();
        for ( ; a <= fib.length(); a += gap )
        {
            Interpolation i = fib.interpolateM(a);
            Vector pos = i.pos();
            Vector dir = i.diff() * beta;
            Vector off = inViewPlane(dir, rad);
            ptr[0] = pos + dir - off;
            ptr[1] = pos - dir;
            ptr[2] = pos - dir;
            ptr[3] = pos + dir + off;
            ptr += 4;
        }
        assert_true( ptr <= flu+cnt );
        gym::unmapBufferVD();
        gym::ref_view();
        gym::disableLighting();
        gym::color(fib.disp->color);
        gym::drawLines(pixwidth(1), 0, ptr-flu);
    }
    else if ( style == 4 )
    {
        gym::color_load(fib.disp->color);
        gym::color_back(dis->back_color);
        // display only middle of fiber:
        drawObject(fib.posMiddle(), pixscale(2*dis->point_size), gle::sphere2);
    }
}

//------------------------------------------------------------------------------
#pragma mark - Lattice

/**
 Draw lattice using one vertex for each site, positionned at the center of the range
 OpenGL will interpolate the colors, and each site will be covered by a gradient.
 */
void Display::drawFiberLattice1(Fiber const& fib, VisibleLattice const& lat, float rad) const
{
    rad = pixwidth(rad);
    const real uni = lat.unit();
    const auto inf = lat.indexM();
    const auto sup = lat.indexP();
    assert_true( inf <= sup );
    
    FiberDisp const*const dis = fib.prop->disp;
    gym_color c, col = dis->color;
    const real fac = 1 / dis->lattice_scale;
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
        if ( dis->lattice_rescale )
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
    gym::ref_view();
    gym::unmapBufferC4VD();
    gym::disableLighting();
    gym::drawLineStrip(rad, 0, ptr-flu);
    gym::cleanupCV();
}


/**
 Draw lattice using two vertices for each site, at both extremities of the range,
 and each site is entirely covered by the color corresponding to the value.
 */
void Display::drawFiberLattice2(Fiber const& fib, VisibleLattice const& lat, float rad) const
{
    rad = pixwidth(rad);
    const real uni = lat.unit();
    const auto inf = lat.indexM();
    const auto sup = lat.indexP();
    assert_true( inf <= sup );
    
    FiberDisp const*const dis = fib.prop->disp;
    gym_color c, col = dis->color;
    const real fac = 1 / dis->lattice_scale;
    size_t cnt = 2 + 2 * ( sup - inf );
    flute4D* flu = gym::mapBufferC4VD(cnt);
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
        if ( dis->lattice_rescale )
        {
            facM = ( lenM > 0.01*uni ? fac*uni/lenM : fac );
            facP = ( lenP > 0.01*uni ? fac*uni/lenP : fac );
        }

        // the terminal site may be truncated
        c = lattice_color(col, facM*lat.data(inf));
        *ptr++ = {c, fib.posEndM()};

        Vector P;
        for ( auto h = inf+1; h <= sup; ++h )
        {
            P = fib.pos(uni*h);
            ptr[0] = {c, P};
            c = lattice_color(col, fac*lat.data(h));
            ptr[1] = {c, P};
            ptr += 2;
        }
        
        // the terminal site may be truncated, so we change its color:
        c = lattice_color(col, facP*lat.data(sup));
        ptr[-1] = {c, P};
        *ptr++ = {c, fib.posEndP()};
    }
    assert_true( ptr <= flu+cnt );
    gym::ref_view();
    gym::unmapBufferC4VD();
    gym::disableLighting();
    gym::drawLines(rad, 0, ptr-flu);
    gym::cleanupCV();
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
    assert_true( ptr <= flu+cnt );
    gym::unmapBufferVD();
    gym::ref_view();
    gym::disableLighting();
    gym::color(fib.disp->color);
    gym::drawPoints(rad, 0, ptr-flu);
    gym::cleanupV();
}


/**
 Display the value of each site in the lattice, in binary form if integers are used
 */
void Display::drawFiberLatticeValues(Fiber const& fib, VisibleLattice const& lat) const
{
    char str[66] = { 0 };
    const real uni = lat.unit();
    const auto inf = lat.indexM();
    const auto sup = lat.indexP();
    
    real abs = (inf+0.25) * uni - fib.abscissaM();
    for ( auto h = inf; h <= sup; ++h, abs += uni )
    {
        snprintf(str, sizeof(str), "%.3f", real(lat.data(h)));
        strokeText(fib.posM(abs), str, pixelSize);
    }
    gym::restoreAlphaTest();
}


/**
 Display the value of each site in the lattice, in binary form if integers are used
 */
void Display::drawFiberLatticeBits(Fiber const& fib, FiberLattice const& lat) const
{
    char str[66] = { 0 };
    const real uni = lat.unit();
    const auto inf = lat.indexM();
    const auto sup = lat.indexP();
    
    real abs = (inf+0.25) * uni - fib.abscissaM();
    for ( auto h = inf; h <= sup; ++h, abs += uni )
    {
#if FIBER_HAS_LATTICE > 0
        //Cymath::binary_representation(str, sizeof(str), 4, lat.data(h));
        if ( lat.data(h) )
        {
            Vector vec = fib.posM(abs);
            gym::translate_ref(vec.XX, vec.y(), vec.z());
            gym::paintSequence(0, 0, 2*pixelSize, 12*pixelSize, lat.data(h), 8*sizeof(FiberLattice::cell_t));
        }
#else
        snprintf(str, sizeof(str), "%.3f", lat.data(h));
        strokeText(fib.posM(abs), str, pixelSize);
#endif
    }
}


void Display::drawFiberLabels(Fiber const& fib, int style) const
{
    char str[32] = { 0 };
    if ( style & 1 )
    {
        // draw fiber identity and vertex indices
        int C = snprintf(str, sizeof(str), " %u:", fib.identity());
        for ( size_t i = 0; i < fib.nbPoints(); ++i )
        {
            snprintf(str+C, sizeof(str)-C, "%lu", i);
            strokeText(fib.posP(i), str, pixelSize);
        }
    } 
    else if ( style & 2 )
    {
        // draw fiber identity and abscissa value at vertices
        int C = snprintf(str, sizeof(str), " %u:", fib.identity());
        for ( size_t i = 0; i < fib.nbPoints(); ++i )
        {
            snprintf(str+C, sizeof(str)-C, "%.3f", fib.abscissaPoint(i));
            strokeText(fib.posP(i), str, pixelSize);
        }
    }
    if ( style & 4 )
    {
        // display integral abscissa along the fiber
        snprintf(str, sizeof(str), "%.3f", fib.abscissaM());
        strokeText(fib.posEndM(), str, pixelSize);
        
        int s = (int)std::ceil(fib.abscissaM());
        int e = (int)std::floor(fib.abscissaP());
        for ( int a = s; a <= e; ++a )
        {
            snprintf(str, sizeof(str), "%i", a);
            strokeText(fib.pos(a), str, pixelSize);
        }
        
        snprintf(str, sizeof(str), "%.3f", fib.abscissaP());
        strokeText(fib.posEndP(), str, pixelSize);
    }
    if ( style & 8 )
    {
        // indicate tension values in the segments
        Vector a = fib.posEndM();
        for ( size_t i = 1; i < fib.nbPoints(); ++i )
        {
            Vector b = fib.posP(i);
            snprintf(str, sizeof(str), "%+4.1f", fib.tension(i-1));
            strokeText(0.5*(a+b), str, pixelSize);
            a = b;
        }
    }
}


/// display forces acting on the fiber's vertices, using lines scaled by 'mag'
void Display::drawFiberForces(Fiber const& fib, real mag, float size) const
{
    gym_color col = fib.prop->disp->force_color;
    gym_color lor = col.alpha_scaled(0.5f);
    unsigned cnt = 2 * fib.nbPoints();
    flute4D* flu = gym::mapBufferC4VD(cnt);
    for ( size_t i = 0; i < fib.nbPoints(); ++i )
    {
        Vector P = fib.posP(i);
        Vector F = mag * fib.netForce(i);
        flu[  2*i] = { col, P };
        flu[1+2*i] = { lor, P+F };
    }
    gym::ref_view();
    gym::unmapBufferC4VD();
    gym::disableLighting();
    gym::drawLines(size, 0, cnt);
    gym::cleanupCV();
}

//------------------------------------------------------------------------------
#pragma mark - Wide and Striped Fiber styles


/** Using triangles to draw a broken line of thickness `2*rad` */
void Display::drawFiberWide(Fiber const& fib, float rad) const
{
    fluteD* flu = gym::mapBufferVD(2*fib.nbPoints()+4);
    fluteD* ptr = flu;
    
    Vector pos = fib.posPoint(0);
    Vector nxt = fib.posPoint(1);
    Vector off = inViewPlane(nxt-pos, rad);
    ptr[0] = pos + off;
    ptr[1] = pos - off;
    ptr += 2;
    
    Vector old;
    for ( size_t i = 1; i < fib.nbSegments(); ++i )
    {
        old = pos;
        pos = nxt;
        nxt = fib.posPoint(i+1);
        off = inViewPlane(nxt-old, rad);
        ptr[0] = pos + off;
        ptr[1] = pos - off;
        ptr += 2;
    }
    off = inViewPlane(nxt-pos, rad);
    ptr[0] = nxt + off;
    ptr[1] = nxt - off;
    ptr += 2;
    gym::ref_view();
    gym::disableLighting();
    gym::unmapBufferVD();
    gym::drawTriangleStrip(0, ptr-flu);
}


/**
Draw segments of length 'inc' with a triangle of length 'off'
*/
void Display::drawFiberArrowed2D(Fiber const& fib, float rad, real inc,
                                 gym_color col, real len, gym_color lor) const
{
    const real sup = fib.length() + len;
    int cnt = 1 + (int)std::floor(fib.abscissaM()/inc);
    // abs in [0, uni] is now relative to minus end
    real abs = inc * cnt - fib.abscissaM();
    // draw segments
    size_t top = 10 * std::floor(sup/inc) + 14;
    flute4D* flu = gym::mapBufferC4VD(top);
    flute4D* ptr = flu;
    Vector pos = fib.posEndM();
    Vector dir = fib.dirEndM();
    Vector off = inViewPlane(dir, rad);
    Vector nxt = pos;
    while ( abs < sup )
    {
        pos = nxt;
        nxt = fib.displayPosM(abs);
        abs += inc;
        flute3 A(pos + off), B(pos - off);
        flute3 M(pos + dir * len);
        //
        dir = fib.dir(abs);
        off = inViewPlane(dir, rad);
        // paint triangle complement:
        ptr[0] = { lor, A };
        ptr[1] = { lor, nxt + off };
        ptr[2] = { lor, M };
        ptr[3] = { lor, nxt - off };
        ptr[4] = { lor, B };
        // paint plus-end facing triangle:
        ptr[5] = { col, B };
        ptr[6] = { col, M };
        ptr[7] = { col, A };
        ptr[8] = { col, A };
        ptr[9] = { col, nxt + off };
        ptr += 10;
    }
    pos = fib.posEndP();
    ptr[0] = { col, nxt + off };
    ptr[1] = { col, pos + off };
    ptr[2] = { col, nxt - off };
    ptr[3] = { col, pos - off };
    ptr += 4;
    assert_true( ptr <= flu+top );

    gym::ref_view();
    gym::unmapBufferC4VD();
    gym::disableLighting();
    gym::drawTriangleStrip(0, ptr-flu);
    gym::cleanupCV();
}


/**
Draw segments of length 'inc' and 'onk' of alternating colors, in register with Fiber's abscissa
*/
void Display::drawFiberStriped2D(Fiber const& fib, float rad, real inc,
                                 gym_color col, real onk, gym_color lor) const
{
    gym_color black(0,0,0,1);
    if ( fib.mark() ) black.set(0.75,0.75,0.75,1);
    const real uni = inc + onk;
    const real len = fib.length();
    int cnt = 1 + (int)std::floor(fib.abscissaM()/uni);
    real abs = uni * cnt - fib.abscissaM();
    // abs in [0, uni] is now relative to minus end
    if ( abs > onk )
    {
        cnt = 2*cnt - 1; // drawing first slice of size `inc`
        abs -= onk;
        std::swap(inc, onk);
    }
    else
    {
        cnt = 2*cnt; // drawing second slice of size `onk`
        std::swap(col, lor);
    }
    gym_color clr = col;
    // draw segments
    size_t top = 8 * fib.length() / uni + 8;
    flute4D* flu = gym::mapBufferC4VD(top);
    flute4D* ptr = flu;
    Vector pos = fib.posEndM();
    Vector dir = fib.dirEndM();
    Vector off = inViewPlane(dir, rad);
    ptr[0] = { clr, pos - off };
    ptr[1] = { clr, pos + off };
    ptr += 2;

    while ( abs < len )
    {
        dir = fib.dir(abs);
        pos = fib.displayPosM(abs);
        abs += inc;
        off = inViewPlane(dir, rad);
        // alternate different tones:
        ptr[0] = { clr, pos - off };
        ptr[1] = { clr, pos + off };
        std::swap(col, lor);
        std::swap(inc, onk);
        if ( (++cnt & 15) == 1 ) clr = black; else clr = col;
        ptr[2] = { clr, pos - off };
        ptr[3] = { clr, pos + off };
        ptr += 4;
    }
    pos = fib.posEndP();
    ptr[0] = { clr, pos - off };
    ptr[1] = { clr, pos + off };
    ptr += 2;
    assert_true( ptr <= flu+top );

    gym::ref_view();
    gym::unmapBufferC4VD();
    gym::disableLighting();
    gym::drawTriangleStrip(0, ptr-flu);
    gym::cleanupCV();
}


/**
 Draw segments of length 'inc, onk' of alternating colors, in register with Fiber's abscissa
 Every 16th stripe is drawn in black, hence with separation of 256 nm for microtubules
*/
void Display::drawFiberStriped(Fiber const& fib, float rad, real inc,
                               gym_color col, real onk, gym_color lor) const
{
    gym_color black(0,0,0,1);
    if ( fib.mark() ) black.set(0.75,0.75,0.75,1);
    const real uni = inc + onk;
    Vector pos, nxt, old = fib.displayPosM(0);
    const real len = fib.length();
    int cnt = 1 + (int)std::floor(fib.abscissaM()/uni);
    real abs = uni * cnt - fib.abscissaM();
    // abs in [0, uni] is now relative to minus end
    if ( abs > onk )
    {
        cnt = 2*cnt - 1; // drawing first slice of size `inc`
        pos = fib.displayPosM(abs-onk);
        nxt = fib.displayPosM(abs);
        gym::color_load(col);
    }
    else
    {
        cnt = 2*cnt; // drawing second slice of size `onk`
        pos = fib.displayPosM(abs);
        nxt = fib.displayPosM(abs+inc);
        gym::color_load(lor);
        abs += inc;
    }
    // draw first segment near the minus end:
    gym::stretchAlignZ(old, pos, rad);
    gle::tube2(); gle::dome();

    // draw middle segments:
    while ( abs < len )
    {
        abs += (cnt&1)?inc:onk;
        old = pos;
        pos = nxt;
        nxt = fib.displayPosM(abs);
        // alternate different tones:
        gym::color_load((++cnt&1)?col:lor);
        if ( (cnt & 15) == 1 ) gym::color_load(black);
        gym::stretchAlignZ(old, pos, rad);
        gle::tube2();
    }

    // draw last segment, which may be truncated:
    gym::color_load((++cnt&1)?col:lor);
    gym::stretchAlignZ(fib.posEndP(), pos, rad);
    gle::shutTube2();
}


/**
 Draw segments of length 'inc, onk' of alternating colors, in register with Fiber's abscissa
 Every 16th stripe is drawn in black, hence with separation of 256 nm for microtubules
*/
void Display::drawFiberStripedClip(Fiber const& fib, float rad, real inc,
                                   gym_color col, real onk, gym_color lor) const
{
    gym_color black(0,0,0,1);
    if ( fib.mark() ) black.set(0.75,0.75,0.75,1);
    const real uni = inc + onk;
    Vector pos, nxt, old = fib.displayPosM(0);
    const real len = fib.length();
    int cnt = 1 + (int)std::floor(fib.abscissaM()/uni);
    real abs = uni * cnt - fib.abscissaM();
    // abs in [0, uni] is now relative to minus end
    if ( abs > onk )
    {
        cnt = 2*cnt - 1; // drawing first slice of size `inc`
        pos = fib.displayPosM(abs-onk);
        nxt = fib.displayPosM(abs);
        gym::color_load(col);
    }
    else
    {
        cnt = 2*cnt; // drawing second slice of size `onk`
        pos = fib.displayPosM(abs);
        nxt = fib.displayPosM(abs+inc);
        gym::color_load(lor);
        abs += inc;
    }
    Vector dir = normalize(nxt-old);
    
    gym::enableClipPlane(4);
    gym::setClipPlane(4, -dir, pos);
    gym::stretchAlignZ1(old, rad, dir, rad);
    gle::capedTube();
    gym::setClipPlane(4, dir, pos);
#if 0
    gym::disableClipPlane(4);
    gym::color_load(gym_color(1,1,1)); drawObject(old, 0.002, gle::octahedron);
    gym::color_load(gym_color(0,0,1)); drawObject(pos, 0.002, gle::octahedron);
    gym::color_load(gym_color(1,1,1)); drawObject(nxt, 0.002, gle::octahedron);
    gym::enableClipPlane(4);
#endif
    gym::enableClipPlane(5);
    // draw segments
    while ( abs < len )
    {
        abs += (cnt&1)?inc:onk;
        old = pos;
        pos = nxt;
        nxt = fib.displayPosM(abs);
        dir = normalize(nxt-old);
        // alternate different tones:
        gym::color_load((++cnt&1)?col:lor);
        if ( (cnt & 15) == 1 ) gym::color_load(black);
        gym::setClipPlane(5, -dir, pos);
        gym::transAlignZ(old, rad, pos-old);
        gle::innerTube();
        gym::setClipPlane(4, dir, pos);
    }
    gym::disableClipPlane(5);

    // draw last segment, which may be truncated:
    gym::color_load((++cnt&1)?col:lor);
    gym::stretchAlignZ1(nxt, -rad, fib.dirEndP(), -rad);
    gle::endedTube();
    gym::disableClipPlane(4);
}

//------------------------------------------------------------------------------
#pragma mark - Specific Fiber styles

/**
 Renders a protofilament by drawing spheres of alternating colors,
 at distance `inc` from each other along the backbone of the `Fiber`.
 */
void Display::drawFilament(Fiber const& fib, const real inc,
                           gym_color col, gym_color lor, gym_color colE) const
{
    // enlarge radius of monomers to make them overlap
    const float rad = 0.65 * inc;
    gym::enableClipPlane(4);
    
    int cnt = (int)std::ceil(fib.abscissaM()/inc);
    real abs = inc * cnt;
    Vector3 p(fib.pos(abs)), q;
    // draw a monomer every 'inc'
    while ( abs < fib.abscissaP() )
    {
        q = p;
        abs += inc;
        p = Vector3(fib.pos(abs));

        // alternate tones:
        gym::color_load((++cnt&1)?col:lor);
        // change color for the last monomer:
        if ( abs + inc > fib.abscissaP() )
        {
            gym::color_front(colE);
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
void Display::drawActin(Fiber const& fib, gym_color col, gym_color lor, gym_color colE) const
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
    // rotate until we reach the minus end
    while ( ab <= fib.abscissaM() )
    {
        ++cnt;
        ab += dab;
        n = d.rotateOrtho(n, cs, sn);
    }
    p = Vector3(fib.pos(ab)) + off * n;
    // draw the monomers until the plus end:
    while ( ab < fib.abscissaP() )
    {
        q = p;
        ab += dab;
        d = Vector3(fib.dir(ab));
        n = d.rotateOrtho(n, cs, sn);
        p = Vector3(fib.pos(ab)) + off * n;
        
        // use different tones to individualize the two strands:
        gym::color_load((++cnt&1)?col:lor);

        // change color for the last monomer:
        if ( ab + dab > fib.abscissaP() )
        {
            gym::color_load(colE);
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
void Display::drawMicrotubule(Fiber const& fib, gym_color colA, gym_color colB, gym_color colE) const
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

        gym::color_load(colA);
        for ( int i = 0; i < 13; ++i )
            drawMonomer(p+dx[i]*d+dy[i]*e+dz[i]*f, rad);

        gym::color_load(colB);
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
                gym::color_load(colA);
                drawMonomer(p+dx[i]*d+dy[i]*e+dz[i]*f, rad);
                if ( ab+5.2*sa+dx[i] < abmax )
                    gym::color_load(colB);
                else
                    gym::color_load(colE);
                drawMonomer(p+(sa+dx[i])*d+dy[i]*e+dz[i]*f, rad);
            }
        }
        ab += dab;
    }
}


//------------------------------------------------------------------------------
#pragma mark - Draw Fiber


void Display::drawFiber(Fiber const& fib)
{
    FiberDisp const*const dis = fib.prop->disp;
    int style = dis->line_style;
    
#if FIBER_HAS_LATTICE || FIBER_HAS_MESH
    if ( dis->lattice_style )
    {
        VisibleLattice const* lat = fib.visibleLattice();
        if ( lat )
        {
            const float rad = dis->line_width;
            switch ( dis->lattice_style )
            {
                case 1: drawFiberLattice1(fib, *lat, rad); break;
                case 2: drawFiberLattice2(fib, *lat, rad); break;
                case 3: drawFiberLattice3(fib, *lat, rad); break;
                case 4: drawFiberLatticeEdges(fib, *lat, rad); break;
            }
            style = 0;
        }
    }
#endif

    if ( dis->style == 1 )
        drawFiberBackbone(fib, fib.disp->color, prop->bone_width);
    // if the Lattice was displayed, do not draw fancy styles:
    else if ( style && dis->style )
    {
        gym_color col1 = fib.disp->color.alpha(1.0);
        gym_color col2 = col1.darken(0.75);
        gym_color colE = fib.disp->end_color[0];
        
        // load backface color:
        if ( dis->coloring )
            gym::color_back(col1);
        else
            gym::color_back(dis->back_color);
        gym::enableLighting();

        switch( dis->style )
        {
            case 2: drawFiberStriped(fib, pixscale(dis->line_width), 0.008, col1, 0.024, col2); break;
            case 3: drawFilament(fib, 0.008, col1, col2, colE); break;
            case 4: drawActin(fib, col1, col2, colE); break;
            case 5: drawMicrotubule(fib, col1, col2, colE); break;
            case 6: drawFiberStriped2D(fib, pixscale(dis->line_width), 0.008, col1, 0.024, col2); break;
            case 7: {
                float w = pixscale(dis->line_width);
                float l = std::max(w, 0.008f);
                drawFiberArrowed2D(fib, w, 4*l, col1, l, col2);
                //drawFiberWide(fib, w);
            } break;
        }
        style = 0;
    }

#if ( DIM >= 3 )
    /*
     Handle styles in 3D that are using transparency to draw fiber's segments
     */
    if (( style==1 && fib.disp->color.transparent()) || ( style==2 || style==3 ))
    {
        for ( size_t i = 0; i < fib.lastPoint(); ++i )
            zObjects.emplace(&fib, i);
        style = 0;
    }
    else if ( style == 6 )
    {
        // color according to the distance from the minus end
        for ( size_t i = 0; i < fib.lastPoint(); ++i )
            if ( color_by_distanceM(fib, i).visible() )
                zObjects.emplace(&fib, i);
        style = 0;
    }
    else if ( style == 7 || ( style == 8 && fib.endStateP() == STATE_GREEN ))
    {
        // color according to the distance from the plus end
        for ( size_t i = 0; i < fib.lastPoint(); ++i )
            if ( color_by_distanceP(fib, i+1).visible() )
                zObjects.emplace(&fib, i);
        style = 0;
    }
#endif

    if ( style && dis->line_width > 0 )
    {
        /*
        if ( style == 1 ) // adding a black outline to filaments
            drawFiberBackbone(fib, gym::background_color, 2*dis->line_width);
         */
        drawFiberLines(fib, style, dis->line_width);
    }
    
    if ( dis->end_style[0] )
    {
        gym::color_load(fib.disp->end_color[0]);
        //gym::color_load(fib.disp->color);
        gym::color_back(dis->back_color);
        drawFiberEndPlus(fib, dis->end_style[0], dis->end_size[0]);
    }

    if ( dis->end_style[1] )
    {
        gym::color_load(fib.disp->end_color[1]);
        //gym::color_load(fib.disp->color);
        gym::color_back(dis->back_color);
        drawFiberEndMinus(fib, dis->end_style[1], dis->end_size[1]);
    }

    if ( dis->growth_style )
        drawFiberGrowth(fib, dis->line_width);
    
    if ( dis->point_style > 0 )
        drawFiberPoints(fib);
    
    if ( dis->speckle_style > 0 )
        drawFiberSpeckles(fib);

    // draw other fiber elements only if fiber is fully visible:
    if ( dis->force_style && fib.disp->visible > 0 )
        drawFiberForces(fib, dis->force_scale, pixwidth(dis->point_size));
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
        if ( fib->disp && fib->disp->visible )
            drawFiber(*fib);
    }
    CHECK_GL_ERROR("in Display::drawFibers()");
}


void Display::drawFiberTexts(FiberSet const& set)
{
    gym::ref_view();
    gym::cancelRotation();
    gym::disableLighting();
    gym::disableAlphaTest();
    // display Fibers in a random (ever changing) order:
    for ( Fiber const* fib = set.first(); fib ; fib=fib->next() )
    {
        if ( fib->disp && fib->disp->visible > 0 )
        {
            FiberDisp const*const dis = fib->prop->disp;
#if FIBER_HAS_LATTICE || FIBER_HAS_MESH
            if ( dis->lattice_style == 5 )
            {
                VisibleLattice const* lat = fib->visibleLattice();
                if ( lat )
                {
                    gym::color(fib->disp->color);
                    drawFiberLatticeValues(*fib, *lat);
                }
            }
#endif
#if FIBER_HAS_LATTICE
            if ( dis->lattice_style == 6 )
            {
                FiberLattice const* lat = fib->lattice();
                if ( lat )
                {
                    gym::color(fib->disp->color);
                    drawFiberLatticeBits(*fib, *lat); break;
                }
            }
#endif
            if ( dis->label_style )
            {
                gym::color(fib->disp->color);
                drawFiberLabels(*fib, dis->label_style);
            }
        }
    }
    CHECK_GL_ERROR("in Display::drawFiberTexts()");
    gym::restoreAlphaTest();
}

//------------------------------------------------------------------------------
#pragma mark - Couples


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
        if ( !cx->fiber1()->disp->visible && !cx->fiber2()->disp->visible )
            continue;
        
        // do not display Couple if both hands are hidden
        if ( !cx->disp1()->visible && !cx->disp2()->visible )
            continue;
        
        drawCoupleB(cx);
    }
}

//------------------------------------------------------------------------------
#pragma mark - Solid

void Display::drawSolid(Solid const& obj)
{
    PointDisp const* dis = obj.prop->disp;
    gym_color col = bodyColorF(obj);
#if NEW_SOLID_HAS_TWIN
    if ( obj.elder_twin() )
        col = bodyColorF(*obj.twin()).tweak(obj.signature());
#endif

    //display points:
    if (( dis->style & 2 ) && dis->perceptible )
    {
        const float rad = pixscale(dis->size);
        gym::color_both(col);
        gym::enableLighting();
        for ( size_t i = 0; i < obj.nbPoints(); ++i )
        {
            if ( 0 < obj.hasTriad(i) )
            {
                Vector P = obj.posP(i);
                gym::transScale(P, rad);
                gle::star();
                for ( int d = 1; d <= DIM; ++d )
                {
                    Vector Q = obj.posP(i+d);
                    gym::transAlignZ(Q, rad, Q-P);
                    gle::coneC();
                }
                i += DIM;
                continue;
            }
            drawObject(obj.posP(i), rad, gle::hedron(obj.radius(i)>0));
        }
    }
    
#if NEW_SOLID_CLAMP
    if ( obj.clampStiffness() > 0 )
    {
        gym::color_both(col);
        gym::enableLighting();
        drawObject(obj.clampPosition(), pixscale(dis->size), gle::star);
    }
#endif
    
#if NEW_SOLID_HAS_TWIN
    // display links between twin solids
    if ( obj.elder_twin() && ( obj.prop->twin_stiffness > 0 ) )
    {
        Solid const* twi = obj.twin();
        gym_color lor = bodyColorF(*twi).mix(col);
        real rad = M_SQRT2 * pixscale(dis->size);
        gym::ref_view();
        gym::enableLighting();
        gym::disableCullFace();
        const unsigned off = DIM+1; // index of point to be linked
        const unsigned sup = std::min(obj.lastPoint(), 2*off+1); // last point to be linked
        for ( unsigned inx = off; inx <= sup; inx += off+1 )
        {
            gym::stretchAlignZ(obj.posPoint(inx), twi->posPoint(inx), rad);
            gym::color_back(lor.darken(0.5));
            gym::color_front(lor, 1);
            gle::helix();
        }
        gym::restoreCullFace();
    }
#endif

#if OLD_SOLID_HAS_TWIN
    // display links between twin solids
    if ( obj.elder_twin() && ( dis->style & 6 ) && dis->perceptible )
    {
        Solid const* twi = obj.twin();
        gym_color lor = bodyColorF(*twi);
        // draw links between 'obj' and its Twin
        real rad = M_SQRT2 * pixscale(dis->size);
        const size_t inx = 0;
        // three points are expected:
        if ( obj.radius(inx) > 0 && obj.nbPoints() > inx+4 )
        {
            gym::enableLighting();
            gym::color_both(col, 1);
            Vector S(obj.posPoint(inx));
            Vector A(obj.posPoint(inx+1));
            Vector B(obj.posPoint(inx+2));
            Vector T(twi->posPoint(inx));
            Vector X(twi->posPoint(inx+1));
            Vector Y(twi->posPoint(inx+2));
            Vector S3(obj.posPoint(inx+DIM+1));
            Vector T3(twi->posPoint(inx+DIM+1));
#if ( DIM >= 3 )
            real R = 3 * pixscale(dis->size) / obj.radius(inx);
            Vector C(obj.posPoint(inx+3));
            Vector dA = R*(A-S), dB = R*(B-S), dC = R*(C-S);
            Vector Z(twi->posPoint(inx+3));
            Vector dX = R*(X-T), dY = R*(Y-T), dZ = R*(Z-T);
#endif
            gym::ref_view();
            if ( obj.prop->twin_stiffness > 0 )
            {
                gym::color_both(col.mix(lor), 1);
                gym::stretchAlignZ(S3, T3, rad);
                gle::hexTube();
#if ( DIM >= 4 )
                gle::paintSpikyPrism(dX, dY, dZ, X, dA, dB, dC, A);
                gle::paintSpikyPrism(dX, dY, dZ, Y, dA, dB, dC, B);
                gle::paintSpikyPrism(dX, dY, dZ, Z, dA, dB, dC, C);
#elif ( DIM == 2 )
                gle::drawSpikyBand(A, X, rad);
                gle::drawSpikyBand(B, Y, rad);
#endif
            }
            else
            {
#if ( DIM >= 3 )
                // A and B are swapped since the triad is inverted
                gym::color_both(col, 1);
                gle::paintTetrahedron(dB, dA, dC, A);
                gle::paintTetrahedron(dB, dA, dC, B);
                gle::paintTetrahedron(dB, dA, dC, C);
                // draw the other guy's handles
                gym::color_both(lor, 1);
                gle::paintTetrahedron(dX, dY, dZ, X);
                gle::paintTetrahedron(dX, dY, dZ, Y);
                gle::paintTetrahedron(dX, dY, dZ, Z);
#elif ( DIM == 2 )
                Vector2 N = 2 * S - ( A + B );
                gym::transAlignZ(A, rad, N); gle::pyramid();
                gym::transAlignZ(B, rad, N); gle::pyramid();
                Vector2 D = 2 * T - ( A + B );
                gym::transAlignZ(X, rad, D); gle::pyramid();
                gym::transAlignZ(Y, rad, D); gle::pyramid();
#endif
            }
        }
    }
#endif

    //display outline of spheres in 2D
    if ( dis->style & 8 )
    {
        gym::color(col);
#if ( DIM == 2 )
        gym::disableLighting();
        for ( size_t i = 0; i < obj.nbPoints(); ++i )
        {
            if ( obj.radius(i) > pixelSize )
            {
                gym::transScale(obj.posP(i), obj.radius(i));
                gle::circle1(dis->width);
            }
        }
        gym::enableLighting();
#elif ( DIM >= 3 )
        //special display for ParM simulations (DYCHE 2006; KINETOCHORES 2019)
        if ( obj.mark()  &&  obj.nbPoints() > 1 )
        {
            gym::enableLighting();
            //drawObject(obj.posP(0), obj.diffPoints(1, 0), obj.radius(0), gle::circle);
            gym::stretchAlignZ(obj.posP(0), obj.posP(1), obj.radius(0));
            gle::cylinder1();
        }
#endif
    }
    
    //print the number for each Solid
    if ( dis->style & 16 )
    {
        char tmp[32];
        gym::color(col);
        snprintf(tmp, sizeof(tmp), "0:%u", obj.identity());
        strokeText(obj.posP(0), tmp, pixelSize);
        for ( size_t i = 1; i < obj.nbPoints(); ++i )
        {
            snprintf(tmp, sizeof(tmp), "%lu", i);
            strokeText(obj.posP(i), tmp, pixelSize);
        }
    }
    
    //draw polygon line joining vertices of Solid
    if ( dis->style & 32 )
    {
        gym::disableLighting();
        gym::color(col);
        gym::ref_view();
        gym::loadPoints(obj.nbPoints(), obj.addrPoints());
        gym::drawLineStrip(dis->widthX, 0, obj.nbPoints());
        gym::enableLighting();
    }
}

/**
 Display a semi-transparent disc / sphere
 */
void Display::drawSolidT(Solid const& obj, unsigned inx) const
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
    gym_color col = bodyColorF(obj);
#if NEW_SOLID_HAS_TWIN
    if ( obj.elder_twin() )
        col = bodyColorF(*obj.twin()).tweak(obj.signature());
#endif
    col.set_alpha(obj.prop->disp->color.alpha());
#if ( DIM > 2 )
    gym::color_both(col);
    drawBallT(X, obj.radius(inx), obj.mark());
#else
    gym::color(col);
    drawDiscT(X, obj.radius(inx));
#endif
    for ( size_t i = 0; i < num; ++i )
        gym::disableClipPlane(5-i);
}


void Display::drawSolids(SolidSet const& set)
{
    for ( Solid * obj = set.first(); obj; obj=obj->next() )
    {
        PointDisp const* dis = obj->prop->disp;
        if ( dis->visible && dis->style )
        {
            drawSolid(*obj);
            if ( dis->style & 1 )
            {
                // draw Solid's balls which can be transparent or not
                size_t sup = obj->nbPoints();
                if ( dis->style & 4 ) sup = 1;
#if ( DIM >= 3 )
#if NEW_SOLID_HAS_TWIN
                Solid const * twi = obj->twin();
                if ( twi )
                {
                    if ( twi != obj ) // only select elder twin
                    for ( size_t inx = 0; inx+DIM < obj->nbPoints(); inx += DIM+2 )
                    {
                        gym::enableLighting();
                        gym_color black(0,0,0,1);
                        gym_color tcol = bodyColorF(*twi);
                        gym_color scol = tcol.tweak(obj->signature());
                        drawFootball(*obj, inx, scol, black, true);
                        drawFootball(*twi, inx, tcol, black, false);
                    }
                } else
#endif
                if ( dis->color.transparent() )
                {
                    for ( size_t i = 0; i < sup; ++i )
                        if ( obj->radius(i) > 0 )
                            zObjects.emplace(obj, i);
                }
                else
#endif
                {
                    for ( size_t i = 0; i < sup; ++i )
                        if ( obj->radius(i) > 0 )
                            drawSolidT(*obj, i);
                }
            }
        }
    }
}

//------------------------------------------------------------------------------
#pragma mark - Beads

void drawVector(Vector pos, Vector vec, real a, real b, gym_color col, float rad)
{
    gym::color_both(col);
    Vector A = pos - a * vec;
    Vector B = pos + b * vec;
    gle::drawArrow(A, B, rad);

}

void Display::drawBead(Bead const& obj)
{
    PointDisp const* dis = obj.prop->disp;
    gym_color col = bodyColorF(obj);
    
    // display center:
    if ( dis->style & 2 )
    {
        gym::color_both(col, 1);
        drawObject(obj.position(), pixscale(dis->size), gle::tetrahedron);
    }
    
#if NEW_SOLID_CLAMP
    if ( obj.clampStiffness() > 0 )
    {
        gym::color_both(col);
        gym::enableLighting();
        drawObject(obj.clampPosition(), pixscale(dis->size), gle::star);
    }
#endif

    if ( dis->style & 4 )
    {
#if ( DIM <= 2 )
        // display outline:
        if ( obj.radius() > pixelSize )
        {
            gym::disableLighting();
            gym::color(col);
            gym::transScale(obj.position(), obj.radius());
            gle::circle1(dis->width);
        }
#else
        if ( obj.mark() )
        {
            col = gym::get_color(obj.mark());
            gym::color_both(col);
            drawObject(obj.position(), pixscale(dis->size), gle::cube);
        }
#endif
    }
}


/**
 Display a semi-transparent disc / sphere
 */
void Display::drawBeadT(Bead const& obj) const
{
    PointDisp const* dis = obj.prop->disp;
    
    assert_true( dis->style & 1 );
    {
        gym_color col = bodyColorF(obj).match_a(dis->color);
#if ( DIM > 2 )
        gym::color_both(col);
        drawBeadS(obj.position(), obj.radius(), obj.mark());
#else
        gym::color(col);
        drawDiscT(obj.position(), obj.radius());
#endif
    }
}


void Display::drawBeads(BeadSet const& set)
{
    gym::enableLighting();
    for ( Bead * obj = set.first(); obj; obj=obj->next() )
    {
        PointDisp const* dis = obj->prop->disp;
        if ( dis->visible && dis->style )
        {
            drawBead(*obj);
            if ( dis->style & 1 )
            {
#if ( DIM >= 3 )
                if ( dis->color.transparent() )
                    zObjects.emplace(obj);
                else
#endif
                    drawBeadT(*obj);
            }
        }
    }
}


//------------------------------------------------------------------------------
#pragma mark - Sphere

void Display::drawSphere(Sphere const& obj)
{
    PointDisp const* dis = obj.prop->disp;
    
    if ( dis->perceptible )
    {
        const float rad = pixscale(dis->size);
        // display surface points
        if ( dis->style & 2 )
        {
            bodyColor(obj);
            for ( size_t i = obj.nbRefPoints; i < obj.nbPoints(); ++i )
                drawObject(obj.posP(i), rad, gle::sphere1);
        }
        
        // display center and reference points
        if ( dis->style & 8 )
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
    PointDisp const* dis = obj.prop->disp;

    if ( dis->style & 7 )
    {
        gym_color col = bodyColorF(obj).match_a(dis->color);
        const Vector C = obj.posP(0);
#if ( DIM < 3 )
        if ( obj.radius() > pixelSize )
        {
            gym::color(col);
            gym::transScale(C, obj.radius());
            if ( dis->style & 1 )
                gle::circle1(dis->widthX);
            if ( dis->style & 2 )
                gle::disc();
            if ( dis->style & 4 )
                drawDiscT(C, obj.radius());
        }
#else
        /* Note: The rotation matrix for the sphere calculated below from the
         reference points, includes scaling by the radius of the sphere.
         We then use a primitive for a sphere of radius 1.
         */
        gym::color_both(col);
        gym::enableLighting();
        gym::transRotate(C, obj.posP(1)-C, obj.posP(2)-C, obj.posP(3)-C);
        if ( dis->style & 1 )
            gle::dualPassSphere2();
        if ( dis->style & 4 )
            gle::threeArrowStrip(0.5, 1);
#endif
    }
}


void Display::drawSpheres(SphereSet const& set)
{
    for ( Sphere * obj=set.first(); obj ; obj=obj->next() )
    {
        PointDisp const* dis = obj->prop->disp;
        if ( dis->visible && dis->style )
        {
            if ( dis->perceptible )
                drawSphere(*obj);
#if ( DIM >= 3 )
            if ( dis->color.transparent() )
                zObjects.emplace(obj);
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
    Solid const* sol = obj.solid();
    Sphere const* sph = obj.sphere();
    if ( !sol && !sph ) return;
    PointDisp const* dis = ( sol ? sol->prop->disp : sph->prop->disp );
    if ( !dis ) return;
    gym_color col;
    if ( sol ) col = bodyColorF(*sol).match_a(dis->color);
    if ( sph ) col = bodyColorF(*sph).match_a(dis->color);

    size_t i = 0, cnt = obj.nbLinks();

    if ( dis->style & 2 )
    {
        // draw links between Solid/Sphere and Mecables
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
        gym::color(col);
        gym::ref_view();
        gym::disableLighting();
        gym::drawLines(dis->widthX, 0, 2*i);
        gym::rebindBufferVD(2);
        gym::drawPoints(dis->sizeX, 0, i);
    }

    /**
     This display the Solid connecting two Aster as a spindle.
     Used for Cleo Kozlowski simulation of C. elegans (2007)
     */
    if ( dis->style & 1 && obj.tag() == Organizer::FAKE_TAG )
    {
        if ( sol && sol->nbPoints() >= 4 )
        {
#if ( DIM >= 3 )
            gym::enableLighting();
            gym::color_front(col);
            gym::color_back(col.darken(0.5));
            Vector3 a = 0.5*(sol->posP(0) + sol->posP(2));
            Vector3 b = 0.5*(sol->posP(1) + sol->posP(3));
            gym::stretchAlignZ(a, b, 1);
            gle::dualPassBarrel();
#else
            gym::disableLighting();
            gym::color(col);
            gym::ref_view();
            gym::loadPoints(sol->nbPoints(), sol->addrPoints());
            gym::drawLines(dis->widthX, 0, sol->nbPoints());
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


/// set depth relative to given axis
void zObject::calculate_depth(Vector const& axis)
{
    Mecable const * mec = point_.mecable();
    if ( mec->tag() == Fiber::TAG && point_.point() < mec->lastPoint() )
        depth_ = dot(mec->midPoint(point_.point()), axis);
    else
        depth_ = dot(point_.pos(), axis);
}


#if ( DIM >= 3 )

/// display sub-part `inx` of object `obj`
void zObject::draw(Display const* dis) const
{
    assert_false( point_.invalid() );
    Mecable const * mec = point_.mecable();
    switch( mec->tag() )
    {
        case Fiber::TAG:
            dis->drawFiberSegmentT(*static_cast<const Fiber*>(mec), point_.point());
            break;
            
        case Solid::TAG:
            dis->drawSolidT(*static_cast<const Solid*>(mec), point_.point());
            break;
            
        case Bead::TAG:
            dis->drawBeadT(*static_cast<const Bead*>(mec));
            break;
            
        case Sphere::TAG:
            dis->drawSphereT(*static_cast<const Sphere*>(mec));
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
        i.calculate_depth(depthAxis);
    
    // depth-sort objects:
    list.quick_sort(compareZObject);
    //std::clog << " depth sorted " << list.size() << " zObjects axis: "<< depthAxis << "\n";

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


