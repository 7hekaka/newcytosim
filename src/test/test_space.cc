// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
/*
 test_space provides a visual test of Cytosim's Space
*/

#include <ctime>
#include "dim.h"
#include "timer.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "glossary.h"
#include "real.h"
#include "vector.h"
#include "random.h"

#include "space_prop.h"
#include "space.h"
#include "space_set.h"
#include "glapp.h"
#include "glut.h"
#include "gle.h"
#include "gym_color.h"
#include "gym_flute.h"
#include "gym_draw.h"
#include "gym_view.h"
#include "gym_check.h"
#include "gym_cap.h"
#include "point_disp.h"

// List of options
Glossary opt;

// property
SpaceProp prop("test_space");

// Space to be tested:
Space * spc = nullptr;

// number of points
const size_t maxpts = 1<<17;
size_t n_pts = 1024;
size_t n_bin = 16;

// INFLATION of the rectangle containing point to be projected
const real INFLATION = 1;

// parameters for planar_distribution
const real grain = 0.1;
real thickness = grain;

// regular or random distribution of the test-points
bool regular_distribution = false;
bool points_around_mouse = false;
bool points_on_edges = false;
int planar_distribution = 0;

//coordinates of the points:
Vector point[maxpts];

//true if inside
bool inside[maxpts];

//coordinates of the projections
Vector project[maxpts];

//coordinates of the projections of the projections
Vector project2[maxpts];

//normals to the projections
Vector normal[maxpts];

// point + normal
Vector upward[maxpts];

//max distance from projection to second projection
real error = 0;

real intercept = 0;
Vector3 origin(0, 0, 0);
Vector3 axis(1, 0, 0);
Vector3 mouse(0,0,0);

//show or hide points in or outside
int showInside    = true;
int showOutside   = true;
int showProjection = true;
int showProjected = true;
int showReproject = true;
int showNormals   = false;
int showSpace     = 3;

//use timer function on or off
int timerOn = false;
int timerDelay = 50;

//display parameter for OpenGL
float LW = 2.0f;
float PS = 3.0f;

//amount of white added to colors
const GLfloat COL = 0.8f;

//------------------------------------------------------------------------------

void generatePoints(real len)
{
    Vector inf, dif;
    spc->boundaries(inf, dif);
    inf -= Vector(len, len, len);
    dif += Vector(len, len, len) - inf;
    
    if ( points_around_mouse )
    {
        Vector m(mouse);
        for ( size_t i = 0; i <= n_pts; ++i )
            point[i] = m + thickness * Vector::randB();
    }
    else if ( points_on_edges )
    {
        for ( size_t i = 0; i <= n_pts; ++i )
            point[i] = spc->placeOnEdge(1);
    }
    else if ( regular_distribution )
    {
        dif /= n_bin;
        size_t kk = 0;
        n_pts = 0;
        //follow a regular lattice:
        for ( size_t ii = 0; ii <= n_bin; ++ii )
        for ( size_t jj = 0; jj <= n_bin; ++jj )
#if ( DIM >= 3 )
        for ( kk = 0; kk <= n_bin; ++kk )
#endif
        {
            point[n_pts++] = inf + dif.e_mul(Vector(ii, jj, kk));
            if ( n_pts >= maxpts )
                return;
        }
    }
    else if ( planar_distribution )
    {
        for ( size_t i = 0; i <= n_pts; ++i )
        {
            point[i] = inf + dif.e_mul(Vector::randP());
            real a = intercept + RNG.shalf() * thickness;
            switch ( planar_distribution )
            {
                case 1: point[i].XX = a; break;
                case 2: point[i].YY = a; break;
#if ( DIM >= 3 )
                case 3: point[i].ZZ = a; break;
#endif
            }
        }
    }
    else
    {
        for ( size_t i = 0; i <= n_pts; ++i )
            point[i] = inf + dif.e_mul(Vector::randP());
        //point[ii] = Vector::randU();
        //point[ii] = spc->placeNearEdge(0.1);
    }
}


void calculateNormals()
{
    for ( size_t i = 0; i < n_pts; ++i )
    {
        normal[i] = spc->normalToEdge(project[i]);
        upward[i] = project[i] + normal[i];
    }
}


void distributePoints(real len = INFLATION)
{
    generatePoints(len);
    error = 0;
    
    for ( size_t i = 0; i < n_pts; ++i )
    {
        normal[i].reset();
        //see if space finds it inside:
        inside[i] = spc->inside(point[i]);
        //calculate the projection:
        project[i] = spc->project(point[i]);
        
        //calculate the projection of the projection:
        project2[i] = spc->project(project[i]);
        
        real d = ( project[i] - project2[i] ).normSqr();
        error = std::max(d, error);
    }
    error = std::sqrt(error);
    if ( showNormals )
        calculateNormals();
    
    char tmp[128];
    snprintf(tmp, sizeof(tmp), "Projection error %.9f", error);
    glApp::setMessage(tmp);
}

//------------------------------------------------------------------------------
void timerFunction(int)
{
    if ( timerOn )
    {
        distributePoints();
        glApp::postRedisplay();
        glutTimerFunc(timerDelay, timerFunction, 0);
    }
}

//------------------------------------------------------------------------------
void checkVolume(size_t CNT)
{
    real vol = spc->volume();
    double avg = 0, dev = 0;
    
    for ( size_t i = 0; i < CNT; ++i )
    {
        real e = spc->estimateVolume(1<<21) - vol;
        //printf("Monte-Carlo estimated volume %.6f\n", e+vol);
        avg += e;
        dev += e * e;
    }
    avg /= CNT;
    dev = ( dev - avg * avg * CNT ) / (CNT-1);
    dev = std::sqrt(max_real(0, dev));
    
    printf("Monte-Carlo estimated volume of `%s` is", spc->prop->shape.c_str());
    printf("  %.3f +/- %.3f;  given volume is %.3f", avg+vol, dev, vol);
    
    if ( abs_real(avg) > 3*dev )
         printf("WARNING: POSSIBLE VOLUME MISMATCH!!!!\n");
}


void setGeometry()
{
    if ( ! prop.disp )
        prop.disp = new PointDisp("space:display", "test_space");
    prop.read(opt);
    
    try {
        delete(spc);
        spc = prop.newSpace(opt);
        if ( spc )
        {
            checkVolume(8);
            Outputter out(stdout, false);
            spc->write(out);
            fprintf(stdout, "\n");
        }
    }
    catch( Exception & e )
    {
        printf("Error: `%s'\n", e.what());
    }
    
    try {
        if ( spc )
            distributePoints(INFLATION);
    }
    catch( Exception & e )
    {
        printf("Error: `%s'\n", e.what());
    }

    glApp::postRedisplay();
}

//------------------------------------------------------------------------------
enum MENUS_ID {
    MENU_QUIT = 102, MENU_RESETVIEW = 103,
    MENU_INSIDE = 104, MENU_OUTSIDE = 105, MENU_PROJECT = 106,
    MENU_XSLICING = 107, MENU_YSLICING = 108, MENU_ZSLICING = 109,
    MENU_PROJECTED = 110, MENU_EDGES = 111
};


void toggleSlicing(int d)
{
    if ( planar_distribution == d )
        planar_distribution = 0;
    else
        planar_distribution = d;
    switch ( planar_distribution )
    {
        case 0: break;
        case 1: axis.set(1,0,0); break;
        case 2: axis.set(0,1,0); break;
        case 3: axis.set(0,0,1); break;
    }
    distributePoints();
}


void processMenu(int item)
{
    switch( item )
    {
        case MENU_QUIT:
            exit(EXIT_SUCCESS);
        case MENU_RESETVIEW:
            glApp::resetView();
            break;
        case MENU_INSIDE:
            showInside = ! showInside;
            break;
        case MENU_OUTSIDE:
            showOutside = ! showOutside;
            break;
        case MENU_EDGES:
            points_on_edges = ! points_on_edges;
            distributePoints();
            break;
        case MENU_PROJECT:
            showProjection = ! showProjection;
            break;
        case MENU_PROJECTED:
            showProjected = ! showProjected;
            break;
        case MENU_XSLICING:
            toggleSlicing(1);
            break;
        case MENU_YSLICING:
            toggleSlicing(2);
            break;
        case MENU_ZSLICING:
            toggleSlicing(3);
            break;
    }
    glApp::postRedisplay();
}


void initMenus()
{
    int gm = glApp::buildMenu();
    glutCreateMenu(processMenu);
    glutAddSubMenu("Control", gm);
    
    glutAddMenuEntry("Reset",                MENU_RESETVIEW);
    glutAddMenuEntry("Quit",                 MENU_QUIT);
    glutAddMenuEntry("-", 0);
    glutAddMenuEntry("Toggle inside  (i)",   MENU_INSIDE);
    glutAddMenuEntry("Toggle outside (o)",   MENU_OUTSIDE);
    glutAddMenuEntry("Toggle edges   (e)",   MENU_EDGES);
    glutAddMenuEntry("Toggle project (p)",   MENU_PROJECT);
    glutAddMenuEntry("Toggle projected (s)", MENU_PROJECTED);

    glutAddMenuEntry("Toggle x-slicing (x)", MENU_XSLICING);
    glutAddMenuEntry("Toggle y-slicing (y)", MENU_YSLICING);
    glutAddMenuEntry("Toggle z-slicing (z)", MENU_ZSLICING);
    
    glutAttachMenu(GLUT_RIGHT_BUTTON);
}

//------------------------------------------------------------------------------


///set callback for shift-click, with unprojected click position
void processMouseClick(int, int, const Vector3 & a, int)
{
    origin = a;
    if ( points_around_mouse )
    {
        distributePoints();
        glApp::postRedisplay();
    }
}

///set callback for shift-drag, with unprojected mouse and click positions
void processMouseDrag(int, int, Vector3 & a, const Vector3 & b, int)
{
    origin = a;
    mouse = b;
    if ( points_around_mouse )
    {
        distributePoints();
        glApp::postRedisplay();
    }
}


void processSpecialKey(int key, int x=0, int y=0)
{
    switch (key)
    {
        case GLUT_KEY_LEFT:
            intercept -= grain;
            break;
        case GLUT_KEY_RIGHT:
            intercept += grain;
            break;
        case GLUT_KEY_UP:
            thickness += grain;
            break;
        case GLUT_KEY_DOWN:
            thickness = std::max(grain, thickness-grain);
            break;
        default:
            break;
    }
    distributePoints();
    glApp::postRedisplay();
}

void processNormalKey(unsigned char c, int x=0, int y=0)
{
    switch (c)
    {
        case 27:
        case 'q':
            exit(EXIT_SUCCESS);
            
        case ' ':
            distributePoints();
            break;
            
        case '0':
            glApp::resetView();
            break;
            
        case ']':
            n_bin *= 2;
            n_pts *= 2;
            if ( n_pts > maxpts )
                n_pts = maxpts;
            distributePoints();
            break;
            
        case '[':
            if ( n_bin > 2 ) n_bin /= 2;
            if ( n_pts > 2 ) n_pts /= 2;
            distributePoints();
            break;
            
        case 'x':
            toggleSlicing(1);
            break;
            
        case 'y':
            toggleSlicing(2);
            break;
            
        case 'z':
            toggleSlicing(3);
            break;

        case 'i':
            showInside = ! showInside;
            break;
            
        case 'm':
            points_around_mouse = ! points_around_mouse;
            break;
            
        case 'o':
            showOutside = ! showOutside;
            break;
            
        case 'r':
            showReproject = ! showReproject;
            break;
            
        case 'p':
            showProjection = ! showProjection;
            break;
            
        case 's':
            showProjected = ! showProjected;
            break;

        case 'e':
            points_on_edges = ! points_on_edges;
            distributePoints();
            break;
            
        case 'n':
            showNormals = ! showNormals;
            if ( showNormals ) calculateNormals();
            break;
            
        case 'R':
            regular_distribution = !regular_distribution;
            distributePoints();
            break;
            
        case 'd':
            showSpace = ( showSpace + 1 ) & 3;
            glApp::flashText("showSpace = %i", showSpace);
            break;

        case 'f':
        {
            real val[] = { -2, -1, 0, 1, 2, 5 };
            opt.define("inflate", val[ RNG.pint32(6) ]);
            setGeometry();
        } break;
            
        case '9': PS += 1; break;
        case '8': PS = std::max(1.f, PS-1); break;
        
        case 't':
            timerOn = ! timerOn;
            if ( timerOn )
                glutTimerFunc(timerDelay, timerFunction, 0);
            break;
            
        default:
            glApp::processNormalKey(c,x,y);
    }
    glApp::postRedisplay();
}

//------------------------------------------------------------------------------

bool visible(size_t i)
{
    if ( inside[i] )
    {
        if ( !showInside ) return false;
    }
    else
    {
        if ( !showOutside ) return false;
    }
    return true;
}

void drawLines(Vector const A[], gym_color const& col, Vector const B[], gym_color const& lor)
{
    flute8* flu = gym::mapBufferC4V4(2*n_pts);
    size_t n = 0;
    for ( size_t i = 0; i < n_pts; ++i )
    {
        if ( visible(i) )
        {
            flu[n++] = { col, A[i] };
            flu[n++] = { lor, B[i] };
        }
    }
    gym::unmapBufferC4V4();
    gym::drawLines(2*LW, 0, n);
    gym::cleanupCV();
}


int display(View& view)
{
    view.back_color.set(0,0,0,1);
    view.openDisplay();
    //gym::printCaps("space");

#if ( DIM >= 3 )
    if ( spc && ( showSpace & 2 ))
    {
        // draw flat back side
        gym::disableLighting();
        gym::color(0, 0, 0, 1);
        gym::enableCullFace(GL_FRONT);
        spc->draw3D();
    }
#endif
    if ( spc && ( showSpace & 1 ))
    {
        gym::ref_view();
#if ( DIM >= 3 )
        // draw transparent front side
        gym::enableLighting();
        gym::color_front(0, 0, 1, 0.2);
        gym::enableCullFace(GL_BACK);
        gym::closeDepthMask();
        spc->draw3D();
        gym::openDepthMask();
        gym::disableCullFace();
#else
        gym::color(1, 1, 1, 1);
        spc->draw2D(1);
#endif
    }
    
    flute8* flu = gym::mapBufferC4V4(2*n_pts);
    gym_color col(0.f, COL, 0.f), lor(0.f, 0.f, COL);
    size_t n = 0;
    for ( size_t i = 0; i < n_pts; ++i )
    {
        if ( visible(i) )
        {
            real d = 1 - std::tanh(distance(point[i], project[i]));
            gym_color c = inside[i] ? col : lor;
            flu[n++] = { c, project[i] };
            flu[n++] = { c.alpha(d), point[i] };
        }
    }
    gym::unmapBufferC4V4();
    gym::ref_view();
    gym::disableLighting();
    gym::drawPoints(PS, 0, 2*n);

    if ( showProjection )
    {
        gym::enableBlending();
        gym::drawLines(LW, 0, n);
    }
    gym::cleanupCV();

    if ( showNormals )
        drawLines(project, gym_color(1.f, 1.f, 1.f), upward, gym_color(1.f, 1.f, 1.f, 0.f));
    
    if ( showReproject )
        drawLines(project, gym_color(COL, 0.f, 0.f), project2, gym_color(COL, 0.5f, 0.5f));

    view.closeDisplay();
    return 0;
}



void speedTest(size_t cnt)
{
    tick();
    for ( size_t i = 0; i < cnt; ++i )
        distributePoints();
    printf(" %lu tests: cpu %5.0f\n", cnt, tock());
}

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
    glutInit(&argc, argv);
    glApp::setDimensionality(DIM);
    glApp::normalKeyFunc(processNormalKey);
    glApp::specialKeyFunc(processSpecialKey);
    glApp::actionFunc(processMouseClick);
    glApp::actionFunc(processMouseDrag);
    glApp::newWindow(display);
    glApp::setScale(20);
    gle::initialize();

    initMenus();
    RNG.seed();

    int mode = 1;
    if ( argc > 1 && isdigit(argv[1][0]) )
        mode = 2;
    if ( argc > mode )
    {
        if ( opt.read_strings(argc-mode, argv+mode) )
            return EXIT_FAILURE;
        setGeometry();
    }
    if ( ! spc )
    {
        printf("A geometry should be given in the command line, for example:\n");
        printf("    test_space shape=ellipse length=2,3,4\n");
        exit(EXIT_SUCCESS);
    }
    if ( mode == 2 )
    {
        speedTest(strtoul(argv[1], nullptr, 10));
    }
    else
        glutMainLoop();
}

