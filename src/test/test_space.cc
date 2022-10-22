// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
/*
 test_space provides a visual test of Cytosim's Space
*/

#include <ctime>
#include "dim.h"
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
size_t nbpts = 1024;
size_t scan  = 16;

// INFLATION of the rectangle containing point to be projected
const real INFLATION = 1;

// regular or random distribution of the test-points
bool regular_distribution = false;


//coordinates of the points:
Vector point[maxpts];

//true if inside
int    inside[maxpts];

//coordinates of the projections
Vector project[maxpts];

//coordinates of the projections of the projections
Vector project2[maxpts];

//normals to the projections
Vector normal[maxpts];

//located on the edges
Vector edge[maxpts];

//max distance from projection to second projection
real  error = 0;

//slicing parameters
int slicing = 0;
const real grain = 0.1;
real thickness = grain;

Vector origin(0, 0, 0);
Vector axis(1, 0, 0);

//show or hide points in or outside
int showInside    = true;
int showOutside   = true;
int showProject   = true;
int showProjected = true;
int showReproject = true;
int showNormals   = false;
int showEdges     = false;
int showSpace     = 3;

//use timer function on or off
int timerOn = false;
int timerDelay = 50;

//display parameter for OpenGL
GLfloat LW = 1.0f;

//amount of white added to colors
const GLfloat COL = 0.8f;

//------------------------------------------------------------------------------

void generatePoints(real len)
{
    Vector inf, dif;
    spc->boundaries(inf, dif);
    inf -= Vector(len, len, len);
    dif += Vector(len, len, len) - inf;
    
    if ( regular_distribution )
    {
        dif /= scan;
        size_t kk = 0;
        nbpts = 0;
        //follow a regular lattice:
        for ( size_t ii = 0; ii <= scan; ++ii )
        for ( size_t jj = 0; jj <= scan; ++jj )
#if ( DIM >= 3 )
        for ( kk = 0; kk <= scan; ++kk )
#endif
        {
            point[nbpts++] = inf + dif.e_mul(Vector(ii, jj, kk));
            if ( nbpts >= maxpts )
                return;
        }
    }
    else
    {
        for ( size_t ii = 0; ii <= nbpts; ++ii )
            point[ii] = inf + dif.e_mul(Vector::randP());
        //point[ii] = Vector::randU();
        //point[ii] = spc->placeNearEdge(0.1);
    }
}


void calculateNormals()
{
    for ( size_t ii = 0; ii < nbpts; ++ii )
        normal[ii] = spc->normalToEdge(project[ii]);
}


void distributePoints(real len = INFLATION)
{
    generatePoints(len);
    error = 0;
    
    for ( size_t ii = 0; ii < nbpts; ++ii )
    {
        normal[ii].reset();
        //see if space finds it inside:
        inside[ii] = spc->inside(point[ii]);
        //calculate the projection:
        project[ii] = spc->project(point[ii]);
        
        //calculate the projection of the projection:
        project2[ii] = spc->project(project[ii]);
        
        edge[ii] = spc->placeOnEdge(1);
        
        real d = ( project[ii] - project2[ii] ).normSqr();
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
    slicing = ( slicing == d ? 0 : d & 3 );
    switch ( d )
    {
        case 0: break;
        case 1: axis.set(1,0,0); break;
        case 2: axis.set(0,1,0); break;
        case 3: axis.set(0,0,1); break;
    }
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
            showEdges = ! showEdges;
            break;
        case MENU_PROJECT:
            showProject = ! showProject;
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
    glutAddMenuEntry("Toggle projected (s)",   MENU_PROJECTED);

    glutAddMenuEntry("Toggle x-slicing (x)", MENU_XSLICING);
    glutAddMenuEntry("Toggle y-slicing (y)", MENU_YSLICING);
    glutAddMenuEntry("Toggle z-slicing (z)", MENU_ZSLICING);
    
    glutAttachMenu(GLUT_RIGHT_BUTTON);
}

//------------------------------------------------------------------------------
void processSpecialKey(int key, int x=0, int y=0)
{
    switch (key)
    {
        case GLUT_KEY_LEFT:
            origin -= grain * axis;
            break;
        case GLUT_KEY_RIGHT:
            origin += grain * axis;
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
            scan *= 2;
            nbpts *= 2;
            if ( nbpts > maxpts )
                nbpts = maxpts;
            distributePoints();
            break;
            
        case '[':
            if ( scan > 2 ) scan /= 2;
            if ( nbpts > 2 ) nbpts /= 2;
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
                
        case 'c':
            toggleSlicing(slicing+1);
            break;

        case 'i':
            showInside = ! showInside;
            break;
            
        case 'o':
            showOutside = ! showOutside;
            break;
            
        case 'r':
            showReproject = ! showReproject;
            break;
            
        case 'p':
            showProject = ! showProject;
            break;
            
        case 's':
            showProjected = ! showProjected;
            break;

        case 'e':
            showEdges = ! showEdges;
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
        if ( ! showInside )
            return false;
    }
    else
    {
        if ( ! showOutside )
            return false;
    }
    
    if ( slicing )
    {
        //return ( abs_real(dot(project[i]-origin, axis)) < thickness );
        return ( abs_real(dot(point[i]-origin, axis)) < thickness );
    }
    return true;
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
        gym::restoreCullFace();
        gym::restoreLighting();
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
        gym::restoreCullFace();
        gym::restoreLighting();
#else
        gym::color(1, 1, 1, 1);
        gym::disableLighting();
        spc->draw2D(1);
#endif
    }
    gym::ref_view();
    if ( 1 )
    {
        //use green for points inside, magenta for point outside:
        flute8* flu = gym::mapBufferC4V4(nbpts);
        gym_color col(0.f, COL, 0.f), lor(0.f, 0.f, COL);
        size_t n = 0;
        for ( size_t i=0; i < nbpts; ++i )
        {
            if ( visible(i) )
                flu[n++] = { inside[i] ? col : lor, point[i] };
         }
        gym::unmapBufferC4V4();
        gym::drawPoints(2, 0, n);
        gym::cleanup();
    }
    
    if ( showProjected )
    {
        //use green for points inside, magenta for point outside:
        flute8* flu = gym::mapBufferC4V4(nbpts);
        gym_color col(COL, COL, 0.f), lor(COL, 0.f, COL);
        for ( size_t i=0; i < nbpts; ++i )
            flu[i] = { inside[i] ? col : lor, project[i] };
        gym::unmapBufferC4V4();
        gym::drawPoints(2, 0, nbpts);
        gym::cleanup();
    }

    if ( showProject )
    {
        flute8* flu = gym::mapBufferC4V4(2*nbpts);
        gym_color col(0.f, COL, 0.f), lor(0.f, 0.f, COL);
        size_t n = 0;
        for ( size_t i = 0; i < nbpts; ++i )
        {
            if ( visible(i) )
            {
                gym_color c = inside[i] ? col : lor;
                flu[n++] = { c, point[i] };
                flu[n++] = { c, project[i] };
            }
        }
        gym::unmapBufferC4V4();
        gym::drawLines(LW, 0, n);
        gym::cleanup();
    }
    
    if ( showNormals )
    {
        flute8* flu = gym::mapBufferC4V4(2*nbpts);
        gym_color col(1.f, 1.f, 1.f), lor(1.f, 1.f, 1.f, 0.f);
        size_t n = 0;
        for ( size_t i = 0; i < nbpts; ++i )
        {
            flu[n++] = { col, project[i] };
            flu[n++] = { lor, project[i]+normal[i] };
        }
        gym::unmapBufferC4V4();
        gym::drawLines(LW, 0, n);
        gym::cleanup();
    }
    
    if ( showReproject )
    {
        flute8* flu = gym::mapBufferC4V4(2*nbpts);
        gym_color col(COL, 0.f, 0.f), lor(COL, 0.f, 0.f, 0.5f);
        size_t n = 0;
        for ( size_t i = 0; i < nbpts; ++i )
        {
            if ( visible(i) )
            {
                flu[n++] = { col, project[i] };
                flu[n++] = { lor, project2[i] };
            }
        }
        gym::unmapBufferC4V4();
        gym::drawLines(2*LW, 0, n);
        gym::cleanup();
    }
    
    if ( showEdges )
    {
        flute8* flu = gym::mapBufferC4V4(2*nbpts);
        gym_color col(0.f, COL, COL), lor(0.f, COL, 0.f);
        size_t n = 0;
        for ( size_t i = 0; i < nbpts; ++i )
        {
            if ( visible(i) )
            {
                flu[n++] = { col, edge[i] };
                flu[n++] = { lor, project[i] };
            }
         }
        gym::unmapBufferC4V4();
        gym::drawPoints(2, 0, n);
        gym::cleanup();
    }
    view.closeDisplay();
    return 0;
}

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
    glutInit(&argc, argv);
    glApp::setDimensionality(DIM);
    glApp::normalKeyFunc(processNormalKey);
    glApp::specialKeyFunc(processSpecialKey);
    glApp::newWindow(display);
    glApp::setScale(20);
    gle::initialize();

    initMenus();
    RNG.seed();

    if ( argc > 1 )
    {
        if ( opt.read_strings(argc-1, argv+1) )
            return EXIT_FAILURE;
        setGeometry();
    }
    
    if ( ! spc )
    {
        printf("A geometry should be given in the command line, for example:\n");
        printf("    test_space shape=ellipse length=2,3,4\n");
        exit(EXIT_SUCCESS);
    }

    glutMainLoop();
}

