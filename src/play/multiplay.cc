/*
 * Multiplay runs simulations in child threads in parallel and display them live
 * This is the second GLFW-based Cytosim player
 * Francois J Nedelec, Cambridge University, 1--24 May 2022
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "simul.h"
#include "glossary.h"
#include "messages.h"
#include "filepath.h"
#include "sim_thread.h"

#include "view.h"
#include "display_prop.h"
#include "display3.h"
#include "gym_draw.h"
#include "gym_view.h"
#include "gym_flat.h"
#include "glfw.h"
#include "gle.h"

#include "random_pcg.h"
#include "random_seed.h"
using namespace PCG32;

int tileX = 2;
int tileY = 3;
constexpr int TOP = 128;

inline int nbBugs() { return std::min(tileX*tileY, TOP); }

// size of each cell (updated if window is resized)
int bugW = 320;
int bugH = 128;

// size of window (updated if window is resized)
int winW = bugW * tileX;
int winH = bugH * tileY;

//------------------------------------------------------------------------------

unsigned period = 16;
View view("multiplay", DIM==3);
DisplayProp disp("multiplay");
Display3 display(&disp);
PropertyList mFiberDisp, mPointDisp;

int selected = 0;
int fate[TOP] = { 0 };

// pointers to hold the Config text, same for all simuls:
char *config_text = nullptr;
size_t code_size = 0;

//------------------------------------------------------------------------------
#pragma mark -

/* returns cell index corresponding to mouse position (mx, my) */
int whichBug(double mx, double my)
{
    mx = std::max(mx, 0.0);
    my = std::max(my, 0.0);
    int x = std::min(tileX-1, int(tileX * mx));
    int y = std::min(tileY-1, int(tileY * my));
    return y + tileY * x;
}

/* respond to mouse cursor movements */
void mouseMotionCallback(GLFWwindow* win, double mx, double my)
{
    int state = glfwGetMouseButton(win, GLFW_MOUSE_BUTTON_RIGHT);
    //printf("mouse @ %8.2f %8.2f (%i)\n", mx, my, state);
    if ( state == GLFW_PRESS )
    {
        int i = whichBug(mx/winW, 1.0-my/winH);
        fate[i] = 1;
    }
}

/* respond to mouse clicks */
void mouseButtonCallback(GLFWwindow* win, int button, int action, int mods)
{
    double mx, my;
    glfwGetCursorPos(win, &mx, &my);
    //printf("click @ %8.2f %8.2f (%i %i)\n", mx, my, button, action);
    if ( button == GLFW_MOUSE_BUTTON_LEFT )
    {
        int i = whichBug(mx/winW, 1.0-my/winH);
        if ( action == GLFW_PRESS )
            selected = i;
        else if ( action == GLFW_RELEASE )
        {
            if ( i == selected )
                fate[i] = 1;
        }
    }
}

/* respond to scroll events (mouse wheel) */
void scrollCallback(GLFWwindow* win, double dx, double dy)
{
    double mx, my;
    glfwGetCursorPos(win, &mx, &my);
    double Z = std::max(0.5, 1.0 - 0.0625 * dy);
    Vector3 pos = view.unproject(mx*bugW/winW, bugH-my*bugH/winH, 0.5);
    //printf("scroll @ %8.2f %8.2f (%8.2f %8.2f) %8.2f %8.2f\n", mx, my, dx, dy, pos.XX, pos.YY);
    view.zoom_out(Z);
    view.move_to((1-Z)*pos+Z*(view.focus+view.focus_shift));
}

/* enter/exit full screen mode */
void toggleFullScreen(GLFWwindow* win)
{
    static int x = 0, y = 0, w = winW, h = winH;
    GLFWmonitor* moni = glfwGetWindowMonitor(win);
    if ( moni )
    {
        glfwSetWindowMonitor(win, nullptr, x, y, w, h, GLFW_DONT_CARE);
    }
    else
    {
        int mx, my, mw, mh;
        glfwGetWindowPos(win, &x, &y);
        glfwGetWindowSize(win, &w, &h);
        int cnt = 0;
        size_t area = 0;
        GLFWmonitor** list = glfwGetMonitors(&cnt);
        // select largest monitor:
        for ( int u = 0; u < cnt; ++u )
        {
            glfwGetMonitorWorkarea(list[u], &mx, &my, &mw, &mh);
            size_t a = mw * mh;
            if ( a > area )
            {
                area = a;
                moni = list[u];
            }
        }
        if ( moni )
        {
            const GLFWvidmode* mode = glfwGetVideoMode(moni);
            glfwSetWindowMonitor(win, moni, 0, 0, mode->width, mode->height, mode->refreshRate);
        }
    }
}

/* respond to keyboard events based on keyboard layout: 'k' is capitalized */
void keyCallback(GLFWwindow* win, int k, int s, int action, int mods)
{
    //printf("keyCallback %c (%i %i)\n", k, action, mods);
    if ( action == GLFW_PRESS ) switch (k)
    {
        case GLFW_KEY_ESCAPE:
            if ( glfwGetWindowMonitor(win) )
                toggleFullScreen(win);
            else
                glfwSetWindowShouldClose(win, GLFW_TRUE);
            break;
#if ( DIM > 2 )
        case GLFW_KEY_UP:    view.rotate_by(Quaternion<real>(0.99,-.1,0,0)); break;
        case GLFW_KEY_DOWN:  view.rotate_by(Quaternion<real>(0.99,0.1,0,0)); break;
        case GLFW_KEY_LEFT:  view.rotate_by(Quaternion<real>(0.99,0,-.1,0)); break;
        case GLFW_KEY_RIGHT: view.rotate_by(Quaternion<real>(0.99,0,0.1,0)); break;
#endif
    }
}

/* respond to keyboard event based on character emitted */
void charCallback(GLFWwindow* win, unsigned int k, int mods)
{
    //printf("charCallback %i (%c)\n", k, k);
    switch (k)
    {
        case 'f': toggleFullScreen(win); break;
        case '-': view.zoom_in(0.91700404320); break;
        case '=': view.zoom_in(1.09050773266); break;
        case 'x': view.axes = ( view.axes ? 0 : 3 ); break;
        case 'z': view.reset(); break;
    }
}

/* change window size, adjust display to maintain isometric axes */
void reshapeCallback(GLFWwindow* win, int W, int H)
{
    glfwGetWindowSize(win, &winW, &winH);
    bugW = W / tileX;
    bugH = H / tileY;
    view.reshape(bugW, bugH);
    //printf("reshape window %ix%i : tile %ix%i\n", winW, winH, bugW, bugH);
}


void glfwError(int error, const char* text)
{
    fprintf(stderr, "GLFW Error: %s\n", text);
}

/* program & OpenGL initialization */
GLFWwindow * initWindow(int W, int H)
{
    if ( !glfwInit() )
    {
        fprintf(stderr, "Failed to initialize GLFW\n");
        exit(EXIT_FAILURE);
    }
    glfwSetErrorCallback(glfwError);

    glfwDefaultWindowHints();
    glfwWindowHint(GLFW_DEPTH_BITS, 24);
    glfwWindowHint(GLFW_DOUBLEBUFFER, GLFW_FALSE);
    glfwWindowHint(GLFW_CLIENT_API, GLFW_OPENGL_API);
    glfwWindowHint(GLFW_CONTEXT_CREATION_API, GLFW_NATIVE_CONTEXT_API);

    GLFWwindow* win = glfwCreateWindow(W, H, "Cytosim", nullptr, nullptr);
    if (!win)
    {
        fprintf(stderr, "Failed to open GLFW window\n");
        glfwTerminate();
        exit(EXIT_FAILURE);
    }
    
    // Set callback functions
    glfwSetKeyCallback(win, keyCallback);
    glfwSetCharModsCallback(win, charCallback);
    glfwSetFramebufferSizeCallback(win, reshapeCallback);
    glfwSetCursorPosCallback(win, mouseMotionCallback);
    glfwSetMouseButtonCallback(win, mouseButtonCallback);
    glfwSetScrollCallback(win, scrollCallback);
    glfwMakeContextCurrent(win);
    //gladLoadGL(glfwGetProcAddress);
    glfwSwapInterval(1);

    unsigned char pixels[16 * 16 * 4];
    memset(pixels, 0xaf, sizeof(pixels));
    GLFWimage image { 16, 16, pixels };
    GLFWcursor * cursor = glfwCreateCursor(&image, 0, 0);
    glfwSetCursor(win, cursor);

    // canvas size in pixels is not necessarily the window size (in screen coordinates)
    glfwGetFramebufferSize(win, &W, &H);
    gle::initialize();
    reshapeCallback(win, W, H);
    view.initGL();
    return win;
}

//------------------------------------------------------------------------------
#pragma mark -

/* get ready to draw System */
void prepareView()
{
    view.setLights();
    view.adjust(bugW, bugH);
    view.loadView();
    //printf("allDisp: %lu\n", allDisp.size());
}


//set pixel size, unit-size and direction:
void prepareDisplay(Simul const& sim)
{
    if ( view.auto_scale > 0 )
        view.set_scale(2*sim.spaces.maxExtension());

    if ( sim.prop.display_fresh )
    {
        Glossary glos(sim.prop.display);
        disp.read(glos);
        view.read(glos);
        sim.prop.display_fresh = false;
        view.setModelView();
    }
    float mag = view.magnify;
    float pix = view.pixelSize() / mag;
    /*
     if `disp.point_value` is set, widths of lines and point sizes are understood to
     be specified in 'real' units, while by default, they are understood in pixels.
     */
    if ( disp.point_value > 0 )
        mag = disp.point_value / pix;
    
    display.setParameters(pix, mag, view.depthAxis());
    display.setStencil(view.stencil);
    
    for ( Property * p : mPointDisp )
        static_cast<PointDisp*>(p)->setPixels(pix, mag, disp.style==2);
    display.prepareDrawing(sim, mFiberDisp, mPointDisp);
}


/* prepare for drawing bug at (x,y) */
inline void selectPanel(int x, int y)
{
    glScissor(bugW*x, bugH*y, bugW, bugH);
    view.setViewport(bugW*x, bugH*y, bugW, bugH);
}


/* draw System */
void drawCytosim(Simul const& sim)
{
    prepareDisplay(sim);
    gym::clearPixels(view.back_color);
    display.drawSimul(sim);
    //gym::enableLighting(); gym::scale(0.2); gle::star();
}

//------------------------------------------------------------------------------

/// read 'X:Y' or 'INTxINT' where the separating character (':') is returned
int readDimensions(int& X, int& Y, const char str[])
{
    char sep = 0;
    if ( !isdigit(*str) )
        return 1;
    int x = X, y = Y;
    char* ptr = nullptr;
    errno = 0;
    x = (int)strtol(str, &ptr, 10);
    if ( x <= 0 )
    {
        if ( errno ) perror(str);
        return 2;
    }
    sep = *ptr++;
    if ( !isprint(sep) )
        return 3;
    y = (int)strtol(ptr, &ptr, 10);
    if ( y <= 0 )
    {
        if ( errno ) perror(str);
        return 4;
    }
    if ( *ptr && !isspace(*ptr) )
        return 5;
    X = x;
    Y = y;
    return sep;
}


/* program entry */
int main(int argc, char *argv[])
{
    Simul simul[TOP];
    SimThread worker[TOP];

    //parse the command line:
    if ( argc > 1 && readDimensions(tileX, tileY, argv[1]) )
        *argv[1] = 0;
    if ( tileX * tileY > TOP )
        printf("Error: the number of cells is limited to %i\n", TOP);
    Glossary arg;
    if ( arg.read_strings(argc-1, argv+1) )
        return 1;
    
    int cnt = std::min(nbBugs(), (int)arg.num_values(".cym"));
    for ( int i = 0; i < cnt; ++i )
        arg.set(simul[i].prop.config_file, ".cym", i);
    if ( cnt < 2 )
    {
        const char * str = simul[0].prop.config_file.c_str();
        if ( !FilePath::read_file(str, config_text, code_size) )
            return EXIT_FAILURE;
        for ( int i = 0; i < nbBugs(); ++i )
            worker[i].config_code = config_text;
    }
    arg.clear(".cym");
    if ( !arg.empty() )
    {
        view.read(arg);
        disp.read(arg);
        arg.set(period, "period");
        arg.print_warnings(stderr, 1, "\n");
    }
    
    GLFWwindow* win = initWindow(bugW*tileX, bugH*tileY);
    Cytosim::silent();

    // seed all RNGs used by sim using a common random generator:
    seed_pcg32();

    for ( int i = 0; i < nbBugs(); ++i )
    {
        simul[i].prop.read(arg);
        simul[i].prop.random_seed = pcg32();
        worker[i].simul(simul+i);
        worker[i].period(period);
        worker[i].start();
    }

    prepareView();
    glEnable(GL_SCISSOR_TEST);
    while ( !glfwWindowShouldClose(win) )
    {
        //usleep(100000);
        int refresh = 0;
        for (int x = 0; x < tileX; ++x )
        for (int y = 0; y < tileY; ++y )
        {
            int i = y + tileY * x;
            if ( worker[i].holding() && 0 == worker[i].trylock() )
            {
                ++refresh;
                selectPanel(x, y);
                drawCytosim(simul[i]);
                if ( worker[i].holding() > 1 )
                    worker[i].restart();
                worker[i].unlock();
                if ( fate[i] > 0 )
                    worker[i].cancel();
                else
                    worker[i].signal();
            }
            else if ( worker[i].dead() )
            {
                //fprintf(stderr, "%2i: dead %i\n", i, worker[i].dead() );
                if ( fate[i] > 0 )
                {
                    ++refresh;
                    selectPanel(x, y);
                    gym::clearPixels(0.5f, 0.5f, 0.5f, 1.0f);
                    fate[i] = -1 - RNG.pint32(5*tileX*tileY);
                }
                else if ( fate[i] < 0 )
                {
                    if ( ++fate[i] == 0 )
                    {
                        //fprintf(stderr, "%2i: start\n", i);
                        worker[i].eraseSimul(1);
                        worker[i].start();
                    }
                }
            }
        }
        if ( refresh )
        {
            glFlush();
            //glfwSwapBuffers(win);
        }
        glfwPollEvents();
    }
    //stopPreconfig();
    glDisable(GL_SCISSOR_TEST);
    for ( int i = 0; i < nbBugs(); ++i )
        worker[i].cancel();
    for ( int i = 0; i < nbBugs(); ++i )
        worker[i].join();
    glfwDestroyWindow(win);
    glfwTerminate();
    free(config_text);
    mFiberDisp.erase();
    mPointDisp.erase();
}

