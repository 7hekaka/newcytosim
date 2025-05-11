/*
 * Monoplay runs one simulation in a child thread and displays it live
 * Francois J Nedelec, Cambridge University, 6 Sep. 2022
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>

#include "simul.h"
#include "glossary.h"
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

// size of framebuffer (updated if window is resized)
int bugW = 1024;
int bugH = 1024;

// size of window (updated if window is resized)
int winW = 1024;
int winH = 1024;

//------------------------------------------------------------------------------

View view("monoplay", DIM==3);
DisplayProp disp("monoplay");
Display3 display(&disp);
PropertyList mFiberDisp, mPointDisp;

//------------------------------------------------------------------------------
#pragma mark -

/* respond to mouse cursor movements */
void mouseMotionCallback(GLFWwindow* win, double mx, double my)
{
    int state = glfwGetMouseButton(win, GLFW_MOUSE_BUTTON_RIGHT);
    //Vector3 pos = view.unproject(mx*bugW/winW, bugH-my*bugH/winH, 0.5);
    //printf("mouse @ %8.2f %8.2f (%i) %8.2f %8.2f \n", mx, my, state, pos.XX, pos.YY);
    if ( state == GLFW_PRESS )
    {
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
        //printf("exit fullscreen\n");
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
            //printf("monitor %i:  %ix%i\n", u, mw, mh);
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
    //printf("keyCallback %i (%c) %i %i\n", k, k, action, mods);
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

/* respond to keyboard events based on character emitted */
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
    bugW = W;
    bugH = H;
    view.reshape(W, H);
    //printf("GLFW window %i:%i buffer %i:%i\n", winW, winH, W, H);
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
    glfwWindowHint(GLFW_DOUBLEBUFFER, GLFW_TRUE);
    glfwWindowHint(GLFW_CLIENT_API, GLFW_OPENGL_API);
    glfwWindowHint(GLFW_CONTEXT_CREATION_API, GLFW_NATIVE_CONTEXT_API);
    glfwWindowHint(GLFW_COCOA_RETINA_FRAMEBUFFER, GLFW_TRUE);
    
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

    GLFWcursor * cursor = glfwCreateStandardCursor(GLFW_CROSSHAIR_CURSOR);
    glfwSetCursor(win, cursor);

    gle::initialize();
    // canvas size in pixels is not necessarily the window size (in screen coordinates)
    glfwGetFramebufferSize(win, &W, &H);
    reshapeCallback(win, W, H);
    view.initGL();
    return win;
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


/* draw System */
void drawCytosim(Simul const& sim)
{
    prepareDisplay(sim);
    view.openDisplay();
    display.drawSimul(sim);
    view.closeDisplay();
}

//------------------------------------------------------------------------------

/* program entry */
int main(int argc, char *argv[])
{
    Simul simul;
    SimThread worker(&simul);
    
    simul.initCytosim();
    seed_pcg32();

    //parse the command line:
    Glossary arg;
    if ( arg.read_strings(argc-1, argv+1) )
        return 1;
    if ( !arg.empty() )
    {
        simul.prop.read(arg);
        view.read(arg);
        disp.read(arg);
        unsigned P = 1;
        arg.set(P, "period");
        worker.period(P);
        arg.set(winW, "size") && ( arg.set(winH, "size", 1) || ( winH = winW ));
        arg.print_warnings(stderr, 1, "\n");
    }

    GLFWwindow* win = initWindow(winW, winH);
    worker.start();

    while ( !glfwWindowShouldClose(win) )
    {
        //usleep(100000);
        int refresh = 0;
        if ( worker.holding() && 0 == worker.trylock() )
        {
            // data is now locked for this process
            ++refresh;
            drawCytosim(simul);
            if ( worker.holding() > 1 )
                worker.restart();
            worker.unlock(); // release lock
            worker.signal(); // permit other thread to resume
        }
        else if ( worker.dead() )
        {
            usleep(100000);
            //fprintf(stderr, "%2i: dead %i\n", i, worker[i].dead() );
            worker.eraseSimul(1);
            worker.start();
        }
        if ( refresh )
        {
            glFlush();
            glfwSwapBuffers(win);
        }
        glfwPollEvents();
    }
    worker.cancel();
    worker.join();
    glfwDestroyWindow(win);
    glfwTerminate();
    mFiberDisp.erase();
    mPointDisp.erase();
}

