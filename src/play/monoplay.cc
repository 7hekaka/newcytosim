/*
 * Multiplay runs one simulation in a child thread and displays it live
 * Francois J Nedelec, Cambridge University, 6 Sep. 2022
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#include "simul.h"
#include "glossary.h"
#include "sim_thread.h"

#include "view.h"
#include "display_prop.h"
#include "display1.h"
#include "gym_draw.h"
#include "gym_view.h"
#include "gym_flat.h"
#include "gle.h"
#include "glfw.h"

// size of window (updated if window is resized)
int winW = 1024;
int winH = 1024;

// size of framebuffer (updated if window is resized)
int bugW = 1024;
int bugH = 1024;

//------------------------------------------------------------------------------

View view("monoplay", DIM==3);
DisplayProp disp("monoplay");
Display1 display(&disp);
PropertyList dispList;

//------------------------------------------------------------------------------
#pragma mark -

/* respond to mouse cursor movements */
void mouseMotionCallback(GLFWwindow* win, double mx, double my)
{
    int state = glfwGetMouseButton(win, GLFW_MOUSE_BUTTON_RIGHT);
    //printf("mouse @ %8.2f %8.2f (%i)\n", mx, my, state);
    if ( state == GLFW_PRESS )
    {
    }
}

/* respond to mouse clicks */
void mouseButtoncallback(GLFWwindow* win, int button, int action, int mods)
{
    double mx, my;
    glfwGetCursorPos(win, &mx, &my);
    //printf("click @ %8.2f %8.2f (%i %i)\n", mx, my, button, action);
    if ( button == GLFW_MOUSE_BUTTON_LEFT )
    {
    }
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
void keysCallback(GLFWwindow* win, int k, int s, int action, int mods)
{
    if ( action != GLFW_PRESS )
        return;
    //printf("keyCallback %c (%i %i)\n", k, action, mods);
    switch (k)
    {
        case GLFW_KEY_ESCAPE:
            if ( glfwGetWindowMonitor(win) )
                toggleFullScreen(win);
            else
                glfwSetWindowShouldClose(win, GLFW_TRUE);
            break;
        case 'F': toggleFullScreen(win); break;
        case '-': view.zoom_in(0.91700404320); break;
        case '=': view.zoom_in(1.09050773266); break;
        case 'X': view.axes = ( view.axes ? 0 : 3 ); break;
#if ( DIM > 1 )
        case GLFW_KEY_UP:    view.rotate_by(Quaternion<real>(0.99,-.1,0,0)); break;
        case GLFW_KEY_DOWN:  view.rotate_by(Quaternion<real>(0.99,0.1,0,0)); break;
        case GLFW_KEY_LEFT:  view.rotate_by(Quaternion<real>(0.99,0,-.1,0)); break;
        case GLFW_KEY_RIGHT: view.rotate_by(Quaternion<real>(0.99,0,0.1,0)); break;
#endif
        default: printf("keysCallback %c (%i)\n", k, mods); return;
    }
}

/* respond to keyboard event based on characters emitted */
void charCallback(GLFWwindow* win, unsigned int k, int mods)
{
    printf("charCallback %c (%i)\n", k, mods);
    switch (k)
    {
        case '-': view.zoom_in(M_SQRT1_2); break;
        case '=': view.zoom_in(M_SQRT2); break;
    }
}

/* change window size, adjust display to maintain isometric axes */
void reshape(GLFWwindow* win, int W, int H)
{
    glfwGetWindowSize(win, &winW, &winH);
    bugW = W;
    bugH = H;
    view.reshape(W, H);
    printf("reshape window %ix%i : framebuffer %ix%i\n", winW, winH, W, H);
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

    GLFWwindow* win = glfwCreateWindow(W, H, "Cytosim", nullptr, nullptr);
    if (!win)
    {
        fprintf(stderr, "Failed to open GLFW window\n");
        glfwTerminate();
        exit(EXIT_FAILURE);
    }
    
    // Set callback functions
    glfwSetKeyCallback(win, keysCallback);
    glfwSetFramebufferSizeCallback(win, reshape);
    glfwSetCursorPosCallback(win, mouseMotionCallback);
    glfwSetMouseButtonCallback(win, mouseButtoncallback);
    //glfwSetCharModsCallback(win, charCallback);
    glfwMakeContextCurrent(win);
    //gladLoadGL(glfwGetProcAddress);
    glfwSwapInterval(1);

    unsigned char pixels[16 * 16 * 4];
    memset(pixels, 0xaf, sizeof(pixels));
    GLFWimage image { 16, 16, pixels };
    GLFWcursor * cursor = glfwCreateCursor(&image, 0, 0);
    glfwSetCursor(win, cursor);

    gle::initialize();
    // canvas size in pixels is not necessarily the window size (in screen coordinates)
    glfwGetFramebufferSize(win, &W, &H);
    reshape(win, W, H);
    view.initGL();
    return win;
}


/* draw System */
void drawBug(Simul const& sim)
{
    view.openDisplay();
    display.setPixelFactors(view.pixelSize(), 1);
    display.prepareForDisplay(sim, dispList, view.depthAxis());
    display.drawSimul(sim);
    view.closeDisplay();
}

//------------------------------------------------------------------------------

/* program entry */
int main(int argc, char *argv[])
{
    Simul simul;
    SimThread worker;

    //parse the command line:
    Glossary arg;
    if ( arg.read_strings(argc-1, argv+1) )
        return 1;
    if ( !arg.empty() )
    {
        view.read(arg);
        disp.read(arg);
        unsigned P = 1;
        arg.set(P, "period");
        worker.period(P);
        simul.initialize(arg);
        arg.print_warnings(std::cerr, 1, "\n");
    }
    
    GLFWwindow* win = initWindow(winW, winH);
    worker.set_simul(&simul);
    worker.start();

    while( !glfwWindowShouldClose(win) )
    {
        //usleep(100000);
        int refresh = 0;
        if ( worker.holding() && 0 == worker.trylock() )
        {
            ++refresh;
            drawBug(simul);
            if ( worker.holding() > 1 )
                worker.restart();
            worker.unlock();
            worker.signal();
        }
        else if ( worker.dead() )
        {
            //fprintf(stderr, "%2i: dead %i\n", i, worker[i].dead() );
            worker.erase_simul(1);
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
}

