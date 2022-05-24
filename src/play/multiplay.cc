/*
 * Multiplay runs simulations in child threads in parallel and display them live
 * This is the second GLFW-based Cytosim player
 * Francois J Nedelec, Cambridge University, 1--23 May 2022
 */

#include <stdio.h>
#include <stdlib.h>

#include "simul.h"
#include "glossary.h"
#include "messages.h"
#include "filepath.h"
#include "sim_thread.h"

#include "view.h"
#include "display_prop.h"
#include "display1.h"
#include "gym_draw.h"
#include "gym_view.h"
#include "gym_flat.h"
#include "gle.h"
#include "glfw.h"

constexpr int TILEX = 5;
constexpr int TILEY = 6;
constexpr int TOP = TILEX * TILEY;
constexpr int RENEW_DELAY = TOP;
constexpr int PERIOD = 16;

// size of each cell (updated if window is resized)
int bugW = 320;
int bugH = 128;

// size of window (updated if window is resized)
int winW = bugW * TILEX;
int winH = bugH * TILEY;

//------------------------------------------------------------------------------

View view("multiplay", DIM==3);
DisplayProp disp("multiplay");
Display1 display(&disp);
PropertyList dispList;

int selected = 0;
int fate[TOP] = { 0 };

std::string config_file = "config.cym";

// pointers to hold the Config text:
char *config_code = nullptr;
size_t code_size = 0;

//------------------------------------------------------------------------------
#pragma mark -

/* returns cell index corresponding to mouse position (mx, my) */
int whichBug(double mx, double my)
{
    mx = std::max(mx, 0.0);
    my = std::max(my, 0.0);
    int x = std::min(TILEX, int(TILEX * mx));
    int y = std::min(TILEY, int(TILEY * my));
    return y + TILEY * x;
}

/* respond to mouse cursor movements */
void mouseMotionCallback(GLFWwindow* win, double mx, double my)
{
    int state = glfwGetMouseButton(win, GLFW_MOUSE_BUTTON_LEFT);
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

/* enter/exit full screen mode */
void toggleFullScreen(GLFWwindow* win)
{
    static int x, y, w, h;
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
            glfwGetMonitorWorkarea(moni, &mx, &my, &mw, &mh);
            glfwSetWindowMonitor(win, moni, mx, my, mw, mh, GLFW_DONT_CARE);
        }
    }
}

/* respond to keyboard events based on keyboard layout */
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
        case '-': view.zoom_in(0.91700404320); break;
        case '=': view.zoom_in(1.09050773266); break;
        case 'F': toggleFullScreen(win); break;
#if ( 0 )
        case GLFW_KEY_UP:
            view.eyeDistance += 0.1; printf("eye %.2f\n", view.eyeDistance); break;
        case GLFW_KEY_DOWN:
            view.eyeDistance -= 0.1; printf("eye %.2f\n", view.eyeDistance); break;
#endif
#if ( DIM > 1 )
        case GLFW_KEY_LEFT:
            view.rotate_by(Quaternion<real>(0.99,0,-.1,0)); break;
        case GLFW_KEY_RIGHT:
            view.rotate_by(Quaternion<real>(0.99,0,0.1,0)); break;
        default:
#endif
            return;
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
    glfwGetFramebufferSize(win, &bugW, &bugH);
    bugW /= TILEX;
    bugH /= TILEY;
    view.reshape(bugW, bugH);
    //printf("Window %ix%i : tile %ix%i\n", winW, winH, bugW, bugH);
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

    // canvas size in pixels is not necessarily the window size (in screen coordinates)
    glfwGetFramebufferSize(win, &W, &H);
    gle::initialize();
    reshape(win, W, H);
    view.initGL();
    return win;
}

//------------------------------------------------------------------------------
#pragma mark -

/* get ready to draw System */
void prepareDraw()
{
    view.setLights();
    view.adjust(bugW, bugH);
    view.loadView();
    display.setPixelFactors(view.pixelSize(), 1);
    //printf("dispList: %lu\n", dispList.size());
}


/* prepare for drawing bug at (x,y) */
inline void selectPanel(int x, int y)
{
    glScissor(bugW*x, bugH*y, bugW, bugH);
    view.setViewport(bugW*x, bugH*y, bugW, bugH);
}


/* draw System */
void drawBug(Simul const& sim)
{
    display.prepareForDisplay(sim, dispList, view.depthAxis());
    gym::clearPixels(view.back_color);
    display.drawSimul(sim);
    //gym::enableLighting(); gym::scale(0.2); gle::star();
}

//------------------------------------------------------------------------------

/* program entry */
int main(int argc, char *argv[])
{
    Simul simul[TOP];
    SimThread worker[TOP];

    //parse the command line:
    Glossary arg;
    if ( arg.read_strings(argc-1, argv+1) )
        return 1;
    if ( !arg.empty() )
    {
        view.read(arg);
        disp.read(arg);
        arg.set(config_file, ".cym");
        arg.print_warnings(std::cerr, 1, "\n");
    }
    
    RNG.seed();
    if ( !FilePath::read_file(config_file.c_str(), config_code, code_size) )
        exit(1);
    GLFWwindow* win = initWindow(bugW*TILEX, bugH*TILEY);
    Cytosim::silent();

    for ( int i = 0; i < TOP; ++i )
    {
        //usleep(100000);
        simul[i].initialize(arg);
        simul[i].prop.config_file = config_file;
        simul[i].prop.random_seed = RNG.pint32();
        worker[i].config_code = config_code;
        worker[i].set_simul(simul+i);
        worker[i].period(PERIOD);
        worker[i].start();
    }

    prepareDraw();
    glEnable(GL_SCISSOR_TEST);
    while( !glfwWindowShouldClose(win) )
    {
        //usleep(100000);
        int refresh = 0;
        for (int x = 0; x < TILEX; ++x )
        for (int y = 0; y < TILEY; ++y )
        {
            int i = y + TILEY * x;
            if ( worker[i].holding() && 0 == worker[i].trylock() )
            {
                ++refresh;
                selectPanel(x, y);
                drawBug(simul[i]);
                if ( worker[i].holding() > 1 )
                    worker[i].restart();
                worker[i].unlock();
                worker[i].signal();
                if ( fate[i] > 0 )
                    worker[i].cancel();
            }
            else if ( worker[i].dead() )
            {
                //fprintf(stderr, "%2i: dead %i\n", i, worker[i].dead() );
                if ( fate[i] > 0 )
                {
                    ++refresh;
                    selectPanel(x, y);
                    gym::clearPixels(0.5f, 0.5f, 0.5f, 1.0f);
                    fate[i] = -1-RNG.pint32(RENEW_DELAY);
                }
                else if ( fate[i] < 0 )
                {
                    if ( ++fate[i] == 0 )
                    {
                        //fprintf(stderr, "%2i: start\n", i);
                        worker[i].erase_simul(1);
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
    for ( int i = 0; i < TOP; ++i )
    {
        worker[i].cancel();
        worker[i].config_code = nullptr;
    }
    for ( int i = 0; i < TOP; ++i )
        worker[i].join();
    glfwDestroyWindow(win);
    glfwTerminate();
    free(config_code);
}

