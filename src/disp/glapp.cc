// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University
// started glApp in Dec 2005

#include "glut.h"
#include "glapp.h"
#include "exceptions.h"
#include "glossary.h"
#include "time_date.h"
#include "gle.h"
#include "gym_check.h"
#include "gym_menu.h"
#include "gym_draw.h"
#include "offscreen.h"
#include "save_image_gl.h"

/** A layer on top of GLUT adding some user interface */
namespace glApp
{
    std::vector<View> views;
    
    int mDIM        = 3;  ///< current dimensionality mode
    int mFullScreen = 0;  ///< flag indicating full-screen mode
    int specialKeys = 0;  ///< state of special keys given by GLUT
    int drawRefresh = 0;  ///< true if display needs to be updated
    
    // --------------- MOUSE
    
    /// actions that can be performed with the mouse
    enum UserMode
    {
        MOUSE_TURN      = 0,
        MOUSE_MOVE      = 1,
        MOUSE_ACTIVE    = 2,
        MOUSE_PUSH      = 3,
        MOUSE_SPIN      = 4,
        MOUSE_MAGNIFIER = 5,
        MOUSE_EDIT_ROI  = 6,
        MOUSE_SELECT    = 7,
        MOUSE_PASSIVE   = 8
    };
    
    /// Specifies in which dimensionality each action is valid
    int actionDimensionality[] = { 3, 1, 1, 3, 2, 1, 1, 4, 4, 4, 0 };
    
    /// the action controlled with the mouse
    UserMode userMode = MOUSE_TURN;

    /// change action
    void switchUserMode(int dir);
    
    View savedView("savedView", 0);

    UserMode mouseAction = MOUSE_MOVE;  ///< the action being performed by the mouse
    int mouseS = 1;     ///< current state of mouse button (up/down)
    int mouseX, mouseY; ///< current position of mouse in pixels
    Vector3  mouseDown(0,0,0);  ///< position where mouse button was pressed down
    Vector3  axle(0,0,0);       ///< vector normal to desired rotation
    Vector3  pole(0,0,0);       ///< axis of rotation for MOUSE_SPIN
    Vector3  viewFocus(0,0,0);  ///< unprojected center of screen
    
    float nearZ = 0;    ///< normalized device Z-coordinate of front-plane
    float midZ  = 0.5;  ///< normalized device Z-coordinate of middle
    float farZ  = 1;    ///< normalized device Z-coordinate of back-plane

    unsigned imageIndex = 0;  ///< index for image name

    Vector3  ROIdup[2];
    bool     moveROI = false;

    /// function pointer for shift-click actions
    void (*mouseClickCallback)(int, int, const Vector3 &, int) = nullptr;
    
    /// function pointer for shift-motion actions
    void (*mouseDragCallback)(int, int, Vector3 &, const Vector3 &, int) = nullptr;
    
    /// function pointer for shift-click actions
    void (*normalKeyCallback)(unsigned char, int, int) = processNormalKey;
    
    /// function pointer for shift-motion actions
    void (*specialKeyCallback)(int, int, int) = processSpecialKey;
    
}


//------------------------------------------------------------------------------
#pragma mark -

/**
 initialize a single View: views[0]
 */
void glApp::initialize()
{
    views.clear();
    views.push_back(View("*", mDIM==3));
}


/**
 This will disable OpenGL depth-test for DIM<3
 */
void glApp::setDimensionality(const int d)
{
    if ( mDIM != d )
    {
        //flashText("dimensionality changed to %i", d);
        mDIM = d;
        userMode = ( d == 3 ) ? MOUSE_TURN : MOUSE_MOVE;
    }
    
    if ( 0 == views.size() )
        initialize();
}

//------------------------------------------------------------------------------

void glApp::toggleFullScreen()
{
    static int  savedPos[4] = { 24, 24, 800, 800 };
    if ( mFullScreen )
    {
        mFullScreen = 0;
        // restore window dimensions:
        if ( savedPos[2] > 8 && savedPos[3] > 8 )
        {
            //std::clog << "saveWindow " << savedPos[2] << " " << savedPos[3] << '\n';
            glutReshapeWindow(savedPos[2], savedPos[3]);
            glutPositionWindow(savedPos[0], savedPos[1]);
        }
        else
            glutReshapeWindow(800, 800);
    }
    else
    {
        mFullScreen = 1;
        savedPos[0] = glutGet(GLUT_WINDOW_X);
        savedPos[1] = glutGet(GLUT_WINDOW_Y);
        savedPos[2] = glutGet(GLUT_WINDOW_WIDTH);
        savedPos[3] = glutGet(GLUT_WINDOW_HEIGHT);
        //invoke full screen from GLUT
        glutFullScreen();
        //std::clog << "Fullscreen window " << glutGetWindow() << '\n';
    }
}

#if !defined(GLUT_WINDOW_SCALE)
#    define GLUT_WINDOW_SCALE 199
#endif

/**
 Adjust the size of window to maximize its dimension within the screen,
 without changing the aspect ratio of the window.
 */
void glApp::toggleWindowSize()
{
#if defined(__APPLE__)    /* For OSX */
    int Y = 53; // size of menubar on macs
#else
    int Y = 0;
#endif
    int W = glutGet(GLUT_WINDOW_WIDTH);
    int H = glutGet(GLUT_WINDOW_HEIGHT);
    int X = glutGet(GLUT_WINDOW_X) + W/2;

    /// using addition by Renaud Blanch to handle Retina display:
    float S = 1; //std::max(1, glutGet(GLUT_WINDOW_SCALE));

    int maxW = glutGet(GLUT_SCREEN_WIDTH);
    int maxH = glutGet(GLUT_SCREEN_HEIGHT) - Y;

    if ( H >= maxH || W >= maxW )
        S = 0.5;
    
    float zW(S*maxW/(float)W);
    float zH(S*maxH/(float)H);
    
    if ( zW < zH )
    {
        W = S * maxW;
        H = H * zW;
    }
    else
    {
        W = W * zH;
        H = S * maxH;
    }
    // keep window within the edges:
    X = std::min(X - W/2, maxW-W);
    X = std::max(X, 0);
    
    glutReshapeWindow(W, H);
    glutPositionWindow(X, Y);
}

//------------------------------------------------------------------------------
#pragma mark -

/**
 @todo: share OpenGL context between new window and main one,
 which is possible with GLFW
 */
int glApp::newWindow(const char name[])
{
    View & view = views[0];
    
    std::ostringstream oss;
    oss << "rgba";
    if ( view.buffered )    oss << " double";
    if ( view.depth_test )  oss << " depth";
    if ( view.stencil )     oss << " stencil~" << view.stencil;
    if ( view.multisample ) oss << " samples~" << view.multisample;
    if ( view.retina )      oss << " hidpi";
    std::string mode = oss.str();
    
    //std::clog << "GLUT string mode " << mode << '\n';
    
    // set GLUT display mode:
    glutInitDisplayString(mode.c_str());

    // set window size:
    //printf("Window: %i %i\n", view.window_size[0], view.window_size[1]);
    glutInitWindowSize(view.window_size[0], view.window_size[1]);
    glutInitWindowPosition(view.window_position[0], view.window_position[1]);
    
    // create window with title containing current working directory:
    int win = glutCreateWindow(name);
    assert_true( win > 0 );
    //std::clog << "new window " << win << '\n';

    // create new View for this window, duplicating the current View:
    if ( win >= (int)views.size() )
        views.resize(win+1, view);
    
    // switch the new view:
    views[win].window(win);
    views[win].initGL();
    
    // set window position:
    views[win].window_position[0] += 20;
    views[win].window_position[1] += 20;

    glutReshapeFunc(resizeWindow);
    glutKeyboardFunc(normalKeyCallback);
    glutSpecialFunc(specialKeyCallback);
    glutMouseFunc(processMouseClick);
    glutMotionFunc(processMouseDrag);
    glutPassiveMotionFunc(processPassiveMouseMotion);

    if ( win <= 1 )
        glutDisplayFunc(displayMain);
    else
        glutDisplayFunc(displayPlain);
    return win;
}

int glApp::newWindow(const char name[], int (*func)(View&))
{
    int win = newWindow(name);
    views[win].setDisplayFunc(func);
    return win;
}

int glApp::newWindow(const char name[], int (*func)(View&), void (*mag)(View&))
{
    int win = newWindow(name);
    views[win].setDisplayFunc(func);
    views[win].setDrawMagFunc(mag);
    return win;
}


/**
 This will not destroy the main window
 */
void glApp::deleteWindow(int win)
{
    if ( win == 0 )
        win = glutGetWindow();
    
    if ( 1 < win  &&  win < (int)views.size()  &&  views[win].window() > 0 )
    {
        //std::clog << "Destroy window " << win << '\n';
        assert_true( views[win].window() == win );
        glutDestroyWindow(views[win].window());
        views[win].window(0);
    }
}


void glApp::resizeWindow(int w, int h)
{
    unsigned win = glutGetWindow();
    flashText("window size %i %i", w, h);
    views[win].reshape(w, h);
    views[win].clearPixels();
    postRedisplay();
}


void glApp::setScale(float s)
{
    if ( 0 == views.size() )
        initialize();
    
    // update all window-associated views:
    for ( size_t n = 0; n < views.size(); ++n )
        views[n].view_scale = s;
}

#ifdef SAVE_IMAGE_H
int glApp::saveImage(char const* name, unsigned mag, unsigned downsample)
{
    View view = currentView();
    int W = mag * view.width();
    int H = mag * view.height();
    GLint vp[] = { 0, 0, W, H };
    if ( mag > 1 )
    {
        view.reshape(W, H);
        if ( OffScreen::openBuffer(W, H, 0) )
        {
            view.initGL();
            displayPlain();
            int res = SaveImage::saveImage(name, "png", vp, downsample);
            OffScreen::releaseBuffer();
            return res;
        }
    }
    return SaveImage::saveImage(name, "png", vp, downsample);
}
#endif

//------------------------------------------------------------------------------
#pragma mark -


/// this works even if no window is open
View& glApp::currentView()
{
    assert_true( views.size() > 0 );
    
    if ( views.size() <= 1 )
        return views[0];
    else
        return views[glutGetWindow()];
}


void glApp::resetView()
{
    assert_true( glutGetWindow() < (int)views.size() );
    glApp::currentView().reset();
}


//------------------------------------------------------------------------------
//------------------------------ keyboard commands -----------------------------
//------------------------------------------------------------------------------
#pragma mark -

void glApp::help(std::ostream& os)
{
    os << "                       Mouse Controls\n";
    os << "\n";
    os << " The display can be manipulated with click-and-drag movements,\n";
    os << " depending on the current `mode' selected by pressing `TAB':\n";
    os << "      Rotate                                     (3D only)\n";
    os << "      Translate in XY     (plane of the camera front lens)\n";
    os << "      Active                        (click to grab fibers)\n";
    os << "      Translate in XZ          (away/closer to the camera)\n";
    os << "      Spin & Zoom                   (Rotation in XY plane)\n";
    os << "      Select/Move region-of-interest\n";
    os << "\n";
    os << "  In the default mode, a SHIFT-CLICK can grab the filaments\n";
    os << "  A Menu is accessed by a right click\n";
    os << "  You might perhaps be able to zoom in/out with the mouse wheel\n";
    os << "\n";
    os << "                       Keyboard Controls\n\n";
    os << " + -         Zoom in and out; hold SHIFT for finer motion\n";
    os << " arrow keys  Translate; hold SHIFT for finer motion; hold ALT for rotation\n";
    os << " z           Reset view and refresh display\n";
    os << " h           Hide/show help\n";
    os << " b x         Show/hide a 10 um scale bar; Show/hide axes\n";
    os << " f ESC       Toggle fullscreen mode; exit full screen mode\n";
    os << " y           Save PNG image\n";
}


//------------------------------------------------------------------------------
void glApp::switchUserMode(int dir)
{
    int u = userMode;
    do {
        u = ( u + dir + MOUSE_PASSIVE ) % MOUSE_PASSIVE;
    } while ( actionDimensionality[u] > mDIM );

    userMode = (UserMode)u;
    switch ( userMode )
    {
        case MOUSE_TURN:      flashText("Mouse: Rotate");       break;
        case MOUSE_MOVE:      flashText("Mouse: Translate");    break;
        case MOUSE_ACTIVE:    flashText("Mouse: Active");       break;
        case MOUSE_PUSH:      flashText("Mouse: Translate-Z");  break;
        case MOUSE_SPIN:      flashText("Mouse: Spin & Zoom");  break;
        case MOUSE_MAGNIFIER: flashText("Mouse: Magnifier");    break;
        case MOUSE_EDIT_ROI:  flashText("Mouse: Edit ROI");     break;
        case MOUSE_SELECT:    flashText("Mouse: Select");       break;
        case MOUSE_PASSIVE:   flashText("Mouse: Passive");      break;
    }
}


///\todo flexible key assignment map to accomodate different keyboard layouts
void glApp::processNormalKey(unsigned char c, int modifiers)
{
    View & view = glApp::currentView();
    
    /* In the switch below:
     - use 'break' if the display needs a refresh
     - use 'return' if redrawing is not necessary.
    */
    switch (c)
    {
        case 17:
            if ( modifiers & GLUT_ACTIVE_CTRL )
                exit(EXIT_SUCCESS);
            break;
        
            
        case 9: // ascii 9 is horizontal TAB
        case 25: // ascii 25 is SHIFT-TAB
            switchUserMode(c==9 ? 1 : -1);
            break;
        
        
        case 27: // ascii 27 is ESCAPE
            if ( mFullScreen )
                toggleFullScreen();
            else
                deleteWindow(glutGetWindow());
            break;
        
            
        case 'f':
            toggleFullScreen();
            break;
        
        case 'F':
            toggleWindowSize();
            break;

        case 'z':
            view.reset();
            postRedisplay();
            break;
        
        
        case 'v':
            if ( mDIM == 3 )
            {
                view.slice = ( view.slice + 1 ) % 8;
                flashText("view:slice = %i", view.slice);
            }
            break;
        
        case 'V':
            view.toggleDepthClamp();
            flashText("view:depth_clamp = %i", view.depth_clamp);
            break;
            
        case 'b':
            view.scalebar = ( view.scalebar + 1 ) & 3;
            flashText("view:pixel_size = %.2f nm", 1000*view.pixelSize());
            break;
        
        case 'B':
            view.back_color.set(0,0,0);
            break;
            
        case 'i':
            view.invertColors();
            break;
            
        case 'I':
            view.blackAndWhite();
            break;

        case 'h':
            view.draw_memo = ( view.draw_memo + 1 ) % 3;
            if ( view.draw_memo == 2 )
            {
                std::ostringstream oss;
                help(oss);
                view.memo = oss.str();
            }
            else
                view.memo = "Please, visit www.cytosim.org";
            break;
        
        
        case 'x':
        case 'X':
            view.axes = ( view.axes ? 0 : mDIM );
            break;

#ifdef SAVE_IMAGE_H
        case 'y':
        if ( modifiers )
        {
            char name[1024] = { 0 };
            GLint vp[] = { 0, 0, view.width(), view.height() };
            snprintf(name, sizeof(name), "image%04i.png", imageIndex++);
            SaveImage::saveDepthBuffer(name, vp);
        }
        else
        {
            char name[1024] = { 0 };
            snprintf(name, sizeof(name), "image%04i.png", imageIndex++);
            saveImage(name, 1, 1);
        } break;
        
        case 'Y': {
            char name[1024] = { 0 };
            snprintf(name, sizeof(name), "image%04i.png", imageIndex++);
            saveImage(name, 4, 2);
        } break;
#endif
        //------------------------- Zoom in and out:
        
        case '_':
        case '-':
            if ( modifiers )
                view.zoom_out(1.071773463f);
            else
                view.zoom_out(1.4142135623f);
            break;
        
        case '=':
        case '+':
            if ( modifiers )
                view.zoom_in(1.071773463);
            else
                view.zoom_in(1.4142135623f);
            break;

        default:
            flashText("ignored key %i [%c]", c, c);
            return;
    }
    
    //if break was used, redisplay is needed:
    postRedisplay();
}


void glApp::processNormalKey(unsigned char c, int, int)
{
    processNormalKey(c, glutGetModifiers());
}

void glApp::normalKeyFunc(void (*func)(unsigned char, int, int))
{
    normalKeyCallback = func;
}

//------------------------------------------------------------------------------

/**
 arrow-keys controls translation, and
 arrow-keys with 'ALT' pressed controls rotation.
 
 motion is reduced by holding down SHIFT.
 */
void glApp::processSpecialKey(int key, int modifiers)
{
    Vector3 vec(0,0,0), dxy(0, 0, 0);
    View & view = glApp::currentView();
    real F = ( modifiers & GLUT_ACTIVE_SHIFT ) ? 0.0625 : 1;

    //std::clog << "special key " << key << ": " << modifiers << "  ";
    switch ( key )
    {
        case GLUT_KEY_HOME:      view.reset();            postRedisplay(); return;
        case GLUT_KEY_PAGE_UP:   view.zoom_in(1.4142f);   postRedisplay(); return;
        case GLUT_KEY_PAGE_DOWN: view.zoom_out(1.4142f);  postRedisplay(); return;
        case GLUT_KEY_LEFT:      dxy.set(-F,0,0); break;
        case GLUT_KEY_RIGHT:     dxy.set(+F,0,0); break;
        case GLUT_KEY_DOWN:      dxy.set(0,-F,0); break;
        case GLUT_KEY_UP:        dxy.set(0,+F,0); break;
    }

    // inverse the rotation of the current view:
    Quaternion<real> rot = view.rotation.conjugated();
    
    if ( modifiers & GLUT_ACTIVE_ALT )
    {
        rot.rotateVector(vec, cross(Vector3(0, 0, 1), dxy));
        view.rotation.rotateVector(dxy, view.focus);
        rot.setFromAxis(vec, F * (M_PI/8));
        view.rotate_by(rot);
        rot = view.rotation.conjugated();
        rot.rotateVector(vec, dxy);
        view.focus = vec;
    }
    else
    {
        // Translate view in the depth direction
        if ( (modifiers & GLUT_ACTIVE_CTRL) ^ (userMode == MOUSE_PUSH) )
            dxy.set(0.5*dxy.XX, 0, -0.5*dxy.YY);
        rot.rotateVector(vec, dxy);
        //std::clog << "vec " << dxy << " >>> " << vec << "\n";
        view.move_by((128*view.pixelSize())*vec);
    }
    postRedisplay();
}

void glApp::processSpecialKey(int key, int, int)
{
    processSpecialKey(key, glutGetModifiers());
}

void glApp::specialKeyFunc(void (*func)(int, int, int))
{
    specialKeyCallback = func;
}

//------------------------------------------------------------------------------
#pragma mark -

int buildFogMenu()
{
    int menu = gym::createMenu(glApp::processMenuEvent);
    gym::addMenuEntry("Disable",         100);
    gym::addMenuEntry("Linear ",         101);
    gym::addMenuEntry("Exponential 16x", 102);
    gym::addMenuEntry("Exponential 8x",  103);
    gym::addMenuEntry("Exponential 4x",  104);
    gym::addMenuEntry("Exponential 2x",  105);
    gym::addMenuEntry("Exponential 1x",  106);
    gym::addMenuEntry("Exponential 1/2x", 107);
    gym::addMenuEntry("Exponential 1/4x", 108);
    gym::addMenuEntry("Exponential 1/8x", 109);
    gym::addMenuEntry("Exponential 1/16x", 110);
    return menu;
}

int buildWindowSizeMenu()
{
    int menu = gym::createMenu(glApp::processMenuEvent);
    gym::addMenuEntry("256x256",   200);
    gym::addMenuEntry("384x384",   201);
    gym::addMenuEntry("512x256",   202);
    gym::addMenuEntry("512x384",   203);
    gym::addMenuEntry("512x512",   204);
    gym::addMenuEntry("768x256",   205);
    gym::addMenuEntry("768x512",   206);
    gym::addMenuEntry("768x768",   207);
    gym::addMenuEntry("1024x128",  208);
    gym::addMenuEntry("1024x256",  209);
    gym::addMenuEntry("1024x512",  210);
    gym::addMenuEntry("1024x768",  211);
    gym::addMenuEntry("1024x1024", 212);
    gym::addMenuEntry("-", 0);
    gym::addMenuEntry("426x240 (240p)",    220);
    gym::addMenuEntry("640x360 (360p)",    221);
    gym::addMenuEntry("854x480 (480p)",    222);
    gym::addMenuEntry("1280x720 (720p)",   223);
    gym::addMenuEntry("1920x1080 (1080p)", 224);
    gym::addMenuEntry("2560x1440 (1440p)", 225);
    return menu;
}


int buildClipMenu()
{
    int menu = gym::createMenu(glApp::processMenuEvent);
    gym::addMenuEntry("Disable",    300);
    
    gym::addMenuEntry(" X > 0",     301);
    gym::addMenuEntry(" X < 0",     302);
    gym::addMenuEntry("-1 < X < 1", 303);
    
    gym::addMenuEntry(" Y > 0",     311);
    gym::addMenuEntry(" Y < 0",     312);
    gym::addMenuEntry("-1 < Y < 1", 313);
    
    gym::addMenuEntry(" 0 < Z",     321);
    gym::addMenuEntry(" Z < 0",     322);
    gym::addMenuEntry(" 0.25 < Z",  323);
    gym::addMenuEntry(" Z < 0.25",  324);
    gym::addMenuEntry("-1 < Z < 1", 325);
    gym::addMenuEntry("-0.5 < Z < 0.5", 326);
    gym::addMenuEntry("-0.25 < Z < 0.25", 327);
    return menu;
}


int glApp::buildMenu()
{
    static int menu = 0;
    static int menu1, menu2, menu3;
    
    //std::clog << "buildMenu" << '\n';
    if ( !menu )
    {
        menu1 = buildFogMenu();
        menu2 = buildWindowSizeMenu();
        menu3 = buildClipMenu();
        menu = gym::createMenu(processMenuEvent);
        gym::addSubMenu("Fog",         menu1);
        gym::addSubMenu("Window Size", menu2);
        gym::addSubMenu("Slice",       menu3);
        gym::addMenuEntry("Reset View",         1);
        gym::addMenuEntry("Match ROI to View",  2);
        gym::addMenuEntry("Match View to ROI",  3);
        gym::addMenuEntry("Show/hide Scalebar", 4);
        gym::addMenuEntry("Show/hide XYZ-axes", 5);
        gym::addMenuEntry("Toggle fullscreen mode", 6);
        gym::addMenuEntry("Use 2D mouse gestures", 8);
        gym::addMenuEntry("Use 3D mouse gestures", 9);
        gym::addMenuEntry("Quit", 20);
    }
    return menu;
}

void glApp::attachMenu(int m)
{
    if ( !m ) buildMenu();
    gym::attachMenu(GLUT_RIGHT_BUTTON);
}

//------------------------------------------------------------------------------

/// shortcut
static void reshapeWindow(int w, int h) { glutReshapeWindow(w, h); }

void glApp::processMenuEvent(int item)
{
    View & view = glApp::currentView();
    switch( item )
    {
        case 0: return;
        case 1: view.reset(); break;
        case 2: view.adjustROI(nearZ); break;
        case 3: view.matchROI(); break;
        case 4: view.scalebar = ! view.scalebar; break;
        case 5: view.axes = (view.axes?0:mDIM ); break;
        case 6: toggleFullScreen(); break;
        case 8: setDimensionality(2); break;
        case 9: setDimensionality(3); break;

        case 20:  exit(EXIT_SUCCESS); break;
        
        case 100: view.enableFog(0, 0);      break;
        case 101: view.enableFog(1, 0);      break;
        case 102: view.enableFog(2, 0.0625); break;
        case 103: view.enableFog(2, 0.125);  break;
        case 104: view.enableFog(2, 0.25);   break;
        case 105: view.enableFog(2, 0.5);    break;
        case 106: view.enableFog(2, 1);      break;
        case 107: view.enableFog(2, 2);      break;
        case 108: view.enableFog(2, 4);      break;
        case 109: view.enableFog(2, 8);      break;
        case 110: view.enableFog(2, 16);     break;
            
        case 200: reshapeWindow(256, 256);   break;
        case 201: reshapeWindow(384, 384);   break;
        case 202: reshapeWindow(512, 256);   break;
        case 203: reshapeWindow(512, 384);   break;
        case 204: reshapeWindow(512, 512);   break;
        case 205: reshapeWindow(768, 256);   break;
        case 206: reshapeWindow(768, 512);   break;
        case 207: reshapeWindow(768, 768);   break;
        case 208: reshapeWindow(1024, 128);  break;
        case 209: reshapeWindow(1024, 256);  break;
        case 210: reshapeWindow(1024, 512);  break;
        case 211: reshapeWindow(1024, 768);  break;
        case 212: reshapeWindow(1024, 1024); break;
        case 213: reshapeWindow(1280, 640);  break;
        case 214: reshapeWindow(1280, 1280); break;
        case 220: reshapeWindow(426, 240);   break;
        case 221: reshapeWindow(640, 360);   break;
        case 222: reshapeWindow(854, 480);   break;
        case 223: reshapeWindow(1280, 720);  break;
        case 224: reshapeWindow(1920, 1080); break;
        case 225: reshapeWindow(2560, 1440); break;
        
        case 300:
            view.disableClipPlane(0);
            view.disableClipPlane(1);
            break;
            
        case 301:
            view.enableClipPlane(0,+1,0,0,0);
            view.disableClipPlane(1);
            break;
            
        case 302:
            view.enableClipPlane(0,-1,0,0,0);
            view.disableClipPlane(1);
            break;
            
        case 303:
            view.enableClipPlane(0,+1,0,0,1);
            view.enableClipPlane(1,-1,0,0,1);
            break;
 
        case 311:
            view.enableClipPlane(0,0,+1,0,0);
            view.disableClipPlane(1);
            break;
            
        case 312:
            view.enableClipPlane(0,0,-1,0,0);
            view.disableClipPlane(1);
            break;
            
        case 313:
            view.enableClipPlane(0,0,+1,0,1);
            view.enableClipPlane(1,0,-1,0,1);
            break;
 
        case 321:
            view.enableClipPlane(0,0,0,+1,0);
            view.disableClipPlane(1);
            break;
            
        case 322:
            view.enableClipPlane(0,0,0,-1,0);
            view.disableClipPlane(1);
            break;
 
        case 323:
            view.enableClipPlane(0,0,0,+1,0.25);
            view.disableClipPlane(1);
            break;
            
        case 324:
            view.enableClipPlane(0,0,0,-1,0.25);
            view.disableClipPlane(1);
            break;
            
        case 325:
            view.enableClipPlane(0,0,0,+1,1);
            view.enableClipPlane(1,0,0,-1,1);
            break;

        case 326:
            view.enableClipPlane(0,0,0,+1,0.5);
            view.enableClipPlane(1,0,0,-1,0.5);
            break;

        case 327:
            view.enableClipPlane(0,0,0,+1,0.25);
            view.enableClipPlane(1,0,0,-1,0.25);
            break;

        default: ABORT_NOW("unknown menu item");
    }
    postRedisplay();
}

//------------------------------------------------------------------------------
//--------------------------------  MOUSE  -------------------------------------
//------------------------------------------------------------------------------
#pragma mark -

void glApp::actionFunc(void (*func)(int, int, const Vector3 &, int))
{
    mouseClickCallback = func;
}

void glApp::actionFunc(void (*func)(int, int, Vector3 &, const Vector3 &, int))
{
    mouseDragCallback = func;
}

#if !defined(GLUT_WHEEL_UP)
#  define GLUT_WHEEL_UP    3
#  define GLUT_WHEEL_DOWN  4
#endif

//------------------------------------------------------------------------------
void glApp::processMouseClick(int button, int state, int mX, int mY)
{
    View & view = glApp::currentView();
    float cenX = 0.5 * view.width();
    float cenY = 0.5 * view.height();

    //printf("mouse button %i (%4i %4i) state %i key %i\n", button, mx, my, state, specialKeys);

    mouseS = state;
    mouseX = mX;
    mouseY = view.height()-mY;
    
    savedView = view; // copy the current Model-View transformation
    
    if ( state == GLUT_UP )
    {
         /*
         Zooming with the mouse-wheel requires an extended GLUT.
         http://iihm.imag.fr/blanch/howtos/MacOSXGLUTMouseWheel.html
         
         The zoom preserves the position pointed at by the mouse.
         */
        float wZ = 1.f;

        if ( button == GLUT_WHEEL_UP )
            wZ = 0.96969696f;
        if ( button == GLUT_WHEEL_DOWN )
            wZ = 1.031250f;

        if ( wZ != 1 )
        {
            /* 
             in 2D, we do not allow any shift in Z,
             and in 3d, we zoom in on the middle Z-plane
             */
            mouseDown = savedView.unproject(mouseX, mouseY, midZ);
            view.zoom_out(wZ);
            view.move_to((1-wZ)*mouseDown+wZ*view.focus);
            postRedisplay();
        }

        glutSetCursor(GLUT_CURSOR_INHERIT);
        mouseAction = MOUSE_PASSIVE;
        return;
    }

    glutSetCursor(GLUT_CURSOR_CROSSHAIR);
    
    // action is primarily decided by current mode
    mouseAction = userMode;
    specialKeys = glutGetModifiers();

    // change the mouse action if the CONTROL is pressed:
    if ( specialKeys & GLUT_ACTIVE_ALT )
    {
        switch ( userMode )
        {
            case MOUSE_MOVE: mouseAction = (mDIM==2)?MOUSE_SPIN:MOUSE_TURN; break;
            case MOUSE_SPIN: mouseAction = (mDIM==2)?MOUSE_MOVE:MOUSE_PUSH; break;
            case MOUSE_EDIT_ROI: mouseAction = MOUSE_MOVE; break;
            case MOUSE_TURN: mouseAction = MOUSE_MOVE; break;
            case MOUSE_PUSH:mouseAction = (mDIM==2)?MOUSE_MOVE:MOUSE_TURN;  break;
            default: break;
        }
    }
    //std::clog << "mouseAction " << mouseAction << " modifiers " << specialKeys << "\n";
    mouseDown = savedView.unproject(mouseX, mouseY, nearZ);
    viewFocus = savedView.unproject(cenX, cenY, nearZ);

    // change the mouse action because the shift key is down:
    if ( specialKeys & GLUT_ACTIVE_SHIFT )
    {
        mouseAction = MOUSE_ACTIVE;
        specialKeys ^= GLUT_ACTIVE_SHIFT;
    }
    
    if ( button == GLUT_RIGHT_BUTTON )
        mouseAction = MOUSE_MAGNIFIER;
    
    switch( mouseAction )
    {
        case MOUSE_MOVE:
            return;
            
        case MOUSE_PUSH:
        {
            axle = normalize(viewFocus - savedView.focus);
            Vector3 U = savedView.unproject(cenX, 2*cenY, nearZ);
            pole = normalize(U - viewFocus);
        } break;
            
        case MOUSE_TURN:
        {
            /* 
            Choose the amplification factor for mouse controlled rotation:
            for a value of one, the rotation exactly follows the mouse pointer 
            */
            const real amplification = 3.0;
            axle  = mouseDown - savedView.focus;
            axle *= amplification / axle.normSqr();
        } break;
            
        case MOUSE_SPIN:
        {
            pole = normalize(viewFocus - savedView.focus);
            axle = mouseDown - viewFocus;
        } break;
        
        case MOUSE_MAGNIFIER:
        {
            mouseDown = savedView.unproject(mouseX, mouseY, midZ);
            viewFocus = savedView.unproject(cenX, cenY, midZ);
            glutSetCursor(GLUT_CURSOR_NONE);
        } break;
            
        case MOUSE_EDIT_ROI:
        {
            mouseDown = savedView.unproject(mouseX, mouseY, midZ);
            ROIdup[0] = view.roi(0);
            ROIdup[1] = view.roi(1);
            moveROI = view.insideROI(mouseDown);
            flashText("click at %.4f %.4f %.4f", mouseDown.XX, mouseDown.YY, mouseDown.ZZ);
        } break;
            
        case MOUSE_ACTIVE:
        {
            if ( mouseClickCallback )
            {
                mouseDown = savedView.unproject(mouseX, mouseY, midZ);
                //std::clog << "Action down at "<<mouseDown<<'\n';
                mouseClickCallback(mouseX, mouseY, mouseDown, specialKeys);
            }
        }
        
        case MOUSE_SELECT:
        case MOUSE_PASSIVE:
            return;
    }
    postRedisplay();
}


//------------------------------------------------------------------------------
void glApp::processMouseDrag(int mX, int mY)
{
    //printf("mouse motion (%i %i) %i\n", mx, my, mouseAction);
    View & view = glApp::currentView();

    mouseX = mX;
    mouseY = view.height()-mY;

    Vector3 mouse = savedView.unproject(mouseX, mouseY, nearZ);

    switch( mouseAction )
    {
        case MOUSE_TURN:
        {
            /* we should rotate after: Q <- dQ * sQ, however dQ is defined in the
            reference frame rotated by sQ already, so dQ = sQ * dM * inv(sQ).
            This leads to the multiplication on the right: Q <- sQ * dM. */
            Quaternion<real> Q;
            Q.setFromAxis(cross(axle, mouse-mouseDown));
            view.rotate_to(savedView.rotation * Q);
        } break;
        
        
        case MOUSE_SPIN:
        {
            real C = dot(axle, mouse - viewFocus);
            real S = dot(pole, cross(axle, mouse - viewFocus));
            Quaternion<real> Q;
            Q.setFromAxis(pole, std::atan2(S, C));
            view.rotate_to(savedView.rotation * Q);
            float Z = norm( mouse - viewFocus ) / norm( mouseDown - viewFocus );
            if ( Z > 0.001 ) view.zoom_to(savedView.zoom * Z);
        } break;

        
        case MOUSE_MOVE:
        {
            view.move_to(savedView.focus - ( mouse - mouseDown ));
        } break;
        
        
        case MOUSE_PUSH:
        {
            real S = dot(mouse - mouseDown, pole);
            Vector3 move = mouse - mouseDown - S * ( axle + pole );
            view.move_to( savedView.focus - move );
        } break;
        
        
        case MOUSE_MAGNIFIER:
        {
            mouseDown = savedView.unproject(mouseX, mouseY, midZ);
        } break;

        
        case MOUSE_EDIT_ROI:
        {
            mouse = savedView.unproject(mouseX, mouseY, midZ);
            if ( moveROI )
            {
                Vector3 d = mouse - mouseDown;
                view.setROI(ROIdup[0]+d, ROIdup[1]+d);
                d = 0.5 * ( view.roi(0) + view.roi(1) );
                flashText("ROI center %.3f %.3f", d.XX, d.YY);
            }
            else
            {
                view.setROI(mouseDown, mouse);
                Vector3 d = view.roi(1) - view.roi(0);
                flashText("ROI %.3fx%.3f diagonal %.3f", d.XX, d.YY, d.norm());
            }
        } break;
            
        
        case MOUSE_ACTIVE:
        {
            if ( mouseDragCallback )
            {
                mouse = savedView.unproject(mouseX, mouseY, midZ);
                //std::clog << "Action move at " << mouse << '\n';
                mouseDragCallback(mouseX, mouseY, mouseDown, mouse, specialKeys);
            }
        } break;
        
        
        case MOUSE_SELECT:
        case MOUSE_PASSIVE: 
            break;
    }
    postRedisplay();
}

//------------------------------------------------------------------------------
void glApp::processPassiveMouseMotion(int mx, int my)
{
    //printf("passive mouse (%i %i)\n", mx, my);
}

//------------------------------------------------------------------------------
#pragma mark -

void glApp::flashText(std::string const& str)
{
    if ( views.size() > 1 )
    {
        View & view = glApp::currentView();
        if ( view.flashText(str) )
            glutPostWindowRedisplay(1);
    }
}


/**
 This is used for any secondary window.
 It does not show the interactive feedback to user.
 */
void glApp::displayPlain()
{
    CHECK_GL_ERROR("before glApp::displayPlain()");
    View & view = glApp::currentView();

    view.callDraw();

    if ( view.buffered )
        glutSwapBuffers();
    else
        glFlush();
}


/**
 This is used for the main window
 */
void glApp::displayMain()
{
    CHECK_GL_ERROR("before glApp::displayMain()");
    View & view = views[1];
    
    if ( 0 == view.callDraw() )
    {
#if 0
        static double millsec = 0;
        std::clog << "--- display " << TimeDate::milliseconds() - millsec << "\n";
        millsec = TimeDate::milliseconds();
#endif
        view.drawInteractiveFeatures();
        
        if ( mouseAction == MOUSE_MAGNIFIER )
        {
            //std::cerr << "mag@ " << mouseX << " " << mouseY << "\n";
            view.drawMagnifier(7, mouseDown, viewFocus, mouseX, mouseY, 256);
        }
        
        CHECK_GL_ERROR("in glApp::displayMain()");
        glFlush();
        if ( view.buffered )
            glutSwapBuffers();
    }
    drawRefresh = 0;
}


void glApp::displayOtherWindows(int (*drawFunc)(View&))
{
    for ( size_t n = 2; n < views.size(); ++n )
    {
        if ( views[n].window() > 0 )
        {
            glutSetWindow(views[n].window());
            drawFunc(views[n]);
            if ( n == 1 )
                views[n].drawInteractiveFeatures();
            if ( views[n].buffered )
                glutSwapBuffers();
            else
                glFlush();
        }
    }
}


void glApp::postRedisplay()
{
    //std::clog << " postRedisplay\n";
    glutPostRedisplay();
    ++drawRefresh;
}


void glApp::setMessage(std::string const& msg)
{
    glApp::currentView().setMessage(msg);
}


void glApp::setSubtitle(std::string const& msg)
{
    glApp::currentView().setSubtitle(msg);
}

