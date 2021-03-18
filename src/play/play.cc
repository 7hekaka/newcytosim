// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University

#include "simul_prop.h"
#include "glossary.h"
#include "messages.h"
#include "offscreen.h"
#include "saveimage.h"
#include "filepath.h"
#include "splash.h"
#include "print_color.h"

#include "opengl.h"
#include "player.h"
#include "view.h"
#include "gle.h"

Player player;

SimThread& thread = player.thread;
Simul&      simul = player.simul;
PlayerProp&  prop = player.prop;
DisplayProp& disp = player.disp;

void drawLive(View& view, int mag);

/// create a player suitable for command-line offscreen rendering only
#define HEADLESS_PLAYER 0

#if HEADLESS_PLAYER
void helpKeys(std::ostream& os) { os << "This is a headless display\n"; }
#else
#  include "glut.h"
#  include "glapp.h"
#  include "fiber_prop.h"
#  include "fiber_disp.h"
#  include "point_disp.h"
using glApp::flashText;
void buildMenus();
#  include "play_keys.cc"
#  include "play_menus.cc"
#  include "play_mouse.cc"
#endif

void goodbye()
{
    //printf("Goodbye...\n");
    player.clear();
}

//------------------------------------------------------------------------------
#pragma mark - Display

/**
 call drawScene() if data can be accessed by current thread
 */
void drawLive(View& view, int mag)
{
    CHECK_GL_ERROR("before drawLive");
    if ( 0 == thread.trylock() )
    {
        thread.executePipedCommands(32);
        player.drawScene(view, mag);
        thread.unlock();
    }
    else
    {
        //thread.debug("display: trylock failed");
        //postRedisplay();
    }
}


// This must be a plain C-function
void drawOffscreen(View & view, int mag)
{
    //std::clog << "drawOffscreen " << glApp::views.size() << '\n';
    player.drawScene(view, mag);
}


/// copy data from multisample to normal buffer
void blitBuffers(GLuint normal, GLuint multi, GLint W, GLint H)
{
    //std::clog << "blitting multisample buffer\n";
    glBindFramebuffer(GL_READ_FRAMEBUFFER, multi);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, normal);
    glBlitFramebuffer(0, 0, W, H, 0, 0, W, H, GL_COLOR_BUFFER_BIT, GL_NEAREST);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, multi);
    glBindFramebuffer(GL_READ_FRAMEBUFFER, normal);
}

//------------------------------------------------------------------------------
#pragma mark - main


void help(std::ostream& os)
{
    os << "play [OPTIONS] [PATH] [FILE]\n"
          "     live                     enter live simulation mode directly\n"
          "     PATH                     change working directory as specified\n"
          "     FILE.cym                 specify input configuration file\n"
          "     FILE.cmo                 specify trajectory file\n"
          "     FILE.cyp                 specify display configuration file\n"
          "     PARAMETER=value          set parameter value (example size=512)\n"
          "     image frame=INT          render specified frame offscreen\n"
          "     image frame=INT,INT,...  render several frames offscreen\n"
          "     image magnify=INT        render frames at higher resolution\n"
          "     movie                    render all frames offscreen\n"
          "     'movie on', 'image on'   use on-screen rendering\n"
          "     movie period=INT         render one frame every INT frames\n"
          " (there should be no whitespace around the equal sign)\n";
}


/// different modes:
enum Style { ONSCREEN, OFFSCREEN };
enum Mode { NORMAL, SAVE_IMAGE, SAVE_MOVIE };


void print_error(Exception const& e)
{
    print_magenta(std::cerr, e.brief());
    std::cerr << e.info() << '\n';
}


int main(int argc, char* argv[])
{
    Style style = OFFSCREEN;
    Mode mode = NORMAL;
    int magnify = 1;
    Glossary arg;
    
    Cytosim::out.silent();
    Cytosim::log.silent();
    
    std::atexit(goodbye);

    if ( arg.read_strings(argc-1, argv+1) )
        return EXIT_FAILURE;
    
    // check for major options:
    if ( arg.use_key("-") )
        Cytosim::warn.silent();

    if ( arg.use_key("help") )
    {
        splash(std::cout);
        help(std::cout);
        return EXIT_SUCCESS;
    }

    if ( arg.use_key("info") || arg.use_key("--version") )
    {
        splash(std::cout);
        print_version(std::cout);
        if ( SaveImage::supported("png") )
            std::cout << "   PNG enabled\n";
        return EXIT_SUCCESS;
    }
    
    if ( arg.use_key("live") || arg.has_key(".cym") )
        prop.goLive = true;

    if ( arg.use_key("image") )
        mode = SAVE_IMAGE;

    if ( arg.use_key("poster") )
    {
        mode = SAVE_IMAGE;
        magnify = 3;
    }
    
    if ( arg.use_key("movie") )
        mode = SAVE_MOVIE;
    
    if ( arg.use_key("on") )
        style = ONSCREEN;
    
    if ( arg.use_key("tga") )
        prop.image_format = "tga";
    
    if ( arg.use_key("png") )
        prop.image_format = "png";

    // get image over-sampling:
    arg.set(magnify, "magnify") || arg.set(magnify, "mag");

    // change working directory if specified:
    if ( arg.has_key("directory") )
    {
        FilePath::change_dir(arg.value("directory"));
        //std::clog << "Cytosim working directory is " << FilePath::get_cwd() << '\n';
    }
    
    // The user can specify a frame index to be loaded:
    size_t frm = 0;
    bool has_frame = false;

    try
    {
        has_frame = arg.set(frm, "frame");
        simul.initialize(arg);
    }
    catch( Exception & e )
    {
        print_error(e);
        return EXIT_FAILURE;
    }
    
#if HEADLESS_PLAYER
    View view("*");
    view.setDisplayFunc(drawOffscreen);
#else
    glApp::setDimensionality(DIM);
    if ( arg.use_key("fullscreen") )
        glApp::setFullScreen(1);
    View& view = glApp::views[0];
    view.setDisplayFunc(drawLive);
#endif
    
    //---------Open trajectory file and load simulation world

    try
    {
        std::string file = simul.prop->config_file;
        
        if ( ! prop.goLive || has_frame )
        {
            // config file is `properties.cmo` in replay mode:
            file = simul.prop->property_file;
            
            // read config to create all properties
            Parser(simul, 1, 1, 0, 0, 0).readConfig(file);
            
            // open trajectory file and load requested frame:
            thread.openFile(simul.prop->trajectory_file);
            
            if ( thread.loadFrame(frm) )
            {
                // load last frame in file:
                if ( thread.loadLastFrame() )
                    std::cerr << "Warning: could only load frame " << thread.currentFrame() << ' ';
            }
            frm = thread.currentFrame();
        }
        else
        {
            // get the name of 'simul' and simul:display from config file
            Parser(simul, 0, 1, 0, 0, 0).readConfig(file);
        }
    }
    catch( Exception & e )
    {
        print_error(e);
        return EXIT_FAILURE;
    }

    //---------Open setup file and read display parameters from command line

    try
    {
        std::string setup;
        // check for play's configuration file specified on the command line:
        if ( arg.set(setup, ".cyp") )
        {
            // extract "simul:display" from setup
            if ( FilePath::is_file(setup) )
                Parser(simul, 0, 1, 0, 0, 0).readConfig(setup);
            else
                std::cerr << " warning: could not read `" << setup << "'\n";
        }
        
        // read settings, but keep anything set on the command-line:
        arg.read(simul.prop->display, 1);
        simul.prop->display_fresh = false;
        
        if ( !arg.empty() )
        {
            view.read(arg);
            disp.read(arg);
            prop.read(arg);
        }
    }
    catch( Exception & e )
    {
        arg.print_warning(std::cerr, 1, "\n");
        print_error(e);
        return EXIT_FAILURE;
    }
    
#ifndef __APPLE__
    // it is necessary under Linux/Windows to initialize GLUT to display fonts
    glutInit(&argc, argv);
#endif
    
    //-------- off-screen (non interactive) rendering -------
    
    if ( mode != NORMAL && style == OFFSCREEN )
    {
        const int W = view.width() * magnify;
        const int H = view.height() * magnify;
        
        if ( !OffScreen::openContext() )
        {
            std::cerr << "Failed to create off-screen context\n";
            return EXIT_FAILURE;
        }
        GLuint fbo = OffScreen::createBuffer(W, H, 0);
        if ( !fbo )
        {
            std::cerr << "Failed to create off-screen pixels\n";
            return EXIT_FAILURE;
        }
        GLuint multi = 0;
        if ( view.multisample > 1 )
        {
            multi = OffScreen::createBuffer(W, H, view.multisample);
        }
        
        gle::initialize();
        player.setStyle(disp.style);
        view.initGL();

        if ( mode == SAVE_IMAGE )
        {
            size_t inx = 0;
            // it is possible to specify multiple frame indices:
            do {
                thread.loadFrame(frm);
                // only save requested frames:
                if ( thread.currentFrame() == frm )
                {
                    player.drawScene(view, magnify);
                    if ( multi )
                        blitBuffers(fbo, multi, W, H);
                    player.saveView("image", frm, prop.downsample);
                }
            } while ( arg.set(frm, "frame", ++inx) );
        }
        else if ( mode == SAVE_MOVIE )
        {
            // save every prop.period
            unsigned s = prop.period;
            do {
                if ( ++s >= prop.period )
                {
                    player.drawScene(view, magnify);
                    if ( multi )
                        blitBuffers(fbo, multi, W, H);
                    player.saveView("movie", frm++, prop.downsample);
                    s = 0;
                }
            } while ( 0 == thread.loadNextFrame() );
        }
        printf("\n");
        OffScreen::closeContext();
        arg.print_warning(std::cerr, 1, "\n");
        return EXIT_SUCCESS;
    }
    
    arg.print_warning(std::cerr, 1, "\n");

    //--------- initialize Window system and create Window
#if ( HEADLESS_PLAYER )
    print_green(std::cout, "This player can only do offscreen rendering.\n");
#else
#ifdef __APPLE__
    glutInit(&argc, argv);
#endif
    
    //register all the GLUT callback functions:
    glApp::actionFunc(processMouseClick);
    glApp::actionFunc(processMouseDrag);
    glApp::normalKeyFunc(processNormalKey);
    glApp::newWindow(drawLive);

    if ( mode == SAVE_IMAGE )
    {
        prop.auto_exit = 2;
        prop.save_images = 1;
        prop.play = 1;
    }
    
    if ( mode == SAVE_MOVIE )
    {
        prop.auto_exit = 1;
        prop.save_images = 9999;
        prop.play = 1;
    }

    //-------- initialize graphical user interface and graphics

    try
    {
        gle::initialize();
        player.setStyle(disp.style);
        buildMenus();
        glutAttachMenu(GLUT_RIGHT_BUTTON);
        glutMenuStatusFunc(menuCallback);
        if ( glApp::isFullScreen() )
            glutFullScreen();
        glutTimerFunc(200, timerCallback, 0);
    }
    catch ( Exception & e )
    {
        print_error(e);
        return EXIT_FAILURE;
    }
    
    if ( prop.goLive )
    {
        thread.period(prop.period);
        try
        {
            if ( has_frame )
                thread.extend();
            else
                thread.start();
        }
        catch( Exception & e )
        {
            print_error(e);
            return EXIT_FAILURE;
        }
    }

    //start the GLUT event handler:
    glutMainLoop();
#endif
}
