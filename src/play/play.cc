// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University

#include "time_date.h"
#include "simul.h"
#include "glossary.h"
#include "messages.h"
#include "offscreen.h"
#include "save_image_gl.h"
#include "filepath.h"
#include "splash.h"
#include "print_color.h"
#include "signal_handlers.h"

#include "gle.h"
#include "player.h"
#include "view.h"

Player player;

#include "random_pcg.h"
#include "random_seed.h"
using namespace PCG32;

SimThread& worker = player.worker;
Simul&      simul = player.simul;
PlayerProp&  prop = player.prop;
DisplayProp& disp = player.disp;

int drawSimulation(View& view);

/// create a player capable of command-line offscreen rendering only
#define HEADLESS_PLAYER 0

#if HEADLESS_PLAYER
void helpKeys(std::ostream& os) { os << "This executable has no display capability\n"; }
void buildMenus() { }
#else
#  include "glut.h"
#  include "glapp.h"
#  include "fiber_prop.h"
#  include "fiber_disp.h"
#  include "point_disp.h"
using glApp::flashText;
#  include "play_keys.cc"
#  include "play_menus.cc"
#  include "play_mouse.cc"
#endif

//------------------------------------------------------------------------------
#pragma mark - Display


/// this does not draw
int drawNot(View& view)
{
    //std::clog << " drawNot(" << std::setprecision(3) << simul.time() << "s)\n";
    //view.clearPixels();
    //player.prepareDisplay(view);
    return 0;
}


/// minimalistic display function
void drawMag(View& view)
{
    //std::clog << " drawMag(" << std::setprecision(3) << simul.time() << "s)\n";
    view.clearPixels();
    view.loadView();
    view.setLights();
    view.setClipping();
    player.setPixelSize(view);
    player.drawCytosim();
    view.endClipping();
}


/**
 calls prepareDisplay() and drawSystem()
 */
int drawSimulation(View& view)
{
    if ( simul.prop.display_fresh )
    {
        player.readDisplayString(view, simul.prop.display);
        simul.prop.display_fresh = false;
    }
    //worker.debug("drawSimulation");
    player.prepareDisplay(view);
    player.drawSystem(view);
    return 0;
}


/**
 call drawSimulation() if data can be accessed by current thread
 */
int drawLive(View& view)
{
    //std::clog << " drawLive(" << std::setprecision(3) << simul.time() << "s)\n";
    if ( 0 == worker.trylock() )
    {
        worker.read_input();
        drawSimulation(view);
        worker.unlock();
        return 0;
    }
    else
    {
#if 0
        static double sec = TimeDate::milliseconds();
        double now = TimeDate::milliseconds();
        fprintf(stderr, "drawLive(failed) %.0f\n", now-sec);
        sec = now;
#endif
        //worker.debug("display: trylock failed");
        //postRedisplay();
    }
    return 1;
}


/// copy color data from 'back' to 'front'
void blitBuffers(unsigned back, unsigned front, int W, int H)
{
    //std::clog << "blitting multisample buffer\n";
    glBindFramebuffer(GL_READ_FRAMEBUFFER, back);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, front);
    glBlitFramebuffer(0, 0, W, H, 0, 0, W, H, GL_COLOR_BUFFER_BIT, GL_NEAREST);
    glBindFramebuffer(GL_READ_FRAMEBUFFER, front);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, back);
}


void byebye()
{
    worker.cancel_join();
}

//------------------------------------------------------------------------------
#pragma mark - main

/// different modes:
enum Style { ONSCREEN, OFFSCREEN };
enum Mode { SAVE_NOTHING = 0, SAVE_IMAGE = 1, SAVE_MOVIE = 2 };


void help(std::ostream& os)
{
    os << "play [OPTIONS] [PATH] [FILE]\n"
          "   live                     enter live simulation mode directly\n"
          "   PATH                     change working directory as specified\n"
          "   FILE.cym                 enter live mode with specified config file\n"
          "   FILE.cmo                 specify trajectory file in replay mode\n"
          "   FILE.cyp                 specify display configuration file\n"
          "   PARAMETER=value          set parameter value (example size=512)\n"
          "   image frame=INT          render specified frame offscreen\n"
          "   image frame=INT,INT,...  render several frames offscreen\n"
          "   image magnify=INT        render frames at higher resolution\n"
          "   movie                    render all frames in trajectory file\n"
          "   'movie on', 'image on'   use on-screen rendering\n"
          "   movie period=INT         render one frame every INT frames\n"
          " (there should be no whitespace on either side of the equal sign)\n";
}


void print_error(Exception const& e)
{
    print_magenta(stderr, e.brief());
    fputs(e.info().c_str(), stderr);
    putc('\n', stderr);
}


int main(int argc, char* argv[])
{
    Style style = OFFSCREEN;
    int mode = SAVE_NOTHING;
    int fullscreen = 0;
    std::string name = "Cytosim     ";
    name += FilePath::get_cwd();
    
    Cytosim::out.silent();
    Cytosim::log.silent();
    register_signal_handlers();
    
#if HEADLESS_PLAYER
    View view("*", DIM==3);
    view.setDisplayFunc(drawSimulation);
#else
    glApp::setDimensionality(DIM);
    View& view = glApp::views[0];
#endif
    if ( DIM == 3 )
        view.back_color = 0x161616FF;
    else
        view.depth_test = 0;

    Glossary arg;
    if ( arg.read_strings(argc-1, argv+1) )
        return 1;
    
    if ( arg.use_key("fullscreen") )
        fullscreen = 1;

    // check for major options:
    if ( arg.use_key("-") )
    {
        simul.prop.verbose = 0;
        Cytosim::out.silent();
        Cytosim::log.silent();
        Cytosim::warn.silent();
    }
    else if ( arg.use_key("+") )
    {
        simul.prop.verbose = 3;
    }
    
    if ( arg.use_key("help") )
    {
        splash(std::cout);
        help(std::cout);
        return 0;
    }

    if ( arg.use_key("info") || arg.use_key("--version") )
    {
        splash(std::cout);
        print_version(std::cout);
        return 0;
    }
    
    int menu = 1;
    arg.set(menu, "menu");

    if ( arg.use_key("live") || arg.has_key(".cym") )
    {
        arg.peek(simul.prop.config_file, ".cym");
        player.goLive = true;
    }
    
    if ( arg.has_key(".cmi") )
        simul.prop.system_file = arg.value(".cmi");
    if ( arg.has_key(".cmo") )
        simul.prop.system_file = arg.value(".cmo");

    if ( arg.use_key("image") )
    {
        prop.image_name = "image%";
        mode = SAVE_IMAGE;
    }
    
    if ( arg.use_key("poster") )
    {
        prop.image_name = "poster%";
        mode = SAVE_IMAGE;
        view.magnify = 3;
    }
    
    if ( arg.use_key("on") )
        style = ONSCREEN;
    
    if ( arg.set(prop.image_name, ".ppm") )
    {
        prop.image_format = "ppm";
        mode |= SAVE_IMAGE;
    }
    else if ( arg.use_key("ppm") )
        prop.image_format = "ppm";
    
    if ( arg.set(prop.image_name, ".png") )
    {
        prop.image_format = "png";
        mode |= SAVE_IMAGE;
    }
    else if ( arg.use_key("png") )
        prop.image_format = "png";
    
    if ( arg.set(prop.image_name, ".tga") )
    {
        prop.image_format = "tga";
        mode |= SAVE_IMAGE;
    }
    else if ( arg.use_key("tga") )
        prop.image_format = "tga";
    
    if ( arg.use_key("movie") )
    {
        if ( mode == SAVE_NOTHING )
            prop.image_name = "movie%";
        mode = SAVE_MOVIE;
    }

    // change working directory if specified:
    if ( arg.set(name, "directory") )
    {
        FilePath::change_dir(name);
        //std::clog << "Cytosim working directory is " << FilePath::get_cwd() << '\n';
    }
    
    // The user can specify a frame index to be loaded:
    size_t frm = 0;
    bool has_frame = false;

    try
    {
        seed_pcg32();
        RNG.seed(pcg32());
        has_frame = arg.set(frm, "frame");
    }
    catch( Exception & e )
    {
        print_error(e);
        return 2;
    }

    //---------Open trajectory file and load simulation world

    try
    {
        if ( ! player.goLive || has_frame )
        {
            // read the `properties` file to import properties
            simul.loadProperties();
            
            // open trajectory file:
            worker.openFile(simul.prop.system_file);
            
            // load requested frame:
            if ( worker.loadFrame(frm) )
            {
                // load last frame in file:
                if ( worker.loadLastFrame() )
                    std::cerr << "Warning: could only load frame " << worker.currentFrame() << ' ';
            }
            frm = worker.currentFrame();
        }
        else
        {
            // get the name of 'simul' object and 'simul:display' from the config file
            Parser(&simul, 1, 0, 0, 0, 0, 0).readConfig();
        }

        // read Simul parameters from command line
        simul.prop.read(arg);
    }
    catch( Exception & e )
    {
        arg.print_warnings(stderr, 1, "\n");
        print_error(e);
        return 3;
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
                Parser(&simul, 1, 0, 0, 0, 0, 0).readConfig(setup);
            else
                std::cerr << " warning: could not read `" << setup << "'\n";
        }
        
        // read settings, but keep anything already set on the command-line:
        arg.read(simul.prop.display, 1);
        simul.prop.display_fresh = false;
        //arg.print(std::cout);
        
        if ( !arg.empty() )
        {
            view.read(arg);
            disp.read(arg);
            prop.read(arg);
        }
    }
    catch( Exception & e )
    {
        arg.print_warnings(stderr, 1, "\n");
        print_error(e);
        return 4;
    }
    
    //-------- off-screen (non interactive) rendering -------
    
    if ( mode != SAVE_NOTHING && style == OFFSCREEN )
    {
        view.resize(view.magnify);
        const int W = view.width();
        const int H = view.height();

        if ( OffScreen::openContext() )
        {
            std::cerr << "Failed to create off-screen context\n";
            return 5;
        }
        unsigned fbo = OffScreen::openBuffer(W, H, 0);
        if ( !fbo )
        {
            std::cerr << "Failed to create off-screen pixels\n";
            return 6;
        }
        unsigned multi = 0;
        if ( view.multisample > 1 )
        {
            multi = OffScreen::openBuffer(W, H, view.multisample);
        }
        gle::initialize();
        player.setStyle(disp.style);
        view.initGL();

        if ( mode & SAVE_MOVIE )
        {
            // save every prop.period
            unsigned s = prop.period;
            do {
                if ( ++s >= prop.period )
                {
                    drawSimulation(view);
                    if ( multi )
                        blitBuffers(multi, fbo, W, H);
                    player.saveView(view, frm++, prop.downsample);
                    s = 0;
                }
            } while ( 0 == worker.loadNextFrame() );
        }
        else if ( mode & SAVE_IMAGE )
        {
            size_t inx = 0;
            // it is possible to specify multiple frame indices:
            do {
                worker.loadFrame(frm);
                // only save requested frames:
                if ( worker.currentFrame() == frm )
                {
                    drawSimulation(view);
                    if ( multi )
                        blitBuffers(multi, fbo, W, H);
                    player.saveView(view, frm, prop.downsample);
                }
            } while ( arg.set(frm, "frame", ++inx) );
        }
        if ( simul.prop.verbose ) printf("\n");
        arg.print_warnings(stderr, 1, "\n");
        OffScreen::closeContext();
        return 0;
    }
    
    arg.print_warnings(stderr, 1, "\n");

    //--------- initialize Window system and create Window
#if ( HEADLESS_PLAYER )
    print_green(std::cout, "This player can only do offscreen rendering.\n");
#else
    glutInit(&argc, argv);
    //register all the GLUT callback functions:
    glApp::actionFunc(processMouseClick);
    glApp::actionFunc(processMouseDrag);
    glApp::normalKeyFunc(processNormalKey);
    glApp::newWindow(name.c_str(), drawLive, drawMag);

    if ( mode & SAVE_MOVIE )
    {
        prop.auto_exit = 1;
        prop.save_images = 9999;
        prop.replay = 1;
    }
    else if ( mode & SAVE_IMAGE )
    {
        prop.auto_exit = 2;
        prop.save_images = 1;
        prop.replay = 1;
    }
    
    //-------- initialize graphical user interface and graphics

    try
    {
        gle::initialize();
        player.setStyle(disp.style);
        if ( menu )
            buildMenus();
        if ( fullscreen )
            glutFullScreen();
        glutTimerFunc(100, timerCallback, 0);
    }
    catch ( Exception & e )
    {
        print_error(e);
        return 7;
    }
    
    if ( player.goLive )
    {
        worker.period(prop.period);
        try
        {
            if ( has_frame )
                worker.prolong();
            else
                worker.start();
        }
        catch( Exception & e )
        {
            print_error(e);
            return 8;
        }
    }

    std::atexit(byebye);
    //start the GLUT event handler:
    glutMainLoop();
#endif
}
