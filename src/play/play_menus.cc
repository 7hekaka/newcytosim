// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University

#include "gym_menu.h"

static void processMenuFiber(int item)
{
    FiberDisp * FD = player.firstFiberDisp();
    
    if ( FD )
    {
        switch (item)
        {
            case 0: break;
            case 1: FD->line_style = FD->line_style?0:1; break;
            case 2: FD->line_style = FD->line_style==2?0:2; break;
                
            case 3: FD->point_style = !FD->point_style; break;
            case 5: FD->point_style = FD->point_style==2?0:2; break;
            case 6: FD->point_style = FD->point_style==3?0:3; break;

            case 7: FD->end_style[1] = 3*!FD->end_style[1]; break;
            case 8: FD->end_style[0] = 2*!FD->end_style[0]; break;
                
            case 9: FD->force_style = !FD->force_style; break;
            case 10: FD->visible = !FD->visible; break;
                
            case 20: FD->coloring = FiberDisp::COLORING_OFF; break;
            case 21: FD->coloring = FiberDisp::COLORING_RANDOM; break;
            case 22: FD->coloring = FiberDisp::COLORING_MARK; break;
            case 23: FD->coloring = FiberDisp::COLORING_FLAG; break;
            case 24: FD->coloring = FiberDisp::COLORING_FAMILY; break;
            case 25: FD->coloring = FiberDisp::COLORING_CLUSTER; break;
            case 26: FD->coloring = FiberDisp::COLORING_DIRECTION; break;
            case 27: FD->coloring = FiberDisp::COLORING_AGE; break;
            case 28: FD->coloring = FiberDisp::COLORING_PSTATE; break;

            case 30: FD->draw_average = 0; break;
            case 31: FD->draw_average = 1; break;
            case 32: FD->draw_average = 2; break;
                
            default:
                std::cerr << "CYTOSIM ERROR: unknown menu code" << item << '\n';
                return;
        }
        player.refresh();
    }
}


static int buildMenuFiber()
{
    static int menuID = 0;
    if ( menuID == 0 )
        menuID = gym::createMenu(processMenuFiber);
    else
        gym::clearMenu(menuID);
    
    FiberDisp * FD = player.firstFiberDisp();
    if ( FD )
    {
        gym::addMenuEntry(FD->visible       ?"Hide"           :"Show",          10);
        gym::addMenuEntry(FD->line_style    ?"Hide Lines"     :"Show Lines",     1);
        gym::addMenuEntry(FD->line_style==2 ?"Hide Tensions"  :"Show Tensions",  2);
        gym::addMenuEntry(FD->point_style   ?"Hide Points"    :"Show Points",    3);
        gym::addMenuEntry(FD->point_style==2?"Hide Arrows"    :"Show Arrows",    5);
        gym::addMenuEntry(FD->point_style==3?"Hide Chevrons"  :"Show Chevrons",  6);
        gym::addMenuEntry(FD->end_style[1]  ?"Hide Minus-ends":"Show Minus-end", 7);
        gym::addMenuEntry(FD->end_style[0]  ?"Hide Plus-ends" :"Show Plus-end",  8);
        gym::addMenuEntry(FD->force_style   ?"Hide Forces"    :"Show Forces", 9);
        gym::addMenuEntry("No coloring",           20);
        gym::addMenuEntry("Coloring by number",    21);
        gym::addMenuEntry("Coloring by mark",      22);
        gym::addMenuEntry("Coloring by flag",      23);
        gym::addMenuEntry("Coloring by family",    24);
        gym::addMenuEntry("Coloring by cluster",   25);
        gym::addMenuEntry("Coloring by direction", 26);
        gym::addMenuEntry("Coloring by age",       27);
        gym::addMenuEntry("Coloring by +end state",28);
        gym::addMenuEntry("draw_average=0", 30);
        gym::addMenuEntry("draw_average=1", 31);
        gym::addMenuEntry("draw_average=2", 32);
    }
    else
        gym::addMenuEntry("no fiber?", 0);

    return menuID;
}

//------------------------------------------------------------------------------
static void processMenuCouple(int item)
{
    switch (item)
    {
        case 0: return;
        case 1: disp.couple_select = 0; break;
        case 2: disp.couple_select = 1; break;
        case 3: disp.couple_select = 2; break;
        case 4: disp.couple_select = 4; break;
        default:
            std::cerr << "CYTOSIM ERROR: unknown menu code" << item << '\n';
            return;
    }
    player.refresh();
}

static int buildMenuCouple()
{
    static int menuID = 0;
    if ( menuID == 0 )
        menuID = gym::createMenu(processMenuCouple);
    else
        gym::clearMenu(menuID);
    
    gym::addMenuEntry("Hide all",    1);
    gym::addMenuEntry("Show free",   2);
    gym::addMenuEntry("Show bound",  3);
    gym::addMenuEntry("Show links",  4);
    return menuID;
}

//------------------------------------------------------------------------------
static void processMenuDisplay(int item)
{
    View & view = glApp::currentView();
    switch (item)
    {
        case 0: return;
        case 1: view.reset(); break;
        case 3: disp.tile = ( disp.tile ? 0 : 8 ); break;
        case 4: glApp::toggleFullScreen(); break;
        case 6: view.track_fibers = !view.track_fibers; break;
        
        case 101: player.setStyle(1); break;
        case 102: player.setStyle(2); break;
        case 103: player.setStyle(3); break;
            
        default:
            std::cerr << "CYTOSIM ERROR: unknown menu code" << item << '\n';
            return;
    }
    player.refresh();
}


static int buildMenuStyle()
{
    static int menuID = 0;
    if ( menuID == 0 )
    {
        menuID = gym::createMenu(processMenuDisplay);
        gym::addMenuEntry("Detailed (style 1)", 101);
        gym::addMenuEntry("Fastest (style 2)", 102);
        gym::addMenuEntry("Best Looking (style 3)", 103);
    }
    return menuID;
}


static int buildMenuDisplay()
{
    static int menuID = 0;
    int m0 = buildMenuStyle();
    int m1 = buildMenuFiber();
    int m2 = buildMenuCouple();
    
    if ( menuID == 0 )
        menuID = gym::createMenu(processMenuDisplay);
    else
        gym::clearMenu(menuID);
    
    gym::addMenuEntry("Reset View",  1);
    gym::addSubMenu("Style",   m0);
    gym::addSubMenu("Fibers",  m1);
    gym::addSubMenu("Couple",  m2);
    
    View const& view = glApp::currentView();
    gym::addMenuEntry("Toggle fullscreen mode (f)", 4);
    gym::addMenuEntry(disp.tile?"Non-tiled Display":"Tiled Display", 3);
    gym::addMenuEntry(view.track_fibers?"stop tracking":"Track Fibers", 6);
    
    return menuID;
}


//------------------------------------------------------------------------------
#pragma mark -

static void processMenuFiberSelect(int item)
{
    FiberDisp * FD = player.firstFiberDisp();
    if ( FD )
    {
        switch (item)
        {
            case 0: return;
            case 1: FD->hide  = 0; break;
            case 2: FD->hide ^= 1; break;
            case 3: FD->hide ^= 2; break;
            default:
                std::cerr << "CYTOSIM ERROR: unknown menu code" << item << '\n';
                return;
        }
        player.refresh();
    }
}

static int buildMenuFiberSelect()
{
    static int menuID = 0;
    if ( menuID == 0 )
        menuID = gym::createMenu(processMenuFiberSelect);
    else
        gym::clearMenu(menuID);
    
    gym::addMenuEntry("Hide All", 1);
    FiberDisp * FD = player.firstFiberDisp();
    if ( FD )
    {
        gym::addMenuEntry(FD->hide&1?"Show right pointing":"Hide right pointing", 2);
        gym::addMenuEntry(FD->hide&2?"Show left pointing":"Hide left pointing", 3);
    }
    return menuID;
}


//------------------------------------------------------------------------------
static void processMenuCoupleSelect(int item)
{
    switch (item)
    {
        case 0: return;
        case 1: disp.couple_select  = 0; break;
        case 2: disp.couple_select ^= 1; break;
        case 3: disp.couple_select ^= 2; break;
        case 4: disp.couple_select ^= 4; break;
        default:
            std::cerr << "CYTOSIM ERROR: unknown menu code" << item << '\n';
            return;
    }
    player.refresh();
}

static int buildMenuCoupleSelect()
{
    static int menuID = 0;
    if ( menuID == 0 )
        menuID = gym::createMenu(processMenuCoupleSelect);
    else
        gym::clearMenu(menuID);
    
    gym::addMenuEntry("Hide All", 1);
    gym::addMenuEntry(disp.couple_select&1?"Hide Free":"Show Free",     2);
    gym::addMenuEntry(disp.couple_select&2?"Hide Bound":"Show Bound",   3);
    gym::addMenuEntry(disp.couple_select&4?"Hide Links":"Show Links", 4);
    return menuID;
}

//------------------------------------------------------------------------------
static void processMenuSingleSelect(int item)
{
    switch (item)
    {
        case 0: return;
        case 1: disp.single_select  = 0; break;
        case 2: disp.single_select ^= 1; break;
        case 3: disp.single_select ^= 2; break;
        
        default:
            std::cerr << "CYTOSIM ERROR: unknown menu code" << item << '\n';
            return;
    }
    player.refresh();
}

static int buildMenuSingleSelect()
{
    static int menuID = 0;
    if ( menuID == 0 )
        menuID = gym::createMenu(processMenuSingleSelect);
    else
        gym::clearMenu(menuID);
    
    gym::addMenuEntry("Hide All",     1);
    gym::addMenuEntry(disp.single_select&1?"Hide Free":"Show Free",     2);
    gym::addMenuEntry(disp.single_select&2?"Hide Bound":"Show Bounds", 3);
    return menuID;
}


//------------------------------------------------------------------------------

static int buildMenuSelect()
{
    static int menuID = 0;
    int m1 = buildMenuFiberSelect();
    int m2 = buildMenuCoupleSelect();
    int m3 = buildMenuSingleSelect();
    
    if ( menuID == 0 )
        menuID = gym::createMenu(nullptr);
    else
        gym::clearMenu(menuID);

    gym::addSubMenu("Fibers",  m1);
    gym::addSubMenu("Couple",  m2);
    gym::addSubMenu("Singles", m3);
    return menuID;
}


//------------------------------------------------------------------------------
#pragma mark -

static void processMenuAnimation(int item)
{
    switch (item)
    {
        case 0: return;
        case 1: processKey('z'); return;
        case 2: processKey('a'); return;
        case 4: processKey('s'); return;
        case 5: processKey('r'); return;
        default:
            std::cerr << "CYTOSIM ERROR: unknown menu code" << item << '\n';
            return;
    }
}

static int buildMenuAnimation()
{
    static int menuID = 0;
    
    if ( menuID == 0 )
    {
        menuID = gym::createMenu(processMenuAnimation);
        gym::addMenuEntry("(z) Reset State",      1);
        gym::addMenuEntry("(a) Start Live",       2);
        gym::addMenuEntry("(s) One Step & Stop",  4);
        gym::addMenuEntry("(r) Read Parameters",  5);
    }
    return menuID;
}


//------------------------------------------------------------------------------

static void processMenuReplay(int item)
{
    switch (item)
    {
        case 0: return;
        case 1: processKey('p'); return;
        case 2: processKey('o'); return;
        case 3: processKey('s'); return;
        case 4: processKey('z'); return;
        case 5: player.previousFrame(); break;
        case 6: player.nextFrame(); break;
        case 7: prop.loop = 0; return;
        case 8: prop.loop = 1; return;
        default:
            std::cerr << "CYTOSIM ERROR: unknown menu code" << item << '\n';
            return;
    }
    player.refresh();
}

static int buildMenuReplay()
{
    static int menuID = 0;
    
    if ( menuID == 0 )
    {
        menuID = gym::createMenu(processMenuReplay);
        gym::addMenuEntry("(p) Play/Faster",    1);
        gym::addMenuEntry("(o) Slower",         2);
        gym::addMenuEntry("(s) Stop",           3);
        gym::addMenuEntry("-",                  0);
        gym::addMenuEntry("(z) First Frame",    4);
        gym::addMenuEntry("(<) Previous Frame", 5);
        gym::addMenuEntry("(>) Next Frame",     6);
        if ( prop.loop )
            gym::addMenuEntry("Do not loop", 7);
        else
            gym::addMenuEntry("Loop movie", 8);
    }
    return menuID;
}


//------------------------------------------------------------------------------
static void processMenuExport(int item)
{
    switch (item)
    {
        case 0: return;
        case 1: player.saveView(glApp::currentView(), prop.image_index++, 1); return;
        case 2: player.saveView(glApp::currentView(), prop.image_index++, 2); return;
        case 3: player.saveScene(3, "image", prop.image_index++, 3); return;
        case 4: player.saveScene(6, "image", prop.image_index++, 3); return;
        case 5: player.saveScene(6, "image", prop.image_index++, 2); return;
        case 6: player.saveScene(4, "poster", prop.poster_index++); return;
        case 7: player.saveScene(8, "poster", prop.poster_index++); return;
        case 8: player.saveScene(16, "poster", prop.poster_index++); return;

        case 9: prop.save_images = 9999; player.startPlayback(); return;
        case 10: prop.image_index = 0; return;
        
        case 20: player.writePlayParameters(std::cout, true); return;
        case 21: player.writeDisplayParameters(std::cout, true); return;
        case 22: worker.writeProperties(std::cout, true); return;
        case 23: worker.exportObjects(false); return;
        case 24: worker.exportObjects(true); return;
            
        default:
            std::cerr << "CYTOSIM ERROR: unknown menu code" << item << '\n';
            return;
    }
}


static int buildMenuExport()
{
    static int menuID = 0;
    if ( menuID == 0 )
        menuID = gym::createMenu(processMenuExport);
    else
        gym::clearMenu(menuID);
    
    gym::addMenuEntry("Save Image (y)",            1);
    gym::addMenuEntry("Save 2x Downsampled Image", 2);
    gym::addMenuEntry("Save Fine Image",           3);
    gym::addMenuEntry("Save 2x Fine Image",        4);
    gym::addMenuEntry("Save 3x Fine Image",        5);
    gym::addMenuEntry("Save 4x Poster",            6);
    gym::addMenuEntry("Save 8x Poster",            7);
    gym::addMenuEntry("Save 16x Poster",           8);
    gym::addMenuEntry("Play & Save Images (Y)",    9);
    gym::addMenuEntry("Reset Image-file Index",   10);
    gym::addMenuEntry("-",                         0);
    gym::addMenuEntry("Write Play Parameters",    20);
    gym::addMenuEntry("Write Display Parameters", 21);
    gym::addMenuEntry("Write Object Properties",  22);
    gym::addMenuEntry("Export Objects",           23);
    gym::addMenuEntry("Export Objects as Binary", 24);
    
    return menuID;
}

//------------------------------------------------------------------------------
//                    MAIN MENU
//------------------------------------------------------------------------------
#pragma mark -

void processTopMenu(int item)
{
    if ( item == 9 )
        exit(EXIT_SUCCESS);
}


int rebuildMenus()
{
    static int menuID = 0;
    int m1 = buildMenuDisplay();
    int m2 = buildMenuSelect();
    int m3 = buildMenuAnimation();
    int m4 = buildMenuReplay();
    int m5 = buildMenuExport();
    int m6 = glApp::buildMenu();

    if ( menuID == 0 )
        menuID = gym::createMenu(processTopMenu);
    else
        gym::clearMenu(menuID);
    
    gym::addSubMenu("Display",           m1);
    gym::addSubMenu("Object-Selection",  m2);
    gym::addSubMenu("Live-Simulation",   m3);
    gym::addSubMenu("File-Replay",       m4);
    gym::addSubMenu("Export",            m5);
    gym::addSubMenu("More",              m6);
    gym::addMenuEntry("Quit",             9);
    return menuID;
}


void menuCallback(int status, int x, int y)
{
    //printf("menu status(%i, %i, %i)\n", status, x, y);
    
    if ( status == GLUT_MENU_NOT_IN_USE )
        rebuildMenus();
}

void buildMenus()
{
    glutMenuStatusFunc(menuCallback);
    glApp::attachMenu(rebuildMenus());
}
