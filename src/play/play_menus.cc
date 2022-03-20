// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University

/// shortcut
static int createMenu(void (*func)(int)) { return glutCreateMenu(func); }

/// shortcut
static void addMenuEntry(char const* str, int val) { glutAddMenuEntry(str, val); }

/// shortcut
static void addSubMenu(char const* str, int val) { glutAddSubMenu(str, val); }



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
                
            case 7: FD->end_style[1] = 3*!FD->end_style[1]; break;
            case 8: FD->end_style[0] = 2*!FD->end_style[0]; break;
                
            case 9: FD->force_style = FD->force_style; break;
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
        glApp::postRedisplay();
    }
}


static int buildMenuFiber()
{
    static int menuID = 0;
    if ( menuID == 0 )
        menuID = createMenu(processMenuFiber);
    else
        glApp::clearMenu(menuID);
    
    FiberDisp * FD = player.firstFiberDisp();
    if ( FD )
    {
        addMenuEntry(FD->visible        ? "Hide"             :"Show",             10);
        addMenuEntry(FD->line_style     ? "Hide Lines"       :"Show Lines",        1);
        addMenuEntry(FD->line_style==2  ? "Hide Tensions"    :"Show Tensions",     2);
        addMenuEntry(FD->point_style    ? "Hide Points"      :"Show Points",       3);
        addMenuEntry(FD->point_style==2 ? "Hide Arrows"      :"Show Arrows",       5);
        addMenuEntry(FD->end_style[1]   ? "Hide Minus-ends"  :"Show Minus-end",    7);
        addMenuEntry(FD->end_style[0]   ? "Hide Plus-ends"   :"Show Plus-end",     8);
        addMenuEntry(FD->force_style    ? "Hide Point-forces":"Show Point-Forces", 9);
        addMenuEntry("No coloring",           20);
        addMenuEntry("Coloring by number",    21);
        addMenuEntry("Coloring by mark",      22);
        addMenuEntry("Coloring by flag",      23);
        addMenuEntry("Coloring by family",    24);
        addMenuEntry("Coloring by cluster",   25);
        addMenuEntry("Coloring by direction", 26);
        addMenuEntry("Coloring by age",       27);
        addMenuEntry("Coloring by +end state",28);
        addMenuEntry("draw_average=0", 30);
        addMenuEntry("draw_average=1", 31);
        addMenuEntry("draw_average=2", 32);
    }
    else
        addMenuEntry("no fiber?", 0);

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
    glApp::postRedisplay();
}

static int buildMenuCouple()
{
    static int menuID = 0;
    if ( menuID == 0 )
        menuID = createMenu(processMenuCouple);
    else
        glApp::clearMenu(menuID);
    
    addMenuEntry("Hide all",    1);
    addMenuEntry("Show free",   2);
    addMenuEntry("Show bound",  3);
    addMenuEntry("Show links",  4);
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
    glApp::postRedisplay();
}


static int buildMenuStyle()
{
    static int menuID = 0;
    if ( menuID == 0 )
    {
        menuID = createMenu(processMenuDisplay);
        addMenuEntry("Detailed (style 1)", 101);
        addMenuEntry("Fastest (style 2)", 102);
        addMenuEntry("Best Looking (style 3)", 103);
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
        menuID = createMenu(processMenuDisplay);
    else
        glApp::clearMenu(menuID);
    
    addMenuEntry("Reset View",  1);
    addSubMenu("Style",   m0);
    addSubMenu("Fibers",  m1);
    addSubMenu("Couple",  m2);
    
    View & view = glApp::currentView();
    addMenuEntry("Toggle fullscreen mode (f)", 4);
    addMenuEntry(disp.tile?"Non-tiled Display":"Tiled Display", 3);
    addMenuEntry(view.track_fibers?"stop tracking":"Track Fibers", 6);
    
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
        glApp::postRedisplay();
    }
}

static int buildMenuFiberSelect()
{
    static int menuID = 0;
    if ( menuID == 0 )
        menuID = createMenu(processMenuFiberSelect);
    else
        glApp::clearMenu(menuID);
    
    addMenuEntry("Hide All", 1);
    FiberDisp * FD = player.firstFiberDisp();
    if ( FD )
    {
        addMenuEntry(FD->hide&1?"Show right pointing":"Hide right pointing", 2);
        addMenuEntry(FD->hide&2?"Show left pointing":"Hide left pointing", 3);
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
    glApp::postRedisplay();
}

static int buildMenuCoupleSelect()
{
    static int menuID = 0;
    if ( menuID == 0 )
        menuID = createMenu(processMenuCoupleSelect);
    else
        glApp::clearMenu(menuID);
    
    addMenuEntry("Hide All", 1);
    addMenuEntry(disp.couple_select&1?"Hide Free":"Show Free",     2);
    addMenuEntry(disp.couple_select&2?"Hide Bound":"Show Bound",   3);
    addMenuEntry(disp.couple_select&4?"Hide Links":"Show Links", 4);
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
    glApp::postRedisplay();
}

static int buildMenuSingleSelect()
{
    static int menuID = 0;
    if ( menuID == 0 )
        menuID = createMenu(processMenuSingleSelect);
    else
        glApp::clearMenu(menuID);
    
    addMenuEntry("Hide All",     1);
    addMenuEntry(disp.single_select&1?"Hide Free":"Show Free",     2);
    addMenuEntry(disp.single_select&2?"Hide Bound":"Show Bounds", 3);
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
        menuID = createMenu(nullptr);
    else
        glApp::clearMenu(menuID);

    addSubMenu("Fibers",  m1);
    addSubMenu("Couple",  m2);
    addSubMenu("Singles", m3);
    
    return menuID;
}


//------------------------------------------------------------------------------
#pragma mark -

static void processMenuAnimation(int item)
{
    switch (item)
    {
        case 0: return;
        case 1: processKey('z'); break;
        case 2: processKey('a'); break;
        case 4: processKey('s'); break;
        case 5: processKey('r'); break;
        default:
            std::cerr << "CYTOSIM ERROR: unknown menu code" << item << '\n';
            return;
    }
    glApp::postRedisplay();
}

static int buildMenuAnimation()
{
    static int menuID = 0;
    
    if ( menuID == 0 )
    {
        menuID = createMenu(processMenuAnimation);
        addMenuEntry("(z) Reset State",      1);
        addMenuEntry("(a) Start Live",       2);
        addMenuEntry("(s) One Step & Stop",  4);
        addMenuEntry("(r) Read Parameters",  5);
    }
    return menuID;
}


//------------------------------------------------------------------------------

static void processMenuReplay(int item)
{
    switch (item)
    {
        case 0: return;
        case 1: processKey('p'); break;
        case 2: processKey('o'); break;
        case 3: processKey('s'); break;
        case 4: processKey('z'); break;
        case 5: player.previousFrame(); break;
        case 6: player.nextFrame(); break;
        case 7: prop.loop = 0; break;
        case 8: prop.loop = 1; break;
        default:
            std::cerr << "CYTOSIM ERROR: unknown menu code" << item << '\n';
            return;
    }
    glApp::postRedisplay();
}

static int buildMenuReplay()
{
    static int menuID = 0;
    
    if ( menuID == 0 )
    {
        menuID = createMenu(processMenuReplay);
        addMenuEntry("(p) Play/Faster",    1);
        addMenuEntry("(o) Slower",         2);
        addMenuEntry("(s) Stop",           3);
        addMenuEntry("-",                  0);
        addMenuEntry("(z) First Frame",    4);
        addMenuEntry("(<) Previous Frame", 5);
        addMenuEntry("(>) Next Frame",     6);
        if ( prop.loop )
            addMenuEntry("Do not loop", 7);
        else
            addMenuEntry("Loop movie", 8);
    }
    return menuID;
}


//------------------------------------------------------------------------------
static void processMenuExport(int item)
{
    switch (item)
    {
        case 0: return;
        case 1: player.saveView(prop.image_index++, 1); return;
        case 2: player.saveView(prop.image_index++, 2); return;
        case 3: player.saveScene(3, "image", prop.image_index++, 3); return;
        case 4: player.saveScene(6, "image", prop.image_index++, 3); return;
        case 5: player.saveScene(9, "image", prop.image_index++, 3); return;
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
    glApp::postRedisplay();
}


static int buildMenuExport()
{
    static int menuID = 0;
    if ( menuID == 0 )
        menuID = createMenu(processMenuExport);
    else
        glApp::clearMenu(menuID);
    
    addMenuEntry("Save Image (y)",            1);
    addMenuEntry("Save Downsampled Image",    2);
    addMenuEntry("Save Fine Image",           3);
    addMenuEntry("Save 2x Fine Image",        4);
    addMenuEntry("Save 3x Fine Image",        5);
    addMenuEntry("Save 4x Poster",            6);
    addMenuEntry("Save 8x Poster",            7);
    addMenuEntry("Save 16x Poster",           8);
    addMenuEntry("Play & Save Images (Y)",    9);
    addMenuEntry("Reset Image-file Index",   10);
    addMenuEntry("-",                         0);
    addMenuEntry("Write Play Parameters",    20);
    addMenuEntry("Write Display Parameters", 21);
    addMenuEntry("Write Object Properties",  22);
    addMenuEntry("Export Objects",           23);
    addMenuEntry("Export Objects as Binary", 24);
    
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


void buildMenus()
{
    static int menuID = 0;
    int m1 = buildMenuDisplay();
    int m2 = buildMenuSelect();
    int m3 = buildMenuAnimation();
    int m4 = buildMenuReplay();
    int m5 = buildMenuExport();
    int m6 = glApp::buildMenu();

    if ( menuID == 0 )
        menuID = createMenu(processTopMenu);
    else
        glApp::clearMenu(menuID);
    
    addSubMenu("Display",           m1);
    addSubMenu("Object-Selection",  m2);
    addSubMenu("Live-Simulation",   m3);
    addSubMenu("File-Replay",       m4);
    addSubMenu("Export",            m5);
    addSubMenu("More",              m6);
    addMenuEntry("Quit",             9);
}


void menuCallback(int status, int x, int y)
{
    //printf("menu status(%i, %i, %i)\n", status, x, y);
    
    if ( GLUT_MENU_NOT_IN_USE )
        buildMenus();
}

