// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University

/**
 It is possible to save images or many images by pressing a key..
 disabling this capacity in public distribution might be safer:
 It can be done by disabling ENABLE_WRITE below:
 */
#define ENABLE_WRITE 1


/// defines the increment for user's modifications (eg. key pressed)
static float grained(float x, int inc)
{
    const float grain = 0.25f;
    float dx = inc * ( 1 + ( x >= 4 ) + 2 * ( x >= 8 ) + 4 * ( x >= 16 ) );
    float nx = std::nearbyint( x / grain + dx );
    float ii = std::abs(inc);
    return grain * std::max(ii, nx);
}

//------------------------------------------------------------------------------

template< typename T >
static void setVisible(T* p, int val)
{
    p->visible = val;
    flashText("%s:visible = %i", p->name_str(), val);
}

static void flipVisible(PointDisp* p, int val)
{
    p->visible = ( p->visible != val ) * val;
    flashText("%s:visible = %i", p->name_str(), p->visible);
}

static void changeStyle(PointDisp * p, int)
{
    p->style = ( p->style + 1 ) % 8;
    flashText("%s:style = %i", p->name_str(), p->style);
}

static void setSize(PointDisp * p, float s)
{
    if ( s >= 0.5 )
    {
        p->size = s;
        flashText("%s:size = %.2f", p->name_str(), s);
    }
}

static void setWidth(PointDisp * p, float s)
{
    if ( s > 0.5 )
    {
        p->width = s;
        flashText("%s:width = %.2f", p->name_str(), s);
    }
}

static void changeSize(PointDisp * p, int inc)
{
    float s = grained(p->size, inc);
    if ( s > 32.f )
        s = 0.5f;
    if ( s > 0 )
    {
        float w = p->width;
        p->size = s;
        p->width *= s / w;
        flashText("%s:size = %.2f", p->name_str(), s);
    }
}

//------------------------------------------------------------------------------
#pragma mark - PointDisp lists

static inline PointDisp* toPointDisp(Property * ptr)
{
    return static_cast<PointDisp*>(ptr);
}


/// apply function to all PointDisp is plist
static void setPointDisp(PropertyList const& plist, void(*func)(PointDisp*, int), int val)
{
    if ( plist.empty() )
        flashText("no relevant object");

    for ( Property * i : plist )
        func(toPointDisp(i), val);
}


static void changePointDispSize(PropertyList const& plist, int inc, bool dos, bool dow)
{
    for ( Property * i : plist )
    {
        PointDisp * p = toPointDisp(i);
        if ( dos ) p->size = grained(p->size, inc*2);
        if ( dow ) p->width = grained(p->width, inc);
    }
    
    if ( disp.style == 1 || plist.size() > 1 )
    {
        if ( dow ) {
            disp.link_width = grained(disp.link_width, inc);
            flashText("simul:link_width %.2f", disp.link_width);
        }
        if ( dos ) {
            disp.point_size = grained(disp.point_size, inc*2);
            flashText("simul:point_size %.2f", disp.point_size);
        }
    }
    else if ( plist.size() == 1 )
    {
        PointDisp * p = toPointDisp(plist.front());
        if ( dow ) flashText("%s:width %.2f", p->name_str(), p->width);
        if ( dos ) flashText("%s:size %.2f", p->name_str(), p->size);
    }
}


static void setPointDispVisible(PropertyList const& plist, int val)
{
    if ( plist.empty() )
        flashText("no relevant object");
    
    for ( Property * i : plist )
        toPointDisp(i)->visible = val;
}


static PointDisp * nextVisiblePointDisp(PropertyList const& plist, size_t& cnt)
{
    PointDisp* one = nullptr;
    cnt = 0;
    // find first one which is visible:
    for ( PropertyList::const_iterator i = plist.begin(); i < plist.end(); ++i )
    {
        PointDisp * dsp = toPointDisp(*i);
        if ( dsp->visible )
        {
            ++cnt;
            // choose follower:
            if ( i+1 < plist.end() )
                one = toPointDisp(*(i+1));
            else
                one = nullptr;
        }
    }
    return one;
}


static void shufflePointDispVisible(const PropertyList& plist, int val)
{
    if ( plist.empty() )
        flashText("no relevant object");

    if ( plist.size() == 1 )
    {
        flipVisible(toPointDisp(plist.front()), val);
    }
    else
    {
        size_t cnt = 0;
        PointDisp * p = nextVisiblePointDisp(plist, cnt);

        if ( cnt > 1 )
        {
            setPointDispVisible(plist, 0);
            p = toPointDisp(plist.front());
            p->visible = val;
            flashText("Only `%s' is visible", p->name_str());
        }
        else if ( p != nullptr )
        {
            setPointDispVisible(plist, 0);
            p->visible = val;
            flashText("Only `%s' is visible", p->name_str());
        }
        else if ( cnt == 1 )
        {
            setPointDispVisible(plist, 0);
            flashText("All hidden");
        }
        else
        {
            setPointDispVisible(plist, val);
            flashText("All visible");
        }
    }
}

//------------------------------------------------------------------------------
#pragma mark - Single Couple


static void changeSingleSelect()
{
    unsigned int & select = disp.single_select;
    switch( select )
    {
        case 3:  select = 0; flashText("single:select=0: hidden");      break;
        case 0:  select = 2; flashText("single:select=2: bound only");  break;
        case 2:  select = 1; flashText("single:select=1: free only");   break;
        default: select = 3; flashText("single:select=3: all");         break;
    }
}


static void changeCoupleSelect()
{
    unsigned int & select = disp.couple_select;
    switch( select )
    {
        case 7:  select = 0; flashText("couple:select=0: hidden");        break;
        case 0:  select = 2; flashText("couple:select=2: bound only");    break;
        case 2:  select = 4; flashText("couple:select=4: bridging only"); break;
        case 4:  select = 1; flashText("couple:select=1: free only");     break;
        default: select = 7; flashText("couple:select=7: all");           break;
    }
}

static void changeCoupleSelect2()
{
    unsigned int & select = disp.couple_select;
    if ( select & 8 )
    {
        select = 16+4;
        flashText("couple: parallel bridging only");
    }
    else if ( select & 16 )
    {
        select = 7;
        flashText("couple: all");
    }
    else
    {
        select = 8+4;
        flashText("couple: antiparallel bridging only");
    }
}

//---------------------------------------------------------------------
#pragma mark - Fibers


static void changeExclude(FiberDisp* p, int val)
{
    if ( val )
        p->hide >>= 2;
    p->hide = ( p->hide + 1 ) & 3;
    if ( val )
        p->hide <<= 2;
    
    switch ( p->hide )
    {
        case 0: flashText("All fibers");                break;
        case 1: flashText("Right-pointing fibers");     break;
        case 2: flashText("Left-pointing fibers");      break;
        case 3: flashText("No fibers");                 break;
        case 4: flashText("Counter-clockwise fibers");  break;
        case 8: flashText("Clockwise fibers");          break;
        case 12: flashText("No fibers");                break;
    }
}


static void flipExplode(FiberDisp* p)
{
    p->explode_style = ! p->explode_style;
    if ( p->explode_style && p->explode_range == 0 )
        p->explode_range = 1;
    flashText("fiber:explode = %i", p->explode_style);
}


static void changeScale(real& scale, int d)
{
    real s = std::log2(std::fabs(scale)) + d * 0.125;
    if ( s < -14 ) s =  10;
    if ( s >  10 ) s = -14;
    scale = std::copysign(std::exp2(s), scale);
}


static void changeScale(FiberDisp* p, int d)
{
    if ( p->lattice_style )
    {
        changeScale(p->lattice_scale, d);
        flashText("fiber:lattice_scale = %.5f", p->lattice_scale);
    }
    else if ( p->line_style == 2 || p->line_style == 3 )
    {
        changeScale(p->tension_scale, d);
        flashText("fiber:tension_scale = %.5f", p->tension_scale);
    }
    else if ( p->force_style )
    {
        changeScale(p->force_scale, d);
        flashText("fiber:force_scale = %.5f", p->force_scale);
    }
    else if ( p->speckle_style )
    {
        changeScale(p->speckle_gap, d);
        flashText("fiber:speckle_gap = %.5f", p->speckle_gap);
    }
    else if ( p->point_style == 3 )
    {
        changeScale(p->point_gap, d);
        flashText("fiber:point_gap = %.5f", p->point_gap);
    }
    else if ( p->line_style == 4 || p->line_style == 6 || p->line_style == 7 || p->line_style == 8 )
    {
        changeScale(p->length_scale, d);
        flashText("fiber:length_scale = %.5f", p->length_scale);
    }
    else if ( disp.style == 2 )
        flipExplode(p);
}


static void invertScale(FiberDisp* p, int)
{
    if ( p->lattice_style )
    {
        p->lattice_scale = -p->lattice_scale;
        flashText("fiber:lattice_scale = %.5f", p->lattice_scale);
    }
    else if ( p->line_style == 2 || p->line_style == 3 )
    {
        p->tension_scale = -p->tension_scale;
        if ( p->tension_scale > 0 )
            flashText("fiber:tension_scale: pulling");
        else
            flashText("fiber:tension_scale: pushing");
    }
    else if ( p->line_style == 4 || p->line_style == 6 || p->line_style == 7 || p->line_style == 8 )
    {
        p->length_scale = -p->length_scale;
        flashText("fiber:length_scale = %.5f", p->length_scale);
    }
}


static void flashColoring(int val)
{
    switch ( val )
    {
        case FiberDisp::COLORING_OFF:       flashText("Fibers: no coloring");         break;
        case FiberDisp::COLORING_RANDOM:    flashText("Fibers randomly colored");     break;
        case FiberDisp::COLORING_DIRECTION: flashText("Fibers colored by direction"); break;
        case FiberDisp::COLORING_MARK:      flashText("Fibers colored by mark");      break;
        case FiberDisp::COLORING_FLAG:      flashText("Fibers colored by flag");      break;
        case FiberDisp::COLORING_FAMILY:    flashText("Fibers colored by family");    break;
        case FiberDisp::COLORING_CLUSTER:   flashText("Fibers colored by cluster");   break;
        case FiberDisp::COLORING_AGE:       flashText("Fibers colored by age");       break;
        case FiberDisp::COLORING_PSTATE:    flashText("Fibers colored by +end state");break;
        default: flashText("unknown fiber:coloring mode"); break;
    }
}


static void setColoring(FiberDisp* p, int val)
{
    p->coloring = ( p->coloring ? 0 : val );
    flashColoring(p->coloring);
}


static void changeColoring(FiberDisp* p, int inc)
{
    p->coloring = ( p->coloring + inc + 9 ) % 9;
    flashColoring(p->coloring);
}


static void setMask(FiberDisp* p, int val)
{
    p->mask = val;
    p->mask_bitfield = RNG.distributed_bits(p->mask);
    flashText("fiber:mask_bitfield=0x%X (%i bits)", p->mask_bitfield, p->mask);
}

static void changeMask(FiberDisp* p, int val)
{
    p->mask = (( p->mask << val ) & 31 ) + ( p->mask == 0 );
    p->mask_bitfield = RNG.distributed_bits(p->mask);
    flashText("fiber:mask_bitfield=0x%X (%i bits)", p->mask_bitfield, p->mask);
}

static void changePointStyle(FiberDisp* p, int arg)
{
    p->point_style = ( p->point_style + 1 ) % arg;
    switch ( p->point_style )
    {
        case 0: flashText("Fibers: no points"); break;
        case 1: flashText("Fibers: vertices"); break;
        case 2: flashText("Fibers: arrowheads"); break;
        case 3: flashText("Fibers: chevrons"); break;
        case 4: flashText("Fibers: center points"); break;
        default: flashText("unknown fiber:point_style"); break;
    }
}


static void flashLineStyle(int val)
{
    switch ( val )
    {
        case 0: flashText("Fibers: no lines"); break;
        case 1: flashText("Fibers: lines"); break;
        case 2: flashText("Fiber color by axial tensions"); break;
        case 3: flashText("Fiber jet color by axial tensions"); break;
        case 4: flashText("Fiber color by curvature"); break;
        case 5: flashText("Fiber color by orientation"); break;
        case 6: flashText("Fiber highlight near minus-end"); break;
        case 7: flashText("Fiber highlight near plus-end"); break;
        case 8: flashText("Fiber color by height"); break;
        case 9: flashText("Fiber color by grid (if style=3)"); break;
        default: flashText("unknown fiber:line style"); break;
    }
}

static void flashFiberStyle(int val)
{
    switch ( val )
    {
        case 0: flashText("Fibers: default style"); break;
        case 1: flashText("Fibers: style=filament"); break;
        case 2: flashText("Fibers: style=actin"); break;
        case 3: flashText("Fibers: style=microtubule"); break;
        case 4: flashText("Fibers: style=backbone"); break;
        default: flashText("unknown fiber:style"); break;
    }
}

static void toggleLineStyle(FiberDisp* p, int val)
{
    p->line_style = ( p->line_style != val ) * val;
    if ( p->style == 4 ) // override the 'backbone' style
        p->style = 0;
    flashLineStyle(p->line_style);
}

static void changeLineStyle(FiberDisp* p, int inc)
{
    p->line_style = ( p->line_style + inc ) % 9;
    flashLineStyle(p->line_style);
}

static void toggleStyle(FiberDisp* p, int val)
{
    p->style = ( p->style != val ) * val;
    flashFiberStyle(p->style);
}


static void changeSpeckleStyle(FiberDisp* p, int)
{
    p->speckle_style = ( p->speckle_style + 1 ) % 3;
    switch ( p->speckle_style )
    {
        case 0: flashText("Fibers: no speckles");       break;
        case 1: flashText("Fibers: random speckles");   break;
        case 2: flashText("Fibers: regular speckles");  break;
    }
}


static void changeSpeckleSize(FiberDisp* p, int inc)
{
    p->speckle_size = grained(p->speckle_size, inc);
    flashText("%s:speckle_size=%0.2f", p->name_str(), p->speckle_size);
}


static void changeLatticeStyle(FiberDisp* p, int)
{
#if FIBER_HAS_LATTICE || FIBER_HAS_MESH
    p->lattice_style = ( 1 + p->lattice_style ) % 5;
    flashText("Fibers: lattice_style=%i", p->lattice_style);
#else
    flashText("Warning: cannot display fiber:lattice");
#endif
}


static void changePointSize(FiberDisp* p, int inc)
{
    if ( p->speckle_style ) changeSpeckleSize(p, inc);
    p->point_style = grained(p->point_style, inc);
    flashText("%s:point_size=%0.2f", p->name_str(), p->point_style);
}


static void changeLineWidth(FiberDisp* p, int inc)
{
    p->line_width = grained(p->line_width, inc);
    flashText("%s:line_width=%0.2f", p->name_str(), p->line_width);
}


static void changeEndStyle(FiberDisp* p, int val)
{
    const int P = 1+val;
    const int M = 1+val*2;
    int * style = p->end_style;
    // showing the plus ends -> the minus ends -> both -> none
    switch( bool(style[1]) + 2*bool(style[0]) )
    {
        case 0:
            style[0] = P;
            style[1] = 0;
            break;
        case 1:
            style[0] = 0;
            style[1] = 0;
            break;
        case 2:
            style[0] = P;
            style[1] = M;
            break;
        case 3:
        default:
            style[0] = 0;
            style[1] = M;
            break;
    }
    
    switch( (style[0]?1:0) + (style[1]?2:0) )
    {
        case 0: flashText("Fibers: no ends");    break;
        case 1: flashText("Fibers: plus-ends");  break;
        case 2: flashText("Fibers: minus-ends"); break;
        case 3: flashText("Fibers: both ends");  break;
    }
}


static void changeEndSize(FiberDisp* p, int inc)
{
    float* size = p->end_size;
    if ( p->end_style[0] && p->end_style[1] )
    {
        size[0] = grained(size[0], inc);
        size[1] = grained(size[1], inc);
        flashText("%s:end_size %.2f %.2f", p->name_str(), size[0], size[1]);
    }
    else if ( p->end_style[0] )
    {
        size[0] = grained(size[0], inc);
        flashText("%s::plus_end %.2f", p->name_str(), size[0]);
    }
    else if ( p->end_style[1] )
    {
        size[1] = grained(size[1], inc);
        flashText("%s::minus_end %.2f", p->name_str(), size[1]);
    }
    else
        changePointSize(p, inc);
}


/// change the size of all features that are visible
static void changeSize(FiberDisp* p, int inc)
{
    if ( p->line_style ) changeLineWidth(p, inc);
    if ( p->point_style ) changePointSize(p, inc);
    if ( p->speckle_style ) changeSpeckleSize(p, inc);
    changeEndSize(p, inc);
}

//---------------------------------------------------------------------
#pragma mark - FiberDisp lists

static inline FiberDisp* toFiberDisp(Property * ptr)
{
    return static_cast<FiberDisp*>(ptr);
}

static void setFiberDisp(PropertyList const& plist, void(*func)(FiberDisp*, int), int val)
{
    for ( Property * i : plist )
        func(toFiberDisp(i), val);
}

static PointDisp * findVisibleFiberDisp(PropertyList const& plist, int& cnt)
{
    PointDisp * one = nullptr;
    cnt = 0;
    // find first one which is visible:
    for ( Property * i : plist )
    {
        if ( toFiberDisp(i)->visible )
        {
            ++cnt;
            if ( !one )
                one = toPointDisp(i);
        }
    }
    return one;
}


static void setFiberDispVisible(PropertyList const& plist, int val)
{
    for ( Property * i : plist )
        toFiberDisp(i)->visible = val;
}


static FiberDisp * nextVisibleFiberDisp(PropertyList const& plist, size_t& cnt)
{
    FiberDisp* one = nullptr;
    cnt = 0;
    // find first one which is visible:
    for ( PropertyList::const_iterator i = plist.begin(); i < plist.end(); ++i )
    {
        FiberDisp * dsp = toFiberDisp(*i);
        if ( dsp->visible )
        {
            ++cnt;
            // choose follower:
            if ( i+1 < plist.end() )
                one = toFiberDisp(*(i+1));
            else
                one = nullptr;
        }
    }
    return one;
}


static void shuffleVisible(FiberDisp* p)
{
    if ( p->visible && p->line_style )
    {
        p->visible = 0;
        flashText("%s:visible = %i", p->name_str(), p->visible);
    }
    else if ( p->visible )
    {
        p->line_style = 1;
        p->speckle_style = 0;
        flashText("%s:line_style = %i", p->name_str(), p->line_style);
    }
    else
    {
        p->visible = 1;
        p->line_style = 0;
        p->speckle_style = 1;
        flashText("%s:speckle_style = %i", p->name_str(), p->speckle_style);
    }
}


static void shuffleFiberDispVisible(const PropertyList& plist, int val)
{
    if ( plist.empty() )
        flashText("no relevant fibers");
    
    if ( plist.size() == 1 )
    {
        shuffleVisible(toFiberDisp(plist.front()));
    }
    else
    {
        size_t cnt = 0;
        FiberDisp * p = nextVisibleFiberDisp(plist, cnt);
        
        if ( cnt > 1 )
        {
            setFiberDispVisible(plist, 0);
            p = toFiberDisp(plist.front());
            p->visible = val;
            flashText("Only `%s' is visible", p->name_str());
        }
        else if ( p != nullptr )
        {
            setFiberDispVisible(plist, 0);
            p->visible = val;
            flashText("Only `%s' is visible", p->name_str());
        }
        else if ( cnt == 1 )
        {
            setFiberDispVisible(plist, 0);
            flashText("All hidden");
        }
        else
        {
            setFiberDispVisible(plist, val);
            flashText("All visible");
        }
    }
}

//------------------------------------------------------------------------------
//---------------------------- keyboard commands -------------------------------
//------------------------------------------------------------------------------
#pragma mark - Keyboard Commands

/// provide minimal on-screen summary of the most important key combinations
void helpKeys(std::ostream& os)
{
    os << "                          Keyboard Commands\n";
    os << "\n";
    os << "   SPACE       Start-stop animation or replay\n";
    os << "   < >         Show previous; show next frame ( , . also works)\n";
    os << "   O s o p     Play reverse; stop; play slower; play faster\n";
    os << "   z           Rewind to first frame / Restart live simulation\n";
    os << "   ALT-SPACE   Reset view (i.e. zoom, translation, rotation)\n";
    os << "   f F         Toggle full-screen mode; maximize window size\n";
    os << "   i v b       Invert colors; toggle slice view; toggle scale bar\n";
    os << "   l k         Read parameter file; Print display parameters\n";
    os << "   r R         Report various informations on display window\n";
#if ENABLE_WRITE
    os << "   y Y         Save current image; Play and save all images\n";
#endif
    os << "\nSimulation\n";
    os << "   a s         Start live mode; Allow one simulation step\n";
    os << "   A a         Double period (num. steps/display); reset period\n";
    os << "   g G         Delete mouse-controlled handles; release handle\n";
    os << "\nFibers\n";
    os << "   `           Address another type of fibers for modifications\n";
    os << "   1           Change display: line / color-coded tension / hide\n";
    os << "   / /         Toggle vertices; toggle backbone display\n";
    os << "   2 3         Decrease; increase line width (ALT: point size)\n";
    os << "   !           Change display of tips: off / plus / both / minus\n";
    os << "   @ #         Decrease; increase fiber_end display size\n";
    os << "   4 $         Change lattice display style; change speckle display\n";
    os << "   c d         Toggle fiber coloring; hide Right/left-pointing\n";
    os << "   m M         Mask a fraction of the fibers; change mask value\n";
    os << "   w e         decrease/increase tension/lattice/explode scale\n";
    os << "   t T         Toggle auto-tracking: 't':nematic; 'T':polar mode\n";
    os << "\nBeads - Solids - Spheres\n";
    os << "   5           Rotate various bead/sphere display style\n";
    os << "   %           Change point size\n";
    os << "\nSingles - Couples\n";
    os << "   6           Change Single visibility based on state\n";
    os << "   7 ALT-7     Change Couple visibility based on state\n";
    os << "   0           Change visibility flags of Hands\n";
    os << "   8 9         Decrease; Increase point size of visible Hands\n";
    os << "   ALT-8 ALT-9 Decrease; Increase line width of visible Hands\n";
    os << "\nSpaces\n";
    os << "   u           Rotate visibility\n";
}


void processKey(unsigned char key, int modifiers = 0)
{
    //std::cerr << "processing key `" << key << "'\n";
    // the view associated with the current window
    View & view = glApp::currentView();
    
    const bool altKeyDown = modifiers & GLUT_ACTIVE_ALT;
    const bool shiftKeyDown = modifiers & GLUT_ACTIVE_SHIFT;
    /*
     In the switch below:
     - use break if the display need to be refreshed,
     - otherwise, use return.
    */
    switch (key)
    {
        case 'h':
            view.draw_memo = ( view.draw_memo + 1 ) % 6;
            break;
        
        case 'N':
            /**Need to share OpenGL context with the main window */
            //glApp::newWindow(drawLive);
            break;

#if ENABLE_WRITE
            
        case 'y':
            // save current image, without decorations
            player.drawScene(glApp::currentView());
            player.saveView(prop.image_index++, 1);
            // with over sampling and downsampling to get super-resolution:
            //player.saveScene(3, "image", prop.image_index++, 3);
            return;
            
        case 'Y':
            // start player to save all images in file
            if ( prop.save_images == 0 )
            {
                if ( player.startPlayback() || worker.alive() )
                    prop.save_images = 9999;
            }
            else
            {
                prop.save_images = 0;
            }
            break;

#endif

        //------------------------- Global controls ----------------------------

        case 'k':
        {
            if ( altKeyDown )
                worker.writeProperties(std::cout, true);
            else
            {
                std::cout << '\n';
                player.writePlayParameters(std::cout, true);
                player.writeDisplayParameters(std::cout, true);
            }
        } break;
            
#if DRAW_MECA_LINKS
        case 'K':
            disp.draw_links = !disp.draw_links;
            flashText("draw_links = %i", disp.draw_links);
            break;
#endif

        case 'l': {
            try {
                std::string file = simul.prop.config_file;
                worker.reloadParameters(file);
                flashText("Reloaded %s", file.c_str());
            }
            catch( Exception & e ) {
                flashText("Error in config: %s", e.what());
            }
        } break;
            
        case 'z':
            if ( worker.goodFile() )
                player.rewind();
            else
                player.restart();
            break;
            
        case 'Z':
            worker.cancel_join();
            player.restart();
            break;
            
        case 'a':
            if ( altKeyDown )
            {
                player.setStyle(1);
                flashText("Style 1");
            }
            else
            {
                prop.period = 1;
                worker.period(prop.period);
                flashText("period = 1");
                player.extendLive();
            }
            break;
            
        case 'A':
            if ( worker.alive() )
            {
                prop.period = 2 * prop.period;
                if ( prop.period > 1024 ) prop.period = 1;
                worker.period(prop.period);
                flashText("period = %i", prop.period);
            }
            else
            {
                // this will initialize the simulation engine without making a step
                worker.prolong_run(0);
            }
            break;
            
        case 's':
            if ( altKeyDown )
            {
                player.setStyle(2);
                flashText("Style 2");
            }
            else
            {
                if ( worker.holding() )
                {
                    glApp::displayAll();
                    worker.signal();
                }
                else if ( worker.alone() )
                    player.nextFrame();
                player.stop();
            }
            break;
            
        case 'S':
            prop.period = 1;
            worker.period(prop.period);
            flashText("period = 1");
            break;
            
        case 'g':
            worker.deleteHandles();
            flashText("Deleted mouse-controled handles");
            break;
            
        case 'G':
            worker.releaseHandle();
            break;
            
        case 'i':
            glApp::currentView().invertColors();
            break;

        case 'r':
            prop.toggleReport(altKeyDown?-1:1);
            break;
            
        case 'R':
            prop.toggleReport(0);
            break;

        //------------------------- play / stop / reverse ----------------------
            
        case '<':
        case ',':
            if ( prop.replay == 1 )
                player.stop();
            else
                player.previousFrame();
            break;
            
        case '>':
        case '.':
            if ( prop.replay == -1 )
                player.stop();
            else
                player.nextFrame();
            break;
            
        case 'o':
            if ( prop.delay < 1 << 13 )
                prop.delay *= 2;
            flashText("Delay %i ms", prop.delay);
            return;
            
        case 'O':
            if ( !player.startBackward() )
                player.accelerate();
            return;
            
        case 'p':
            if ( !player.startPlayback() )
                player.accelerate();
            return;
        
        case ' ':
            if ( altKeyDown )
            {
                view.clearPixels();
                view.reset();
                flashText("");
            }
            else
                player.startstop();
            return;
            
        //------------------------------ Fibers --------------------------------
           
        case '`':
            shuffleFiberDispVisible(player.allFiberDisp(), 1);
            break;
            
        case '~':
            shuffleFiberDispVisible(player.allFiberDisp(), 1);
            break;

        case 't':
            view.track_fibers ^= 1;
            flashText("view.track_fibers = %i (translation)", view.track_fibers);
            break;
            
        case 'T':
            view.track_fibers ^= 3;
            flashText("view.track_fibers = %i (rotation)", view.track_fibers);
            break;
            
        case 'd':
            if ( altKeyDown )
            {
                player.setStyle(3);
                flashText("Style 3");
            }
            else
            {
                setFiberDisp(player.allVisibleFiberDisp(), changeExclude, 0);
            }
            break;
            
        case 'D':
            setFiberDisp(player.allVisibleFiberDisp(), changeExclude, 1);
            break;
                
        case 'q':
            setFiberDisp(player.allVisibleFiberDisp(), invertScale, 0);
            break;

        case 'w':
            setFiberDisp(player.allVisibleFiberDisp(), changeScale, -8);
            break;

        case 'e':
            setFiberDisp(player.allVisibleFiberDisp(), changeScale, 8);
            break;
            
        case 'W':
            setFiberDisp(player.allVisibleFiberDisp(), changeScale, -1);
            break;

        case 'E':
            setFiberDisp(player.allVisibleFiberDisp(), changeScale, 1);
            break;

        case 'm':
            if ( altKeyDown )
                setFiberDisp(player.allVisibleFiberDisp(), setMask, 0);
            else
                setFiberDisp(player.allVisibleFiberDisp(), changeMask, 1);
            break;

        case 'M':
            setFiberDisp(player.allVisibleFiberDisp(), changeMask, 0);
            break;
            
        case 'c':
            setFiberDisp(player.allVisibleFiberDisp(), changeColoring, 1);
            break;
                
        case 'C':
            setFiberDisp(player.allVisibleFiberDisp(), setColoring, 1);
            break;

        case 167:
            setFiberDisp(player.allVisibleFiberDisp(), toggleLineStyle, 1);
            break;
            
        case '?':
        case 177:
            setFiberDisp(player.allVisibleFiberDisp(), toggleStyle, 4);
            break;

        case '1':
            if ( altKeyDown )
                setFiberDisp(player.allVisibleFiberDisp(), changePointStyle, 5);
            else
                setFiberDisp(player.allVisibleFiberDisp(), changeLineStyle, 1);
            break;
            
        case '!':
            setFiberDisp(player.allVisibleFiberDisp(), changeEndStyle, !altKeyDown);
            break;
            
        case '2':
            if ( altKeyDown)
                setFiberDisp(player.allVisibleFiberDisp(), changePointSize, -1);
            else if ( shiftKeyDown )
                setFiberDisp(player.allVisibleFiberDisp(), changeLineWidth, -1);
            else
                setFiberDisp(player.allVisibleFiberDisp(), changeSize, -1);
            break;

        case '@':
            setFiberDisp(player.allVisibleFiberDisp(), changeEndSize, -1);
            break;

        case '3':
            if ( altKeyDown)
                setFiberDisp(player.allVisibleFiberDisp(), changePointSize, 1);
            else if ( shiftKeyDown )
                setFiberDisp(player.allVisibleFiberDisp(), changeLineWidth, 1);
            else
                setFiberDisp(player.allVisibleFiberDisp(), changeSize, 1);
            break;
            
        case '#':
            setFiberDisp(player.allVisibleFiberDisp(), changeEndSize, 1);
            break;
            
        case '$':
            setFiberDisp(player.allVisibleFiberDisp(), changeSpeckleStyle, 0);
            break;
            
        case '4':
            setFiberDisp(player.allVisibleFiberDisp(), changeLatticeStyle, 0);
            break;
            
        case '/':
            setFiberDisp(player.allVisibleFiberDisp(), changePointStyle, 2);
            break;

        case '*':
            setFiberDisp(player.allVisibleFiberDisp(), changeLatticeStyle, 0);
            break;

        //------------------------ Solid, Bead & Sphere ------------------------
  
        case '5':
            if ( altKeyDown )
                shufflePointDispVisible(player.allSphereDisp(), 1);
            else
                setPointDisp(player.allVisibleSphereDisp(), changeStyle, 0);
            break;
        
        case '%':
            setPointDisp(player.allVisibleSphereDisp(), changeSize, 2);
            break;
            
        //------------------------ Single/Couple + Hands -----------------------
           
        case '6':
            changeSingleSelect();
            break;
            
        case '&':
            shufflePointDispVisible(player.allSpaceDisp(), 3);
            break;
            
        case 'u':
            shufflePointDispVisible(player.allSpaceDisp(), 3);
            break;

        case 'U':
            shufflePointDispVisible(player.allSpaceDisp(), 1);
            break;

        case '7':
            if ( altKeyDown )
                changeCoupleSelect2();
            else
                changeCoupleSelect();
            break;
            
        case '8':
            changePointDispSize(player.allVisibleHandDisp(), -1, !altKeyDown, !shiftKeyDown);
            break;
            
        case '9':
            changePointDispSize(player.allVisibleHandDisp(), +1, !altKeyDown, !shiftKeyDown);
            break;

        case '0':
            if ( altKeyDown )
                setPointDispVisible(player.allHandDisp(), 1);
            else
                shufflePointDispVisible(player.allHandDisp(), 1);
            break;
        
#if 0
        case 185: //that is the key left of '=' on the numpad
        case '=':
            break;
        case '-':
            break;
        case '+':
            break;
#endif

        default:
            // other keys are handles by glApp
            glApp::processNormalKey(key, 0, 0);
            return;
    }
    
    // if break was called, redraw the scene:
    glApp::displayMain();
}


void processNormalKey(const unsigned char key, const int x, const int y)
{
    // check for user-defined `magic_key`
    for ( int k = 0; k < PlayerProp::NB_MAGIC_KEYS; ++k )
    {
        if ( key == prop.magic_key[k] )
        {
            worker.evaluate(prop.magic_code[k]);
            flashText("%s", prop.magic_code[k].c_str());
            return;
        }
    }
    
    processKey(key, glutGetModifiers());
}
