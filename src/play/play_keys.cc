// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University

/**
 It is possible to save images or many images by pressing a key..
 disabling this capacity in public distribution might be safer:
 It can be done by disabling ENABLE_WRITE below:
 */
#define PLAY_CAN_WRITE 1


/// change `x` by adding `inc` while keeping discrete values
static float grained(float x, int inc)
{
    const float grain = 0.125f;
    // for large values, the increment is larger:
    float g = 1 + ( x >= 2 ) + 2 * ( x >= 4 ) + 4 * ( x >= 8 ) + 8 * ( x >= 16 );
    float n = std::nearbyint( x / grain + inc * g );
    float i = std::abs(inc);
    return grain * std::max(i, n);
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

static void changeColoring(PointDisp* p, int)
{
    p->coloring = ( p->coloring + 1 ) % 4;
    flashText("%s:coloring = %i", p->name_str(), p->coloring);
}

[[maybe_unused]]
static void setSize(PointDisp * p, float s)
{
    if ( s >= 0.5 )
    {
        p->size = s;
        flashText("%s:size = %.2f", p->name_str(), s);
    }
}

[[maybe_unused]]
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
    
    // change global values for style 1:
    if ( disp.style == 1 || plist.size() > 1 )
    {
        if ( dos ) disp.point_size = grained(disp.point_size, inc*2);
        if ( dow ) disp.link_width = grained(disp.link_width, inc);
        if ( dos & dow ) flashText("simul:point_size %.2f link_width %.2f", disp.point_size, disp.link_width);
        else if ( dos ) flashText("simul:point_size %.2f", disp.point_size);
        else if ( dow ) flashText("simul:link_width %.2f", disp.link_width);
    }
    else if ( plist.size() > 0 )
    {
        PointDisp * p = toPointDisp(plist.front());
        char const* n = p->name_str();
        if ( dos & dow ) flashText("%s:size %.2f width %.2f", n, p->size, p->width);
        else if ( dow ) flashText("%s:width %.2f", n, p->width);
        else if ( dos ) flashText("%s:size %.2f", n, p->size);
    }
}


static void setPointDispVisible(PropertyList const& plist, int val)
{
    if ( plist.empty() )
    {
        flashText("no relevant object");
        return;
    }
    
    for ( Property * i : plist )
        toPointDisp(i)->visible = val;

    std::string cat = plist.front()->category();

    if ( val )
        flashText("All "+cat+" visible");
    else
        flashText("No "+cat+" visible");
}


/* this set `cnt` number of PointDisp for which visible != 0. */
static PointDisp * nextVisiblePointDisp(PropertyList const& plist, size_t& cnt)
{
    cnt = 0;
    PointDisp * pick = nullptr;
    for ( Property * i : plist )
    {
        PointDisp * d = toPointDisp(i);
        //std::clog << i->name() << ".visible=" << d->visible << " (" << cnt << ")\n";
        if ( cnt == 1 && !pick )
            pick = d;
        cnt += ( 0 != d->visible );
    }
    //std::clog << " --> pick " << (pick?pick->name():"null") << "\n";
    return pick;
}


static void shufflePointDispVisible(const PropertyList& plist, int val)
{
    if ( plist.empty() )
        return flashText("no relevant object");

    if ( plist.size() == 1 )
    {
        flipVisible(toPointDisp(plist.front()), val);
    }
    else
    {
        size_t cnt = 0;
        PointDisp * p = nextVisiblePointDisp(plist, cnt);
        if ( cnt == 0 )
        {
            setPointDispVisible(plist, val);
        }
        else
        {
            setPointDispVisible(plist, 0);
            if ( cnt > 1 )
                p = toPointDisp(plist.front());
            if ( p )
            {
                p->visible = val;
                flashText("Only `%s' is visible", p->name_str());
            }
        }
    }
}

//------------------------------------------------------------------------------
#pragma mark - Single Couple


static void shuffleSingleSelect()
{
    unsigned int & select = disp.single_select;
    switch( select )
    {
        case 3:  select = 2; flashText("single:select=2: bound only");  break;
        case 2:  select = 1; flashText("single:select=1: free only");   break;
        case 1:  select = 0; flashText("single:select=0: hidden");      break;
        default: select = 3; flashText("single:select=3: all");         break;
    }
}


static void shuffleCoupleSelect()
{
    unsigned int & select = disp.couple_select;
    switch( select )
    {
        case 7:  select = 4; flashText("couple:select=4: bridging only"); break;
        case 4:  select = 2; flashText("couple:select=2: bound only");    break;
        case 2:  select = 1; flashText("couple:select=1: free only");     break;
        case 1:  select = 0; flashText("couple:select=0: hidden");        break;
        default: select = 7; flashText("couple:select=7: all");           break;
    }
}

static void shuffleCoupleSelect2()
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
#pragma mark - Fibers flash

static void flashColoring(FiberDisp const* p)
{
    char const* n = p->name_str();
    switch ( p->coloring )
    {
        case FiberDisp::COLORING_OFF:       flashText("%s: no coloring", n);      break;
        case FiberDisp::COLORING_RANDOM:    flashText("%s: randomly colored", n);     break;
        case FiberDisp::COLORING_DIRECTION: flashText("%s: colored by direction", n); break;
        case FiberDisp::COLORING_MARK:      flashText("%s: colored by mark", n);      break;
        case FiberDisp::COLORING_FLAG:      flashText("%s: colored by flag", n);      break;
        case FiberDisp::COLORING_FAMILY:    flashText("%s: colored by family", n);    break;
        case FiberDisp::COLORING_CLUSTER:   flashText("%s: colored by cluster", n);   break;
        case FiberDisp::COLORING_AGE:       flashText("%s: colored by age", n);       break;
        case FiberDisp::COLORING_PSTATE:    flashText("%s: colored by +end state", n);break;
        default: flashText("unknown %s:coloring mode", n); break;
    }
}


static void flashPointStyle(FiberDisp const* p)
{
    char const* n = p->name_str();
    switch ( p->point_style )
    {
        case 0: flashText("%s: no points", n); break;
        case 1: flashText("%s: vertices", n); break;
        case 2: flashText("%s: arrowheads", n); break;
        case 3: flashText("%s: chevrons", n); break;
        case 4: flashText("%s: center points", n); break;
        default: flashText("unknown %s:point_style", n); break;
    }
}

static void flashFiberStyle(FiberDisp const* p)
{
    char const* n = p->name_str();
    switch ( p->style )
    {
        case 0: flashText("%s: style=0 (default)", n); break;
        case 1: flashText("%s: style=1 (backbone)", n); break;
        case 2: flashText("%s: style=2 (striped)", n); break;
        case 3: flashText("%s: style=3 (filament)", n); break;
        case 4: flashText("%s: style=4 (actin)", n); break;
        case 5: flashText("%s: style=5 (microtubule)", n); break;
        default: flashText("unknown %s:style", n); break;
    }
}

static void flashLineStyle(FiberDisp const* p)
{
    char const* n = p->name_str();
    if ( p->style )
        flashFiberStyle(p);
    else switch ( p->line_style )
    {
        case 0: flashText("%s: no lines", n); break;
        case 1: flashText("%s: lines", n); break;
        case 2:
            if ( p->tension_scale > 0 )
                flashText("%s: color by axial tensions (pulling)", n);
            else
                flashText("%s: color by axial tensions (pushing)", n);
            break;
        case 3:
            if ( p->tension_scale > 0 )
                flashText("%s: jet color by axial tensions (pulling)", n);
            else
                flashText("%s: jet color by axial tensions (pushing)", n);
            break;
        case 4: flashText("%s: color by orientation", n); break;
        case 5: flashText("%s: color by curvature", n); break;
        case 6: flashText("%s: gradient from minus-end", n); break;
        case 7: flashText("%s: gradient from plus-end", n); break;
        case 8: flashText("%s: gradient from growing plus-end", n); break;
        case 9: flashText("%s: color by height", n); break;
        case 10: flashText("%s: color by grid (if style=3)", n); break;
        default: flashText("unknown %s:line style", n); break;
    }
}

static void flashHide(FiberDisp const* p)
{
    char const* n = p->name_str();
    switch ( p->hide )
    {
        case 0: flashText("%s: Right and left pointing", n); break;
        case 1: flashText("%s: Right-pointing only", n);     break;
        case 2: flashText("%s: Left-pointing only", n);      break;
        case 3: flashText("%s: all hidden", n);              break;
        case 4: flashText("%s: Counter-clockwise only", n);  break;
        case 8: flashText("%s: Clockwise fibers only", n);   break;
        case 12: flashText("%s: all hidden", n);             break;
    }
}

static void flashTracking(unsigned mode)
{
    switch ( mode )
    {
        case 0: flashText("no tracking"); break;
        case 1: flashText("tracking fiber position"); break;
        case 2: flashText("tracking fiber nematic direction"); break;
        case 3: flashText("tracking fiber position & nematic"); break;
        case 4: flashText("tracking fiber position spread"); break;
        case 5: flashText("tracking fiber position & spread"); break;
        case 8: flashText("tracking solids position"); break;
        case 10: flashText("tracking solids position & nematic direction"); break;
    }
}

//---------------------------------------------------------------------
#pragma mark - Fibers

static void changeHide(FiberDisp* p, int val)
{
    if ( val )
        p->hide >>= 2;
    p->hide = ( p->hide + 1 ) & 3;
    if ( val )
        p->hide <<= 2;
    flashHide(p);
}


static void changeMarked(FiberDisp* p, int val)
{
    unsigned &n = p->show_marked;
    ++n;
    if ( n == 4 )
    {
        n = ~0U;
        flashText("%s: showing all marks", p->name_str());
    }
    else
        flashText("%s: showing only mark=%i", p->name_str(), n);
}


static void flipExplode(DisplayProp& p)
{
    p.explode_style = ! p.explode_style;
    if ( p.explode_style && p.explode_range == 0 )
        p.explode_range = 1;
    flashText("display:explode = %i", p.explode_style);
}


template < typename REAL >
static void changeScale(REAL& scale, int d)
{
    REAL s = std::log2(std::fabs(scale)) + d * 0.125;
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
    else if ( p->line_style > 4 )
    {
        changeScale(p->length_scale, d);
        flashText("fiber:length_scale = %.5f", p->length_scale);
    }
    else if ( disp.style == 1 && disp.explode_style )
    {
        changeScale(disp.explode_range, d);
        flashText("fiber:explode_range = %.5f", disp.explode_range);
    }
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
            flashText("fiber:tension_scale > 0: pulling");
        else
            flashText("fiber:tension_scale < 0: pushing");
    }
    else if ( p->line_style > 4 )
    {
        p->length_scale = -p->length_scale;
        flashText("fiber:length_scale = %.5f", p->length_scale);
    }
}

[[maybe_unused]]
static void setColoring(FiberDisp* p, int val)
{
    p->coloring = ( p->coloring ? 0 : val );
    flashColoring(p);
}

static void changeColoring(FiberDisp* p, int inc)
{
    p->coloring = ( p->coloring + inc + 9 ) % 9;
    flashColoring(p);
}

static void setMask(FiberDisp* p, int val)
{
    p->mask = (unsigned)val;
    p->mask_bitfield = distribute_bits(p->mask, pcg32_state);
    flashText("fiber:mask_bitfield=0x%X (%i bits)", p->mask_bitfield, p->mask);
}

static void changeMask(FiberDisp* p, int val)
{
    p->mask = (( p->mask << val ) & 31 ) + ( p->mask == 0 );
    p->mask_bitfield = distribute_bits(p->mask, pcg32_state);
    flashText("fiber:mask_bitfield=0x%X (%i bits)", p->mask_bitfield, p->mask);
}

static void changePointStyle(FiberDisp* p, int arg)
{
    p->point_style = ( p->point_style + 1 ) % arg;
    flashPointStyle(p);
}

static void toggleLineStyle(FiberDisp* p, int val)
{
    if ( p->style )
    {
        p->style = 0;
        flashFiberStyle(p);
    }
    else
    {
        p->line_style = ( p->line_style != val ) * val;
        flashLineStyle(p);
    }
}

static void changeLineStyle(FiberDisp* p, int inc)
{
    p->line_style = ( p->line_style + inc ) % 10;
    flashLineStyle(p);
}

static void toggleFiberStyle(FiberDisp* p, int inc)
{
    p->style = ( p->style + inc ) % 6;
    flashFiberStyle(p);
}


static void changeSpeckleStyle(FiberDisp* p, int)
{
    p->speckle_style = ( p->speckle_style + 1 ) % 3;
    char const* n = p->name_str();
    switch ( p->speckle_style )
    {
        case 0: flashText("%s: no speckles", n);       break;
        case 1: flashText("%s: random speckles", n);   break;
        case 2: flashText("%s: regular speckles", n);  break;
    }
}


static void changeSpeckleSize(FiberDisp* p, int inc)
{
    p->speckle_size = grained(p->speckle_size, inc);
    flashText("%s:speckle_size=%0.2f", p->name_str(), p->speckle_size);
}


static void changeLatticeStyle(FiberDisp* p, int)
{
#if FIBER_HAS_LATTICE || FIBER_HAS_DENSITY
    p->lattice_style = ( 1 + p->lattice_style ) % 5;
    flashText("%s: lattice_style=%i", p->name_str(), p->lattice_style);
#else
    flashText("Warning: no fiber:lattice support");
#endif
}


static void changePointSize(FiberDisp* p, int inc)
{
    if ( p->speckle_style )
        changeSpeckleSize(p, inc);
    else if ( p->point_style )
    {
        p->point_size = grained(p->point_size, inc);
        flashText("%s:point_size=%0.2f", p->name_str(), p->point_size);
    }
}


static void changeLineWidth(FiberDisp* p, int inc)
{
    p->line_width = grained(p->line_width, inc);
    flashText("%s:line_width=%0.2f", p->name_str(), p->line_width);
}


static void changeEndStyle(FiberDisp* d, int val)
{
    const int P = 1+val;
    const int M = 1+val*2;
    int * style = d->end_style;
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
    
    char const* n = d->name_str();
    switch( (style[0]?1:0) + (style[1]?2:0) )
    {
        case 0: flashText("%s: no ends", n);    break;
        case 1: flashText("%s: plus-ends", n);  break;
        case 2: flashText("%s: minus-ends", n); break;
        case 3: flashText("%s: both ends", n);  break;
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
}


/// change the size of all features that are visible
static void changeSize(FiberDisp* p, int inc)
{
    if ( p->line_style ) changeLineWidth(p, inc);
    if ( p->point_style ) changePointSize(p, inc);
    if ( p->speckle_style ) changeSpeckleSize(p, inc);
    if ( p->style == 1 )
    {
        p->bone_width = grained(p->bone_width, inc);
        flashText("display:bone_width=%0.2f", p->bone_width);
    }
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

[[maybe_unused]]
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
    if ( plist.size() < 1 )
        flashText("no fiber visible!");
    for ( Property * i : plist )
    {
        FiberDisp * d = toFiberDisp(i);
        if ( d )
        {
            d->visible = val;
            if ( val && !d->speckle_style )
                d->line_style = 1;
        }
    }
}


static void shuffleFiberStyle(FiberDisp* p, int val)
{
    char const* n = p->name_str();
    if ( val && p->line_style )
    {
        p->line_style = 0;
        p->speckle_style = 1;
        flashText("%s:speckle_style = %i", n, p->speckle_style);
    }
    else if ( p->speckle_style )
    {
        p->speckle_style = 0;
        p->style = 1;
        flashText("%s:style = %i, width = %.2f", n, p->style, p->bone_width);
    }
    else
    {
        p->style = 0;
        p->line_style = 1;
        flashText("%s:line_style = %i, width = %.2f", n, p->line_style, p->line_width);
    }
}


static FiberDisp * nextVisibleFiberDisp(PropertyList const& plist, size_t& cnt)
{
    cnt = 0;
    FiberDisp* pick = nullptr;
    for ( Property * i : plist )
    {
        FiberDisp * d = toFiberDisp(i);
        //std::clog << i->name() << ".visible=" << d->visible << " (" << cnt << ")\n";
        if ( cnt == 1 && !pick )
            pick = d;
        cnt += ( 0 != d->visible );
    }
    //std::clog << " --> pick " << (pick?pick->name():"null") << "\n";
    return pick;
}


static void shuffleFiberDispVisible(const PropertyList& plist, int val)
{
    if ( plist.size() == 1 )
    {
        for ( Property * i : plist )
            shuffleFiberStyle(toFiberDisp(i), 1);
    }
    else
    {
        size_t cnt = 0;
        FiberDisp * p = nextVisibleFiberDisp(plist, cnt);
        if ( cnt == 0 )
        {
            setFiberDispVisible(plist, val);
            if ( val )
                flashText("All fibers visible");
            else
                flashText("No fiber visible");
        }
        else
        {
            setFiberDispVisible(plist, 0);
            if ( cnt > 1 )
                p = toFiberDisp(plist.front());
            if ( p )
            {
                p->visible = val;
                p->line_style = val;
                flashText("Only `%s' is visible", p->name_str());
            }
            else
                flashText("No fiber visible");
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
    os << "   l L         Read parameter file; Print display parameters\n";
    os << "   r R         Report various informations on display window\n";
#if PLAY_CAN_WRITE
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
    os << "   * (         Decrease; Increase line width of visible Hands\n";
    os << "\nSpaces\n";
    os << "   u           Rotate visibility\n";
}


void processKey(unsigned char key, int modifiers = 0)
{
    // the view associated with the current window
    View & view = glApp::currentView();
    
    const bool altKeyDown = modifiers & GLUT_ACTIVE_ALT;
    const bool shiftKeyDown = modifiers & GLUT_ACTIVE_SHIFT;
    //std::cerr<<"processKey("<<key<<") SHIFT "<<shiftKeyDown<<" ALT "<<altKeyDown<<"\n";
    
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
            //glApp::newWindow(drawSimulation);
            break;

#if PLAY_CAN_WRITE
        case 'y': {
            // save current image, without decorations
            player.drawSystem(view);
            player.saveView(view, prop.image_index++, 1);
            // with over sampling and downsampling to get super-resolution:
            //player.saveScene(3, "image", prop.image_index++, 3);
        } return;
            
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

        case 'L':
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
        case 'k':
            setFiberDisp(player.allVisibleFiberDisp(), changeMarked, 0);
            break;

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
       
#if ENABLE_EXPLODED_DISPLAY
        case 'X':
            flipExplode(disp);
            break;
#endif
            
        case 'z':
            if ( worker.goodFile() )
                player.rewind();
            else
                player.restart(0);
            break;
            
        case 'Z':
            worker.cancel_join();
            player.restart(1);
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
                    worker.signal();
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
            
        case 'P':
            player.startPlayback();
            player.setTimelapse(1);
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
            setFiberDisp(player.allVisibleFiberDisp(), shuffleFiberStyle, 0);
            break;

        case 't':
            if ( altKeyDown )
                view.track_fibers ^= 3;
            else
                view.track_fibers ^= 1;
            flashTracking(view.track_fibers);
            break;
            
        case 'T':
            view.track_fibers = ( view.track_fibers ? 0 : 26 );
            flashTracking(view.track_fibers);
            break;
            
        case 'd':
            if ( altKeyDown )
            {
                player.setStyle(3);
                flashText("Style 3");
            }
            else
            {
                setFiberDisp(player.allVisibleFiberDisp(), changeHide, 0);
            }
            break;
            
        case 'D':
            setFiberDisp(player.allVisibleFiberDisp(), changeHide, 1);
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
            setPointDisp(player.allVisibleSphereDisp(), changeColoring, 1);
            break;

        case 167:
            setFiberDisp(player.allVisibleFiberDisp(), toggleLineStyle, 1);
            break;
            
        case '?':
        case 177:
            setFiberDisp(player.allVisibleFiberDisp(), toggleFiberStyle, 1);
            break;

        case '1':
            if ( altKeyDown || shiftKeyDown )
                setFiberDisp(player.allVisibleFiberDisp(), changePointStyle, 5);
            else
                setFiberDisp(player.allVisibleFiberDisp(), changeLineStyle, 1);
            break;
            
        case '!':
            setFiberDisp(player.allVisibleFiberDisp(), changeEndStyle, !altKeyDown);
            break;
            
        case '2':
            if ( altKeyDown || shiftKeyDown )
                setFiberDisp(player.allVisibleFiberDisp(), changePointSize, -1);
            else
                setFiberDisp(player.allVisibleFiberDisp(), changeSize, -1);
            break;

        case '@':
            setFiberDisp(player.allVisibleFiberDisp(), changeEndSize, -1);
            break;

        case '3':
            if ( altKeyDown || shiftKeyDown )
                setFiberDisp(player.allVisibleFiberDisp(), changePointSize, 1);
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
            
        //------------------------ Solid, Bead & Sphere ------------------------
  
        case '5':
            if ( altKeyDown || shiftKeyDown )
                shufflePointDispVisible(player.allSphereDisp(), shiftKeyDown);
            else
                setPointDisp(player.allVisibleSphereDisp(), changeStyle, 0);
            break;
        
        case '%':
            setPointDisp(player.allVisibleSphereDisp(), changeSize, 2);
            break;
            
        //------------------------ Single/Couple + Hands -----------------------
           
        case '6':
            shuffleSingleSelect();
            break;
            
        case 'u':
            if ( altKeyDown )
                shufflePointDispVisible(player.allSpaceDisp(), 1);
            else
                shufflePointDispVisible(player.allSpaceDisp(), 3);
            break;

        case 'U':
            shufflePointDispVisible(player.allSpaceDisp(), 2);
            break;

        case '7':
            if ( altKeyDown )
                shuffleCoupleSelect2();
            else
                shuffleCoupleSelect();
            break;
            
        case '8':
            if ( shiftKeyDown )
                changePointDispSize(player.allVisibleHandDisp(), -1, 0, 1);
            else
                changePointDispSize(player.allVisibleHandDisp(), -1, 1, 1);
            break;
            
        case '*':
            changePointDispSize(player.allVisibleHandDisp(), -1, 0, 1);
            break;
            
        case '9':
            if ( shiftKeyDown )
                changePointDispSize(player.allVisibleHandDisp(), +1, 0, 1);
            else
                changePointDispSize(player.allVisibleHandDisp(), +1, 1, 1);
            break;

        case '(':
            changePointDispSize(player.allVisibleHandDisp(), +1, 0, 1);
            break;

        case '0':
            shufflePointDispVisible(player.allHandDisp(), 1);
            break;
        
        case ')': // flip between hiding/showing all hands
            setPointDispVisible(player.allHandDisp(), player.allVisibleHandDisp().empty());
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
            // other keys are passed-on to glApp
            glApp::processNormalKey(key, 0, 0);
            return;
    }
    
    // if break was called, redraw the scene:
    simul.fresh_ = 1;
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
