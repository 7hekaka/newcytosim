// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University

//------------------------------------------------------------------------------
#pragma mark - Color schemes for Fibers

inline gle_color color_fiber(Fiber const& fib, size_t)
{
    return fib.disp->color;
}

inline gle_color color_by_tension(Fiber const& fib, size_t seg)
{
    real x = fib.prop->disp->tension_alpha * fib.tension(seg);
    return fib.disp->color.alpha(x);
}

inline gle_color color_by_tension_jet(Fiber const& fib, size_t seg)
{
    real x = fib.prop->disp->tension_alpha * fib.tension(seg);
    // use jet coloring, where Lagrange multipliers are negative under compression
    return gle_color::jet_color_alpha(x);
    //return fib.disp->color.alpha(x-1);
}

inline gle_color color_by_curvature(Fiber const& fib, size_t pti)
{
    real beta = fib.prop->disp->length_scale;
    return gle_color::jet_color(beta*fib.curvature(pti));
}

inline gle_color color_seg_curvature(Fiber const& fib, size_t seg)
{
    real beta = 0.5 * fib.prop->disp->length_scale;
    real c = fib.curvature(seg);
    real d = fib.curvature(seg+1);
    return gle_color::jet_color(beta*(c+d));
}

inline gle_color color_by_direction(Fiber const& fib, size_t seg)
{
    return gle::radial_color(fib.dirSegment(seg));
}

/// using distance from the minus end to the start of segment `seg`
inline gle_color color_by_distanceM(Fiber const& fib, real pti)
{
    real beta = fib.segmentation() / fib.prop->disp->length_scale;
    real x = std::min(pti*beta, (real)32.0);
    return fib.disp->color.alpha(std::exp(-x));
}

/// using the distance from the plus end to vertex `pti`
inline gle_color color_by_distanceP(Fiber const& fib, real pti)
{
    real beta = fib.segmentation() / fib.prop->disp->length_scale;
    // using the distance at the vertex
    real x = std::min((fib.lastPoint()-pti)*beta, (real)32.0);
    return fib.disp->color.alpha(std::exp(-x));
}

/// using distance from the plus end to the end of segment `seg`
inline gle_color color_by_abscissaM(Fiber const& fib, size_t seg)
{
    real beta = fib.segmentation() / fib.prop->disp->length_scale;
    real x = std::min(seg*beta, (real)32.0);
    return fib.disp->color.alpha(std::exp(-x));
}

/// using distance from the plus end to the end of segment `seg`
inline gle_color color_by_abscissaP(Fiber const& fib, size_t seg)
{
    real beta = fib.segmentation() / fib.prop->disp->length_scale;
    real x = std::min((fib.lastSegment()-seg)*beta, (real)32.0);
    return fib.disp->color.alpha(std::exp(-x));
}

/// color set according to distance to the confining Space
inline gle_color color_by_height(Fiber const& fib, size_t pti)
{
    real beta = 1 / fib.prop->disp->length_scale;
    real Z = 0;
    Space const* spc = fib.prop->confine_space_ptr;
    if ( spc )
        Z = -spc->signedDistanceToEdge(fib.posPoint(pti));
#if ( DIM > 2 )
    else
        Z = fib.posPoint(pti).ZZ;
#endif
    return gle_color::jet_color(Z*beta);
}

/// color set according to steric grid
inline gle_color color_by_grid(Fiber const& fib, size_t seg)
{
    Map<DIM> const& map = fib.simul().locusGrid.map();
    Vector w = fib.posPoint(seg, 0.5);
    size_t i = map.index(w);
    return gle::alt_color(i);
}


//------------------------------------------------------------------------------
#pragma mark - Color schemes for Fiber Lattice

inline gle_color color_alternate(Fiber const& fib, long ix, real)
{
    if ( ix & 1 )
        return fib.disp->color;
    else
        return fib.disp->color.darken(0.75);
}

inline gle_color color_by_lattice(Fiber const& fib, long ix, real beta)
{
    return fib.disp->color.darken(beta*fib.visibleLattice()->data(ix));
}

inline gle_color color_by_lattice_jet(Fiber const& fib, long ix, real beta)
{
    return gle_color::jet_color(beta*fib.visibleLattice()->data(ix));
}

inline gle_color color_by_lattice_striped(Fiber const& fib, long ix, real beta)
{
    real x = beta * fib.visibleLattice()->data(ix);
    return gle_color::jet_color(x).lighten(0.9375 + ( ix & 1 ) * 0.125);
}

inline gle_color lattice_color(gle_color const& col, real val)
{
    return col.darken(val);
    //return gle_color::jet_color(val);
}

//------------------------------------------------------------------------------
#pragma mark - Color schemes for Sphere / Solid / Bead

/**
 returns the front color used to display Object
 if `coloring` is enabled, this loads the N-th bright color,
 otherwise load the object's display color
 */
template < typename T >
gle_color bodyColorF(T const& obj)
{
    const PointDisp * disp = obj.prop->disp;
    if ( disp->coloring )
    {
        size_t i = ( disp->coloring == 2 ? obj.mark() : obj.signature());
        return gle::bright_color(i).match_a(disp->color);
    }
    return disp->color;
}

/**
 returns the front color used to display Object
 if `coloring` is enabled, this loads the N-th bright color,
 otherwise load the object's display color
 */
inline gle_color bodyColorF(PointDisp const* disp, size_t s)
{
    if ( disp->coloring )
        return gle::bright_color(s).match_a(disp->color);
    return disp->color;
}

/**
Sets color material for lighting mode
if `coloring` is enabled, this loads a bright color,
otherwise load the object's display color
*/
template < typename T >
void bodyColor(T const& obj)
{
    const PointDisp * disp = obj.prop->disp;
    if ( disp->coloring )
    {
        size_t i = ( disp->coloring == 2 ? obj.mark() : obj.signature());
        gle_color col = gle::bright_color(i);
        col.load_front();
        col.darken(0.5).load_back();
    }
    else
    {
        disp->color.load_front(1.0);
        disp->color2.load_back();
    }
}

/**
 Sets color material for lighting mode
 if `coloring` is enabled, this loads the N-th bright color,
 otherwise load the object's display color
 */
inline void bodyColor(PointDisp const* disp, size_t s)
{
    if ( disp->coloring )
    {
        gle_color col = gle::bright_color(s);
        col.load_front();
        col.darken(0.5).load_back();
    }
    else
    {
        disp->color.load_front(1.0);
        disp->color2.load_back();
    }
}

