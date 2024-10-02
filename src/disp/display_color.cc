// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University

// hide definitions in anonymous namespace
namespace
{
    // colors that vary with the direction of a vector:
    static gym_color radial_color(const Vector3& d, gym_color::COLOF a = 1.f)
    { return gym_color::radial_color((gym_color::COLOF)d.XX, (gym_color::COLOF)d.YY, (gym_color::COLOF)d.ZZ, a); }
    
    gym_color radial_color(const Vector2& d, gym_color::COLOF a = 1.f)
    { return gym_color::radial_colorXY((gym_color::COLOF)d.XX, (gym_color::COLOF)d.YY, a); }
    
    gym_color radial_color(const Vector1& d, gym_color::COLOF a = 1.f)
    { if ( d.XX > 0 ) return gym_color(1,1,1,a); else return gym_color(0,1,0,a); }
    
    //------------------------------------------------------------------------------
#pragma mark - Color schemes to draw Fibers
    
    real color_scale(Fiber const* fib, int style)
    {
        switch ( style )
        {
            case 1: return 1; // color_by_fiber
            case 2: // color_by_tension
            case 3: return 1 / fib->prop->disp->tension_scale;
            case 4: return 1; // color_by_direction
            case 5: return fib->prop->disp->length_scale; // color_by_curvature
            case 6: case 7: // 6, 7 and 8: color_by_distanceP
            case 8: return fib->segmentation() / fib->prop->disp->length_scale;
            case 9: return 1 / fib->prop->disp->length_scale; // color_by_height
            case 10: return 1; // color_by_grid
        }
        return 1;
    }
    
    gym_color color_fiber(Fiber const& fib, index_t)
    {
        return fib.disp->color;
    }
    
    gym_color color_alternate(Fiber const& fib, index_t seg)
    {
        if ( seg & 1 )
            return fib.disp->color.darken(0.75);
        else
            return fib.disp->color;
    }
    
    gym_color color_by_tension(Fiber const& fib, index_t seg)
    {
        real x = fib.disp->color_scale * fib.tension(seg);
        return fib.disp->color.alpha(x);
    }
    
    gym_color color_by_tension_jet(Fiber const& fib, index_t seg)
    {
        real x = fib.disp->color_scale * fib.tension(seg);
        // use jet coloring, where Lagrange multipliers are negative under compression
        return gym_color::jet_color_alpha(x);
        //return fib.disp->color.alpha(x-1);
    }
    
    gym_color color_by_curvature(Fiber const& fib, index_t pti)
    {
        const real beta = fib.disp->color_scale;
        return gym_color::jet_color(beta*fib.curvature(pti));
    }
    
    gym_color color_seg_curvature(Fiber const& fib, index_t seg)
    {
        const real beta = fib.disp->color_scale;
        real c = 0.5 * (fib.curvature(seg) + fib.curvature(seg+1));
        return gym_color::jet_color(beta*c);
    }
    
    gym_color color_by_direction(Fiber const& fib, index_t seg)
    {
        return radial_color(fib.dirSegment(seg));
    }
    
    /// using distance from the minus end to the start of segment `seg`
    gym_color color_by_distanceM(Fiber const& fib, real pti)
    {
        const real beta = fib.disp->color_scale;
        real x = min_real(pti*beta, 32);
        return fib.disp->color.alpha(std::exp(-x));
    }
    
    /// using the distance from the plus end to vertex `pti`
    gym_color color_by_distanceP(Fiber const& fib, real pti)
    {
        const real beta = fib.disp->color_scale;
        // using the distance at the vertex
        real x = min_real((fib.lastPoint()-pti)*beta, 32);
        return fib.disp->color.alpha(std::exp(-x));
    }
    
    /// using distance from the plus end to the end of segment `seg`
    gym_color color_by_abscissaM(Fiber const& fib, index_t seg)
    {
        const real beta = fib.disp->color_scale;
        real x = min_real(seg*beta, 32);
        return fib.disp->color.alpha(std::exp(-x));
    }
    
    /// using distance from the plus end to the end of segment `seg`
    gym_color color_by_abscissaP(Fiber const& fib, index_t seg)
    {
        const real beta = fib.disp->color_scale;
        real x = min_real((fib.lastSegment()-seg)*beta, 32);
        return fib.disp->color.alpha(std::exp(-x));
    }
    
    /// color set according to distance to the confining Space
    gym_color color_by_height(Fiber const& fib, index_t pti)
    {
        const real beta = fib.disp->color_scale;
        real Z = 0;
        Space const* spc = fib.prop->confine_space;
        if ( spc && fib.prop->confine )
            Z = -spc->signedDistanceToEdge(fib.posPoint(pti));
#if ( DIM > 2 )
        else
            Z = fib.posPoint(pti).ZZ;
#endif
        return gym_color::jet_color(Z*beta);
    }
    
    
    /// color set according to some Map
    gym_color color_by_grid(Fiber const& fib, index_t seg)
    {
        Map<DIM> const& map = fib.simul().visibleMap();
        Vector w = fib.midPoint(seg);
        index_t i = map.index(w);
        return gym::alt_color(i);
    }
    
    
    //------------------------------------------------------------------------------
#pragma mark - Color schemes for Fiber Lattice
    
    gym_color color_alternate(Fiber const& fib, int ix, real)
    {
        if ( ix & 1 )
            return fib.disp->color.darken(0.75);
        else
            return fib.disp->color;
    }
    
    gym_color color_by_lattice(Fiber const& fib, int ix, real beta)
    {
        return fib.disp->color.darken(beta*fib.visibleLattice()->data(ix));
    }
    
    gym_color color_by_lattice_jet(Fiber const& fib, int ix, real beta)
    {
        return gym_color::jet_color(beta*fib.visibleLattice()->data(ix));
    }
    
    gym_color color_by_lattice_white(Fiber const& fib, int ix, real beta)
    {
        real x = beta * fib.visibleLattice()->data(ix);
        return gym_color(x, x, x);
    }
    
    gym_color lattice_color(gym_color const& col, real val)
    {
        return col.darken(val);
        //return gym_color::jet_color(val);
    }
    
    //------------------------------------------------------------------------------
#pragma mark - Color schemes for Sphere / Solid / Bead
    
    /**
     returns the front color to be used to display the object.
     
     If `coloring` is enabled, the color depends on the object's `mark` and `signature`.
     `coloring==1` will randomly color according to signature
     `coloring==2` will color according to mark
     `coloring==3` will use different tones of colors around the normal display color
     otherwise load the object's display color
     */
    template < typename T >
    gym_color bodyColorF(T const& obj)
    {
        PointDisp const* dis = obj.prop->disp;
        if ( dis->coloring && obj.mark() )
            return gym::bright_color(obj.mark());
        switch ( dis->coloring )
        {
            case 1: return gym::bright_color(obj.signature());
            case 2: return gym::bright_color(obj.mark());
            case 3: return dis->color.tweak(obj.signature());
            default: return dis->color;
        }
    }
    
    /**
     Sets color material for lighting mode
     if `coloring` is enabled, this loads a bright color,
     otherwise load the object's display color
     */
    template < typename T >
    void bodyColor(T const& obj)
    {
        PointDisp const* dis = obj.prop->disp;
        if ( dis->coloring )
        {
            index_t i = ( dis->coloring == 2 ? obj.mark() : obj.signature());
            gym_color col = gym::bright_color(i);
            gym::color_load(col);
            gym::color_back(col.darken(0.5));
        }
        else
        {
            gym::color_front(dis->color, 1.0);
            gym::color_back(dis->color2);
        }
    }
    
}
