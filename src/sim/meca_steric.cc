// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University.

template <typename T>
static bool has_steric(T const& obj, int sup)
{
    int S = obj->prop->steric_key;
    if ( S < 0 || S > sup )
        throw InvalidParameter(obj->prop->name()+":steric is out-of-range");
    return S;
}


template < typename GRID >
static void setStericGrid(GRID& grid, Vector const& inf, Vector const& sup, real& range)
{
    real res = range;

    // adjust grid size to avoid excessive memory footprint:
    const size_t top = 1 << 17;
    while ( grid.setGrid(inf, sup, res) > top )
        res *= M_SQRT2;

    if ( res != range )
    {
        Cytosim::log("simul:steric_max_range <-- ", res, "\n");
        range = res;
    }
    
    /** If the Space dimensions have not changed, we would only need
     to allocate the Grid once, but this is not garanteed in general */
    grid.createCells();
}


void Meca::selectStericEngine(Simul const& sim, SimulProp const& prop)
{
    steric_ = 0;
    if ( prop.steric_mode )
    {
        Space const* spc = sim.spaces.master();
        if ( ! spc )
            throw InvalidParameter("simul:steric could not determine its Space");

        /* initialize the Grid to cover 'spc' entirely, increasing the
         cell size until we get acceptable memory requirements */
        Vector inf, sup;
        spc->boundaries(inf, sup);

        // the steric area can be restricted by parameters:
        inf = inf.e_max(prop.steric_region[0]);
        sup = sup.e_min(prop.steric_region[1]);
        
        // without attractive forces, LocusGrid will be used:
        steric_ = 1 + ( prop.steric_stiff_pull[0] <= 0 );
        
        // the grid size can be specified, but otherwise will be computed automatically:
        if ( prop.steric_max_range <= REAL_EPSILON )
            prop.steric_max_range = sim.minimumStericRange();
        
        if ( prop.steric_max_range <= REAL_EPSILON )
            throw InvalidParameter("simul:steric_max_range must be defined");

        switch ( steric_ )
        {
            case 1:
                pointGrid.stiffness(prop.steric_stiff_push[0], prop.steric_stiff_pull[0]);
                //if ( !pointGrid.hasGrid() )
                setStericGrid(pointGrid, inf, sup, prop.steric_max_range);
                break;
                
            case 2:
                locusGrid.stiffness(prop.steric_stiff_push[0]);
                //if ( !locusGrid.hasGrid() )
                setStericGrid(locusGrid, inf, sup, prop.steric_max_range);
                break;
        }
    }
}


//------------------------------------------------------------------------------
#pragma mark - PointGrid with attraction

/**
 The `prop->steric_key` is a bit-field specifying one or more 'pane' where the
 object is present. Each panes is then treated consecutively and independently,
 and only objects in the same pane may interact.
 
     for ( int pane = 1; pane <= Simul::steric_mode; ++pane )
     {
         if ( obj->prop->steric_key & pane )
         ...
     }
 
 With this mechanism, the user can flexibly configure which objects
 may see each other and thus control the steric interactions.
 
 At present, we only support 1 pane (Simul property steric).
 This can be extended if necessary, but the steric_stiffness[]
 properties should be extended as well.
 */
void Meca::addStericInteractions(PointGrid& grid, Simul const& sim)
{
    grid.clear();
    
    // distribute Fiber-points on the grid
    for ( Fiber const* F=sim.fibers.first(); F; F=F->next() )
    {
        if ( has_steric(F, grid.nbPanes()) )
        {
            const real rad = F->prop->steric_radius; // equilibrium radius
            const real rge = rad + F->prop->steric_range; // range of interaction
            const real sup = rge + 0.5 * F->segmentation();

            // include segments, in the cell associated with their center
            for ( index_t i = 0; i < F->nbSegments(); ++i )
                grid.add(F, i, rad, rge, sup);
        }
    }
    
    // include Spheres
    for ( Sphere const* O=sim.spheres.first(); O; O=O->next() )
    {
        if ( has_steric(O, grid.nbPanes()) )
        {
            real R = O->radius();
            grid.add(O, 0, R, R+O->prop->steric_range);
        }
    }
    
    // include Beads
    for ( Bead const* B=sim.beads.first(); B; B=B->next() )
    {
        if ( has_steric(B, grid.nbPanes()) )
        {
            real R = B->radius();
            grid.add(B, 0, R, R+B->prop->steric_range);
        }
    }
    
    // from Solids, include Points with radius > 0
    for ( Solid const* S=sim.solids.first(); S; S=S->next() )
    {
        if ( has_steric(S, grid.nbPanes()) )
        {
            for ( index_t i = 0; i < S->nbPoints(); ++i )
            {
                if ( S->radius(i) > REAL_EPSILON )
                    grid.add(S, i, S->radius(i), S->radius(i)+S->prop->steric_range);
            }
        }
    }
    
#if ( NUM_STERIC_PANES == 1 )
    
    grid.setSterics();

#else
    if ( sim.prop.steric_mode == 1 )
    {
        // add steric interactions inside pane 1:
        grid.setSterics(1);
        // add steric interactions between panes 1 and 2:
        grid.setSterics(1, 2);
        // add steric interactions in pane 2:
    }
    else
    {
        // add steric interactions within each pane:
        for ( size_t p = 1; p <= NUM_STERIC_PANES; ++p )
            grid.setSterics(p);
    }
#endif
}


//------------------------------------------------------------------------------
#pragma mark - streamlined LocusGrid


/// distribute Mecables with steric enabled to `grid`
static void distributeStericMecables(LocusGrid& grid, Simul const& sim)
{
    // distribute Fiber-points on the grid
    for ( Fiber const* F=sim.fibers.first(); F; F=F->next() )
    {
        if ( has_steric(F, grid.nbPanes()) )
        {
            const real rad = F->prop->steric_radius;
            const real rge = rad + 0.5 * F->segmentation();
            // include segments, in the cell associated with their center
            for ( index_t i = 0; i < F->nbSegments(); ++i )
                grid.add(F, i, rad, rge);
        }
    }

    // include Spheres
    for ( Sphere const* O=sim.spheres.first(); O; O=O->next() )
    {
        if ( has_steric(O, grid.nbPanes()) )
            grid.add(O, 0, O->radius());
    }
    
    // include Beads
    for ( Bead const* B=sim.beads.first(); B; B=B->next() )
    {
        if ( has_steric(B, grid.nbPanes()) )
            grid.add(B, 0, B->radius());
    }
    
    // from Solids, include Points with radius > 0
    for ( Solid const* S=sim.solids.first(); S; S=S->next() )
    {
        if ( has_steric(S, grid.nbPanes()) )
        {
            for ( index_t i = 0; i < S->nbPoints(); ++i )
            {
                if ( S->radius(i) > REAL_EPSILON )
                    grid.add(S, i, S->radius(i));
            }
        }
    }
}

#if ENABLE_PERIODIC_BOUNDARIES
/// distribute Mecables with steric enabled to `grid`, given periodic boundaries
static void distributeStericMecablesModulo(LocusGrid& grid, Simul const& sim)
{
    // distribute Fiber-points on the grid
    for ( Fiber const* F=sim.fibers.first(); F; F=F->next() )
    {
        if ( has_steric(F, grid.nbPanes()) )
        {
            const real rad = F->prop->steric_radius;
            const real rge = rad + 0.5 * F->segmentation();
            // include segments, in the cell associated with their center
            for ( index_t i = 0; i < F->nbSegments(); ++i )
                grid.add_modulo(F, i, rad, rge);
        }
    }
    
    // include Spheres
    for ( Sphere const* O=sim.spheres.first(); O; O=O->next() )
    {
        if ( has_steric(O, grid.nbPanes()) )
            grid.add_modulo(O, 0, O->radius());
    }
    
    // include Beads
    for ( Bead const* B=sim.beads.first(); B; B=B->next() )
    {
        if ( has_steric(B, grid.nbPanes()) )
            grid.add_modulo(B, 0, B->radius());
    }
    
    // from Solids, include Points with radius > 0
    for ( Solid const* S=sim.solids.first(); S; S=S->next() )
    {
        if ( has_steric(S, grid.nbPanes()) )
        {
            for ( index_t i = 0; i < S->nbPoints(); ++i )
            {
                if ( S->radius(i) > REAL_EPSILON )
                    grid.add_modulo(S, i, S->radius(i));
            }
        }
    }
}
#endif

/**
 The `prop->steric_key` is a bit-field specifying one or more 'pane' where the
 object is present. Each panes is then treated consecutively and independently,
 and only objects in the same pane may interact.

     for ( int pane = 1; pane <= Simul::steric_mode; ++pane )
     {
         if ( obj->prop->steric_key & pane )
         ...
     }

 With this mechanism, the user can flexibly configure which objects
 may see each other and thus control the steric interactions.
 
 At present, we only support 1 pane (Simul property steric).
 This can be extended if necessary, but the steric_stiffness[]
 properties should be extended as well.
 */
void Meca::addStericInteractions(LocusGrid& grid, Simul const& sim)
{
    grid.clear();

#if ENABLE_PERIODIC_BOUNDARIES
    if ( modulo )
        distributeStericMecablesModulo(grid, sim);
    else
#endif
        distributeStericMecables(grid, sim);
    
#if ( MAX_STERIC_PANES == 1 )
    
    grid.setSterics();

#else
    if ( sim.prop.steric_mode == 1 )
    {
        // add steric interactions inside pane 1:
        grid.setSterics(1);
        // add steric interactions between panes 1 and 2:
        grid.setSterics(1, 2);
    }
    else
    {
        // add steric interactions within each pane:
        for ( size_t p = 1; p <= MAX_STERIC_PANES; ++p )
            grid.setSterics(p);
    }
#endif
    //std::clog << "LocusGrid has capacity " << locusGrid.capacity() << "\n";
}

