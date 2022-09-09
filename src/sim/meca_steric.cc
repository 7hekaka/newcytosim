// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University.



template < typename GRID >
static void setStericGrid(GRID& grid, Space const* spc, real& range, real inf)
{
    assert_true(spc);
    real res = range;
    
    res = std::max(res, inf);

    if ( res <= 0 )
        throw InvalidParameter("simul:steric_max_range must be defined");

    const size_t sup = 1 << 17;
    while ( grid.setGrid(spc, res) > sup )
        res *= M_SQRT2;

    if ( res != range )
    {
        Cytosim::log("simul:steric_max_range <-- %.3f\n", res);
        range = res;
    }
    
    grid.createCells();
}


void Meca::selectStericEngine(Simul const& sim)
{
    Space const* spc = sim.spaces.master();

    if ( !spc || !sim.prop.steric_mode )
    {
        steric_ = 0;
        return;
    }
    
    // without pulling, use `locusGrid` which is simpler:
    steric_ = 1 + ( sim.prop.steric_stiff_pull[0] <= 0 );
    
    switch ( steric_ )
    {
        case 1:
            pointGrid.stiffness(sim.prop.steric_stiff_push[0], sim.prop.steric_stiff_pull[0]);
            if ( !pointGrid.hasGrid() )
            {
                real est = sim.minimumStericRange();
                setStericGrid(pointGrid, spc, sim.prop.steric_max_range, est);
            }

        case 2:
            locusGrid.stiffness(sim.prop.steric_stiff_push[0]);
            if ( !locusGrid.hasGrid() )
            {
                real est = sim.minimumStericRange();
                setStericGrid(locusGrid, spc, sim.prop.steric_max_range, est);
            }
    }
}


/**
 The prop->steric of each object is a bit-field that
 specify one or more 'pane' where the object is present.
 The different panes are then treated consecutively and independently,
 and only objects in the same pane may interact.
 
     for ( int pane=1; pane<=2 && pane<=prop->steric; ++pane )
     {
         if ( obj->prop->steric & pane )
         ...
     }
 
 With this mechanism, the user can flexibly configure which objects
 may see each other and thus control the steric interactions.
 
 At present, we only support 1 pane (Simul property steric).
 This can be extended if necessary, but the steric_stiffness[]
 properties should be extended as well.
 */
void Meca::addStericInteractionsP(Simul const& sim)
{
    // clear grid
    pointGrid.clear();
    
    // distribute Fiber-points on the grid
    for ( Fiber const* F=sim.fibers.first(); F; F=F->next() )
    {
        if ( F->prop->steric )
        {
            const real rad = F->prop->steric_radius;        // equilibrium radius
            const real rge = rad + F->prop->steric_range;   // extended range of interaction
            const real sup = rge + 0.5 * F->segmentation();

            // include segments, in the cell associated with their center
            for ( size_t i = 0; i < F->nbSegments(); ++i )
                pointGrid.add(F, i, rad, rge, sup);
        }
    }
    
    // include Spheres
    for ( Sphere const* O=sim.spheres.first(); O; O=O->next() )
    {
        if ( O->prop->steric )
            pointGrid.add(O, 0, O->radius(), O->radius()+O->prop->steric_range);
    }
    
    // include Beads
    for ( Bead const* B=sim.beads.first(); B; B=B->next() )
    {
        if ( B->prop->steric )
            pointGrid.add(B, 0, B->radius(), B->radius()+B->prop->steric_range);
    }
        
    // from Solids, include Points with radius > 0
    for ( Solid const* S=sim.solids.first(); S; S=S->next() )
    {
        if ( S->prop->steric )
        {
            for ( size_t i = 0; i < S->nbPoints(); ++i )
            {
                if ( S->radius(i) > REAL_EPSILON )
                    pointGrid.add(S, i, S->radius(i), S->radius(i)+S->prop->steric_range);
            }
        }
    }
    
#if ( NUM_STERIC_PANES == 1 )
    
    pointGrid.setSterics();

#elif ( NUM_STERIC_PANES == 2 )
    
    // add steric interactions inside pane 1:
    pointGrid.setSterics(1);
    // add steric interactions between panes 1 and 2:
    pointGrid.setSterics(1, 2);
    //pointGrid.setSterics(*this, pam, 2, 1);

#else
    
    // add steric interactions within each pane:
    for ( size_t p = 1; p <= NUM_STERIC_PANES; ++p )
        pointGrid.setSterics(p);

#endif
}




/// add Mecables with steric enabled to the steric grid
static void addStericMecables(LocusGrid& grid, Simul const& sim)
{
    // distribute Fiber-points on the grid
    for ( Fiber const* F=sim.fibers.first(); F; F=F->next() )
    {
        if ( F->prop->steric )
        {
            const real rad = F->prop->steric_radius;
            const real rge = rad + 0.5 * F->segmentation();
            // include segments, in the cell associated with their center
            for ( size_t i = 0; i < F->nbSegments(); ++i )
                grid.add(F, i, rad, rge);
        }
    }
    
    // mark edge between Fiber segments and other type of elements
    grid.delimit();
    
    // include Spheres
    for ( Sphere const* O=sim.spheres.first(); O; O=O->next() )
    {
        if ( O->prop->steric )
            grid.add(O, 0, O->radius());
    }
    
    // include Beads
    for ( Bead const* B=sim.beads.first(); B; B=B->next() )
    {
        if ( B->prop->steric )
            grid.add(B, 0, B->radius());
    }
        
    // from Solids, include Points with radius > 0
    for ( Solid const* S=sim.solids.first(); S; S=S->next() )
    {
        if ( S->prop->steric )
        {
            for ( size_t i = 0; i < S->nbPoints(); ++i )
            {
                if ( S->radius(i) > REAL_EPSILON )
                    grid.add(S, i, S->radius(i));
            }
        }
    }
}

#if GRID_HAS_PERIODIC
/// add Mecables with steric enabled to the steric grid, for periodic boundary conditions
static void addStericMecablesModulo(LocusGrid& grid, Simul const& sim)
{
    // distribute Fiber-points on the grid
    for ( Fiber const* F=sim.fibers.first(); F; F=F->next() )
    {
        if ( F->prop->steric )
        {
            const real rad = F->prop->steric_radius;
            const real rge = rad + 0.5 * F->segmentation();
            // include segments, in the cell associated with their center
            for ( size_t i = 0; i < F->nbSegments(); ++i )
                grid.add_modulo(F, i, rad, rge);
        }
    }
    
    // mark edge between Fiber segments and other type of elements
    grid.delimit();
    
    // include Spheres
    for ( Sphere const* O=sim.spheres.first(); O; O=O->next() )
    {
        if ( O->prop->steric )
            grid.add_modulo(O, 0, O->radius());
    }
    
    // include Beads
    for ( Bead const* B=sim.beads.first(); B; B=B->next() )
    {
        if ( B->prop->steric )
            grid.add_modulo(B, 0, B->radius());
    }
        
    // from Solids, include Points with radius > 0
    for ( Solid const* S=sim.solids.first(); S; S=S->next() )
    {
        if ( S->prop->steric )
        {
            for ( size_t i = 0; i < S->nbPoints(); ++i )
            {
                if ( S->radius(i) > REAL_EPSILON )
                    grid.add_modulo(S, i, S->radius(i));
            }
        }
    }
}
#endif

/**
 The prop->steric of each object is a bit-field that
 specify one or more 'pane' where the object is present.
 The different panes are then treated consecutively and independently,
 and only objects in the same pane may interact.
 
     for ( int pane=1; pane<=2 && pane<=prop->steric; ++pane )
     {
         if ( obj->prop->steric & pane )
         ...
     }
 
 With this mechanism, the user can flexibly configure which objects
 may see each other and thus control the steric interactions.
 
 At present, we only support 1 pane (Simul property steric).
 This can be extended if necessary, but the steric_stiffness[]
 properties should be extended as well.
 */
void Meca::addStericInteractionsL(Simul const& sim)
{
    // clear grid
    locusGrid.clear();

#if GRID_HAS_PERIODIC
    if ( modulo )
        addStericMecablesModulo(locusGrid, sim);
    else
#endif
        addStericMecables(locusGrid, sim);
    
#if ( MAX_STERIC_PANES == 1 )
        
    locusGrid.setSterics();

#elif ( MAX_STERIC_PANES == 2 )
    
    // add steric interactions inside pane 1:
    locusGrid.setSterics(1);
    // add steric interactions between panes 1 and 2:
    locusGrid.setSterics(1, 2);
    //pointGrid.setSterics(*this, pam, 2, 1);

#else
    
    // add steric interactions within each pane:
    for ( size_t p = 1; p <= MAX_STERIC_PANES; ++p )
        locusGrid.setSterics(p);

#endif
    //std::clog << "LocusGrid has capacity " << locusGrid.capacity() << "\n";
}


void Meca::addSomeStericInteractions()
{
    locusGrid.clear();
    
    // distribute Fiber-points on the grid
    for ( Mecable * mec : mecables )
    {
        Fiber * F = static_cast<Fiber*>(mec);
        if ( F->prop->steric )
        {
            const real rad = F->prop->steric_radius;
            const real rge = rad + 0.5 * F->segmentation();
            // include segments, in the cell associated with their center
            for ( size_t i = 0; i < F->nbSegments(); ++i )
                locusGrid.add(F, i, rad, rge);
        }
    }
    
    // mark edge between Fiber segments and other type of elements
    locusGrid.delimit();

#if ( MAX_STERIC_PANES == 1 )
    
    locusGrid.setSterics();

#elif ( MAX_STERIC_PANES == 2 )
    
    // add steric interactions inside pane 1:
    locusGrid.setSterics(1);
    // add steric interactions between panes 1 and 2:
    locusGrid.setSterics(1, 2);
    //pointGrid.setSterics(*this, pam, 2, 1);

#else
    
    // add steric interactions within each pane:
    for ( size_t p = 1; p <= MAX_STERIC_PANES; ++p )
        locusGrid.setSterics(p);

#endif
    //std::clog << "LocusGrid has capacity " << locusGrid.capacity() << "\n";
}

