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
        case 1: if ( !pointGrid.hasGrid() )
        {
            real est = sim.estimateStericRange();
            setStericGrid(pointGrid, spc, sim.prop.steric_max_range, est);
        }

        case 2: if ( !locusGrid.hasGrid() )
        {
            real est = sim.estimateStericRange();
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
void Meca::addStericInteractions1(Simul const& sim)
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
#if ( NUM_STERIC_PANES == 1 )
                pointGrid.add(F, i, rad, rge, sup);
#else
                pointGrid.add(F->prop->steric, F, i, rad, rge, sup);
#endif
        }
    }
    
    // include Spheres
    for ( Sphere const* O=sim.spheres.first(); O; O=O->next() )
    {
        if ( O->prop->steric )
#if ( NUM_STERIC_PANES == 1 )
            pointGrid.add(O, 0, O->radius(), O->radius()+O->prop->steric_range);
#else
            pointGrid.add(O->prop->steric, O, 0, O->radius(), O->radius()+O->prop->steric_range);
#endif
    }
    
    // include Beads
    for ( Bead const* B=sim.beads.first(); B; B=B->next() )
    {
        if ( B->prop->steric )
#if ( NUM_STERIC_PANES == 1 )
            pointGrid.add(B, 0, B->radius(), B->radius()+B->prop->steric_range);
#else
            pointGrid.add(B->prop->steric, B, 0, B->radius(), B->radius()+B->prop->steric_range);
#endif
    }
        
    // include Points that have a radius from Solids
    for ( Solid const* S=sim.solids.first(); S; S=S->next() )
    {
        if ( S->prop->steric )
        {
            for ( size_t i = 0; i < S->nbPoints(); ++i )
            {
                if ( S->radius(i) > REAL_EPSILON )
#if ( NUM_STERIC_PANES == 1 )
                    pointGrid.add(S, i, S->radius(i), S->radius(i)+S->prop->steric_range);
#else
                    pointGrid.add(S->prop->steric, S, i, S->radius(i), S->radius(i)+S->prop->steric_range);
#endif
            }
        }
    }
    
    /// create parameters
    Stiffness pam(sim.prop.steric_stiff_push[0], sim.prop.steric_stiff_pull[0]);
    
#if ( NUM_STERIC_PANES == 1 )
    
    pointGrid.setSterics(*this, pam);

#elif ( NUM_STERIC_PANES == 2 )
    
    // add steric interactions inside pane 1:
    pointGrid.setSterics(*this, pam, 1);
    // add steric interactions between panes 1 and 2:
    pointGrid.setSterics(*this, pam, 1, 2);
    //pointGrid.setSterics(*this, pam, 2, 1);

#else
    
    // add steric interactions within each pane:
    for ( size_t p = 1; p <= NUM_STERIC_PANES; ++p )
        pointGrid.setSterics(*this, pam, p);

#endif
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
void Meca::addStericInteractions2(Simul const& sim)
{
    // clear grid
    locusGrid.clear();
    
    // distribute Fiber-points on the grid
    for ( Fiber const* F=sim.fibers.first(); F; F=F->next() )
    {
        if ( F->prop->steric )
        {
            const real rad = F->prop->steric_radius;
            const real rge = rad + 0.5 * F->segmentation();
            // include segments, in the cell associated with their center
            for ( size_t i = 0; i < F->nbSegments(); ++i )
#if ( MAX_STERIC_PANES == 1 )
                locusGrid.add(F, i, rad, rge);
#else
                locusGrid.add(F->prop->steric, F, i, rad, rge);
#endif
        }
    }
    
    // mark edge between Fiber segments and other type of elements
    locusGrid.delimit();
    
    // include Spheres
    for ( Sphere const* O=sim.spheres.first(); O; O=O->next() )
    {
        if ( O->prop->steric )
#if ( MAX_STERIC_PANES == 1 )
            locusGrid.add(O, 0, O->radius());
#else
            locusGrid.add(O->prop->steric, O, 0, O->radius());
#endif
    }
    
    // include Beads
    for ( Bead const* B=sim.beads.first(); B; B=B->next() )
    {
        if ( B->prop->steric )
#if ( MAX_STERIC_PANES == 1 )
            locusGrid.add(B, 0, B->radius());
#else
            locusGrid.add(B->prop->steric, B, 0, B->radius());
#endif
    }
        
    // include Points that have a radius from Solids
    for ( Solid const* S=sim.solids.first(); S; S=S->next() )
    {
        if ( S->prop->steric )
        {
            for ( size_t i = 0; i < S->nbPoints(); ++i )
            {
                if ( S->radius(i) > REAL_EPSILON )
#if ( MAX_STERIC_PANES == 1 )
                    locusGrid.add(S, i, S->radius(i));
#else
                    locusGrid.add(S->prop->steric, S, i, S->radius(i));
#endif
            }
        }
    }
    
    /// stiffness
    real stiff = sim.prop.steric_stiff_push[0];

#if ( MAX_STERIC_PANES == 1 )
        
    locusGrid.setSterics(*this, stiff);

#elif ( MAX_STERIC_PANES == 2 )
    
    // add steric interactions inside pane 1:
    locusGrid.setSterics(*this, stiff, 1);
    // add steric interactions between panes 1 and 2:
    locusGrid.setSterics(*this, stiff, 1, 2);
    //pointGrid.setSterics(*this, pam, 2, 1);

#else
    
    // add steric interactions within each pane:
    for ( size_t p = 1; p <= MAX_STERIC_PANES; ++p )
        locusGrid.setSterics(*this, stiff, p);

#endif
    //std::clog << "LocusGrid has capacity " << locusGrid.capacity() << "\n";
}


void Meca::addSomeStericInteractions(real stiff)
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
#if ( MAX_STERIC_PANES == 1 )
                locusGrid.add(F, i, rad, rge);
#else
                locusGrid.add(F->prop->steric, F, i, rad, rge);
#endif
        }
    }
    
    // mark edge between Fiber segments and other type of elements
    locusGrid.delimit();

#if ( MAX_STERIC_PANES == 1 )
    
    locusGrid.setSterics(*this, stiff);

#elif ( MAX_STERIC_PANES == 2 )
    
    // add steric interactions inside pane 1:
    locusGrid.setSterics(*this, stiff, 1);
    // add steric interactions between panes 1 and 2:
    locusGrid.setSterics(*this, stiff, 1, 2);
    //pointGrid.setSterics(*this, pam, 2, 1);

#else
    
    // add steric interactions within each pane:
    for ( size_t p = 1; p <= MAX_STERIC_PANES; ++p )
        locusGrid.setSterics(*this, stiff, p);

#endif
    //std::clog << "LocusGrid has capacity " << locusGrid.capacity() << "\n";
}

