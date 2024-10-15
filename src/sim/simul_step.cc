// Cytosim was created by Francois Nedelec. Copyright 2024 Cambridge University

/**
 This is the master Monte-Carlo step function.
 
 Lists are mixed such that objects are considered in a different
 and random order at each step, to avoid biais in the simulation

 steps() is called for every list, calling step() for every Object
 */
void Simul::steps()
{
    //auto rdt = timer();
    // increment time:
    prop.time += prop.time_step;
    //fprintf(stderr, "\n----------------------------------- time is %8.3f\n", prop.time);

    // Monte-Carlo step for all objects
    events.steps();
    organizers.steps();
    fields.steps();
    spaces.steps();
    spheres.steps();
    beads.steps();
    solids.steps();

    //printf("     ::steps    %16llu\n", (timer()-rdt)>>5); rdt = timer();

#if POOL_UNATTACHED > 1
    doAttachCounter = ( doAttachCounter + 1 ) % POOL_UNATTACHED;
    if ( doAttachCounter )
    {
        couples.stepsSkippingUnattached();
        singles.stepsSkippingUnattached();
        //printf("     ::noattach %16llu\n", (timer()-rdt)>>3);
    }
    else
#endif
    {
        // distribute Fibers over a grid preparing for binding of Hands:
        const real range = maxBindingRange();
        fiberGrid.paintGrid(fibers.first(), nullptr, range);
        
        //printf("     ::paint    %16llu\n", (timer()-rdt)>>5); rdt = timer();
        
        if ( 0 )
        {
            // This code continuously tests the binding algorithm.
            HandProp hp("test_binding");
            hp.binding_rate  = INFINITY;
            hp.binding_range = RNG.preal() * range;
            hp.bind_also_end = BOTH_ENDS;
            hp.complete(*this);
            
            Space const* spc = spaces.master();
            for ( size_t i = 0; i < 16; ++i )
            {
                Vector pos = spc->place();
                fiberGrid.testAttach(stdout, pos, fibers, &hp);
            }
        }

        // step Hand-containing objects, giving them a possibility to attach Fibers:
        couples.steps();
        singles.steps();
        //printf("     ::attach   %16llu\n", (timer()-rdt)>>3);
    }

    // This will also update all the attached Hands
    fibers.steps();

    fresh_ = 1;
}


//------------------------------------------------------------------------------
#pragma mark -


// calculate grid range from Hand's binding range:
real Simul::maxBindingRange() const
{
    real res = 0.0;
    for ( Property const* i : properties.find_all("hand") )
        res = std::max(res, static_cast<HandProp const*>(i)->binding_range);
    return res;
}


/**
 @returns the largest segmentation of all known FiberProp, or -1.0 if none is defined
 */
real Simul::estimateFiberGridStep() const
{
    real res = -1.0;
    
    for ( Property const* i : properties.find_all("fiber") )
    {
        FiberProp const* fp = static_cast<FiberProp const*>(i);
        res = std::max(res, fp->segmentation);
    }
    
    return res;
}


/**
 The FiberGrid is used to quickly find the fibers that are close to any point.
 In brief:
 1. if `binding_grid_step` is not set, attempt to find a suitable value for it,
 2. if the number of cells is superior to 1e5, double the step size,
 3. initialize the grid with the estimated step size.
 
 Note: if no FiberProp is defined, there cannot be any fiber,
 and the grid is not needed. In this case, the grid is initialized with step=1
 */
void Simul::setFiberGrid(FiberGrid& grid, Space const* spc, real& grid_step)
{
    assert_true(spc);
    real res = grid_step;

    /* initialize the Grid to cover 'spc' entirely, increasing the
     cell size until we get acceptable memory requirements */
    Vector inf, sup;
    spc->boundaries(inf, sup);

    const size_t too_much = 1 << 17;
    while ( grid.setGrid(inf, sup, abs_real(res)) > too_much )
        res *= M_SQRT2;

    if ( res > 0 && res != grid_step )
    {
        Cytosim::log("simul:binding_grid_step <-- ", res, "\n");
        grid_step = res;
    }

    // create the grid voxels:
    grid.createCells();
}


/**
 Will pepare the simulation engine to make it able to execute step():
 - set FiberGrid used for attachment of Hands,
 - set StericGrid
 - call complete() for all registered Property
 .
 The simulated objects should not be changed.
 
 */
void Simul::prepare()
{
    //fprintf(stderr, "Simul:%p:prepare()\n", this);
    if ( !spaces.master() )
        throw InvalidSyntax("A space must be defined first!");

    primed_ = 1;

    // make sure properties are ready for simulations:
    prop.complete(*this);
    
    // try to determine cell size from the filaments characteristics
    if ( prop.binding_grid_step <= 0 )
        prop.binding_grid_step = estimateFiberGridStep();

    // prepare grid for attachments:
    setFiberGrid(fiberGrid, spaces.master(), prop.binding_grid_step);
    
    // this will allocate Fiber::Lattice
    fibers.updateFibers();
    // this is necessary for diffusion in Field:
    fields.prepare();
    
    // this prepares for 'fast_diffusion':
    singles.prepare();
    couples.prepare();
    
    primed_ = 2;
}


void Simul::relax()
{
    singles.relax();
    couples.relax();
    primed_ = 0;
}


void Simul::drawLinks() const
{
    prop.complete(*this);
    sMeca.getReady(*this);
    sMeca.drawLinks = 1;
    setAllInteractions(sMeca);
    sMeca.drawLinks = 0;
}
