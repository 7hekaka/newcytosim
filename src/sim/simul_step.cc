// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

/**
 @returns the largest segmentation of all existing FiberProp, multiplied by 0.5
 @returns -1.0 if no FiberProp is defined
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
void Simul::setFiberGrid(Space const* spc, real& grid_step) const
{
    assert_true(spc);
    real res = grid_step;

    // try to find cell size from the filaments characteristics
    if ( res <= 0 )
        res = estimateFiberGridStep();

    /* initialize the Grid to cover 'spc' entirely, increasing the
     cell size until we get acceptable memory requirements */
    const size_t sup = 1 << 17;
    while ( fiberGrid.setGrid(spc, abs_real(res)) > sup )
        res *= M_SQRT2;

    if ( res > 0 && res != grid_step )
    {
        Cytosim::log("simul:binding_grid_step <-- %.3f\n", res);
        grid_step = res;
    }

    // create the grid cells:
    fiberGrid.createCells();
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
    
    // prepare grid for attachments:
    setFiberGrid(spaces.master(), prop.binding_grid_step);
    
    // this will allocate Fiber::Lattice
    fibers.updateFibers();
    // this is necessary for diffusion in Field:
    fields.prepare();
    
    // this prepares for 'fast_diffusion':
    singles.prepare(properties);
    couples.prepare(properties);
    
    primed_ = 2;
}


// calculate grid range from Hand's binding range:
real Simul::maxBindingRange() const
{
    real res = 0.0;
    for ( Property const* i : properties.find_all("hand") )
        res = std::max(res, static_cast<HandProp const*>(i)->binding_range);
    return res;
}


/**
 This is the master Monte-Carlo step function.
 
 Lists are mixed such that objects are considered in a different
 and random order at each step, to avoid biais in the simulation

 step() is called for every list, i.e. for every Object
 */
void Simul::step()
{
    //auto rdt = timer();
    // increment time:
    prop.time += prop.time_step;
    //fprintf(stderr, "\n----------------------------------- time is %8.3f\n", prop.time);

    // mix object lists
    if ( events.size() > 1 ) events.shuffle();
    if ( organizers.size() > 1 ) organizers.shuffle();
    if ( tubules.size() > 1 ) tubules.shuffle();
    if ( beads.size() > 1 ) beads.shuffle();
    if ( solids.size() > 1 ) solids.shuffle();
    if ( fibers.size() > 1 ) fibers.shuffle();
    if ( spheres.size() > 1 ) spheres.shuffle();
    if ( couples.size() > 1 ) couples.shuffle();
    if ( singles.size() > 1 ) singles.shuffle();
    if ( spaces.size() > 1 ) spaces.shuffle();
    if ( fields.size() > 1 ) fields.shuffle();

    //printf("Simul::shuffles %16llu\n", (timer()-rdt)>>5); rdt = timer();

    // Monte-Carlo step for all objects
    events.step();
    organizers.step();
    tubules.step();
    fields.step();
    spaces.step();
    spheres.step();
    beads.step();
    solids.step();
    
    //printf("     ::steps    %16llu\n", (timer()-rdt)>>5); rdt = timer();

#if POOL_UNATTACHED > 1
    doAttachCounter = ( doAttachCounter + 1 ) % POOL_UNATTACHED;
    if ( doAttachCounter )
    {
        couples.stepSkipUnattached();
        singles.stepSkipUnattached();
        //printf("     ::noattach %16llu\n", (timer()-rdt)>>3);
    }
    else
#endif
    {
        // distribute Fibers over a grid for binding of Hands:
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
        couples.step();
        singles.step();
        //printf("     ::attach   %16llu\n", (timer()-rdt)>>3);
    }
    
    // This will also update all the attached Hands
    fibers.step();
    fresh_ = 1;
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
