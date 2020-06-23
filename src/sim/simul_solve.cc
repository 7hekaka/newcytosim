// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

/**
 Calculate the minimum grid cell size, given the segmentation of Fibers that
 is within a certain range because their length may vary, and taking into
 account the radius of Bead, Sphere and Solid.
 
 This function can be used to set SimulProp::steric_max_range.
 
 It assumes that Fiber::adjustSegmentation() is used, such that at any time:
     actual_segmentation <  4/3 * FiberProp::segmentation
 */
real Simul::estimateStericRange() const
{
    real ran = 0;
    real len = 0;
    
    // check all FiberProp with enabled steric:
    for ( Property const* i : properties.find_all("fiber") )
    {
        FiberProp const* fp = static_cast<FiberProp const*>(i);
        if ( fp->steric )
        {
            // The maximum length of a segment is 4/3 * segmentation
            len = std::max(len, (real)(1.34) * fp->segmentation);
            
            // check extended range of interaction
            ran = std::max(ran, fp->steric_radius + fp->steric_range);
        }
    }
    
    // verify against the actual segmentations of the Fibers:
    for ( Fiber const* fib=fibers.first(); fib; fib=fib->next() )
    {
        if ( fib->prop->steric )
            len = std::max(len, fib->segmentation());
    }

    /*
     The interaction can be aligned with the fiber, and we must add the distances:
     2 * range if two fibers of radius 'range' interact.
     + 2 * ( len / 2 ) since len/2 is the distance between the center of the segment
     and its most distal point.
     */
    ran = len + 2*ran;
    
    
    for ( Sphere const* sp=spheres.first(); sp; sp=sp->next() )
    {
        if ( sp->prop->steric )
            ran = std::max(ran, 2 * sp->radius() + sp->prop->steric_range);
    }
    
    for ( Bead const* bd=beads.first(); bd; bd=bd->next() )
    {
        if ( bd->prop->steric )
            ran = std::max(ran, 2 * bd->radius() + bd->prop->steric_range);
    }
    
    for ( Solid const* so=solids.first(); so; so=so->next() )
    {
        if ( so->prop->steric )
        {
            for ( size_t p = 0; p < so->nbPoints(); ++p )
                ran = std::max(ran, 2 * so->radius(p) + so->prop->steric_range);
        }
    }
    
    if ( ran < REAL_EPSILON )
        Cytosim::warn << "could not estimate simul:steric_max_range automatically!\n";
    
    return ran;
}


template < typename GRID >
void Simul::setStericGrid(GRID& grid, Space const* spc) const
{
    assert_true(spc);
    real res = prop->steric_max_range;
    real inf = estimateStericRange();
    
    res = std::max(res, inf);

    if ( res <= 0 )
        throw InvalidParameter("simul:steric_max_range must be defined");

    const size_t sup = 1 << 17;
    while ( grid.setGrid(spc, res) > sup )
        res *= M_SQRT2;

    if ( res != prop->steric_max_range )
    {
        Cytosim::log("adjusting simul:steric_max_range = %.3f\n", res);
        prop->steric_max_range = res;
    }
    
    grid.createCells();
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
void Simul::setStericInteractions(Meca& meca) const
{
    if ( !pointGrid.hasGrid() )
    {
        if (!spaces.master())
            return;
        setStericGrid(pointGrid, spaces.master());
    }

    // clear grid
    pointGrid.clear();
    
    // distribute Fiber-points on the grid
    for ( Fiber const* fib=fibers.first(); fib; fib=fib->next() )
    {
        if ( fib->prop->steric )
        {
            const real rad = fib->prop->steric_radius;        // equilibrium radius
            const real ran = rad + fib->prop->steric_range;   // extended range of interaction
        
            // include segments, in the cell associated with their center
            for ( size_t r = 0; r < fib->nbSegments(); ++r )
#if ( N_STERIC_PANES == 1 )
                pointGrid.add(FiberSegment(fib, r), rad, ran);
#else
                pointGrid.add(fib->prop->steric, FiberSegment(fib, r), rad, ran);
#endif
        }
    }
    
    // include Spheres
    for ( Sphere const* sp=spheres.first(); sp; sp=sp->next() )
    {
        if ( sp->prop->steric )
#if ( N_STERIC_PANES == 1 )
            pointGrid.add(Mecapoint(sp, 0), sp->radius(), sp->radius()+sp->prop->steric_range);
#else
            pointGrid.add(sp->prop->steric, Mecapoint(sp, 0), sp->radius(), sp->radius()+sp->prop->steric_range);
#endif
    }
    
    // include Beads
    for ( Bead const* bd=beads.first(); bd; bd=bd->next() )
    {
        if ( bd->prop->steric )
#if ( N_STERIC_PANES == 1 )
            pointGrid.add(Mecapoint(bd, 0), bd->radius(), bd->radius()+bd->prop->steric_range);
#else
            pointGrid.add(bd->prop->steric, Mecapoint(bd, 0), bd->radius(), bd->radius()+bd->prop->steric_range);
#endif
    }
        
    // include Points that have a radius from Solids
    for ( Solid const* so=solids.first(); so; so=so->next() )
    {
        if ( so->prop->steric )
        {
            for ( size_t i = 0; i < so->nbPoints(); ++i )
            {
                if ( so->radius(i) > REAL_EPSILON )
#if ( N_STERIC_PANES == 1 )
                    pointGrid.add(Mecapoint(so, i), so->radius(i), so->radius(i)+so->prop->steric_range);
#else
                    pointGrid.add(so->prop->steric, Mecapoint(so, i), so->radius(i), so->radius(i)+so->prop->steric_range);
#endif
            }
        }
    }
    
    /// create parameters
    StericParam pam(prop->steric_stiffness_push[0], prop->steric_stiffness_pull[0]);
    
    assert_true(prop->steric_stiffness_push[0] >= 0);
    assert_true(prop->steric_stiffness_pull[0] >= 0);
    
#if ( N_STERIC_PANES == 1 )
    
    pointGrid.setInteractions(meca, pam);

#elif ( N_STERIC_PANES == 2 )
    
    // add steric interactions inside pane 1:
    pointGrid.setInteractions(meca, pam, 1);
    // add steric interactions between panes 1 and 2:
    pointGrid.setInteractions(meca, pam, 1, 2);
    //pointGrid.setInteractions(meca, pam, 2, 1);

#else
    
    // add steric interactions between different panes:
    for ( size_t p = 1; p <= N_STERIC_PANES; ++p )
        pointGrid.setInteractions(meca, pam, p);

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
void Simul::setStericInteractionsF(Meca& meca) const
{
    if ( !pointGridF.hasGrid() )
    {
        if (!spaces.master())
            return;
        setStericGrid(pointGridF, spaces.master());
    }

    // clear grid
    pointGridF.clear();
    
    // distribute Fiber-points on the grid
    for ( Fiber const* fib=fibers.first(); fib; fib=fib->next() )
    {
        if ( fib->prop->steric )
        {
            const real rad = fib->prop->steric_radius;
            // include segments, in the cell associated with their center
            for ( size_t r = 0; r < fib->nbSegments(); ++r )
#if ( MAX_STERIC_PANES == 1 )
                pointGridF.add(FiberSegment(fib, r), rad);
#else
                pointGridF.add(fib->prop->steric, FiberSegment(fib, r), rad);
#endif
        }
    }
    
    // include Spheres
    for ( Sphere const* sp=spheres.first(); sp; sp=sp->next() )
    {
        if ( sp->prop->steric )
#if ( MAX_STERIC_PANES == 1 )
            pointGridF.add(Mecapoint(sp, 0), sp->radius());
#else
            pointGridF.add(sp->prop->steric, Mecapoint(sp, 0), sp->radius());
#endif
    }
    
    // include Beads
    for ( Bead const* bd=beads.first(); bd; bd=bd->next() )
    {
        if ( bd->prop->steric )
#if ( MAX_STERIC_PANES == 1 )
            pointGridF.add(Mecapoint(bd, 0), bd->radius());
#else
            pointGridF.add(bd->prop->steric, Mecapoint(bd, 0), bd->radius());
#endif
    }
        
    // include Points that have a radius from Solids
    for ( Solid const* so=solids.first(); so; so=so->next() )
    {
        if ( so->prop->steric )
        {
            for ( size_t i = 0; i < so->nbPoints(); ++i )
            {
                if ( so->radius(i) > REAL_EPSILON )
#if ( MAX_STERIC_PANES == 1 )
                    pointGridF.add(Mecapoint(so, i), so->radius(i));
#else
                    pointGridF.add(so->prop->steric, Mecapoint(so, i), so->radius(i));
#endif
            }
        }
    }
    
    /// create parameters
    real stiff = prop->steric_stiffness_push[0];
    
    assert_true(prop->steric_stiffness_push[0] >= 0);

#if ( MAX_STERIC_PANES == 1 )
        
    pointGridF.setInteractions(meca, stiff);

#elif ( MAX_STERIC_PANES == 2 )
    
    // add steric interactions inside pane 1:
    pointGridF.setInteractions(meca, stiff, 1);
    // add steric interactions between panes 1 and 2:
    pointGridF.setInteractions(meca, stiff, 1, 2);
    //pointGrid.setInteractions(meca, pam, 2, 1);

#else
    
    // add steric interactions between different panes:
    for ( size_t p = 1; p <= MAX_STERIC_PANES; ++p )
        pointGridF.setInteractions(meca, stiff, p);

#endif
}


//------------------------------------------------------------------------------
/**
 This will:
 - call setInteractions() for all objects in the system,
 - call setStericInteractions() if prop->steric is true.
 .
 */
void Simul::setAllInteractions(Meca& meca) const
{
    for ( Space const* s=spaces.first(); s; s=s->next() )
        s->setInteractions(meca);
    
    for ( Fiber const* f=fibers.first(); f ; f=f->next() )
        f->setInteractions(meca);
    
    for ( Solid const* s=solids.first(); s ; s=s->next() )
        s->setInteractions(meca);
    
    for ( Sphere const* o=spheres.first(); o ; o=o->next() )
        o->setInteractions(meca);
    
    for ( Bead const* b=beads.first(); b ; b=b->next() )
        b->setInteractions(meca);

    for ( Single const* i=singles.firstA(); i ; i=i->next() )
        i->setInteractions(meca);

    for ( Couple const* c=couples.firstAA(); c ; c=c->next() )
        c->setInteractions(meca);

#ifdef NEW_DANGEROUS_CONFINEMENTS
    for ( Couple const* c=couples.firstAF(); c ; c=c->next() )
        c->setInteractionsAF(meca);

    for ( Couple const* c=couples.firstFA(); c ; c=c->next() )
        c->setInteractionsFA(meca);
#endif
    
    for ( Organizer const* x = organizers.first(); x; x=x->next() )
        x->setInteractions(meca);
    
    for ( Tubule const* t = tubules.first(); t; t=t->next() )
        t->setInteractions(meca);

    //for ( Event const* e = events.first(); e; e=e->next() )
    //    e->setInteractions(meca);

    // add steric interactions
    if ( prop->steric )
    {
        if ( prop->steric_stiffness_pull[0] > 0 )
            setStericInteractions(meca);
        else
            setStericInteractionsF(meca);
    }
    //addExperimentalInteractions(meca);

#if ( 0 )
    // add steric interaction between a sphere and all fibers
    Sphere * sol = spheres.first();
    if ( sol )
    {
        const Vector cen = sol->posPoint(0);
        const real rad = sol->radius();
        const real rad2 = square(rad);
        const real stiff = prop->steric_stiffness_push[0];

        for ( Fiber const* fib = fibers.first(); fib; fib = fib->next() )
        {
            for ( size_t n = 0; n < fib->nbSegments(); ++n )
            {
                FiberSegment seg(fib, n);
                real dis = INFINITY;
                real abs = seg.projectPoint(cen, dis);
                if ( dis < rad2 )
                    meca.addSideLink(Interpolation(seg, abs), Mecapoint(sol, 0), rad, stiff);
            }
        }
    }
#endif
}

#pragma mark -

void Simul::solve()
{
    sMeca.prepare(this);
    //auto rdt = __rdtsc();
    setAllInteractions(sMeca);
    //printf("     ::set      %16llu\n", (__rdtsc()-rdt)>>5); rdt = __rdtsc();
    sMeca.solve(prop, prop->precondition);
    //printf("     ::solve    %16llu\n", (__rdtsc()-rdt)>>5); rdt = __rdtsc();
    sMeca.apply();
    //printf("     ::apply    %16llu\n", (__rdtsc()-rdt)>>5);
}


void Simul::solve_half()
{
    sMeca.prepare(this);
    setAllInteractions(sMeca);
    sMeca.solve(prop, prop->precondition);
}



/**
 Solve the system, and automatically select the fastest preconditionning method
 */
void Simul::solve_auto()
{
    sMeca.prepare(this);
    setAllInteractions(sMeca);
    
    // solve the system, recording time:
    //double cpu = TicToc::milliseconds();
    size_t cnt = sMeca.solve(prop, autoPrecond);
    //cpu = TicToc::milliseconds() - cpu;
    double cpu = sMeca.cycles_ >> 10;
    
    sMeca.apply();

    // Automatic selection of preconditionning method:
    constexpr size_t N_METHODS = 4;
    constexpr size_t N_TESTS = 8;
    constexpr size_t PERIOD = 128;
    
    //automatically select the preconditionning mode:
    //by trying each methods N_STEP steps, adding CPU time and use fastest.
    
    if ( ++autoCounter <= N_TESTS*N_METHODS )
    {
        assert_true(autoPrecond < 6);
        autoCPU[autoPrecond] += cpu;
        autoCNT[autoPrecond] += cnt;

        //std::clog << " precond " << autoPrecond << " cnt " << cnt << " CPU " << cpu << "\n";
        
        if ( autoCounter == N_TESTS*N_METHODS )
        {
            /*
             Compare the performance of all methods, and select the fastest.
             Only adopt a more complicated method if the gain is significant,
             as the simpler one uses less memory.
             */
            autoPrecond = 0;
            for ( size_t m : { 0, 1, 2, 4 } )
            {
                if ( autoCPU[m] < autoCPU[autoPrecond] * 0.95 )
                    autoPrecond = m;
            }
            if ( 1 )
            {
                char str[256], *ptr = str;
                char*const end = str+sizeof(str);
                ptr += snprintf(ptr, end-ptr, " precond selection %lu | method cnt cpu", N_TESTS);
                for ( size_t u : { 0, 1, 2, 4 } )
                    ptr += snprintf(ptr, end-ptr, " | %lu %6.1f %6.0f", u, (real)autoCNT[u]/N_TESTS, autoCPU[u]/N_TESTS);
                ptr += snprintf(ptr, end-ptr, " |  -----> %i", autoPrecond);
                Cytosim::log << str << std::endl;
                if ( prop->verbose )
                    std::clog << str << std::endl;
            }
            for ( size_t u = 0; u < 6; ++u )
            {
                autoCPU[u] = 0;
                autoCNT[u] = 0;
            }
        }
        else
        {
            //alternate betwen { 0, 1, 2, 4 }
            autoPrecond = ( 1 + autoPrecond + (autoPrecond==2) ) % 5;
        }
    }
    else
    {
        if ( autoCounter > PERIOD )
            autoCounter = 0;
    }
}


void Simul::computeForces() const
{
    // we could use here an accessory Meca mec;
    try {
        if ( !ready() )
        {
            prop->complete(*this);
            sMeca.prepare(this);
            setAllInteractions(sMeca);
            sMeca.computeForces();
        }
        else
        {
#if ( 0 )
            /* if the simulation is running live, the force are already available.
            This code is used to check that recalculating gives similar results */
            sMeca.prepare(this);
            setAllInteractions(sMeca);
            fibers.firstID()->printTensions(std::clog);
            sMeca.computeForces();
            fibers.firstID()->printTensions(std::clog);
            std::clog<<"\n";
#endif
        }
    }
    catch ( Exception & e )
    {
        std::cerr << "Error, cytosim could not compute forces:\n";
        std::cerr << "   " << e.message() << '\n';
    }
}


void Simul::flagMecaClusters() const
{
    prop->complete(*this);
    sMeca.prepare(this);
    setAllInteractions(sMeca);
    sMeca.flagClusters();
}


//==============================================================================
//                           EXPERIMENTAL-DEBUG
//==============================================================================
#pragma mark - Experimental


void Simul::addExperimentalInteractions(Meca& meca) const
{
    // ALL THE FORCES BELOW ARE FOR DEVELOPMENT/TESTING PURPOSES:
#if ( 0 )
    LOG_ONCE("AD-HOC FUNKY REPULSIVE FORCE ENABLED\n");
    // add pairwise repulsive force:
    for ( Bead const* i=beads.first(); i ; i=i->next() )
        for ( Bead * j=i->next()    ; j ; j=j->next() )
            meca.addCoulomb(Mecapoint(i,0), Mecapoint(j,0), 0.1);
#endif
#if ( 0 )
    if ( beads.size() > 1 )
    {
        LOG_ONCE("AD-HOC BEAD-STRING FORCES ENABLED\n");
        // attach beads together into a open string:
        Bead * p = beads.firstID();
        if ( p )
        {
            Bead * n = beads.nextID(p);
            while ( p )
            {
                meca.addLongLink(Mecapoint(n,0), Mecapoint(p,0), 1, 1000);
                n = p;
                p = beads.nextID(p);
            }
        }
    }
#endif
#if ( 0 )
    if ( beads.size() > 2 )
    {
        LOG_ONCE("AD-HOC BEAD TORQUES ENABLED\n");
        const real sti = 10000;
        const real ang = 2 * M_PI / 12;
        real co = cos(ang), si = sin(ang);
        // attach beads together in a closed loop:
        Bead * a = beads.firstID();
        Bead * b = beads.nextID(a);
        Bead * c = beads.nextID(b);
        const real len = 2 * a->radius();
        Torque dir = normalize(cross(b->pos()-a->pos(), c->pos()-b->pos()));
        MatrixBlock mat = Meca::torqueMatrix(sti, dir, co, si);
        meca.addTorqueLong(Mecapoint(a,0), Mecapoint(b,0), Mecapoint(c,0), mat, sti, len, sti);
    }
#endif
#if ( 0 )
    LOG_ONCE("AD-HOC FUNKY RADIAL FORCES ENABLED\n");
    // attach beads together in a string:
    for( Bead const* b=beads.first(); b; b=b->next() )
    {
        real x = ( b->identity() - 2 ) * ang;
        Vector pos(5*cos(x), 5*sin(x), 0);
        meca.addPointClamp(Mecapoint(b, 0), pos, 1);
    }
#endif
#if ( 0 )
    LOG_ONCE("AD-HOC CALIBRATED FORCE ENABLED\n");
    // add calibrated forces, for testing rotation
    for ( Fiber const* fib = fibers.first(); fib; fib = fib->next() )
        meca.addTorqueClamp(fib->interpolateCenter(), Vector(0,1,0), 1);
#endif
#if ( 0 )
    LOG_ONCE("AD-HOC CALIBRATED FORCE ENABLED\n");
    // add calibrated force to test rotation of spheres:
    Vector force(0,1,0);
    for ( Sphere const* sph = spheres.first(); sph; sph = sph->next() )
    {
        meca.addForce(Mecapoint(sph, 1), -force);
        meca.addForce(Mecapoint(sph, 2), +force);
    }
#endif
}

//==============================================================================
//                              SOLVE-X 1D
//==============================================================================
#pragma mark -

#include "meca1d.h"

void Simul::solveX()
{
    if ( !pMeca1D )
        pMeca1D = new Meca1D();

    //-----initialize-----
    
    pMeca1D->clear();
    
    for(Fiber * fib = fibers.first(); fib; fib=fib->next())
        pMeca1D->add(fib);

    pMeca1D->prepare(prop->time_step, prop->kT);
    
    //-----set matrix-----

    for ( Couple * c = couples.firstAA(); c ; c=c->next() )
    {
        Hand const* h1 = c->hand1();
        Hand const* h2 = c->hand2();
        
        const size_t i1 = h1->fiber()->matIndex();
        const size_t i2 = h2->fiber()->matIndex();
        assert_true( i1 != i2 );
        
        pMeca1D->addLink(i1, i2, c->prop->stiffness, h2->pos().XX - h1->pos().XX);
    }
    
    for ( Single * s = singles.firstA(); s ; s=s->next() )
    {
        Hand const* h = s->hand();
        const size_t i = h->fiber()->matIndex();
        
        pMeca1D->addClamp(i, s->prop->stiffness, s->position().XX - h->pos().XX);
    }
    
    //-----resolution-----

    real noise = pMeca1D->setRightHandSide(prop->kT);
    
    pMeca1D->solve(prop->tolerance * noise);
    pMeca1D->apply();
}

