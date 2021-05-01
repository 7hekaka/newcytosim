// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

/**
 Calculate the minimum grid cell size, given the segmentation of Fibers that
 is within a certain range because their length may vary, and taking into
 account the radius of Bead, Sphere and Solid.
 
 This function is used to estimate SimulProp::steric_max_range, when it is
 not specified by the user.
 
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
    for ( Fiber const* F=fibers.first(); F; F=F->next() )
    {
        if ( F->prop->steric )
            len = std::max(len, F->segmentation());
    }

    /*
     The interaction can be aligned with the fiber, and we must add the distances:
     2 * range if two fibers of radius 'range' interact.
     + 2 * ( len / 2 ) since len/2 is the distance between the center of the segment
     and its most distal point.
     */
    ran = len + 2 * ran;
    
    
    for ( Sphere const* O=spheres.first(); O; O=O->next() )
    {
        if ( O->prop->steric )
            ran = std::max(ran, 2 * O->radius() + O->prop->steric_range);
    }
    
    for ( Bead const* B=beads.first(); B; B=B->next() )
    {
        if ( B->prop->steric )
            ran = std::max(ran, 2 * B->radius() + B->prop->steric_range);
    }
    
    for ( Solid const* S=solids.first(); S; S=S->next() )
    {
        if ( S->prop->steric )
        {
            for ( size_t p = 0; p < S->nbPoints(); ++p )
                ran = std::max(ran, 2 * S->radius(p) + S->prop->steric_range);
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
    for ( Fiber const* F=fibers.first(); F; F=F->next() )
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
    for ( Sphere const* O=spheres.first(); O; O=O->next() )
    {
        if ( O->prop->steric )
#if ( NUM_STERIC_PANES == 1 )
            pointGrid.add(O, 0, O->radius(), O->radius()+O->prop->steric_range);
#else
            pointGrid.add(O->prop->steric, O, 0, O->radius(), O->radius()+O->prop->steric_range);
#endif
    }
    
    // include Beads
    for ( Bead const* B=beads.first(); B; B=B->next() )
    {
        if ( B->prop->steric )
#if ( NUM_STERIC_PANES == 1 )
            pointGrid.add(B, 0, B->radius(), B->radius()+B->prop->steric_range);
#else
            pointGrid.add(B->prop->steric, B, 0, B->radius(), B->radius()+B->prop->steric_range);
#endif
    }
        
    // include Points that have a radius from Solids
    for ( Solid const* S=solids.first(); S; S=S->next() )
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
    Stiffness pam(prop->steric_stiff_push[0], prop->steric_stiff_pull[0]);
    
#if ( NUM_STERIC_PANES == 1 )
    
    pointGrid.setInteractions(meca, pam);

#elif ( NUM_STERIC_PANES == 2 )
    
    // add steric interactions inside pane 1:
    pointGrid.setInteractions(meca, pam, 1);
    // add steric interactions between panes 1 and 2:
    pointGrid.setInteractions(meca, pam, 1, 2);
    //pointGrid.setInteractions(meca, pam, 2, 1);

#else
    
    // add steric interactions within each pane:
    for ( size_t p = 1; p <= NUM_STERIC_PANES; ++p )
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
void Simul::setStericInteractionsAlt(Meca& meca) const
{
    if ( !locusGrid.hasGrid() )
    {
        if (!spaces.master())
            return;
        setStericGrid(locusGrid, spaces.master());
    }

    // clear grid
    locusGrid.clear();
    
    // distribute Fiber-points on the grid
    for ( Fiber const* F=fibers.first(); F; F=F->next() )
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
    locusGrid.mark();
    
    // include Spheres
    for ( Sphere const* O=spheres.first(); O; O=O->next() )
    {
        if ( O->prop->steric )
#if ( MAX_STERIC_PANES == 1 )
            locusGrid.add(O, 0, O->radius());
#else
            locusGrid.add(O->prop->steric, O, 0, O->radius());
#endif
    }
    
    // include Beads
    for ( Bead const* B=beads.first(); B; B=B->next() )
    {
        if ( B->prop->steric )
#if ( MAX_STERIC_PANES == 1 )
            locusGrid.add(B, 0, B->radius());
#else
            locusGrid.add(B->prop->steric, B, 0, B->radius());
#endif
    }
        
    // include Points that have a radius from Solids
    for ( Solid const* S=solids.first(); S; S=S->next() )
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
    real stiff = prop->steric_stiff_push[0];

#if ( MAX_STERIC_PANES == 1 )
        
    locusGrid.setInteractions(meca, stiff);

#elif ( MAX_STERIC_PANES == 2 )
    
    // add steric interactions inside pane 1:
    locusGrid.setInteractions(meca, stiff, 1);
    // add steric interactions between panes 1 and 2:
    locusGrid.setInteractions(meca, stiff, 1, 2);
    //pointGrid.setInteractions(meca, pam, 2, 1);

#else
    
    // add steric interactions within each pane:
    for ( size_t p = 1; p <= MAX_STERIC_PANES; ++p )
        locusGrid.setInteractions(meca, stiff, p);

#endif
    //std::clog << "LocusGrid has capacity " << locusGrid.capacity() << "\n";
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
    if ( prop->steric_mode )
    {
        // in the presence of pulling, use the most complete steric engine:
        if ( prop->steric_stiff_pull[0] > 0 )
            setStericInteractions(meca);
        else
            setStericInteractionsAlt(meca);
    }
    //addExperimentalInteractions(meca);

#if ( 0 )
    // add steric interaction between the first Sphere and all Fibers
    Sphere * sol = spheres.firstID();
    if ( sol && sol->prop->steric )
    {
        const Vector cen = sol->posPoint(0);
        const real rad = sol->radius();
        const real rad2 = square(rad);
        const real stiff = prop->steric_stiff_push[0];

        for ( Fiber const* F = fibers.first(); F; F = F->next() )
        {
            for ( size_t n = 0; n < F->nbSegments(); ++n )
            {
                FiberSegment seg(F, n);
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
#if ( 0 )
    // check that recalculating gives similar forces
    fibers.firstID()->printTensions(stderr, 47);
    sMeca.computeForces();
    fibers.firstID()->printTensions(stderr, 92);
    putc('\n', stderr);
#endif
}


void Simul::solve_force()
{
    sMeca.prepare(this);
    setAllInteractions(sMeca);
    sMeca.computeForces();
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
    float cpu = sMeca.cycles_ >> 10;
    
    sMeca.apply();

    // Automatic selection of preconditionning method:
    constexpr size_t N_METHODS = 3;
    constexpr size_t N_TESTS = 8;
    constexpr size_t PERIOD = 32;
    
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
             Only adopt a more complicated method if the gain is significant.
             */
            autoPrecond = 0;
            for ( unsigned m : { 0, 1, 4 } )
            {
                if ( autoCPU[m] < autoCPU[autoPrecond] * 0.95 )
                    autoPrecond = m;
            }
            if ( 1 )
            {
                char str[256], *ptr = str;
                char*const end = str+sizeof(str);
                ptr += snprintf(ptr, end-ptr, " precond selection %lu | method count cpu", N_TESTS);
                for ( size_t u : { 0, 1, 4 } )
                    ptr += snprintf(ptr, end-ptr, " | %lu %6.1f %6.0f", u, (real)autoCNT[u]/N_TESTS, autoCPU[u]/N_TESTS);
                snprintf(ptr, end-ptr, " |  -----> %i", autoPrecond);
                Cytosim::log << str << '\n';
            }
            for ( size_t u = 0; u < 8; ++u )
            {
                autoCPU[u] = 0;
                autoCNT[u] = 0;
            }
        }
        else
        {
            //alternate betwen { 0, 1, 4 }
            autoPrecond = ( 1 + autoPrecond + 2*(autoPrecond==1) ) % 5;
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
    try {
        // if the simulation is running live, the force should be available.
        if ( !primed() )
        {
            // we could use here a different Meca for safety
            prop->complete(*this);
            sMeca.prepare(this);
            setAllInteractions(sMeca);
            sMeca.computeForces();
        }
    }
    catch ( Exception & e )
    {
        std::cerr << "Error, Cytosim could not compute forces:\n";
        std::cerr << "   " << e.message() << '\n';
    }
}


void Simul::flagClustersMeca() const
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
    LOG_ONCE("AD-HOC BEAD-STRING FORCES ENABLED\n");
    // attach beads together into an open/closed string:
    const real stiff = 1000;
    Bead * B = beads.firstID();
    if ( B )
    {
        Bead * N = beads.nextID(B);
        while ( N )
        {
            meca.addLongLink(Mecapoint(N,0), Mecapoint(B,0), 1, stiff);
            B = N;
            N = beads.nextID(N);
        }
        N = beads.firstID();
        meca.addLongLink(Mecapoint(N,0), Mecapoint(B,0), 1, stiff);
    }
#endif
#if ( 0 )
    if ( beads.size() > 2 )
    {
        LOG_ONCE("AD-HOC BEAD TORQUES ENABLED\n");
        const real sti = 10000;
        const real angle = 2 * M_PI / 12;
        Vector2 ang(std::cos(angle), std::sin(angle));
        // attach beads together in a closed loop:
        Bead * a = beads.firstID();
        Bead * b = beads.nextID(a);
        Bead * c = beads.nextID(b);
        const real len = 2 * a->radius();
        Torque dir = normalize(cross(b->pos()-a->pos(), c->pos()-b->pos()));
        MatrixBlock mat = Meca::torqueMatrix(sti, dir, ang);
        meca.addTorqueLong(Mecapoint(a,0), Mecapoint(b,0), Mecapoint(c,0), mat, sti, len, sti);
    }
#endif
#if ( 0 )
    LOG_ONCE("AD-HOC BEAD CLAMPS ENABLED\n");
    // attach beads to fixed positions on a circle:
    const real ang = M_PI / beads.size();
    const real rad = 5;
    for( Bead const* B=beads.first(); B; B=B->next() )
    {
        real x = B->identity() * ang;
        Vector pos(rad*std::cos(x), rad*std::sin(x), 0);
        meca.addPointClamp(Mecapoint(B, 0), pos, 1);
    }
#endif
#if ( 0 )
    LOG_ONCE("AD-HOC CALIBRATED FORCE ENABLED\n");
    // add calibrated forces, to test rotation under known torque
    for ( Fiber const* F=fibers.first(); F; F=F->next() )
        meca.addTorqueClamp(F->interpolateCenter(), Vector(0,1,0), 1);
#endif
#if ( 0 )
    LOG_ONCE("AD-HOC CALIBRATED FORCE ENABLED\n");
    // add calibrated force to test rotation of spheres:
    Vector force(0,1,0);
    for ( Sphere const* O=spheres.first(); O; O=O->next() )
    {
        meca.addForce(Mecapoint(O, 1), -force);
        meca.addForce(Mecapoint(O, 2), +force);
    }
#endif
}

//==============================================================================
//                              SOLVE-X 1D
//==============================================================================
#pragma mark -

#include "meca1d.h"

void Simul::solve_onlyX()
{
    if ( !pMeca1D )
        pMeca1D = new Meca1D();

    //-----initialize-----

    pMeca1D->prepare(this, prop->time_step, prop->kT);
    
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

