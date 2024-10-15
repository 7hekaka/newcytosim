// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University.

/**
 Estimate the minimum value of the `grid cell size`, from the segmentation
 parameter of fibers, and the radius of Bead, Sphere and Solid.
 
 This function is used to set SimulProp::steric_max_range, when this parameter
 is not specified by the user.
 
 It assumes that Fiber::adjustSegmentation() is used, such that at any time:
     actual_segmentation <  4/3 * FiberProp::segmentation
 
 In the particular case where the fiber have a fixed length, the real segmentation
 can be used, and it is recommended in that case to set steric_max_range in the config.
 */
real Simul::minimumStericRange() const
{
    /*
     For a spherical object, the steric range must be at least 2 * radius
     */
    real dia = 0;
    
    for ( Sphere const* O=spheres.first(); O; O=O->next() )
    {
        if ( O->prop->steric_key )
            dia = max_real(dia, 2 * O->radius() + O->prop->steric_range);
    }
    
    for ( Bead const* B=beads.first(); B; B=B->next() )
    {
        if ( B->prop->steric_key )
            dia = max_real(dia, 2 * B->radius() + B->prop->steric_range);
    }
    
    for ( Solid const* S=solids.first(); S; S=S->next() )
    {
        if ( S->prop->steric_key )
        {
            for ( index_t p = 0; p < S->nbPoints(); ++p )
                dia = max_real(dia, 2 * S->radius(p) + S->prop->steric_range);
        }
    }
    
    // check all FiberProp with enabled steric:
    /*
     Given that two fibers can be aligned end-on, we must add the distances:
         `2 * steric_radius + steric_range + 2 * ( len / 2 )`
     since `len/2` is the distance between the center of a segment and
     any point on the segment that may be considered for interaction.
     */
    real ran = 0;
    real len = 0;
    for ( Property const* i : properties.find_all("fiber") )
    {
        FiberProp const* fp = static_cast<FiberProp const*>(i);
        if ( fp->steric_key )
        {
            // The maximum length of a segment is 4/3 * segmentation
            len = max_real(len, (real)(1.34) * fp->segmentation);
            
            // check extended range of interaction
            ran = max_real(ran, 2 * fp->steric_radius + fp->steric_range);
        }
    }
    
    // for safety, check the actual segmentation of all fibers:
    for ( Fiber const* F=fibers.first(); F; F=F->next() )
    {
        if ( F->prop->steric_key )
            len = max_real(len, F->segmentation());
    }
    
    // finally compute the minimum range, all objects considered:
    ran = std::max(dia, ran+len);
    
    if ( ran < REAL_EPSILON )
        Cytosim::warn("could not estimate simul:steric_max_range automatically!\n");
    
    //std::clog << "simul:minimumStericRange = " << ran << "\n";
    return ran;
}

//------------------------------------------------------------------------------
/**
 This will:
 - call setInteractions() for all objects in the system,
 - call addStericInteractions() if simul:steric is true.
 .
 */
void Simul::setAllInteractions(Meca& meca) const
{
#if 1
    for ( Mecable const* mec : meca.mecables )
        mec->setInteractions(meca);
#else
    for ( Fiber const* f=fibers.first(); f ; f=f->next() )
        f->setInteractions(meca);
    
    for ( Solid const* d=solids.first(); d ; d=d->next() )
        d->setInteractions(meca);
    
    for ( Sphere const* o=spheres.first(); o ; o=o->next() )
        o->setInteractions(meca);
    
    for ( Bead const* b=beads.first(); b ; b=b->next() )
        b->setInteractions(meca);
#endif

    for ( Space const* e=spaces.first(); e; e=e->next() )
        e->setInteractions(meca, *this);

    for ( Single const* s=singles.firstA(); s ; s=s->next() )
        s->setInteractions(meca);

    for ( Couple const* c=couples.firstAA(); c ; c=c->next() )
        c->setInteractions(meca);
    
    for ( Organizer const* x = organizers.first(); x; x=x->next() )
        x->setInteractions(meca);

    
    // add steric interactions
    if ( meca.steric_ == 2 )
        Meca::addStericInteractions(meca.locusGrid, *this);
    else if ( meca.steric_ == 1 )
        Meca::addStericInteractions(meca.pointGrid, *this);

    //addExperimentalInteractions(meca);

#if ( 0 )
    /*
     Add simplified steric interactions between the first Sphere and all Fibers
     This is not necessarily equivalent to the steric engine, since we do not add
     the 'radius' of the fiber, but it can be faster eg. if there is only one
     sphere in the system. The code can easily be adapted to handle Beads
     */
    Sphere * S = spheres.firstID();
    if ( S && S->prop->steric )
    {
        LOG_ONCE("Limited steric interactions with first Sphere enabled!");
        const real stiff = prop->steric_stiff_push[0];
        const Vector cen = S->posPoint(0);
        const real rad = S->radius();
        const real rad2 = square(rad);

        for ( Fiber const* F = fibers.first(); F; F = F->next() )
        {
            for ( index_t n = 0; n < F->nbSegments(); ++n )
            {
                FiberSegment seg(F, n);
                real dis = INFINITY;
                real abs = seg.projectPoint(cen, dis);
                if ( dis < rad2 )
                    meca.addSideSlidingLink(seg, abs, Mecapoint(S, 0), rad, stiff);
            }
        }
    }
#endif
}

//------------------------------------------------------------------------------
#pragma mark - Solve

void Simul::solve_meca()
{
#if ( 0 )
    ObjectFlag sup = fibers.inventory_.highest();
    Object ** table = new Object*[sup+2]{nullptr};
    ObjectFlag cnt = orderClustersCouple(table, sup);
    std::clog << "Ordered " << cnt << " clusters (" << sup << ")\n";
    delete[] table;
#endif
    sMeca.getReady(*this);
    //auto rdt = timer();
    setAllInteractions(sMeca);
    //printf("     ::set      %16llu\n", (timer()-rdt)>>5); rdt = timer();
    sMeca.solve();
    //printf("     ::solve    %16llu\n", (timer()-rdt)>>5); rdt = timer()
    sMeca.apply();
    //printf("     ::apply    %16llu\n", (timer()-rdt)>>5);
#if ( 0 )
    Fiber * fib = fibers.firstID();
    if ( fib ) {
        // check that recalculating gives similar forces
        fib->printTensions(stderr, 47);
        sMeca.calculateForces();
        fib->printTensions(stderr, 92);
        putc('\n', stderr);
    }
#endif
}


void Simul::prepared_solve()
{
    setAllInteractions(sMeca);
    sMeca.solve();
    sMeca.apply();
}


void Simul::solve_force()
{
    sMeca.getReady(*this);
    setAllInteractions(sMeca);
    sMeca.calculateForces();
}


void Simul::solve_half()
{
    sMeca.getReady(*this);
    setAllInteractions(sMeca);
    sMeca.solve();
}


/**
 Solve the system, and automatically select the fastest preconditionning method
 */
void Simul::solve_auto()
{
    sMeca.getReady(*this);
    sMeca.precond_ = autoPrecond;
    setAllInteractions(sMeca);
    
    // solve the system, recording time:
    //double cpu = TimeDate::milliseconds();
    unsigned cnt = sMeca.solve();
    //float cpu = TimeDate::milliseconds() - cpu;
    // use Meca::cycles_ that only includes preconditionning parts!
    float cpu = sMeca.cycles_;
    
    sMeca.apply();
    
    // list of preconditionning method to be tried:
    std::initializer_list<unsigned> methods { 0, 1, 6 };
    
    // how many timestep accumulated for each method:
    constexpr index_t N_TESTS = 8;
    // number of timestep for next trial series:
    constexpr index_t PERIOD = 32;
    
    //automatically select the preconditionning mode:
    //by trying each methods N_STEP steps, adding CPU time and use fastest.
    
    if ( ++autoCounter <= N_TESTS * methods.size() )
    {
        assert_true(autoPrecond < 8);
        autoCPU[autoPrecond] += cpu;
        autoCNT[autoPrecond] += cnt;

        //std::clog << " precond " << autoPrecond << " cnt " << cnt << " CPU " << cpu << "\n";
        
        if ( autoCounter == N_TESTS * methods.size() )
        {
            unsigned z = * methods.begin();
            for ( unsigned m : methods )
            {
                /*
                 Compare the performance of some methods, and select the fastest.
                 The number of iterations should be decreased, with some CPU gain.
                 Only adopt a more complicated method if the gain is significant.
                 */
                if ( autoCNT[m] < autoCNT[z] && autoCPU[m] < autoCPU[z] )
                    z = m;
            }
            autoPrecond = z;
            if ( 1 )
            {
                char str[256], *ptr = str;
                char*const end = str+sizeof(str);
                ptr += snprintf(ptr, end-ptr, " precond selection %u | method count cpu", N_TESTS);
                for ( unsigned u : methods )
                {
                    real N = (real)autoCNT[u] / N_TESTS;
                    real T = autoCPU[u] / N_TESTS;
                    ptr += snprintf(ptr, end-ptr, " | %u %6.1f %6.0f", u, N, T);
                }
                snprintf(ptr, end-ptr, " |  -----> %i", autoPrecond);
                Cytosim::log(str, '\n');
            }
            for ( unsigned u = 0; u < 8; ++u )
            {
                autoCPU[u] = 0;
                autoCNT[u] = 0;
            }
        }
        else
        {
            // rotate betwen methods
            auto i = std::find(methods.begin(), methods.end(), autoPrecond);
            if ( i == methods.end() || ++i == methods.end() )
                i = methods.begin();
            autoPrecond = *i;
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
    // if the simulation is running live, the force should be available.
    if ( !primed() )
    {
        primed_ = 1;
        try {
            // we could use here a different Meca for safety
            prop.complete(*this);
            sMeca.getReady(*this);
            setAllInteractions(sMeca);
            sMeca.calculateForces();
        }
        catch ( Exception & e )
        {
            std::cerr << "Error, Cytosim could not compute forces: ";
            std::cerr << e.message() << '\n';
        }
        primed_ = 0;
    }
}


//==============================================================================
//                              UNIAXIAL SOLVE
//==============================================================================

#include "meca1d.h"

void Simul::solve_uniaxial()
{
    if ( !pMeca1D )
        pMeca1D = new Meca1D();

    //-----initialize-----

    pMeca1D->getReady(*this);

    //-----set matrix-----

    for ( Couple * c = couples.firstAA(); c ; c=c->next() )
    {
        Hand const* h1 = c->hand1();
        Hand const* h2 = c->hand2();
        
        const index_t i1 = h1->fiber()->matIndex();
        const index_t i2 = h2->fiber()->matIndex();
        assert_true( i1 != i2 );
        
        pMeca1D->addLink(i1, i2, c->prop->stiffness, h2->pos().XX - h1->pos().XX);
    }
    
    for ( Single * s = singles.firstA(); s ; s=s->next() )
    {
        Hand const* h = s->hand();
        const index_t i = h->fiber()->matIndex();
        
        pMeca1D->addClamp(i, s->prop->stiffness, s->position().XX - h->pos().XX);
    }
    
    //-----resolution-----

    real noise = pMeca1D->setRightHandSide(prop.kT);
    pMeca1D->solve(prop.tolerance * noise);
    pMeca1D->apply();
}


//------------------------------------------------------------------------------
#pragma mark - Analysis

void Simul::flagClustersMeca() const
{
    prop.complete(*this);
    sMeca.getReady(*this);
    setAllInteractions(sMeca);
    sMeca.flagClusters();
}


// Join lists: f -> f + g;  g -> null
static void join(Object ** table, ObjectFlag f, ObjectFlag g)
{
    assert_true( f < g );
    assert_true( table[g] );
    Object * F = nullptr;
    Object * G = table[g];
    do {
        F = G;
        G->flag(f);
        G = G->next();
    } while ( G );
    F->next(table[f]);
    table[f] = table[g];
    table[g] = nullptr;
}


/// number of points in list starting at 'ptr' and defined by 'next()'
static index_t depth(Object * ptr)
{
    index_t cnt = 0;
    while ( ptr )
    {
        cnt += static_cast<Mecable*>(ptr)->nbPoints();
        ptr = ptr->next();
    }
    return cnt;
}


ObjectFlag Simul::orderClustersCouple(Object ** table, ObjectFlag sup)
{
    // attribute unique flag to all Fibers:
    Object * F = fibers.first();
    ObjectFlag flg = 0;
    while ( F )
    {
        F->flag(++flg);
        table[flg] = F;
        Object * X = F->next();
        F->next(nullptr);
        F = X;
    }
    assert_true( flg <= sup );
    // join subsets that are connected by a Couple:
    for ( Couple const* C=couples.firstAA(); C ; C=C->next() )
    {
        ObjectFlag f = C->fiber1()->flag();
        ObjectFlag g = C->fiber2()->flag();
        if ( f < g )
            join(table, f, g);
        else if ( g < f )
            join(table, g, f);
    }
#if 1
    // pool smaller clusters together to reach 'target'
    const index_t target = 256;
    ObjectFlag s = 1;
    index_t cnt = depth(table[s]);
    for ( ObjectFlag f = s+1; f <= sup; ++f )
    {
        if ( table[f] )
        {
            index_t d = depth(table[f]);
            if ( cnt + d < target )
            {
                join(table, s, f);
                cnt += d;
            }
            else if ( d < target )
            {
                s = f;
                cnt = d;
            }
        }
    }
#endif
    // put all objects back in list
    Object * G = nullptr;
    flg = 0;
    for ( ObjectFlag f = 0; f <= sup; ++f )
    {
        F = table[f];
        if ( F )
        {
            ++flg;
            F->flag(flg);
            F->prev(G);
            if ( G )
                G->next(F);
            else
                fibers.pool_.front(F);
            G = F;
            F = F->next();
            while ( F )
            {
                F->flag(flg);
                F->prev(G);
                G->next(F);
                G = F;
                F = F->next();
            }
        }
    }
    fibers.pool_.back(G);
    return flg;
}

//==============================================================================
//                           EXPERIMENTAL-DEBUG
//==============================================================================
#pragma mark - Experimental


void Simul::addExperimentalInteractions(Meca& meca) const
{
    // ALL THE FORCES BELOW ARE FOR DEVELOPMENT/TESTING PURPOSES:
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
        meca.addForce(O, 1, -force);
        meca.addForce(O, 2, +force);
    }
#endif
}
