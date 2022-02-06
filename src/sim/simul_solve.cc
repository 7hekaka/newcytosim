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

//------------------------------------------------------------------------------
/**
 This will:
 - call setInteractions() for all objects in the system,
 - call addStericInteractions() if prop->steric is true.
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

#if 0
    for ( Event const* e = events.first(); e; e=e->next() )
        e->setInteractions(meca);
#endif
    
    // add steric interactions
    if ( prop->steric_mode && spaces.master() )
    {
        // in the presence of pulling, use the most complete steric engine:
        if ( prop->steric_stiff_pull[0] > 0 )
        {
            if ( !meca.pointGrid.hasGrid() )
                setStericGrid(meca.pointGrid, spaces.master(), prop->steric_max_range, estimateStericRange());
            meca.addStericInteractions(*this);
        }
        else
        {
            if ( !meca.locusGrid.hasGrid() )
                setStericGrid(meca.locusGrid, spaces.master(), prop->steric_max_range, estimateStericRange());
            meca.addStericInteractionsAlt(*this);
        }
    }
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
            for ( size_t n = 0; n < F->nbSegments(); ++n )
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

#pragma mark -

void Simul::solve()
{
    sMeca.pickMecables(*this);
    sMeca.getReady();
    //auto rdt = timer();
    setAllInteractions(sMeca);
    //printf("     ::set      %16llu\n", (timer()-rdt)>>5); rdt = timer();
    sMeca.solve(prop, prop->precondition);
    //printf("     ::solve    %16llu\n", (timer()-rdt)>>5); rdt = timer()
    sMeca.apply();
    //printf("     ::apply    %16llu\n", (timer()-rdt)>>5);
#if ( 0 )
    // check that recalculating gives similar forces
    fibers.firstID()->printTensions(stderr, 47);
    sMeca.calculateForces();
    fibers.firstID()->printTensions(stderr, 92);
    putc('\n', stderr);
#endif
}


/**
 This is attempting to separate the system into subclusters
 that can be solved independently -- should lead to parallelization
 */
void Simul::solve_separate()
{
    size_t cnt = fibers.size();
    ObjectFlag sup = fibers.inventory_.highest();
#if 0
    resetFlags(fibers);
    flagClustersCouples();
#else
    Object ** table = new Object*[sup+2]{nullptr};
    sup = orderClustersCouple(table, sup);
    std::clog << "Ordered " << sup << " clusters\n";
    delete[] table;
#endif
    
    // ready steric interactions
    if ( prop->steric_mode && spaces.master() )
    {
        if ( !sMeca.locusGrid.hasGrid() )
            setStericGrid(sMeca.locusGrid, spaces.master(), prop->steric_max_range, estimateStericRange());
    }
    //std::clog << "Separating " << cnt << " fibers\n";
    for ( ObjectFlag f = 0; f <= sup; ++f )
    {
        sMeca.mecables.clear();
        for ( Fiber * F=fibers.first(); F; F=F->next() )
            if ( F->flag() == f )
                sMeca.addMecable(F);
        size_t num = sMeca.mecables.size();
        if ( num > 0 )
        {
            std::clog << "  cluster " << f << " has " << num << " fibers\n";
            cnt -= num;
            sMeca.getReady();
            sMeca.setSomeInteractions();
            if ( prop->steric_mode && spaces.master() )
                sMeca.addSomeStericInteractions(prop->steric_stiff_push[0]);
            sMeca.solve(prop, prop->precondition);
            sMeca.apply();
        }
    }
    assert_true(cnt == 0);
}


void Simul::solve_force()
{
    sMeca.pickMecables(*this);
    sMeca.getReady();
    setAllInteractions(sMeca);
    sMeca.calculateForces();
}


void Simul::solve_half()
{
    sMeca.pickMecables(*this);
    sMeca.getReady();
    setAllInteractions(sMeca);
    sMeca.solve(prop, prop->precondition);
}



/**
 Solve the system, and automatically select the fastest preconditionning method
 */
void Simul::solve_auto()
{
    sMeca.pickMecables(*this);
    sMeca.getReady();
    setAllInteractions(sMeca);
    
    // solve the system, recording time:
    //double cpu = TimeDate::milliseconds();
    size_t cnt = sMeca.solve(prop, autoPrecond);
    //float cpu = TimeDate::milliseconds() - cpu;
    // use Meca::cycles_ that only includes preconditionning parts!
    float cpu = sMeca.cycles_;
    
    sMeca.apply();

    // Automatic selection of preconditionning method:
    constexpr size_t N_METHODS = 4;
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
             Compare the performance of some methods, and select the fastest.
             Only adopt a more complicated method if the gain is significant.
             */
            autoPrecond = 0;
            for ( unsigned m : { 0, 1, 4, 6 } )
            {
                if ( autoCPU[m] < autoCPU[autoPrecond] * 0.95 )
                    autoPrecond = m;
            }
            if ( 1 )
            {
                char str[256], *ptr = str;
                char*const end = str+sizeof(str);
                ptr += snprintf(ptr, end-ptr, " precond selection %lu | method count cpu", N_TESTS);
                for ( size_t u : { 0, 1, 4, 6 } )
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
            //rotate betwen { 0, 1, 4, 6 }
            autoPrecond = ( 1 + autoPrecond + 2*(autoPrecond==1) + (autoPrecond==4) ) % 7;
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
            sMeca.pickMecables(*this);
            sMeca.getReady();
            setAllInteractions(sMeca);
            sMeca.calculateForces();
        }
    }
    catch ( Exception & e )
    {
        std::cerr << "Error, Cytosim could not compute forces:\n";
        std::cerr << "   " << e.message() << '\n';
    }
}

//==============================================================================
//                              SOLVE-X 1D
//==============================================================================

#include "meca1d.h"

void Simul::solve_onlyX()
{
    if ( !pMeca1D )
        pMeca1D = new Meca1D();

    //-----initialize-----

    pMeca1D->pickMecables(*this);
    pMeca1D->getReady(prop->time_step, prop->kT);

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

//------------------------------------------------------------------------------
#pragma mark - Analysis


void Simul::flagClustersMeca() const
{
    prop->complete(*this);
    sMeca.pickMecables(*this);
    sMeca.getReady();
    setAllInteractions(sMeca);
    sMeca.flagClusters();
}


// Join lists: f -> f + g;  g -> null
static void join(Object ** table, ObjectFlag f, ObjectFlag g)
{
    assert_true( f < g );
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


static size_t depth(Object * ptr)
{
    size_t cnt = 0;
    while ( ptr )
    {
        cnt += static_cast<Mecable*>(ptr)->nbPoints();
        ptr = ptr->next();
    }
    return cnt;
}


ObjectFlag Simul::orderClustersCouple(Object ** table, ObjectFlag sup)
{
    Object * F = fibers.first();
    Object * G;
    while ( F )
    {
        Object * X = F->next();
        F->matchFlagIdentity();
        table[F->flag()] = F;
        F->next(nullptr);
        F = X;
    }
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
    const size_t target = 256;
    // pool smaller clusters together
    ObjectFlag s = 1;
    size_t cnt = depth(table[s]);
    for ( ObjectFlag f = s+1; f <= sup; ++f )
    {
        if ( table[f] )
        {
            size_t d = depth(table[f]);
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
    G = nullptr;
    ObjectFlag num = 0;
    for ( ObjectFlag f = 0; f <= sup; ++f )
    {
        F = table[f];
        if ( F )
        {
            ++num;
            F->flag(num);
            F->prev(G);
            if ( G )
                G->next(F);
            else
                fibers.pool_.front(F);
            G = F;
            F = F->next();
            while ( F )
            {
                F->flag(num);
                F->prev(G);
                G->next(F);
                G = F;
                F = F->next();
            }
        }
    }
    fibers.pool_.back(G);
    return num;
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

