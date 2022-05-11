// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University

#include <cstdio>
#include <time.h>
#include <unistd.h>

#include "sim_thread.h"
#include "exceptions.h"
#include "print_color.h"
#include "picket.h"

//------------------------------------------------------------------------------

/**
 This uses a Parser that cannot write to disc.
 */
SimThread::SimThread(Simul* sim)
: Parser(sim, 1, 1, 1, 1, 0)
{
    status_ = -2;
    order_  = 0;
    hold_   = 0;
    holding_= 0;
    period_ = 1;
    pthread_mutex_init(&mutex_, nullptr);
    pthread_cond_init(&condition_, nullptr);
}

/**
 Possible issue:
 When quitting the application, this destructor might be called after the
 destructor of Simul(), in which case it will access non-existent data,
 most likely causing a crash().
 */
SimThread::~SimThread()
{
    //std::cerr << "~SimThread()\n";
    stop();
    pthread_cond_destroy(&condition_);
    pthread_mutex_destroy(&mutex_);
}

//------------------------------------------------------------------------------
#pragma mark - Process control


void SimThread::debug(const char* msg) const
{
    if ( isWorker() )
        fprintf(stdout, "\n- - -  %-16s", msg);
    else
        fprintf(stdout, "\n* * *  %-16s", msg);
}


void SimThread::gubed(const char* msg) const
{
    if ( isWorker() )
        fprintf(stdout, "  - - %12s ", msg);
    else
        fprintf(stdout, "  * * %12s ", msg);
}


void SimThread::hold()
{
    assert_true( isWorker() );

    if ( order_ < 0 )
        pthread_exit(nullptr);
    
    if ( ++hold_ >= period_ )
    {
        hold_ = 0;
        //debug("holding");
        holding_ = 1;
        cond_wait();  // this also unlocks and locks the mutex
    }
}


//------------------------------------------------------------------------------
#pragma mark - Lauching threads


void SimThread::run()
{
    assert_true( isWorker() );
    try {
        do {
            order_ = 0;
            Parser::readConfig();
            do {
                usleep(25000);
                holding_ = 2;
                cond_wait();  // this also unlocks and locks the mutex
            } while ( !order_ );
            erase_simul(1);
        } while ( order_ >= 0 );
    }
    catch( Exception & e ) {
        std::cerr << e.brief() << e.info() << '\n';
        //flashText("Error: the simulation died");
    }
}


/** C-style function to cleanup after thread has terminated */
void child_cleanup(void * arg)
{
    SimThread * st = static_cast<SimThread*>(arg);
    assert_true( st->isWorker() );
    st->status_ = -1;
    st->holding_ = 0;
    st->unlock();
    //st->debug("ended");
}


/** C-style function to start a new thread */
void* run_launcher(void * arg)
{
    //std::clog << "slave  " << pthread_self() << '\n';
    SimThread * st = static_cast<SimThread*>(arg);
    st->lock();
    // let the system cleanup upon normal termination:
    pthread_cleanup_push(child_cleanup, arg);
    st->run();
    pthread_cleanup_pop(1);
    return nullptr;
}


/**
 This attempts to start the live simulation by
 calling `run()` in a slave thread
 */
void SimThread::start()
{
    assert_false( isWorker() );
    if ( status_ )
    {
        order_ = 0;
        assert_true(holding_ == 0);
        Parser::restart_ = 0;
        //std::clog << "master " << pthread_self() << '\n';
        status_ = pthread_create(&child_, nullptr, run_launcher, this);
        // let the system cleanup upon normal child termination:
        if ( status_ == 0 )
            pthread_detach(child_);
        else
            printf("%p failed to create thread", this);
    }
}

//------------------------------------------------------------------------------

void SimThread::extend_run(size_t n_steps)
{
    assert_true( isWorker() );
    try {
        sim_->parser(this);
        Parser::execute_run(n_steps);
    }
    catch( Exception & e ) {
        std::cerr << e.brief() << e.info() << '\n';
        //flashText("Error: %s", e.what());
    }
    sim_->parser(nullptr);
}


/** C-style function to start a new thread */
void* extend_launcher(void * arg)
{
    //std::clog << "slave  " << pthread_self() << '\n';
    SimThread * st = static_cast<SimThread*>(arg);
    st->lock();
    pthread_cleanup_push(child_cleanup, arg);
    st->extend_run(1<<20);
    pthread_cleanup_pop(1);
    return nullptr;
}


/// call `extend_run()` in a slave thread
int SimThread::extend()
{
    assert_false( isWorker() );
    if ( status_ )
    {
        order_ = 0;
        holding_ = 0;
        //std::clog << "master " << pthread_self() << '\n';
        status_ = pthread_create(&child_, nullptr, extend_launcher, this);
    }
    return status_;
}

//------------------------------------------------------------------------------
#pragma mark - Thread control & termination

/**
 ask the slave thread to exit at the next spontaneous halt
*/
void SimThread::stop()
{
    assert_false( isWorker() );
    if ( status_ == 0 )
    {
        // request clean termination:
        order_ = -1;
        signal();
        // wait for termination:
        //debug("join...");
        if ( 0 == pthread_join(child_, nullptr) )
            status_ = -2;
    }
}

/**
 kill the slave thread immediately
 */
void SimThread::cancel_join()
{
    assert_false( isWorker() );
    if ( status_ == 0 )
    {
        //debug("cancel...");
        // force termination:
        if ( 0 == pthread_cancel(child_) )
        {
            // wait for termination to reclaim resources:
            if ( 0 == pthread_join(child_, nullptr) )
                status_ = -2;
        }
    }
}


//------------------------------------------------------------------------------
#pragma mark - Loading

int SimThread::loadFrame(size_t f)
{
    int r = 7;
    lock();
    try {
        r = reader_.loadFrame(*sim_, f);
    }
    catch( Exception & e )
    {
        print_blue(std::cerr, e.brief());
        std::cerr << e.info() << " (loading frame)\n";
    }
    unlock();
    return r;
}

int SimThread::loadNextFrame()
{
    int r = 7;
    lock();
    try {
        r = reader_.loadNextFrame(*sim_);
    }
    catch( Exception & e )
    {
        print_blue(std::cerr, e.brief());
        std::cerr << e.info() << " (loading next frame)\n";
    }
    unlock();
    return r;
}

int SimThread::loadLastFrame()
{
    int r = 7;
    lock();
    try {
        r = reader_.loadLastFrame(*sim_);
    }
    catch( Exception & e )
    {
        print_blue(std::cerr, e.brief());
        std::cerr << e.info() << " (loading last frame)\n";
    }
    unlock();
    return r;
}

//------------------------------------------------------------------------------
#pragma mark - Mouse-controlled Single


SingleProp * SimThread::getHandleProperty() const
{
    Property * p = sim_->properties.find("single", "user_single");
    return static_cast<SingleProp*>(p);
}


SingleProp * SimThread::makeHandleProperty(real range)
{
    // Create a Hand that attaches fast and never detach:
    HandProp * hap = new HandProp("user_hand");
    hap->binding_range   = range;
    hap->binding_rate    = 1 / sim_->time_step();
    hap->unbinding_rate  = 0;
    hap->unbinding_force = INFINITY;
    hap->complete(*sim_);
    sim_->properties.deposit(hap);

    SingleProp * sip = new SingleProp("user_single");
    sip->hand = "user_hand";
    sip->stiffness = 2000;
    sip->complete(*sim_);
    sim_->properties.deposit(sip);
    
    return sip;
}


Single * SimThread::createHandle(Vector const& pos, real range)
{
    SingleProp * sip = getHandleProperty();
    if ( !sip )
        sip = makeHandleProperty(range);
    Single * res = new Picket(sip, pos);
    sim_->singles.add(res);
    handle_ = res;
    return res;
}


ObjectList SimThread::allHandles(SingleProp const* sip) const
{
    return sim_->singles.collect(match_property, sip);
}


bool SimThread::selectClosestHandle(Vector const& pos, real range)
{
    SingleProp * sip = getHandleProperty();
    
    if ( sip )
    {
        real dsm = 0;
        Single * res = nullptr;
        for ( Object * i : allHandles(sip) )
        {
            Single * s = static_cast<Single*>(i);
            real d = ( s->posFoot() - pos ).normSqr();
            if ( !res || d < dsm )
            {
                res = s;
                dsm = d;
            }
        }
        if ( res && dsm < range )
        {
            handle_ = res;
            return 1;
        }
    }
    return 0;
}


Single const* SimThread::handle() const
{
    SingleProp * sip = getHandleProperty();
    if ( sip && handle_ )
    {
        for ( Object * i : allHandles(sip) )
            if ( i == handle_ )
                return handle_;
    }
    handle_ = nullptr;
    return nullptr;
}


void SimThread::detachHandle()
{
    if ( handle_ )
    {
        if ( handle_->attached() )
            handle_->detach();
    }
}

void SimThread::moveHandle(Vector const& pos)
{
    if ( handle_ )
    {
        handle_->setPosition(pos);
    }
}


void SimThread::moveHandles(Vector const& vec)
{
    SingleProp * sip = getHandleProperty();
    if ( sip )
        ObjectSet::translateObjects(allHandles(sip), vec);
}


void SimThread::deleteHandles()
{
    lock();
    SingleProp * sip = getHandleProperty();
    if ( sip )
        sim_->erase(allHandles(sip));
    handle_ = nullptr;
    unlock();
}

void SimThread::erase_simul(bool arg) const
{
    Interface::erase_simul(arg);
    handle_ = nullptr;
}

//------------------------------------------------------------------------------
#pragma mark - Parameter modifications

/**
 Read config file from the start, allowing parameters to be changed, while 
 simulation objects remain as they are. This will pause a running simulation 
 is running live, read the config file, and allow it to proceed.
 */
void SimThread::reloadParameters(std::string const& file)
{
    lock();
    // set a parser that can only change properties:
    Parser(sim_, 1, 0, 0, 0, 0).readConfig(file);
    //std::cerr << "reloaded " << file << '\n';
    unlock();
}


/**
 This will execute the given code, with full rights to modify Simul.
 
 A simulation running live will be paused; the code executed in another Parser,
 and the simulation then allowed to proceed.
 
 This can be executed by the parent thread who does not own the data
 */
void SimThread::evaluate(std::string const& code)
{
    lock();
    try {
        Parser::evaluate(code);
    }
    catch( Exception & e ) {
        std::cerr << e.brief() << e.info() << '\n';
    }
    unlock();
}


/**
 Save current state in two files
 */
void SimThread::exportObjects(bool binary)
{
    lock();
    try {
        char str[64] = { '\0' };
        
        snprintf(str, sizeof(str), "properties%04li.cmp", reader_.currentFrame());
        sim_->writeProperties(str, true);
        
        snprintf(str, sizeof(str), "objects%04li.cmo", reader_.currentFrame());
        sim_->writeObjects(str, false, binary);
    }
    catch( Exception & e )
    {
        print_blue(std::cerr, e.brief());
        std::cerr << e.info() << " (export objects)\n";
    }
    unlock();
}


void SimThread::writeProperties(std::ostream& os, bool prune)
{
    lock();
    try {
        sim_->writeProperties(os, prune);
    }
    catch( Exception & e )
    {
        print_blue(std::cerr, e.brief());
        std::cerr << e.info() << " (write properties)\n";
    }
    unlock();
}


