// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University

#include <cstdio>
#include <time.h>
#include <unistd.h>
#include <signal.h>
#include <fstream>

#include "vector.h"
#include "single.h"
#include "simul_prop.h"
#include "simul.h"
#include "sim_thread.h"
#include "exceptions.h"
#include "print_color.h"

static void thread_signal_handler(int sig)
{
    psignal(sig, "SimThread");
    _exit(sig);
}

//------------------------------------------------------------------------------

/**
 This uses a Parser that cannot write to disc.
 */
SimThread::SimThread(Simul* sim)
: Parser(sim, 1, 1, 1, 1, 0, 0)
{
    status_ = -2;
    repeat_ = 0;
    hold_   = 0;
    holding_= 0;
    cycle_ = 1;
    pthread_mutex_init(&mutex_, nullptr);
    pthread_cond_init(&condition_, nullptr);
    config_code = nullptr;
    code_size = 0;
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
    if ( status_ == -1 )
        pthread_join(child_, nullptr);
    pthread_cond_destroy(&condition_);
    pthread_mutex_destroy(&mutex_);
    // release memory only if it was allocated by us:
    if ( code_size > 0 ) free(config_code);
}

//------------------------------------------------------------------------------
#pragma mark - Process control


void SimThread::debug(const char* msg) const
{
    if ( isWorker() )
        fprintf(stdout, "\n- - - %-16s", msg);
    else
        fprintf(stdout, "\n* * * %-16s", msg);
}


void SimThread::gubed(const char* msg) const
{
    if ( isWorker() )
        fprintf(stdout, "  - - %-16s ", msg);
    else
        fprintf(stdout, "  * * %-16s ", msg);
}


void SimThread::hold()
{
    assert_true( isWorker() );

    if ( repeat_ == -1 )
        pthread_exit(nullptr);
    
    if ( ++hold_ >= cycle_ )
    {
        hold_ = 0;
        //debug("holding");
        holding_ = 1;
        cond_wait();  // unlocks mutex, wait for 'signal', lock again
        holding_ = 0;
    }
}


//------------------------------------------------------------------------------
#pragma mark - Launching threads


/** called before thread is started */
void SimThread::getReady()
{
    assert_false( isWorker() );
    repeat_ = 0;
    Parser::restart_ = 0;
    //fprintf(stderr, "ready %i\n", sim_->prop.flag);
    assert_true(holding_ == 0);
    if ( status_ == -1 )
        pthread_join(child_, nullptr);
}


void SimThread::run()
{
    assert_true( isWorker() );
    while ( 1 )
    {
        sim_->initCytosim();
        // read configuration file:
        std::string filename = sim_->prop.config_file;
        if ( config_code )
        {
            // code was provided as string:
            std::stringstream iss(config_code);
            readConfig(iss, filename);
        }
        else
        {
            // read code from a file:
            std::ifstream ifs(filename.c_str(), std::ifstream::in);
            readConfig(ifs, filename);
        }
        // the simulation is now terminated, we shall wait for `repeat_`
        //fprintf(stderr, "Completed simulation %i\n", sim_->prop.flag);
        while ( !repeat_ )
        {
            holding_ = 2;
            cond_wait();
        }
        if ( --repeat_ < 0 )
            break;
        eraseSimul(1);
    }
}


/** C-style function called automatically after thread has terminated */
void child_cleanup(void * arg)
{
    SimThread * st = static_cast<SimThread*>(arg);
    assert_true( st->isWorker() );
    //fprintf(stderr, "cleanup dead %i\n", st->sim_->prop.flag);
    st->status_ = -1;
    st->holding_ = 0;
    st->unlock();
    //st->debug("ended");
}


/** C-style function to start a new thread */
void* child_entry(void * arg)
{
    //std::clog << "slave  " << pthread_self() << '\n';sign
    SimThread * st = static_cast<SimThread*>(arg);
    st->lock();
    // let the system cleanup upon normal termination:
    pthread_cleanup_push(child_cleanup, arg);
    try {
        st->run();
    }
    catch( Exception & e ) {
        std::cerr << e.brief() << e.info() << '\n';
        //flashText("Error: the simulation died");
    }
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
        getReady();
        //std::clog << "master " << pthread_self() << '\n';
        status_ = pthread_create(&child_, nullptr, child_entry, this);
        // let the system cleanup upon normal child termination:
        if ( status_ )
            printf("%p failed to create thread", this);
    }
}

//------------------------------------------------------------------------------

void SimThread::prolong_run(double sec)
{
    assert_true( isWorker() );
    Parser * back = sim_->parser();
    sim_->initCytosim();
    sim_->parser(this);
    try {
        Parser::execute_run(sec);
    }
    catch( Exception & e ) {
        std::cerr << e.brief() << e.info() << '\n';
        //flashText("Error: %s", e.what());
    }
    sim_->parser(back);
}


/** C-style function to start a new thread */
void* prolong_launcher(void * arg)
{
    //std::clog << "slave  " << pthread_self() << '\n';
    SimThread * st = static_cast<SimThread*>(arg);
    st->lock();
    pthread_cleanup_push(child_cleanup, arg);
    st->prolong_run(3600);
    pthread_cleanup_pop(1);
    return nullptr;
}


/// call `prolong_run()` in a slave thread
int SimThread::prolong()
{
    assert_false( isWorker() );
    if ( status_ )
    {
        getReady();
        //std::clog << "master " << pthread_self() << '\n';
        status_ = pthread_create(&child_, nullptr, prolong_launcher, this);
    }
    return status_;
}

//------------------------------------------------------------------------------
#pragma mark - Thread control & termination

/**
 signal the slave thread to exit at the next convenient pause
*/
void SimThread::stop()
{
    if ( status_ == 0 )
    {
        assert_false( isWorker() );
        // request clean termination:
        repeat_ = -1;
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
    if ( status_ == 0 )
    {
        assert_false( isWorker() );
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
        print_blue(stderr, e.brief());
        std::cerr << e.info() << " (loading frame " << f << ")\n";
    }
    unlock();
    return r;
}


int SimThread::loadPreviousFrame()
{
    if ( reader_.currentFrame() > 0 )
        return loadFrame(reader_.currentFrame()-1);
    return 7;
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
        print_blue(stderr, e.brief());
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
        print_blue(stderr, e.brief());
        std::cerr << e.info() << " (loading last frame)\n";
    }
    unlock();
    return r;
}

//------------------------------------------------------------------------------
#pragma mark - User-controlled Single


SingleProp * SimThread::getHandleProperty() const
{
    Property * p = sim_->properties.find("single", "live_single");
    return static_cast<SingleProp*>(p);
}


SingleProp * SimThread::makeHandleProperty(real range)
{
    // Create a Hand that attaches fast and never detach:
    HandProp * hap = new HandProp("live_hand");
    hap->binding_range   = range;
    hap->binding_rate    = 1.0 / sim_->time_step();
    hap->unbinding_rate  = 0;
    hap->unbinding_force = INFINITY;
    hap->complete(*sim_);
    sim_->properties.deposit(hap);

    SingleProp * sip = new SingleProp("live_single");
    sip->hand = "live_hand";
    sip->activity = "fixed";
    sip->confine = CONFINE_OFF;
    sip->stiffness = 2000;
    sip->complete(*sim_);
    sim_->properties.deposit(sip);
    
    return sip;
}


Single * SimThread::createHandle(Vector3 const& pos3, real range)
{
    SingleProp * sip = getHandleProperty();
    if ( !sip )
        sip = makeHandleProperty(range);
    handle_ = sim_->singles.addFreeSingle(sip, Vector(pos3));
    return handle_;
}


ObjectList SimThread::allHandles(SingleProp const* sip) const
{
    return sim_->singles.collect(match_property, sip);
}


bool SimThread::selectClosestHandle(Vector3 const& pos3, real range)
{
    SingleProp * sip = getHandleProperty();
    
    if ( sip )
    {
        real dsm = 0;
        Vector pos(pos3);
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

void SimThread::moveHandle(Vector3 const& vec)
{
    if ( handle_ )
    {
        handle_->setPosition(Vector(vec));
    }
}


void SimThread::moveHandles(Vector3 const& vec)
{
    SingleProp * sip = getHandleProperty();
    if ( sip )
        ObjectSet::translateObjects(allHandles(sip), Vector(vec));
}


void SimThread::deleteHandles()
{
    lock();
    SingleProp * sip = getHandleProperty();
    if ( sip )
        sim_->singles.eraseObjects(allHandles(sip));
    handle_ = nullptr;
    unlock();
}

void SimThread::eraseSimul(bool arg) const
{
    Interface::eraseSimul(arg);
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
    Parser(sim_, 1, 0, 0, 0, 0, 0).readConfig(file);
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
void SimThread::exportObjects(const int binary)
{
    lock();
    try {
        char str[64] = { '\0' };
        snprintf(str, sizeof(str), "properties%04li.cmp", reader_.currentFrame());
        std::ofstream ofs(str);
        sim_->writeProperties(ofs, true);
        
        snprintf(str, sizeof(str), "objects%04li.cmo", reader_.currentFrame());
        sim_->writeObjects(str, false, binary);
    }
    catch( Exception & e )
    {
        print_blue(stderr, e.brief());
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
        print_blue(stderr, e.brief());
        std::cerr << e.info() << " (write properties)\n";
    }
    unlock();
}


