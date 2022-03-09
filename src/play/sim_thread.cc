// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University

#include <sys/stat.h>
#include <sys/time.h>
#include <unistd.h>
#include <cstdio>
#include <time.h>
#include "sim_thread.h"
#include "exceptions.h"
#include "print_color.h"
#include "picket.h"
#include "glapp.h"

//------------------------------------------------------------------------------

/**
 This uses a Parser that cannot write to disc.
 The function `callback` is called when Parser::hold() is reached.
 */
SimThread::SimThread(Simul& sim, void (*callback)(void))
: Parser(sim, 1, 1, 1, 1, 0), hold_callback(callback)
{
    alone_  = 1;
    flag_   = 0;
    hold_   = 0;
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

    if ( flag_ )
        pthread_exit(nullptr);
    
    if ( ++hold_ >= period_ )
    {
        hold_ = 0;
        //debug("holding");
        hold_callback();
        if ( flag_ )
            pthread_exit(nullptr);
        wait();  // this also unlocks and locks the mutex
    }
}


//------------------------------------------------------------------------------
#pragma mark - Lauching threads


void SimThread::run()
{
    assert_true( isWorker() );
    try {
        Parser::readConfig();
    }
    catch( Exception & e ) {
        simul_.relax();
        std::cerr << e.brief() << e.info() << '\n';
        //flashText("Error: the simulation died");
    }
    hold_callback();
}


/** C-style function to cleanup after thread has terminated */
void child_cleanup(void * arg)
{
    SimThread * st = static_cast<SimThread*>(arg);
    st->alone_ = 1;
    st->unlock();
    //st->debug("ended");
}


/** C-style function to start a new thread */
void* run_launcher(void * arg)
{
    //std::clog << "slave  " << pthread_self() << '\n';
    SimThread * st = static_cast<SimThread*>(arg);
    st->lock();
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
    if ( alone_ )
    {
        flag_ = 0;
        //std::clog << "master " << pthread_self() << '\n';
        if ( pthread_create(&worker_, nullptr, run_launcher, this) )
            throw Exception("failed to create thread");
        alone_ = 0;
    }
}


//------------------------------------------------------------------------------


void SimThread::extend_run(size_t n_steps)
{
    assert_true( isWorker() );
    try {
        simul_.parser_ = this;
        Parser::execute_run(n_steps);
    }
    catch( Exception & e ) {
        std::cerr << e.brief() << e.info() << '\n';
        simul_.relax();
        //flashText("Error: %s", e.what());
    }
    simul_.parser_ = nullptr;
    hold_callback();
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
    if ( alone_ )
    {
        flag_ = 0;
        //std::clog << "master " << pthread_self() << '\n';
        if ( pthread_create(&worker_, nullptr, extend_launcher, this) )
            throw Exception("failed to create thread");
        alone_ = 0;
        return 0;
    }
    return 1;
}


//------------------------------------------------------------------------------
#pragma mark - Loading

int SimThread::loadFrame(size_t f)
{
    int r = 7;
    lock();
    try {
        r = reader_.loadFrame(simul_, f);
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
        r = reader_.loadNextFrame(simul_);
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
        r = reader_.loadLastFrame(simul_);
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
#pragma mark - Thread control & termination


void SimThread::step()
{
    assert_false( isWorker() );
    if ( !alone_ )
        signal();
}


/**
 ask the slave thread to exit at the next spontaneous halt
*/ 
void SimThread::stop()
{
    assert_false( isWorker() );
    if ( !alone_ )
    {
        // request clean termination:
        flag_ = 1;
        signal();
        // wait for termination:
        if ( !alone_ )
        {
            //debug("join...");
            // wait for termination:
            pthread_join(worker_, nullptr);
            alone_ = 1;
        }
    }
}

/**
 kill the slave thread immediately
 */
void SimThread::cancel()
{
    assert_false( isWorker() );
    if ( !alone_ )
    {
        flag_ = 2;
        //debug("cancel...");
        // force termination:
        if ( 0 == pthread_cancel(worker_) )
        {
            // wait for termination:
            pthread_join(worker_, nullptr);
            alone_ = 1;
            unlock();
        }
    }
}


void SimThread::restart()
{
    assert_false( isWorker() );
    stop();
    clear();
    start();
}

//------------------------------------------------------------------------------
#pragma mark - Mouse-controlled Single


SingleProp * SimThread::getHandleProperty() const
{
    Property * p = simul_.properties.find("single", "user_single");
    return static_cast<SingleProp*>(p);
}


SingleProp * SimThread::makeHandleProperty(real range)
{
    // Create a Hand that attaches fast and never detach:
    HandProp * hap = new HandProp("user_hand");
    hap->binding_range   = range;
    hap->binding_rate    = 1 / simul_.time_step();
    hap->unbinding_rate  = 0;
    hap->unbinding_force = INFINITY;
    hap->complete(simul_);
    simul_.properties.deposit(hap);

    SingleProp * sip = new SingleProp("user_single");
    sip->hand = "user_hand";
    sip->stiffness = 2000;
    sip->complete(simul_);
    simul_.properties.deposit(sip);
    
    return sip;
}


Single * SimThread::createHandle(Vector const& pos, real range)
{
    SingleProp * sip = getHandleProperty();
    if ( !sip )
        sip = makeHandleProperty(range);
    Single * res = new Picket(sip, pos);
    simul_.singles.add(res);
    handle_ = res;
    return res;
}


ObjectList SimThread::allHandles(SingleProp const* sip) const
{
    return simul_.singles.collect(match_property, sip);
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
        simul_.erase(allHandles(sip));
    handle_ = nullptr;
    unlock();
}

void SimThread::clear()
{
    assert_false( isWorker() );
    simul_.erase_all(1);
    handle_ = nullptr;
}

//------------------------------------------------------------------------------
#pragma mark - Parameter modifications

#if ( 0 )

#include <fcntl.h>

/// set file to not block on read() even if data is not available:
void set_nonblocking(int fd)
{
    fcntl(fd, F_SETFL, O_NONBLOCK);
}

#endif


/// check if file has data for input
inline int has_input(int fd)
{
    fd_set fds;
    FD_ZERO(&fds);
    FD_SET(fd, &fds);
    struct timeval tv = {0, 10};   // seconds, microseconds
    return select(1, &fds, nullptr, nullptr, &tv);
}


/**
 Read standard input and executes incoming commands.
 This should be executed by a process who already owns the lock on the data
 */
size_t SimThread::executePipedCommands(size_t max_nb_lines)
{
    const size_t LINESIZE = 2048;
    clearerr(stdin);
    
    if ( has_input(STDIN_FILENO) > 0 )
    {
        size_t cnt = 0;
        // some input is available, process line-by-line:
        char str[LINESIZE];
        
        // read one line from standard input (including terminating \n):
        while ( fgets(str, LINESIZE, stdin) )
        {
            //write(STDOUT_FILENO, ">>>> ", 5); write(STDOUT_FILENO, str, strlen(str));
            try {
                Parser::evaluate(str);
            }
            catch ( Exception & e ) {
                print_green(std::cerr, e.brief());
                std::cerr << " in: " << str;
            }
            if ( ++cnt >= max_nb_lines )
                break;
            // check if more input is available:
            if ( has_input(STDIN_FILENO) < 1 )
                break;
        }
        glApp::flashText(str);
        //printf("executed %lu lines from standard input\n", cnt);
        return cnt;
    }
    return 0;
}

/**
 Read config file from the start, allowing parameters to be changed, while 
 simulation objects remain as they are. This will pause a running simulation 
 is running live, read the config file, and allow it to proceed.
 */
void SimThread::reloadParameters(std::string const& file)
{
    lock();
    // set a parser that can only change properties:
    Parser(simul_, 1, 0, 0, 0, 0).readConfig(file);
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
        simul_.writeProperties(str, true);
        
        snprintf(str, sizeof(str), "objects%04li.cmo", reader_.currentFrame());
        simul_.writeObjects(str, false, binary);
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
        simul_.writeProperties(os, prune);
    }
    catch( Exception & e )
    {
        print_blue(std::cerr, e.brief());
        std::cerr << e.info() << " (write properties)\n";
    }
    unlock();
}


