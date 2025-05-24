// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University
#ifndef SIM_THREAD_H
#define SIM_THREAD_H

#include <pthread.h>
#include <csignal>

#include "parser.h"
#include "frame_reader.h"

class SingleProp;
class Single;

/// SimThread is used to run a simulation in a dedicated thread
/**
 The SimThread needs to derive from Parser, for overwritting 'hold()'
 */
class SimThread : public Parser
{
    /// cleanup callback
    friend void child_cleanup(void*);
    
private:
    
    /// thread used to run the simulation
    pthread_t child_;
    
    /// mutex protecting write access to simulation state
    pthread_mutex_t mutex_;

    /// condition variable used to control the thread execution
    pthread_cond_t condition_;
    
    /// signal for hold()
    int holding_;

    /// error state of the child thread: 0 = running
    int status_;
    
    /// a flag to indicate if child thread should terminate or restart
    int repeat_;

    /// counter for hold()
    unsigned int hold_;
    
    /// period for hold()
    unsigned int cycle_;
    
    /// Reader used to access frames in a trajectory file
    FrameReader reader_;

    /// the current Single being controlled with the mouse
    mutable Single * handle_;
    
    /// return the SingleProp used for the handles
    SingleProp * getHandleProperty() const;

    /// make a new SingleProp for the handles with given attachment range
    SingleProp * makeHandleProperty(real range);
    
    /// return list of Handles
    ObjectList allHandles(SingleProp const*) const;

    /// True if current thread is the worker thread
    bool isWorker() const { return pthread_equal(pthread_self(), child_); }
    
    /// reset state variable before starting
    void getReady();
    
public:
    
    /// the size allocated for `config_code` below, or zero if externally managed
    size_t code_size;
    
    /// the text of the config file, allocated externally
    char * config_code;

    /// run the simulation live
    void run();
    
    /// continue to run for specified time beyond its normal termination
    void prolong_run(double);

    /// redefining Interface::hold(), which is called repeatedly at each timestep
    void hold();
    
    /// return true if new data is available
    int holding() const { return holding_; }

    /// print message to identify thread
    void debug(const char *) const;
    
    /// print message to identify thread
    void gubed(const char *) const;

    /// create a SimThread with given Simul
    SimThread(Simul*);
    
    /// create a SimThread to be initialized later
    SimThread() : SimThread(nullptr) {}
    
    /// destructor
    ~SimThread();

#if ( 1 )

    /// lock access to the Simulation data
    void lock()   { pthread_mutex_lock(&mutex_); }
    
    /// unlock access to the Simulation data
    void unlock() { pthread_mutex_unlock(&mutex_);}
    
    /// try to lock access to the Simulation data
    int trylock() { return pthread_mutex_trylock(&mutex_); }

    /// unlock access to data and wait for the condition
    int cond_wait() { return pthread_cond_wait(&condition_, &mutex_); }
    
    /// send signal to child thread. Will reset holding state if successful
    void signal() { pthread_cond_signal(&condition_); }

#else
    
    /// lock access to the Simulation data
    void lock()   { debug("  lock..."); pthread_mutex_lock(&mutex_); debug("  locked!"); }
    
    /// unlock access to the Simulation data
    void unlock() { pthread_mutex_unlock(&mutex_); gubed("unlock"); }
    
    /// try to lock access to the Simulation data
    int trylock() { int R=pthread_mutex_trylock(&mutex_); debug(R?"  failed trylock":"  trylock"); return R; }
    
    /// wait for the condition
    int cond_wait() { debug("unlock, wait"); int R=pthread_cond_wait(&condition_, &mutex_); debug("wake, lock"); return R; }
    
    /// signal child thread to continue, called by the parent
    void signal() { debug("signal"); pthread_cond_signal(&condition_); }
    
#endif
    
    /// set how many 'hold()' are necessary to halt the thread
    void period(unsigned int c) { cycle_ = c; }
    
    /// true if child thread is not running
    bool alone() const { return status_ != 0; }

    /// true if child thread is (apparently) running
    bool alive() const { return status_ == 0; }
    
    /// return status of child process, which is 0 if alive and healthy
    int dead() { return status_; }

    /// start the thread that will run a simulation
    void start();
    
    /// continue to run the simulation after its normal termination
    int prolong();
    
    /// gently stop the simulation
    void stop();

    /// ask the simulation to stop
    void cancel() { pthread_cancel(child_); }

    /// kill the simulation
    void terminate() { pthread_kill(child_, SIGTERM); }

    /// wait for child to terminate
    void join() { if ( 0 == pthread_join(child_, nullptr) ) status_ = -2; }
    
    /// stop the simulation and wait for cleaning operations
    void cancel_join();
    
    /// gently restart the simulation engine
    void restart() { repeat_ = 1; signal(); };

    /// clear the simulation world
    void eraseSimul(bool) const;
    
    /// halt the live simulation, read the config file and change the object parameters
    void reloadParameters(std::string const& file);
    
    /// execute given code
    void evaluate(std::string const&);
    
    /// export simulation Propertes and Objects to file
    void exportObjects(int binary);
    
    /// export properties to file
    void writeProperties(std::ostream&, bool prune);

    
    /// open trajectory file for input
    void openFile(std::string const& name) { reader_.openFile(name); }
    
    /// returns true if ready to read from file
    bool goodFile() { return reader_.good(); }
    
    /// end-of-file flag of input file
    int eof() const { return reader_.eof(); }
    
    /// rewind file
    void rewindFile() { lock(); reader_.rewind(); unlock(); }
    
    /// attempt to load specified frame from file (0 = first frame; -1 = last frame)
    int loadFrame(size_t f);
    
    /// load previous frame in file
    int loadPreviousFrame();

    /// load next frame in file
    int loadNextFrame();
    
    /// attempt to load last frame from file
    int loadLastFrame();

    /// index of current frame (0 is lowest valid value)
    size_t currentFrame() const { return reader_.currentFrame(); }

    
    /// return the Single that is manipulated by the User
    Single const* handle() const;

    /// make a new Single that can be controlled by the user
    Single * createHandle(Vector3 const&, real range);
    
    /// switch current handle
    bool selectClosestHandle(Vector3 const&, real range);
    
    /// detach current handle
    void detachHandle();
    
    /// move the current handle
    void moveHandle(Vector3 const&);
    
    /// move all handles
    void moveHandles(Vector3 const&);
    
    /// delete all handles
    void deleteHandles();
    
    /// detach current handle from mouse control
    void releaseHandle() { handle_ = nullptr; }
    
};


#endif

