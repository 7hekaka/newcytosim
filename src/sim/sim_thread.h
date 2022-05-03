// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University
#ifndef SIM_THREAD_H
#define SIM_THREAD_H

#include <pthread.h>
#include "simul.h"
#include "parser.h"
#include "frame_reader.h"


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
    bool holding_;

    /// a flag reflecting if the child thread is running or not
    bool alone_;
    
    /// a flag to indicate that child thread should terminate or restart
    int flag_;

    /// counter for hold()
    unsigned int hold_;
    
    /// period for hold()
    unsigned int period_;

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
    
public:
    
    /// run the simulation live
    void run();
    
    /// continue to run a simulation beyond its normal termination
    void extend_run(size_t n_steps);

    /// redefining Interface::hold(), which is called repeatedly at each timestep
    void hold();
    
    /// return true if new data is available
    int holding() { return holding_; }

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
    
    /// send signal to child thread
    void signal() { holding_ = 0; if ( !alone_ ) pthread_cond_signal(&condition_); }

#else
    
    /// lock access to the Simulation data
    void lock()   {  debug("  lock..."); pthread_mutex_lock(&mutex_); debug("  locked!"); }
    
    /// unlock access to the Simulation data
    void unlock() { pthread_mutex_unlock(&mutex_); gubed("  unlock"); }
    
    /// try to lock access to the Simulation data
    int trylock() { int R=pthread_mutex_trylock(&mutex_); debug(R?"  failed trylock":"  trylock"); return R; }
    
    /// wait for the condition
    int cond_wait() { debug("unlock, wait"); int R=pthread_cond_wait(&condition_, &mutex_); debug("wake, lock"); return R; }
    
    /// signal child thread to continue
    void signal() { debug("signal"); holding_ = 0; pthread_cond_signal(&condition_); }
    
#endif
    
    /// set how many 'hold()' are necessary to halt the thread
    void period(unsigned int c) { period_ = c; }
    
    /// true if child thread is running
    bool alone() const { return alone_; }

    /// true if child thread is running
    bool alive() const { return !alone_; }
    
    /// start the thread that will run a simulation
    void start();
    
    /// continue to run the simulation after its normal termination
    int extend();
    
    /// gently stop the simulation
    void stop();

    /// stop the simulation
    void cancel();
    
    /// restart engine
    void restart();

    /// clear the simulation world
    void erase_simul(bool) const;
    
    /// halt the live simulation, read the config file and change the object parameters
    void reloadParameters(std::string const& file);
    
    /// execute given code
    void evaluate(std::string const&);
    
    /// export simulation Propertes and Objects to file
    void exportObjects(bool binary);
    
    /// export properties to file
    void writeProperties(std::ostream&, bool prune);

    
    /// open trajectory file for input
    void openFile(std::string const& name) { reader_.openFile(name); }
    
    /// true if ready to read from file
    bool goodFile() const { return reader_.good(); }
    
    /// status of file
    int eof() const { return reader_.eof(); }
    
    /// rewind file
    void rewindFile() { lock(); reader_.rewind(); unlock(); }
    
    /// attempt to load specified frame from file (0 = first frame; -1 = last frame)
    int loadFrame(size_t f);

    /// load next frame in file
    int loadNextFrame();
    
    /// attempt to load last frame from file
    int loadLastFrame();

    /// index of current frame (0 is lowest valid value)
    size_t currentFrame() const { return reader_.currentFrame(); }

    
    /// return the Single that is manipulated by the User
    Single const* handle() const;

    /// make a new Single that can be controlled by the user
    Single * createHandle(Vector const&, real range);
    
    /// switch current handle
    bool selectClosestHandle(Vector const&, real range);
    
    /// detach current handle
    void detachHandle();
    
    /// move the current handle
    void moveHandle(Vector const&);
    
    /// move all handles
    void moveHandles(Vector const&);
    
    /// delete all handles
    void deleteHandles();
    
    /// detach current handle from mouse control
    void releaseHandle() { handle_ = nullptr; }
    
};


#endif

