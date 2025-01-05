// Cytosim was created by Francois Nedelec. Copyright 2023 Cambridge University
// Created by Francois Nedelec on 04/04/2012. Updated 01/02/2023, 18/2/2023


#include <list>
#include <cmath>
#include <iostream>

#include "random.h"


real sim_time = 0;

/// Node of a doubly-linked list
class Linkable
{
public:
    
    Linkable * prevL;
    Linkable * nextL;
    
public:
    
    Linkable() : prevL(0), nextL(0) {}
    
    void push_front(Linkable *& head, Linkable *& tail)
    {
        prevL = nullptr;
        nextL = head;
        if ( head )
            head->prevL = this;
        else
            tail = this;
        head = this;
    }

    void push_back(Linkable *& head, Linkable *& tail)
    {
        prevL = tail;
        nextL = nullptr;
        if ( tail )
            tail->nextL = this;
        else
            head = this;
        tail = this;
    }
    
    void pop(Linkable *& head, Linkable *& tail)
    {
        if ( prevL )
            prevL->nextL = nextL;
        else
            head = nextL;
        
        if ( nextL )
            nextL->prevL = prevL;
        else
            tail = prevL;
    }
};


/// Stochastic event
template < typename ENGINE >
class Reaction : public Linkable
{
public:
    
    /// pointer to a member function
    typedef int (ENGINE::*MFPR)(Reaction<ENGINE>*);

private:
    
    real mRand;
    real mTime;
    MFPR mFunc;
    
public:
    
    Reaction(real rate, MFPR mfp)
    {
        mRand = RNG.exponential();
        mTime = mRand / rate;
        mFunc = mfp;
    }
    
    /// use this if the rate is constant
    void step(real interval)
    {
        mTime -= interval;
    }
    
    /// use this if the rate changes with time
    void step(real interval, real rate)
    {
        mRand -= rate * interval;
        if ( mRand < 0 )
            mTime = mRand / rate;
    }
    
    /// call engine's member function:
    int act(ENGINE & obj)
    {
        return (obj.*mFunc)(this);
    }
    
    /// increment Gillespie time
    void refresh(real rate)
    {
        if ( mRand < 0 )
        {
            mRand += RNG.exponential();
            mTime = mRand / rate;
        }
        else
        {
            mTime += RNG.exponential() / rate;
        }
    }
    
    real time()
    {
        return mTime;
    }
    
    void print(std::ostream& os)
    {
        os << mTime << "  " << mFunc << '\n';
    }    
};


/// Stochastic engine
class Engine
{
public:
    
    typedef Reaction<Engine> Event;
    typedef Linkable * iterator;
    
protected:
    
    /// reaction rates
    real rate1, rate2, rate3;

    /// list heads to hold the events:
    Linkable * head_, * tail_;

public:
    
    Engine(real r1, real r2, real r3) : head_(nullptr), tail_(nullptr)
    {
        rate1 = r1;
        rate2 = r2;
        rate3 = r3;
        initialize();
    }
    
    void initialize()
    {
    }

    void add(real r, Event::MFPR mfp)
    {
        (new Event(r, mfp))->push_back(head_, tail_);
    }

    void step(real interval)
    {
        Linkable * H = nullptr;
        Linkable * T = nullptr;
        
        /* Update events; move the ones that will fire to [H,T] list */
        iterator stop = head_;
        iterator i = head_;
        while ( i )
        {
            Event * e = static_cast<Event*>(i);
            i = i->nextL;
            e->step(interval);
            if ( e->time() < 0 )
            {
                e->pop(head_, tail_);
                e->push_back(H, T);
            }
        }
        
        /* Fire events from [H,T] list in the order of time */
        while ( H )
        {
            Event * evt = static_cast<Event*>(H);
            real next_time = evt->time();
            // find earliest event in [H,T] list:
            for ( iterator i = H->nextL; i; i = i->nextL )
            {
                Event * e = static_cast<Event*>(i);
                real t = e->time();
                if ( t < next_time )
                {
                    next_time = t;
                    evt = e;
                }
            }
            // fire event, and proceed depending on outcome:
            if ( 0 == evt->act(*this) )
            {
                // case 1: the event only fires once and is deleted:
                evt->pop(H, T);
                delete(evt);
            }
            else if ( evt->time() >= 0 )
            {
                // case 2: event can fire more, but will only do in the future:
                evt->pop(H, T);
                evt->push_back(head_, tail_);
            }
            // case 3: event fires again in the same interval ( time() < 0 )
        }
    }

    //----------------------Event's action callbacks---------------------------
    
    int event1(Event* e)
    {
        std::cout << "event 1 @ " << sim_time+e->time() << '\n';
        add(rate1, &Engine::event2);
        add(rate1, &Engine::event3);
        return 0;
    }
    
    int event2(Event* e)
    {
        std::cout << "event 2 @ " << sim_time+e->time() << '\n';
        e->refresh(rate2);
        return 1;
    }
    
    int event3(Event* e)
    {
        std::cout << "event 3 @ " << sim_time+e->time() << '\n';
        e->refresh(rate3);
        return 1;
    }
    
};


int main(int argc, char * argv[])
{
    RNG.seed();
    Engine engine(2,1,1);
    engine.add(1, &Engine::event1);
    engine.add(1, &Engine::event1);

    real time_step = 1;
    for ( int i = 0; i < 10; ++i )
    {
        sim_time += time_step;
        engine.step(time_step);
        std::cout << "------  " << sim_time << std::endl;
    }
}

