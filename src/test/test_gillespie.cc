// Cytosim was created by Francois Nedelec. Copyright 2023 Cambridge University
// Created by Francois Nedelec on 04/04/2012. Updated 01/02/2023


#include <list>
#include <cmath>
#include <iostream>

#include "random.h"


real sim_time = 0;

/// Node of a doubly-linked list
class Linkable
{
public:
    
    Linkable * prev_;
    Linkable * next_;
    
public:
    
    Linkable() : prev_(0), next_(0) {}
    
    void push_front(Linkable *& head, Linkable *& tail)
    {
        prev_ = nullptr;
        next_ = head;
        if ( head )
            head->prev_ = this;
        else
            tail = this;
        head = this;
    }

    void push_back(Linkable *& head, Linkable *& tail)
    {
        prev_ = tail;
        next_ = nullptr;
        if ( tail )
            tail->next_ = this;
        else
            head = this;
        tail = this;
    }
    
    void pop(Linkable *& head, Linkable *& tail)
    {
        if ( prev_ )
            prev_->next_ = next_;
        else
            head = next_;
        
        if ( next_ )
            next_->prev_ = prev_;
        else
            tail = prev_;
    }
};


/// Stochastic event
template < typename ENGINE >
class Reaction : public Linkable
{
public:
    
    /// pointer to a member function
    typedef void (ENGINE::*MFPR)(Reaction<ENGINE>*);

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
    void act(ENGINE & obj)
    {
        (obj.*mFunc)(this);
        std::cout << " @ " << sim_time + mTime << '\n';
    }
    
    /// increment Gillespie time
    void renew(real rate)
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
    
    void renew(Event * e, real rate)
    {
        e->push_back(head_, tail_);
        e->renew(rate);
    }

public:
    
    Engine(real r1, real r2, real r3) : head_(nullptr), tail_(nullptr)
    {
        rate1 = r1;
        rate2 = r2;
        rate3 = r3;
        initialize();
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
            i = i->next_;
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
            Event * earliest = static_cast<Event*>(H);
            real next_time = earliest->time();
            // find earliest event in [H,T] list:
            for ( iterator i = H->next_; i; i = i->next_ )
            {
                Event * e = static_cast<Event*>(i);
                real t = e->time();
                if ( t < next_time )
                {
                    next_time = t;
                    earliest = e;
                }
            }
            earliest->pop(H, T);
            earliest->act(*this);
        }
    }
    
    void initialize()
    {
    }

    //----------------------Event's action callbacks---------------------------
    
    void event1(Event* e)
    {
        std::cout << "event 1";
        add(rate1, &Engine::event2);
        add(rate1, &Engine::event3);
        delete(e);
    }
    
    void event2(Event* e)
    {
        std::cout << "event 2";
        renew(e, rate2);
    }
    
    void event3(Event* e)
    {
        std::cout << "event 3";
        renew(e, rate3);
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
        std::cout << "------" << std::endl;
    }
}

