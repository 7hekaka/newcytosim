// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University

#include "player.h"

#include "gle.h"
#include "gle_color.h"
#include "gle_color_list.h"
#include "glapp.h"

using glApp::flashText;

#include "fiber_disp.h"
#include "line_disp.h"
#include "point_disp.h"
#include "player_disp.cc"


Player::Player()
: disp("*"), prop("*"), worker(&simul), mDisplay(nullptr)
{
}

Player::~Player()
{
    clear();
}

void Player::clear()
{
    worker.stop();
    worker.erase_simul(1);
    dispList.erase();
    if ( mDisplay )
        delete(mDisplay);
    mDisplay = nullptr;
}

//------------------------------------------------------------------------------
#pragma mark - I/O


void Player::previousFrame()
{
    if ( worker.currentFrame() > 0 )
        worker.loadFrame(worker.currentFrame()-1);
    else {
        if ( prop.loop )
            worker.loadLastFrame();
        else
            stop();
    }
}

/**
 Reads the next frame from the current file position.
 */
void Player::nextFrame()
{    
    int res = worker.loadNextFrame();
    if ( res == 1 )
    {
        // this means end of file is reached:
        if ( prop.auto_exit )
            exit(EXIT_SUCCESS);
        if ( prop.loop )
            worker.loadFrame(0);
        else
        {
            flashText("end-of-file\n");
            stop();
        }
    }
    else if ( res )
    {
        flashText("Error "+std::to_string(res)+" occured while loading frame\n");
        stop();
    }
    //std::clog << simul.nbObjects() << " objects @ " << std::fixed << simul.time() << "s\n";
}

//------------------------------------------------------------------------------
#pragma mark - Commands


void Player::rewind()
{
    if ( worker.goodFile() )
    {
        stop();
        worker.rewindFile();
        worker.loadFrame(0);
        glApp::postRedisplay();
    }
}


bool Player::startPlayback()
{
    if ( worker.goodFile()  &&  prop.play != 1  && !prop.goLive )
    {
        //rewind file if its end was reached:
        if ( worker.eof() )
            worker.rewindFile();
        prop.play = 1;
        return true;
    }
    return false;
}


bool Player::startBackward()
{
    if ( prop.play != -1 )
    {
        if ( worker.currentFrame() == 0 )
            worker.loadLastFrame();
        else
            flashText("Play reverse");
        prop.play = -1;
        return true;
    }
    return false;
}


void Player::accelerate()
{
    prop.delay /= 2;
    //the delay should be compatible with graphic refresh rates:
    const unsigned int min_delay = 1;
    if ( prop.delay < min_delay )
    {
        prop.delay = min_delay;
        if ( prop.goLive )
            flashText("Delay is %i ms! use 'A' to jump frames", prop.delay);
        else
            flashText("Delay is %i ms!", prop.delay);
    }
    else {
        flashText("Delay %i ms", prop.delay);
    }
}


void Player::stop()
{
    prop.goLive = 0;
    prop.play = 0;
    prop.save_images = 0;
}


void Player::startstop()
{
    if ( worker.alive() )
        prop.goLive = !prop.goLive;
    else if ( worker.goodFile() )
    {
        if ( !prop.play )
            startPlayback();
        else
            stop();
    }
}


void Player::extendLive()
{
    if ( 0 == worker.extend() )
        flashText("Extend simulation...");
    prop.goLive = 1;
}


void Player::restart()
{
    try
    {
        worker.stop();
        worker.erase_simul(1);
        dispList.erase();
        worker.start();
    }
    catch( Exception & e ) {
        flashText("Error: %s", e.what());
    }
}


//------------------------------------------------------------------------------
#pragma mark - Display selection routines

inline FiberDisp* toFiberDisp(Property * ptr)
{
    return static_cast<FiberDisp*>(ptr);
}

inline PointDisp* toPointDisp(Property * ptr)
{
    return static_cast<PointDisp*>(ptr);
}


PropertyList Player::allFiberDisp() const
{
    return dispList.find_all("fiber:display");
}

PropertyList Player::allVisibleFiberDisp() const
{
    PropertyList res, plist = dispList.find_all("fiber:display");
    
    for ( Property * i : plist )
    {
        if ( toFiberDisp(i)->visible )
            res.push_back(i);
    }
    return res;
}

PropertyList Player::allHandDisp() const
{
    return dispList.find_all("hand:display");
}

PropertyList Player::allVisibleHandDisp() const
{
    PropertyList res, plist = dispList.find_all("hand:display");
    
    for ( Property * i : plist )
    {
        if ( toPointDisp(i)->visible )
            res.push_back(i);
    }
    return res;
}

PropertyList Player::allSphereDisp() const
{
    return dispList.find_all("bead:display", "solid:display", "sphere:display");
}

PropertyList Player::allVisibleSphereDisp() const
{
    PropertyList res;
    PropertyList all = dispList.find_all("bead:display", "solid:display", "sphere:display");
    for ( Property * i : all )
    {
        if ( toPointDisp(i)->visible )
            res.push_back(i);
    }
    return res;
}

PropertyList Player::allSpaceDisp() const
{
    return dispList.find_all("space:display");
}

FiberDisp * Player::firstFiberDisp() const
{
    PropertyList plist = allVisibleFiberDisp();
    if ( plist.size() )
        return toFiberDisp(plist.front());
    return nullptr;
}

/**
 Write global parameters that control the display:
 - GlappProp
 - DisplayProp
 - PlayerProp
 .
 */
void Player::writePlayParameters(std::ostream& os, bool prune) const
{
    os << "set " << simul.prop.name() << " display\n{\n";
    if ( glApp::views.size() > 0 )
    {
        View& view = glApp::currentView();
        view.write_values_diff(os, prune);
    }
    disp.write_values_diff(os, prune);
    //output parameters for the main view:
    prop.write_values_diff(os, prune);
    os << "}\n";
}

/**
 Write all the parameters that control the display:
 - GlappProp
 - DisplayProp
 - PlayerProp
 - ObjectDisp
 .
 */
void Player::writeDisplayParameters(std::ostream& os, bool prune) const
{
    dispList.write(os, prune);
}
