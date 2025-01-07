// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University

/**
 Callback for mouse clicks
 */
void processMouseClick(int, int, const Vector3& pos, int)
{
    // distance in pixels where mouse-Hand binds:
    const real pixrad = 10;
    const real range = pixrad * glApp::currentView().pixelSize();

    worker.lock();
    if ( worker.selectClosestHandle(pos, range) )
        worker.moveHandle(pos);
    else
    {
        if ( worker.handle() )
        {
            worker.detachHandle();
            worker.moveHandle(pos);
        }
        else
        {
            Single * s = worker.createHandle(pos, range);
            PointDisp *& pd = s->prop->hand_prop->disp;
            pd = static_cast<PointDisp*>(player.mPointDisp.find("hand:display", "live_hand"));
            if ( !pd )
            {
                pd = new PointDisp("hand:display", "live_hand");
                pd->size   = pixrad;
                pd->color  = glApp::currentView().front_color;
                pd->color2 = pd->color;
                player.mPointDisp.deposit(pd);
            }
        }
    }
    worker.unlock();
}


/**
 Processes mouse motion
 */
void processMouseDrag(int, int, Vector3& ori, const Vector3& pos, int mode)
{
    worker.lock();
    if ( mode )
    {
        worker.moveHandles(pos-ori);
        ori = pos;
    }
    else
        worker.moveHandle(pos);
    worker.unlock();
}


//------------------------------------------------------------------------------
#pragma mark - Timer callback


void timerCallback(const int value)
{
    unsigned millisec = prop.delay;

    //std::clog << "timerCallback(" << value << ") @ " << std::fixed << simul.time() << "s\n";
    if ( player.goLive && worker.alive() )
    {
        //worker.debug("timerCallback");
        if ( worker.holding() )
        {
            /* if the worker is on 'hold', the simulation data was renewed,
             and it should be accessible for this thread to display */
            glApp::displayMain();
            worker.signal();
        }
    }
    else
    {
        // in replay mode, no need to lock the simulation data
        bool draw = simul.fresh_ || glApp::drawRefresh;
        if ( worker.read_input() )
            draw = true;
        
        if ( prop.replay > 0 )
        {
            // skip `prop.period` frames, and at least one if prop.period > 0
            for ( unsigned s = 0; s < prop.period; ++s )
                player.nextFrame();
            draw = true;
        }
        else if ( prop.replay < 0 )
        {
            player.previousFrame();
            draw = true;
        }
        
        if ( draw )
        {
            glApp::displayMain();
            //glApp::displayOtherWindows(drawSimulation);
        }
        else
            millisec = 100;
    }
    
    glutTimerFunc(millisec, timerCallback, 2);
}

