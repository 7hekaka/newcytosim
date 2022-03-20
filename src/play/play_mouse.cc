// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University

/**
 Callback for mouse clicks
 */
void processMouseClick(int, int, const Vector3& pos3, int)
{
    // distance in pixels where mouse-Hand binds:
    const real pixrad = 6;
    const real range = pixrad * glApp::currentView().pixelSize();
    Vector pos(pos3);

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
            pd = static_cast<PointDisp*>(player.dispList.find("hand:display", "user_hand"));
            if ( !pd )
            {
                pd = new PointDisp("hand:display", "user_hand");
                pd->size   = 2 * pixrad;
                pd->color  = glApp::currentView().front_color;
                pd->color2 = pd->color;
                player.dispList.deposit(pd);
            }
        }
    }
    worker.unlock();
}


/**
 Processes mouse motion
 */
void processMouseDrag(int, int, Vector3& ori3, const Vector3& pos3, int mode)
{
    Vector pos(pos3);
    Vector ori(ori3);
    
    worker.lock();
    if ( mode )
    {
        worker.moveHandles(pos-ori);
        ori3 = pos3;
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

    //std::clog << "timerCallback @ " << std::fixed << simul.time() << "s\n";
    if ( prop.goLive && worker.alive() )
    {
        //worker.debug("timerCallback");
        worker.signal();
    }
    else
    {
        // in replay mode, no need to lock the simulation state
        if ( worker.executePipedCommands(32) )
            glApp::postRedisplay();
        
        if ( prop.play > 0 )
        {
            // skip prop.period frames, and at least one
            for ( unsigned s = 1; s < prop.period; ++s )
                player.nextFrame();
            player.nextFrame();
            glApp::displayAll();
        }
        else if ( prop.play < 0 )
        {
            player.previousFrame();
            glApp::displayAll();
        }
        else
            millisec = 100;
    }
    
    glutTimerFunc(millisec, timerCallback, 2);
}

