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

    thread.lock();
    if ( thread.selectClosestHandle(pos, range) )
        thread.moveHandle(pos);
    else
    {
        if ( thread.handle() )
        {
            thread.detachHandle();
            thread.moveHandle(pos);
        }
        else
        {
            Single * s = thread.createHandle(pos, range);
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
    thread.unlock();
}


/**
 Processes mouse motion
 */
void processMouseDrag(int, int, Vector3& ori3, const Vector3& pos3, int mode)
{
    Vector pos(pos3);
    Vector ori(ori3);
    
    thread.lock();
    if ( mode )
    {
        thread.moveHandles(pos-ori);
        ori3 = pos3;
    }
    else
        thread.moveHandle(pos);
    thread.unlock();
}


//------------------------------------------------------------------------------
#pragma mark - Timer callback


void timerCallback(const int value)
{
    unsigned millisec = prop.delay;

    //std::clog << "timerCallback @ " << std::fixed << simul.time() << "s\n";
    if ( prop.goLive && thread.alive() )
    {
        //thread.debug("timerCallback");
        thread.signal();
    }
    else
    {
        // in replay mode, no need to lock the simulation state
        if ( thread.executePipedCommands(32) )
            glApp::postRedisplay();
        
        if ( prop.play > 0 )
        {
            // skip prop.period frames, and at least one
            for ( unsigned s = 1; s < prop.period; ++s )
                player.nextFrame();
            player.nextFrame();
            glApp::postRedisplayAll();
        }
        else if ( prop.play < 0 )
        {
            player.previousFrame();
            glApp::postRedisplayAll();
        }
        else
            millisec = 200;
    }
    
    glutTimerFunc(millisec, timerCallback, 2);
}

