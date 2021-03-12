// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

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
#pragma mark - Input / Timer callback


void timerCallback(const int value)
{
    unsigned millisec = prop.delay;
    
    //thread.debug("timerCallback");
    if ( player.goLive && thread.alive() )
    {
        if ( prop.save_images && 0 == thread.trylock() )
        {
            player.drawScene(glApp::views[1]);
            player.saveView("image", prop.image_index++, prop.downsample);
            thread.unlock();
        }
        
        thread.signal();
    }
    else if ( prop.play )
    {
        if ( prop.save_images )
        {
            player.drawScene(glApp::views[1]);
            player.saveView("movie", thread.currentFrame(), prop.downsample);
        }

        if ( prop.play == 1 )
        {
            // skip prop.period frames, and at least one
            for ( unsigned s = 1; s < prop.period; ++s )
                player.nextFrame();
            player.nextFrame();
        }
        else if ( prop.play == -1 )
            player.previousFrame();
        
        if ( !prop.save_images )
            glApp::displayMain(); //glApp::postRedisplay();
    }
    else
        millisec = 200;
    
    glutTimerFunc(millisec, timerCallback, 2);
}

