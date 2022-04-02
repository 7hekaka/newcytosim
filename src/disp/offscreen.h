// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

#ifndef OFFSCREEN_H
#define OFFSCREEN_H


/// functions to open/close an OpenGL off-screen display
namespace OffScreen
{
    
    /// create OpenGL context suitable for offscreen rendering, returns 0 if success
    int openContext();

    /// allocate display buffer, returns buffer name or 0 if failure
    int createBuffer(int width, int height, int multisample);

    /// copy color data from 'back' to 'front'
    void blitBuffers(unsigned back, unsigned front, int W, int H);

    /// release display buffer
    void releaseBuffer();
    
    /// close the OpenGL context created by openContext()
    void closeContext();
    
}


#endif

