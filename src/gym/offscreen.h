// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

#ifndef OFFSCREEN_H
#define OFFSCREEN_H


/// functions to open/close an OpenGL off-screen display
namespace OffScreen
{
    
    /// create OpenGL context suitable for offscreen rendering, returns 0 if success
    unsigned openContext();

    /// allocate display buffer, returns buffer name or 0 if failure
    unsigned makeBuffer(int width, int height, int multisample);

    /// create display buffer, and make it current
    unsigned openBuffer(int width, int height, int multisample);
    
    /// activate buffer
    void bindBuffer(unsigned);

    /// release display buffer
    void releaseBuffer();
    
    /// close the OpenGL context created by openContext()
    void closeContext();
    
}


#endif

