// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

#include "offscreen.h"
#include "opengl.h"


/// copy color data from 'back' to 'front'
void OffScreen::blitBuffers(unsigned back, unsigned front, int W, int H)
{
    //std::clog << "blitting multisample buffer\n";
    glBindFramebuffer(GL_READ_FRAMEBUFFER, back);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, front);
    glBlitFramebuffer(0, 0, W, H, 0, 0, W, H, GL_COLOR_BUFFER_BIT, GL_NEAREST);
    glBindFramebuffer(GL_READ_FRAMEBUFFER, front);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, back);
}


#if defined(__APPLE__) && defined(GL_VERSION_2_1)

// OpenGL Frame Buffer Objects
#include "offscreen_fbo.cc"

#elif defined(__linux)

// X-windows offscreen rendering routines (Linux)
#include "offscreen_glx.cc"

#else

// dummy routines
#include <cstdio>

int OffScreen::openContext()
{
    //fprintf(stderr,"This program cannot render off-screen\n");
    return 1;
}

int OffScreen::createBuffer(int, int, int)
{
    //fprintf(stderr,"This program cannot render off-screen\n");
    return 0;
}

void OffScreen::releaseBuffer()
{
}

void OffScreen::closeContext()
{
}

#endif
