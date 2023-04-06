// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

#include "offscreen.h"
#include "opengl.h"


#if defined(__APPLE__) && defined(GL_VERSION_2_1)

// OpenGL Frame Buffer Objects
#include "offscreen_fbo.cc"

#elif defined(__linux)

// X-windows offscreen rendering routines (Linux)
#include "offscreen_glx.cc"

#else

// dummy routines
#include <cstdio>

unsigned OffScreen::openContext()
{
    //fprintf(stderr,"This program cannot render off-screen\n");
    return 1;
}

unsigned OffScreen::openBuffer(int, int, int)
{
    //fprintf(stderr,"This program cannot render off-screen\n");
    return 0;
}

void OffScreen::bindBuffer(unsigned) {}

void OffScreen::releaseBuffer() {}

void OffScreen::closeContext() {}

#endif
