// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

#include <cstdio>
#include "assert_macro.h"
#include "opengl.h"


namespace gym
{
    
    /// convert OpenGL error code to string
    const char* errorString(unsigned code);
    
    /// check and print OpenGL error(s)
    void reportErrors(FILE*, const char* msg);
    
    /// print current color properties of OpenGL context
    void printColors(FILE*);
    
    /// print some info for debugging purpose
    void printCaps(const char[]);
    
    /// print OpenGL matrices for debugging purpose
    void printMatrices(FILE*);
    
    /// print depth-buffer range
    void printDepthRange(int X, int Y, int W, int H);
    
};

#ifdef NDEBUG
#  define CHECK_GL_ERROR(ARG) ((void) 0)
#  define assert_enabled(ARG) ((void) 0)
#  define assertLighting() {}
#  define assertCullFace() {}
#  define assertVertexArray() {}
#else
#  define CHECK_GL_ERROR(ARG) gym::reportErrors(stderr, ARG)
#  define assert_enabled(CAP) { if (!glIsEnabled(CAP))\
 { fprintf(stderr, "%s is not enabled in `%s`,  %s:%d\n", #CAP, SFUNC, SFILE, __LINE__); }}
#  define assertLighting() { assert_enabled(GL_LIGHTING); }
#  define assertCullFace() { assert_enabled(GL_CULL_FACE); }
#  define assertVertexArray() { assert_enabled(GL_VERTEX_ARRAY); }
#endif
