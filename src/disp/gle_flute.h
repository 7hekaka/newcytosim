// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#ifndef GLE_FLUTE_H
#define GLE_FLUTE_H

#include "flute.h"

/*
 Attention:
 Since the types depend on the dimensionality, this unit must be compiled
 separately for each different value of DIM.
 */

#if ( DIM >= 3 )
typedef flute3 fluteD;
typedef flute8 fluteD4;
#else
typedef flute2 fluteD;
typedef flute6 fluteD4;
#endif

namespace gle
{
    /// current buffer
    GLuint currStream();
    
    /// switch to next buffer
    GLuint nextStream();
    
    /// init buffers for data streaming
    void initStreams();
    
    /// release buffers
    void releaseStreams();
    
    /// map GPU buffer
    float* mapFloatBuffer(size_t cnt);
    /// unmap GPU buffer
    void unmapFloatBuffer(size_t vertex, size_t normals, size_t colors);
    /// bind GPU buffer
    void bindFloatBuffer(size_t vertex, size_t normals, size_t colors, size_t skip);
    
    /// map / unmap GPU buffer for 2D vertex + 4 color data
    inline flute2* mapBuffer200(size_t n) { return (flute2*)mapFloatBuffer(2*n); }
    inline void  unmapBuffer200() { unmapFloatBuffer(2, 0, 0); }

    /// map / unmap GPU buffer for 2D vertex + 4 color data
    inline flute6* mapBuffer204(size_t n) { return (flute6*)mapFloatBuffer(6*n); }
    inline void  unmapBuffer204() { unmapFloatBuffer(2, 0, 4); }
    
    /// map / unmap GPU buffer for 4D vertex + 4 color data
    inline flute4* mapBuffer400(size_t n) { return (flute4*)mapFloatBuffer(8*n); }
    inline void  unmapBuffer400() { unmapFloatBuffer(4, 0, 0); }

    /// map / unmap GPU buffer for 4D vertex + 4 color data
    inline flute8* mapBuffer404(size_t n) { return (flute8*)mapFloatBuffer(8*n); }
    inline void  unmapBuffer404() { unmapFloatBuffer(4, 0, 4); }
    
    /// map / unmap GPU buffer for 3D vertex + 3 normal data
    inline flute6* mapBuffer330(size_t n) { return (flute6*)mapFloatBuffer(6*n); }
    inline void unmapBuffer330() { unmapFloatBuffer(3, 3, 0); }
    
    /// map / unmap GPU buffer for vertex data only
    fluteD* mapBufferD00(size_t cnt);
    void  unmapBufferD00();
    void   bindBufferD00(size_t skip);

    /// map / unmap GPU buffer for vertex + color data
    fluteD4* mapBufferD04(size_t cnt);
    void  unmapBufferD04();
    void   bindVertexD04(size_t skip);

    unsigned* mapIndexBuffer(size_t cnt);
    void  unmapIndexBuffer();
    void   bindIndexBuffer(size_t skip);
};

#endif
