// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#ifndef GLE_FLUTE_H
#define GLE_FLUTE_H

#include "flute.h"
#include "dim.h"

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
    /// define buffer layout
    void setBufferV(size_t vertex, size_t skip = 1);

    /// define buffer layout
    void setBufferVN(size_t normals, size_t vertex);

    /// define buffer layout
    void setBufferCNV(size_t colors, size_t normals, size_t vertex, size_t skip = 1);

    
    /// current buffer
    GLuint currStream();
    
    /// used in debug mode
    GLuint boundBuffer();

    /// switch to next buffer
    GLuint nextStream();
    
    /// init buffers for data streaming
    void initStreams();
    
    /// release buffers
    void releaseStreams();
    
    /// map GPU buffer
    float* mapFloatBuffer(size_t cnt);
    /// private
    inline void unmap() { glUnmapBuffer(GL_ARRAY_BUFFER); }
    /// private
    inline void bind() { glBindBuffer(GL_ARRAY_BUFFER, currStream()); }
    
    /// map / unmap GPU buffer for 2D vertex
    inline flute2* mapBufferV2(size_t n) { return (flute2*)mapFloatBuffer(2*n); }
    inline void  unmapBufferV2() { unmap(); setBufferV(2); }
    
    /// map / unmap GPU buffer for 3D vertex
    inline flute3* mapBufferV3(size_t n) { return (flute3*)mapFloatBuffer(3*n); }
    inline void  unmapBufferV3() { unmap(); setBufferV(3); }

    /// map / unmap GPU buffer for 4 color data + 2D vertex
    inline flute6* mapBufferC4V2(size_t n) { return (flute6*)mapFloatBuffer(6*n); }
    inline void  unmapBufferC4V2() { unmap(); setBufferCNV(4, 0, 2); }

    /// map / unmap GPU buffer for 4 color data + 4D vertex
    inline flute8* mapBufferC4V4(size_t n) { return (flute8*)mapFloatBuffer(8*n); }
    inline void  unmapBufferC4V4() { unmap(); setBufferCNV(4, 0, 4); }
    
    /// map / unmap GPU buffer for 3D vertex + 3 normal data
    inline flute6* mapBufferV3N3(size_t n) { return (flute6*)mapFloatBuffer(6*n); }
    inline void unmapBufferV3N3() { unmap(); setBufferVN(3, 3); }
    
    /// map / unmap GPU buffer for vertex data only
    inline fluteD* mapBufferVD(size_t n) { return (fluteD*)mapFloatBuffer((DIM>2?3:2)*n); }
    inline void  unmapBufferVD() { unmap(); setBufferV((DIM>2?3:2)); }
    inline void   bindBufferVD(size_t skip) { bind(); setBufferV((DIM>2?3:2), skip); }

    /// map / unmap GPU buffer for color data + vertex
    inline fluteD4* mapBufferC4VD(size_t n) { return (fluteD4*)mapFloatBuffer((DIM>2?8:6)*n); }
    inline void   unmapBufferC4VD() { unmap(); setBufferCNV(4, 0, (DIM>2?4:2)); }
    inline void    bindBufferC4VD(size_t skip) { bind(); setBufferCNV(4, 0, (DIM>2?4:2), skip); }

    unsigned* mapIndexBuffer(size_t n);
    void  unmapIndexBuffer();
    void   bindIndexBuffer(size_t skip);
};

#endif
