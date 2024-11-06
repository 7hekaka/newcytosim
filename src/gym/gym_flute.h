// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#ifndef GYM_FLUTE_H
#define GYM_FLUTE_H

#include "opengl.h"
#include "flute.h"

namespace gym
{
    /// define buffer layout
    void setBufferV(size_t size, size_t stride = 1, size_t off = 0);

    /// define buffer layout
    void setBufferVN(size_t normals, size_t vertex);

    /// define buffer layout
    void setBufferCV(size_t colors, size_t vertex, size_t stride = 1);

    /// define buffer layout
    void setBufferCNV(size_t colors, size_t normals, size_t vertex, size_t stride);
    
    
    /// define buffer layout for a Device buffer
    void bindBufferV2(GLuint, size_t off = 0);
    
    /// define buffer layout for a Device buffer
    void bindBufferV3(GLuint, size_t stride = 1, size_t off = 0);

    /// bind buffer with 3-position coordinates, 3 normal coordinates
    void bindBufferV3N3(GLuint);
    
    /// define buffer layout for a Device buffer, with position == normal
    void setBufferV3N0(GLsizei first);
    
    /// define buffer layout for a Device buffer, with position == normal
    void bindBufferV3N0(GLuint, GLsizei first);

    
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
    inline void rebind() { glBindBuffer(GL_ARRAY_BUFFER, currStream()); }

    /// it is not necessary to unbind buffer, but might be useful for debugging
    inline void unbind1() { } //glBindBuffer(GL_ARRAY_BUFFER, 0); }
    
    /// it is not necessary to unbind buffer, but might be useful for debugging
    inline void unbind2() { } //glBindBuffer(GL_ARRAY_BUFFER, 0); glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0); }

    ///
    inline void subdata(size_t size, float const* ptr) { glBufferSubData(currStream(), 0, size, ptr); }
    
    /// map / unmap GPU buffer for 2D vertex
    inline flute2* mapBufferV2(size_t n) { return (flute2*)mapFloatBuffer(2*n); }
    inline void  unmapBufferV2() { unmap(); setBufferV(2); }
    inline void rebindBufferV2() { rebind(); setBufferV(2); }

    void setBufferVxT2(size_t, size_t);
    
    inline void enableTexture(GLuint tex) { glEnable(GL_TEXTURE_2D); glBindTexture(GL_TEXTURE_2D, tex); }
    inline void cleanupTexture() { glDisableClientState(GL_TEXTURE_COORD_ARRAY); glDisable(GL_TEXTURE_2D); }

    inline flute4* mapBufferV2T2(size_t n) { return (flute4*)mapFloatBuffer(4*n); }
    inline void  unmapBufferV2T2() { unmap(); setBufferVxT2(2, 4); }
    
    inline flute6* mapBufferV4T2(size_t n) { return (flute6*)mapFloatBuffer(6*n); }
    inline void  unmapBufferV4T2() { unmap(); setBufferVxT2(4, 6); }

    /// map / unmap GPU buffer for 3D vertex
    inline flute3* mapBufferV3(size_t n) { return (flute3*)mapFloatBuffer(3*n); }
    inline void  unmapBufferV3() { unmap(); setBufferV(3); }
    inline void rebindBufferV3(size_t stride, size_t off) { rebind(); setBufferV(3, stride, 3*off); }

    inline void unmapBufferV3N0() { unmap(); setBufferV3N0(0); }
    
    /// map / unmap GPU buffer for 4 color data + 2D vertex
    inline flute6* mapBufferC4V2(size_t n) { return (flute6*)mapFloatBuffer(6*n); }
    inline void  unmapBufferC4V2() { unmap(); setBufferCV(4, 2); }

    /// map / unmap GPU buffer for 4 color data + 4D vertex
    inline flute8* mapBufferC4V4(size_t n) { return (flute8*)mapFloatBuffer(8*n); }
    inline void  unmapBufferC4V4() { unmap(); setBufferCV(4, 4); }
    
    /// map / unmap GPU buffer for 3D vertex + 3 normal data
    inline flute6* mapBufferV3N3(size_t n) { return (flute6*)mapFloatBuffer(6*n); }
    inline void unmapBufferV3N3() { unmap(); setBufferVN(3, 3); }

    inline void cleanupV() { }
    inline void cleanupVN() { glDisableClientState(GL_NORMAL_ARRAY); }
    inline void cleanupCV() { glDisableClientState(GL_COLOR_ARRAY); glDisableClientState(GL_NORMAL_ARRAY); }
    
    unsigned short* mapIndexBuffer(size_t n);
    void unmapIndexBuffer();
    void bindIndexBuffer(size_t stride);
};

#endif
