// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include "gle.h"
#include "gle_flute.h"
#include "gle_color.h"

namespace gle
{
    /// number of buffers used to stream data to GPU
    const unsigned N_STREAMS = 32;
    
    /// index of current stream
    unsigned stream_indx = 0;
    
    /// OpenGL buffers objects for streaming
    GLuint stream_[N_STREAMS] = { 0 };

    void initStreams()
    {
        glGenBuffers(N_STREAMS, stream_);
        CHECK_GL_ERROR("glGenBuffers(-, streams_)");
    }
    
    void releaseStreams()
    {
        glDeleteBuffers(N_STREAMS, stream_);
        for (unsigned i=0; i<N_STREAMS; ++i) stream_[i] = 0;
    }
    
    GLuint currStream()
    {
        return stream_[stream_indx];
    }
    
    GLuint nextStream()
    {
        stream_indx = ( stream_indx + 1 ) % N_STREAMS;
        return stream_[stream_indx];
    }
    
    GLuint boundBuffer()
    {
        GLint i = 0;
        glGetIntegerv(GL_ARRAY_BUFFER_BINDING, &i);
        return i;
    }

    float* mapFloatBuffer(size_t cnt)
    {
        glBindBuffer(GL_ARRAY_BUFFER, nextStream());
        glBufferData(GL_ARRAY_BUFFER, cnt*sizeof(float), nullptr, GL_STREAM_DRAW);
        return (float*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
    }
    
    void setBufferV(size_t pts, size_t skip)
    {
        glVertexPointer(pts, GL_FLOAT, skip * pts * sizeof(float), nullptr);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

    void setBufferVN(size_t pts, size_t nor)
    {
        size_t tot = ( pts + nor ) * sizeof(float);
        glVertexPointer(pts, GL_FLOAT, tot, nullptr);
        if ( nor > 1 )
        {
            assert_true(!glIsEnabled(GL_NORMAL_ARRAY));
            glEnableClientState(GL_NORMAL_ARRAY);
            glNormalPointer(GL_FLOAT, tot, (void*)(pts*sizeof(float)));
        }
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }
    
    void setBufferCNV(size_t col, size_t nor, size_t pts, size_t skip)
    {
        assert_true(currStream() == boundBuffer());
        size_t tot = skip * ( pts + nor + col ) * sizeof(float);
        if ( col > 0 )
        {
            glEnableClientState(GL_COLOR_ARRAY);
            glColorPointer(4, GL_FLOAT, tot, nullptr);
        }
        else
        {
            assert_true(!glIsEnabled(GL_COLOR_ARRAY));
        }
        if ( nor > 1 )
        {
            glEnableClientState(GL_NORMAL_ARRAY);
            glNormalPointer(GL_FLOAT, tot, (void*)(col*sizeof(float)));
        }
        else
        {
            assert_true(!glIsEnabled(GL_NORMAL_ARRAY));
        }
        assert_true(glIsEnabled(GL_VERTEX_ARRAY));
        glVertexPointer(pts, GL_FLOAT, tot, (void*)((col+nor)*sizeof(float)));
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }
    
    unsigned* mapIndexBuffer(size_t cnt)
    {
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, nextStream());
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, cnt*sizeof(unsigned), nullptr, GL_STREAM_DRAW);
        return (unsigned*)glMapBuffer(GL_ELEMENT_ARRAY_BUFFER, GL_WRITE_ONLY);
    }
    
    void unmapIndexBuffer()
    {
        glUnmapBuffer(GL_ELEMENT_ARRAY_BUFFER);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    }

};
