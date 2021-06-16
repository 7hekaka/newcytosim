// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include "gle.h"
#include "gle_flute.h"
#include "gle_color.h"

namespace gle
{
    /// number of buffers used to stream data to GPU
    const unsigned N_STREAMS = 128;
    
    /// index of current stream
    unsigned stream_indx = 0;
    
    /// OpenGL buffers objects for streaming
    GLuint streams_[N_STREAMS] = { 0 };

    void initStreams()
    {
        if ( !glIsBuffer(streams_[0]) )
            glGenBuffers(N_STREAMS, streams_);
        CHECK_GL_ERROR("glGenBuffers(-, streams_)");
    }
    
    void releaseStreams()
    {
        glDeleteBuffers(N_STREAMS, streams_);
        streams_[0] = 0;
    }
    
    GLuint currStream()
    {
        return streams_[stream_indx];
    }
    
    GLuint nextStream()
    {
        stream_indx = ( stream_indx + 1 ) % N_STREAMS;
        return streams_[stream_indx];
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
            glEnableClientState(GL_NORMAL_ARRAY);
            glNormalPointer(GL_FLOAT, tot, (void*)(pts*sizeof(float)));
        }
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }
    
    void setBufferCNV(size_t col, size_t nor, size_t pts, size_t skip)
    {
        size_t tot = skip * ( pts + nor + col) * sizeof(float);
        if ( col > 0 )
        {
            glEnableClientState(GL_COLOR_ARRAY);
            glColorPointer(4, GL_FLOAT, tot, nullptr);
        }
        if ( nor > 1 )
        {
            glEnableClientState(GL_NORMAL_ARRAY);
            glNormalPointer(GL_FLOAT, tot, (void*)(pts*sizeof(float)));
        }
        glVertexPointer(pts, GL_FLOAT, tot, (void*)((pts+nor)*sizeof(float)));
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
