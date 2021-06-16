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
    
    void bindFloatBuffer(size_t pts, size_t nor, size_t col, size_t skip)
    {
        glBindBuffer(GL_ARRAY_BUFFER, currStream());
        size_t tot = skip * (pts+nor+col) * sizeof(float);
        glVertexPointer(pts, GL_FLOAT, tot, nullptr);
        if ( nor > 1 )
        {
            glEnableClientState(GL_NORMAL_ARRAY);
            glNormalPointer(GL_FLOAT, tot, (void*)(pts*sizeof(float)));
        }
        if ( col > 0 )
        {
            glEnableClientState(GL_COLOR_ARRAY);
            glColorPointer(4, GL_FLOAT, tot, (void*)((pts+nor)*sizeof(float)));
        }
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

    void unmapFloatBuffer(size_t pts, size_t nor, size_t col)
    {
        glUnmapBuffer(GL_ARRAY_BUFFER);
        size_t tot = (pts+nor+col) * sizeof(float);
        glVertexPointer(pts, GL_FLOAT, tot, nullptr);
        if ( nor > 1 )
        {
            glEnableClientState(GL_NORMAL_ARRAY);
            glNormalPointer(GL_FLOAT, tot, (void*)(pts*sizeof(float)));
        }
        if ( col > 0 )
        {
            glEnableClientState(GL_COLOR_ARRAY);
            glColorPointer(4, GL_FLOAT, tot, (void*)((pts+nor)*sizeof(float)));
        }
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

    
    fluteD* mapBufferD00(size_t cnt)
    {
        glBindBuffer(GL_ARRAY_BUFFER, nextStream());
        glBufferData(GL_ARRAY_BUFFER, cnt*sizeof(fluteD), nullptr, GL_STREAM_DRAW);
        return (fluteD*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
    }
    
    void unmapBufferD00()
    {
        glUnmapBuffer(GL_ARRAY_BUFFER);
        glVertexPointer((DIM>2?3:2), GL_FLOAT, 0, nullptr);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }
    
    void bindBufferD00(size_t skip)
    {
        glBindBuffer(GL_ARRAY_BUFFER, currStream());
        glVertexPointer((DIM>2?3:2), GL_FLOAT, skip*sizeof(fluteD), nullptr);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }
    
    fluteD4* mapBufferD04(size_t cnt)
    {
        glBindBuffer(GL_ARRAY_BUFFER, nextStream());
        glBufferData(GL_ARRAY_BUFFER, cnt*sizeof(fluteD4), nullptr, GL_STREAM_DRAW);
        return (fluteD4*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
    }
    
    void unmapBufferD04()
    {
        glUnmapBuffer(GL_ARRAY_BUFFER);
        glVertexPointer((DIM>2?3:2), GL_FLOAT, sizeof(fluteD4), nullptr);
        glColorPointer(4, GL_FLOAT, sizeof(fluteD4), (void*)((DIM>2?4:2)*sizeof(float)));
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }
    
    void bindVertexD04(size_t skip)
    {
        glBindBuffer(GL_ARRAY_BUFFER, currStream());
        glVertexPointer((DIM>2?3:2), GL_FLOAT, skip*sizeof(fluteD4), nullptr);
        glColorPointer(4, GL_FLOAT, skip*sizeof(fluteD4), (void*)((DIM>2?4:2)*sizeof(float)));
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
