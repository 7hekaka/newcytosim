// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include "gle.h"
#include "gym_flute.h"
#include "gym_check.h"
#include "gle_color.h"
#include "assert_macro.h"

namespace gym
{
    constexpr size_t Q = sizeof(float);

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
        glBufferData(GL_ARRAY_BUFFER, cnt*Q, nullptr, GL_STREAM_DRAW);
        void* ptr = glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
        assert_true(ptr);
        return (float*)ptr;
    }
    
    void setBufferV(size_t pts, size_t gap, size_t off)
    {
        glVertexPointer(pts, GL_FLOAT, pts*gap*Q, (void*)(off*Q));
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

    void setBufferV2T2()
    {
        assert_true(!glIsEnabled(GL_TEXTURE_COORD_ARRAY));
        glEnableClientState(GL_TEXTURE_COORD_ARRAY);
        glVertexPointer(2, GL_FLOAT, 4*Q, nullptr);
        glTexCoordPointer(2, GL_FLOAT, 4*Q, (void*)(2*Q));
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

    void setBufferVN(size_t pts, size_t nor)
    {
        assert_enabled(GL_VERTEX_ARRAY);
        size_t tot = ( pts + nor );
        glVertexPointer(std::min(pts, 3UL), GL_FLOAT, tot*Q, nullptr);
        if ( nor > 1 )
        {
            assert_true(!glIsEnabled(GL_NORMAL_ARRAY));
            glEnableClientState(GL_NORMAL_ARRAY);
            glNormalPointer(GL_FLOAT, tot*Q, (void*)(pts*Q));
        }
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }
    
    void bindBufferV2(GLuint buf, size_t off)
    {
        glBindBuffer(GL_ARRAY_BUFFER, buf);
        setBufferV(2, 1, off);
    }
    
    void bindBufferV3(GLuint buf, size_t off)
    {
        glBindBuffer(GL_ARRAY_BUFFER, buf);
        setBufferV(3, 1, off);
    }

    void bindBufferV3N3(GLuint buf)
    {
        glBindBuffer(GL_ARRAY_BUFFER, buf);
        setBufferVN(3, 3);
    }
    
    // normal = position
    void setBufferV3N0(GLsizei first)
    {
        assert_enabled(GL_VERTEX_ARRAY);
        glVertexPointer(3, GL_FLOAT, 3*Q, (void*)(first*Q));
        glEnableClientState(GL_NORMAL_ARRAY);
        glNormalPointer(GL_FLOAT, 3*Q, (void*)(first*Q));
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

    // normal = position
    void bindBufferV3N0(GLuint buf, GLsizei first)
    {
        glBindBuffer(GL_ARRAY_BUFFER, buf);
        setBufferV3N0(first);
    }

    void setBufferCV(size_t col, size_t pts, size_t gap)
    {
        assert_true(currStream() == boundBuffer());
        size_t tot = ( col + pts ) * gap;
        if ( col > 0 )
        {
            glEnableClientState(GL_COLOR_ARRAY);
            glColorPointer(4, GL_FLOAT, tot*Q, nullptr);
        }
        else
        {
            assert_true(!glIsEnabled(GL_COLOR_ARRAY));
        }
        assert_enabled(GL_VERTEX_ARRAY);
        glVertexPointer(std::min(pts, 3UL), GL_FLOAT, tot*Q, (void*)(col*Q));
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

    void setBufferCNV(size_t col, size_t nor, size_t pts, size_t gap)
    {
        assert_true(currStream() == boundBuffer());
        size_t tot = ( col + pts + nor ) * gap;
        if ( col > 0 )
        {
            glEnableClientState(GL_COLOR_ARRAY);
            glColorPointer(4, GL_FLOAT, tot*Q, nullptr);
        }
        else
        {
            assert_true(!glIsEnabled(GL_COLOR_ARRAY));
        }
        if ( nor > 1 )
        {
            glEnableClientState(GL_NORMAL_ARRAY);
            glNormalPointer(GL_FLOAT, tot*Q, (void*)(col*Q));
        }
        else
        {
            assert_true(!glIsEnabled(GL_NORMAL_ARRAY));
        }
        assert_enabled(GL_VERTEX_ARRAY);
        glVertexPointer(pts, GL_FLOAT, tot*Q, (void*)((col+nor)*Q));
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }
    
    unsigned short* mapIndexBuffer(size_t cnt)
    {
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, nextStream());
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, cnt*sizeof(short), nullptr, GL_STREAM_DRAW);
        return (unsigned short*)glMapBuffer(GL_ELEMENT_ARRAY_BUFFER, GL_WRITE_ONLY);
    }
    
    void unmapIndexBuffer()
    {
        glUnmapBuffer(GL_ELEMENT_ARRAY_BUFFER);
        //glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

};
