// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include "gle.h"
#include "gle_flute.h"
#include "opengl.h"
#include "gle_color.h"

namespace gle
{
    /// defined in gle.cc
    extern GLuint stream_[4];
    
    static GLuint boundBuffer()
    {
        GLint i = 0;
        glGetIntegerv(GL_ARRAY_BUFFER_BINDING, &i);
        return (GLuint)i;
    }
    
    fluteV* mapVertexBuffer(size_t cnt)
    {
        glBindBuffer(GL_ARRAY_BUFFER, stream_[0]);
        glBufferData(GL_ARRAY_BUFFER, cnt*sizeof(fluteV), nullptr, GL_STREAM_DRAW);
        return (fluteV*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
    }
    
    void unmapVertexBuffer()
    {
        assert_true(stream_[0] == boundBuffer());
        //glBindBuffer(GL_ARRAY_BUFFER, stream_[0]);
        glUnmapBuffer(GL_ARRAY_BUFFER);
        glVertexPointer((DIM>2?3:2), GL_FLOAT, 0, nullptr);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }
    
    void bindVertexBuffer(size_t skip)
    {
        glBindBuffer(GL_ARRAY_BUFFER, stream_[0]);
        glVertexPointer((DIM>2?3:2), GL_FLOAT, skip*sizeof(fluteV), nullptr);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }
    
    fluteVN* mapVertexNormalBuffer(size_t cnt)
    {
        glBindBuffer(GL_ARRAY_BUFFER, stream_[1]);
        glBufferData(GL_ARRAY_BUFFER, cnt*sizeof(fluteVN), nullptr, GL_STREAM_DRAW);
        return (fluteVN*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
    }
    
    void unmapVertexNormalBuffer()
    {
        assert_true(stream_[1] == boundBuffer());
        //glBindBuffer(GL_ARRAY_BUFFER, stream_[1]);
        glUnmapBuffer(GL_ARRAY_BUFFER);
        glNormalPointer(GL_FLOAT, 0, nullptr);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }
    
    void bindVertexNormalBuffer(size_t skip)
    {
        glBindBuffer(GL_ARRAY_BUFFER, stream_[1]);
        glNormalPointer(GL_FLOAT, skip*sizeof(fluteVN), nullptr);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }
    
    fluteVC* mapVertexColorBuffer(size_t cnt)
    {
        glBindBuffer(GL_ARRAY_BUFFER, stream_[2]);
        glBufferData(GL_ARRAY_BUFFER, cnt*sizeof(fluteVC), nullptr, GL_STREAM_DRAW);
        return (fluteVC*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
    }
    
    void unmapVertexColorBuffer()
    {
        assert_true(stream_[2] == boundBuffer());
        //glBindBuffer(GL_ARRAY_BUFFER, stream_[2]);
        glUnmapBuffer(GL_ARRAY_BUFFER);
        glVertexPointer((DIM>2?3:2), GL_FLOAT, sizeof(fluteVC), nullptr);
        glColorPointer(4, GL_FLOAT, sizeof(fluteVC), (void*)(DIM>2?0x10:0x8));
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }
    
    void bindVertexColorBuffer(size_t skip)
    {
        glBindBuffer(GL_ARRAY_BUFFER, stream_[2]);
        glVertexPointer((DIM>2?3:2), GL_FLOAT, skip*sizeof(fluteVC), nullptr);
        glColorPointer(4, GL_FLOAT, skip*sizeof(fluteVC), (void*)(DIM>2?0x10:0x8));
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }
    
};
