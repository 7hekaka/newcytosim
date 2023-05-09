// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University
/** 
 off-screen rendering using OpenGL Frame Buffer Objects
 
 http://en.wikipedia.org/wiki/Framebuffer_Object
 
 https://developer.apple.com/library/mac/#documentation/graphicsimaging/Conceptual/OpenGL-MacProgGuide/opengl_offscreen/opengl_offscreen.html
*/

#include <vector>
#include <cstdio>

#include <OpenGL/OpenGL.h>
#include <OpenGL/CGLTypes.h>


CGLContextObj preContext = nullptr;
CGLContextObj offContext = nullptr;

/// Holds the pointers of a set of OpenGL buffers.
class GLBuffer
{
public:
    
    GLuint width;
    GLuint height;
    GLuint frameBuffer;
    GLuint colorBuffer;
    GLuint depthBuffer;
    
    GLBuffer(GLuint W, GLuint H)
    {
        width = W;
        height = H;
        frameBuffer = 0;
        colorBuffer = 0;
        depthBuffer = 0;
    }
    
    /// allocate graphical memory
    void create()
    {
        glGenFramebuffers(1, &frameBuffer);
        glGenRenderbuffers(1, &colorBuffer);
        glGenRenderbuffers(1, &depthBuffer);
        //printf("offscreen- created framebuffer %u\n", frameBuffer);
    }
    
    /// delete graphical memory
    void release()
    {
        glDeleteRenderbuffers(1, &colorBuffer);
        glDeleteRenderbuffers(1, &depthBuffer);
        glDeleteFramebuffers(1, &frameBuffer);
        //printf("offscreen- released framebuffer %u\n", frameBuffer);
        colorBuffer = 0;
        depthBuffer = 0;
        frameBuffer = 0;
        //glPopAttrib();
    }
    
    /// enable buffer for both reading and writing operations
    void bind() const
    {
        glBindFramebuffer(GL_FRAMEBUFFER, frameBuffer);
        //printf("offscreen- bind framebuffer %u : %ix%i\n", frameBuffer, width, height);
        glViewport(0, 0, width, height);
#if ( 0 )
        GLint readbuf = 0, drawbuf = 0;
        glGetIntegerv(GL_READ_BUFFER, &readbuf);
        glGetIntegerv(GL_DRAW_BUFFER, &drawbuf);
        printf("framebuffers: read %i draw %i\n", readbuf, drawbuf);
#endif
   }
};

/// list of known buffers
std::vector<GLBuffer> all_buffers;

//------------------------------------------------------------------------------
char const* glErrorString(GLenum code)
{
    switch ( code )
    {
        case GL_NO_ERROR:          return "GL_NO_ERROR";
        case GL_INVALID_ENUM:      return "GL_INVALID_ENUM";
        case GL_INVALID_VALUE:     return "GL_INVALID_VALUE";
        case GL_INVALID_OPERATION: return "GL_INVALID_OPERATION";
        case GL_STACK_OVERFLOW:    return "GL_STACK_OVERFLOW";
        case GL_STACK_UNDERFLOW:   return "GL_STACK_UNDERFLOW";
        case GL_OUT_OF_MEMORY:     return "GL_OUT_OF_MEMORY";
        case GL_TABLE_TOO_LARGE:   return "GL_TABLE_TOO_LARGE";
        default:                   return "GL_UNKNOWN_ERROR";
    }
}


static void checkGLError(const char msg[])
{
    GLenum e = glGetError();
    while ( e != GL_NO_ERROR )
    {
        fprintf(stderr, "OpenGL error `%s' %s\n", glErrorString(e), msg);
        e = glGetError();
    }
}


void describePixelFormat(CGLPixelFormatObj const& p, GLint n)
{
    GLint cap = 0;
    printf("pixel format %i:\n", n);
    CGLDescribePixelFormat(p, n, kCGLPFAColorSize, &cap);
    printf("    colors      %i\n", cap);
    CGLDescribePixelFormat(p, n, kCGLPFADepthSize, &cap);
    printf("    depth       %i\n", cap);
    CGLDescribePixelFormat(p, n, kCGLPFAMultisample, &cap);
    printf("    multi       %i\n", cap);
    if (cap)
    {
        CGLDescribePixelFormat(p, n, kCGLPFASampleBuffers, &cap);
        printf("       buffers  %i\n", cap);
        CGLDescribePixelFormat(p, n, kCGLPFASamples, &cap);
        printf("       samples  %i\n", cap);
    }
}

//------------------------------------------------------------------------------

CGLContextObj createOffscreenContext()
{
    CGLPixelFormatAttribute attribs[] =
    {
        kCGLPFAMinimumPolicy,
        //kCGLPFAPBuffer,
        //kCGLPFARendererID,
        kCGLPFAColorSize,
        (CGLPixelFormatAttribute)32,
        kCGLPFADepthSize,
        (CGLPixelFormatAttribute)16,
        kCGLPFAMultisample,
        kCGLPFASamples,
        (CGLPixelFormatAttribute)4,
        (CGLPixelFormatAttribute)0
    };

    GLint npix = 0;
    CGLError err;
    CGLPixelFormatObj pix;
    
    err = CGLChoosePixelFormat(attribs, &pix, &npix);
    
    if ( err || npix == 0 )
    {
        fprintf(stderr, "could not find suitable pixel format\n");
        return nullptr;
    }
    
#if ( 0 )
    // display supported pixel formats:
    for ( int n = 0; n < npix; ++n )
    {
        describePixelFormat(pix, n);
        GLint cap = 0;
        CGLDescribePixelFormat(pix, n, kCGLPFAMultisample, &cap);
    }
#endif
    
    CGLContextObj newContext;
    CGLCreateContext(pix, preContext, &newContext);
    CGLReleasePixelFormat(pix);
    
    if ( !newContext )
    {
        fprintf(stderr, "Could not create OpenGL context\n");
        return nullptr;
    }
    
    checkGLError("createOffscreenContext()");
    return newContext;
}


/**
 Set up an OpenGL context for offscreen rendering.
 */
unsigned OffScreen::openContext()
{
    if ( offContext && offContext == CGLGetCurrentContext() )
    {
        fprintf(stderr, "Offscreen context is already active\n");
        return 1;
    }
    
    preContext = CGLGetCurrentContext();
    //printf("offscreen- current context %p\n", preContext);

    if ( !offContext )
        offContext = createOffscreenContext();

    if ( 0 == CGLSetCurrentContext(offContext) )
    {
        //printf("offscreen- switched context %p\n", offContext);
    }
    else
    {
        CGLDestroyContext(offContext);
        preContext = nullptr;
        offContext = nullptr;
        fprintf(stderr, "Could not switch OpenGL context\n");
        return 2;
    }

    return 0;
}


void OffScreen::closeContext()
{
    if ( offContext == CGLGetCurrentContext() )
    {
        while ( !all_buffers.empty() )
        {
            all_buffers.back().release();
            all_buffers.pop_back();
        }
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        CGLSetCurrentContext(preContext);
        CGLDestroyContext(offContext);
        //printf("offscreen- restored context %p\n", preContext);
        preContext = nullptr;
        offContext = nullptr;
    }
    checkGLError("OffScreen::closeContext()");
}


//------------------------------------------------------------------------------

/**
 Set up a Frame Buffer object with two Render buffers attached,
 for color and depth data.
 */
static unsigned makeBuffer(const int width, const int height, int multisample)
{
    GLint s = 0;
    glGetIntegerv(GL_MAX_SAMPLES, &s);
    if ( multisample > s )
    {
        multisample = s;
        //fprintf(stderr, "warning: OpenGL multisamples limited to %i\n", s);
    }
    
    //if ( multisample > 1 )
    //    fprintf(stderr, "OpenGL buffer is multisampled %i\n", multisample);

    //Set up a FBO with two renderBuffer attachment
    GLBuffer buffer(width, height);
    buffer.create();
    
    glBindFramebuffer(GL_FRAMEBUFFER, buffer.frameBuffer);
    //checkGLError("OffScreen:glBindFramebuffer()");
    
    glBindRenderbuffer(GL_RENDERBUFFER, buffer.colorBuffer);
    if ( multisample > 1 )
        glRenderbufferStorageMultisample(GL_RENDERBUFFER, multisample, GL_RGBA8, width, height);
    else
        glRenderbufferStorage(GL_RENDERBUFFER, GL_RGBA8, width, height);
    
    glBindRenderbuffer(GL_RENDERBUFFER, buffer.depthBuffer);
    if ( multisample > 1 )
        glRenderbufferStorageMultisample(GL_RENDERBUFFER, multisample, GL_DEPTH24_STENCIL8, width, height);
    else
        glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, width, height);
    
    GLenum err = glGetError();
    if ( err != GL_NO_ERROR )
    {
        fprintf(stderr, "Offscreen OpenGL error `%s'\n", glErrorString(err));
        fprintf(stderr, " (max. buffer size is %i)\n", GL_MAX_RENDERBUFFER_SIZE);
        return 0;
    }
    
    glBindRenderbuffer(GL_RENDERBUFFER, 0);
    
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, buffer.colorBuffer);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, buffer.depthBuffer);
    //checkGLError("OffScreen::glFramebufferRenderbuffer()");

    GLenum status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
    
    if ( status == GL_FRAMEBUFFER_UNSUPPORTED )
        return 0;
    if ( status != GL_FRAMEBUFFER_COMPLETE )
        return 0;
    
    //printf("offscreen- allocated framebuffer %u  %ux%u\n", buffer.frameBuffer, width, height);
    checkGLError("OffScreen:makeBuffer()");
    
    all_buffers.push_back(buffer);

    return buffer.frameBuffer;
}


unsigned OffScreen::openBuffer(const int width, const int height, int multisample)
{
    GLuint buf = makeBuffer(width, height, multisample);
    glBindFramebuffer(GL_FRAMEBUFFER, buf);
    glViewport(0, 0, width, height);
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    return buf;
}


void OffScreen::bindBuffer(unsigned id)
{
    for ( GLBuffer const& buf : all_buffers )
    {
        if ( buf.frameBuffer == id )
        {
            buf.bind();
            glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
            break;
        }
    }
}


/**
 Release the last buffer that was created, and promote the one before it
 */
void OffScreen::releaseBuffer()
{
    if ( !all_buffers.empty() )
    {
        all_buffers.back().release();
        all_buffers.pop_back();
        if ( all_buffers.empty() )
        {
            glBindFramebuffer(GL_FRAMEBUFFER, 0);
        }
        else
            all_buffers.back().bind();
    }
    checkGLError("OffScreen::close()");
}

