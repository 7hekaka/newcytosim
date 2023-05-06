// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University

#include "save_image_gl.h"
#include "gym_view.h"
#include <cstdlib>
#include <cstring>
#include <new>

#include "opengl.h"


/// destination of error messages (set to zero to suppress output)
static FILE * ERF = stderr;


static uint8_t* new_pixels(size_t s, size_t col)
{
    void * ptr = nullptr;
    if ( posix_memalign(&ptr, 32, s*col) )
        throw std::bad_alloc();
    return (uint8_t*)ptr;
}

static void free_pixels(uint8_t *& ptr)
{
    free(ptr);
    ptr = nullptr;
}

//------------------------------------------------------------------------------
#pragma mark - Interface

/**
 saveImage(...) will read pixels from the current OpenGL read buffer,
 and save them in a file with the requested format
 */
int SaveImage::saveEntireImage(const char * filename,
                               const char format[],
                               int downsample)
{
    GLint vp[4];
    glGetIntegerv(GL_VIEWPORT, vp);
    //printf("saveImage viewport %i %i %i %i\n", vp[0], vp[1], vp[2], vp[3]);
    
    return saveImage(filename, format, vp, downsample);
}

/**
 saveImage(...) will read pixels from the current OpenGL read buffer,
 and save them in a file with the requested format
 */
int SaveImage::saveImage(const char * filename,
                         const char format[],
                         const int vp[4],
                         int downsample)
{
    int res = FAILED_READ;

    //allocate memory to hold image:
    uint8_t* pix = new_pixels(vp[2]*vp[3], 3);

    if ( 0 == readPixels(vp[0], vp[1], vp[2], vp[3], pix) )
        res = savePixels(filename, format, pix, vp[2], vp[3], downsample);
    
    free_pixels(pix);
    return res;
}

//------------------------------------------------------------------------------
#pragma mark - Export

/**
 After setting a higher resolution, this will translate the ModelView to produce several
 images that will be stiched together in memory, into an image with higher resolution.
 This works even if the image is larger than the maximum OpenGL viewPort,
 but there can be artifacts caused by objects in the stitched zones.
 */
int SaveImage::saveCompositeImage(const int mag,
                                  const char * filename,
                                  const char format[],
                                  const uint32_t width, const uint32_t height,
                                  const double pixel_size,
                                  void (*drawFunc)(int, void *), void * arg,
                                  int downsample)
{
    if ( ! supported(format) )
        return UNKNOWN_FORMAT;
    int res = OPENGL_ERROR;
    int mW = mag * width;
    int mH = mag * height;
    
    const int PIX = 3;  //number of bytes for each pixel
    uint8_t* pix = new_pixels(mW*mH, PIX);
    uint8_t* sub = new_pixels(width*height, PIX);
    
    const double cc = ( mag - 1 ) * 0.5;
    const double dx = width * pixel_size;
    const double dy = height * pixel_size;
    
    float mat[16];
    gym::get_view(mat);
    for ( int iy = 0; iy < mag; ++iy )
    for ( int ix = 0; ix < mag; ++ix )
    {
        gym::transScale((cc-ix)*dx, (cc-iy)*dy, 0, mag);
        glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
        drawFunc(mag, arg);
        if ( 0 == readPixels(0, 0, width, height, sub) )
        {
            uint8_t * dst = &pix[width*PIX*(ix+mH*iy)];
            for ( uint32_t u = 0; u < height; ++u )
                memcpy(&dst[u*mW*PIX], &sub[u*width*PIX], width*PIX);
        }
    }
    gym::set_view(mat);
    res = savePixels(filename, format, pix, mW, mH, downsample);
    free_pixels(pix);
    free_pixels(sub);
    return res;
}


/**
 This adjusts the Viewport to produce an image with higher resolution.
 The result should be better than saveCompositeImage, but uses more
 memory on the graphic card.
 */
int SaveImage::saveMagnifiedImage(const int mag,
                                  const char * filename,
                                  const char format[],
                                  const uint32_t width, const uint32_t height,
                                  void (*drawFunc)(int, void *), void * arg,
                                  int downsample)
{
    if ( ! supported(format) )
        return UNKNOWN_FORMAT;
    
    int mW = mag * width;
    int mH = mag * height;
    
    GLint dim[2] = { 0 };
    glGetIntegerv(GL_MAX_VIEWPORT_DIMS, dim);
    if ( mW > dim[0] || mH > dim[1] )
    {
        fprintf(ERF, "SaveImage:: exceeding maximum supported size (%ix%i)\n", (int)dim[0], (int)dim[1]);
        return FAILED_ALLOCATION;
    }
    
    const int PIX = 3;  //number of bytes for each pixel
    //allocate memory to hold the full image:
    uint8_t* pix = new_pixels(mW*mH, PIX);
    uint8_t* sub = new_pixels(width*height, PIX);

    GLint svp[4];
    glGetIntegerv(GL_VIEWPORT, svp);
    for ( int iy = 0; iy < mag; ++iy )
    for ( int ix = 0; ix < mag; ++ix )
    {
        glViewport(-ix*width, -iy*height, mW, mH);
        glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
        drawFunc(mag, arg);
        if ( 0 == readPixels(0, 0, width, height, sub) )
        {
            uint8_t * dst = &pix[width*PIX*(ix+mH*iy)];
            for ( uint32_t h = 0; h < height; ++h )
                memcpy(&dst[h*mW*PIX], &sub[h*width*PIX], width*PIX);
        }
    }
    int res = savePixels(filename, format, pix, mW, mH, downsample);
    free_pixels(pix);
    //restore original viewport:
    glViewport(svp[0], svp[1], svp[2], svp[3]);
    free_pixels(sub);
    return res;
}

//------------------------------------------------------------------------------
#pragma mark - Pixels

int SaveImage::readPixels(int32_t X, int32_t Y, uint32_t W, uint32_t H, GLvoid *pixels)
{
    // set the alignment to double-words
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    
#if ( 0 )
    GLint readbuf = 0, drawbuf = 0;
    glGetIntegerv(GL_READ_BUFFER, &readbuf);
    glGetIntegerv(GL_DRAW_BUFFER, &drawbuf);
    printf("framebuffers: read %i draw %i\n", readbuf, drawbuf);
#endif
    
    //read the pixel values, from top-left corner:
    glReadPixels(X, Y, W, H, GL_RGB, GL_UNSIGNED_BYTE, pixels);
    //printf(" read OpenGL pixels %ux%u\n", W, H);

    GLenum err = glGetError();
    
    if ( err != GL_NO_ERROR )
    {
        fprintf(ERF, "OpenGL error: could not read pixels (error %u)\n", err);
        return OPENGL_ERROR;
    }
    return NO_ERROR;
}


/// used to convert depth buffer image from float to uint16
static void convertPixels(uint32_t S, void * buf)
{
    float * src = reinterpret_cast<float*>(buf);
    uint16_t * dst = reinterpret_cast<uint16_t*>(buf);
    for ( size_t i = 0; i < S; ++i )
        dst[i] = uint16_t( src[i] * 65535 );
}

//------------------------------------------------------------------------------
#pragma mark - Depth buffer export in PNG format


int SaveImage::readDepthPixels(int32_t X, int32_t Y, uint32_t W, uint32_t H, GLvoid *pixels)
{
    // set the alignment to double-words
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    
    //read the pixel values, from top-left corner:
    glReadPixels(X, Y, W, H, GL_DEPTH_COMPONENT, GL_FLOAT, pixels);
    convertPixels(W*H, pixels);

    GLenum err = glGetError();
    
    if ( err != GL_NO_ERROR )
    {
        fprintf(ERF, "OpenGL error: could not read pixels (error %u)\n", err);
        return OPENGL_ERROR;
    }
    return NO_ERROR;
}


int SaveImage::saveDepthBuffer(const char * filename, const int vp[4])
{
    int res = FAILED_READ;
    uint8_t* pix = new_pixels(vp[2]*vp[3], 4);

    if ( 0 == readDepthPixels(vp[0], vp[1], vp[2], vp[3], pix) )
    {
        FILE * file = SaveImage::openFile(filename);
        if ( file )
        {
            res = saveGrayPNG(file, (uint16_t*)pix, vp[2], vp[3]);
            fclose(file);
            if ( res )
            {
                fprintf(ERF, " error %i while saving %s\n", res, filename);
                remove(filename);
            }
            return res;
        }
        return FILE_ERROR;
    }
    free_pixels(pix);
    return res;
}
