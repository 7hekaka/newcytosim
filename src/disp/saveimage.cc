// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University
#include "saveimage.h"
#include <cstdlib>
#include <cstring>

#ifndef NO_OPENGL
#  include "opengl.h"
#endif

bool SaveImage::supported(const char format[])
{
#ifdef HAS_PNG
    if ( 0 == strcasecmp(format, "png") )
        return true;
#endif
    if ( 0 == strcasecmp(format, "ppm") )
        return true;
    if ( 0 == strcasecmp(format, "tga") )
        return true;

    return false;
}


static uint8_t* new_pixels(size_t s)
{
    return (uint8_t*)calloc(s, sizeof(uint8_t));
}

static void free_pixels(uint8_t *& ptr)
{
    free(ptr);
    ptr = nullptr;
}


//------------------------------------------------------------------------------
#pragma mark Images

#ifndef NO_OPENGL

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
    int res = FAILED_ALLOCATION;

    //allocate memory to hold image:
    uint8_t* pixels = new_pixels(3*vp[2]*vp[3]);

    if ( pixels )
    {
        if ( 0 == readPixels(vp[0], vp[1], vp[2], vp[3], pixels) )
            res = savePixels(filename, format, pixels, vp[2], vp[3], downsample);
        free_pixels(pixels);
    }
    else
        fprintf(err, "Error: SaveImage failed to allocate memory!\n");

    return res;
}


/**
 After setting a higher resolution, this will translate the ModelView to produce several
 images that will be stiched together in memory, into an image with higher resolution.
 This works even if the image is larger than the maximum OpenGL viewPort,
 but there can be artifacts caused by objects in the stitched zones.
 */
int SaveImage::saveCompositeImage(const int mag,
                                  const char * filename,
                                  const char format[],
                                  const int width, const int height,
                                  const double pixel_size,
                                  void (*display)(int, void *), void * arg,
                                  int downsample)
{
    if ( ! supported(format) )
        return UNKNOWN_FORMAT;
    
    int res = 1;
    int mW = mag * width;
    int mH = mag * height;
    
    const int PIX = 3;  //number of bytes for each pixel
    uint8_t* pixels = new_pixels(mW*mH*PIX);
    uint8_t* sub    = new_pixels(width*height*PIX);
    
    if ( pixels )
    {
        const double cc = ( mag - 1 ) * 0.5;
        const double dx = width * pixel_size / mag;
        const double dy = height * pixel_size / mag;
        
        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();
        glScalef(mag, mag, mag);
        
        for ( int iy = 0; iy < mag ; ++iy )
        {
            for ( int ix = 0; ix < mag; ++ix )
            {
                glPushMatrix();
                glTranslated((cc-ix)*dx, (cc-iy)*dy, 0);
                display(mag, arg);
                glPopMatrix();
                if ( 0 == readPixels(0, 0, width, height, sub) )
                {
                    uint8_t * dst = &pixels[width*PIX*(ix+mH*iy)];
                    for ( int ii=0; ii<height; ++ii )
                        memcpy(&dst[ii*mW*PIX], &sub[ii*width*PIX], width*PIX);
                }
            }
        }
        glPopMatrix();
        
        res = savePixels(filename, format, pixels, mW, mH, downsample);
        
        free_pixels(pixels);
    }
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
                                  const int width, const int height,
                                  void (*display)(int, void *), void * arg,
                                  int downsample)
{
    if ( ! supported(format) )
        return UNKNOWN_FORMAT;
    
    int res = 1;
    int mW = mag * width;
    int mH = mag * height;
    
    GLint dim[2] = { 0 };
    glGetIntegerv(GL_MAX_VIEWPORT_DIMS, dim);
    if ( mW > dim[0] || mH > dim[1] )
    {
        fprintf(err, "SaveImage:: exceeding maximum supported size (%ix%i)\n", (int)dim[0], (int)dim[1]);
        return FAILED_ALLOCATION;
    }
    
    const int PIX = 3;  //number of bytes for each pixel
    //allocate memory to hold the full image:
    uint8_t* pixels = new_pixels(mW*mH*PIX);
    uint8_t* sub = new_pixels(width*height*PIX);
    if ( pixels )
    {
        GLint svp[4];
        glGetIntegerv(GL_VIEWPORT, svp);
        for ( int iy = 0; iy < mag; ++iy )
        for ( int ix = 0; ix < mag; ++ix )
        {
            glViewport(-ix*width, -iy*height, mW, mH);
            display(mag, arg);
            if ( 0 == readPixels(0, 0, width, height, sub) )
            {
                uint8_t * dst = &pixels[width*PIX*(ix+mH*iy)];
                for ( int h = 0; h < height; ++h )
                    memcpy(&dst[h*mW*PIX], &sub[h*width*PIX], width*PIX);
            }
        }
        res = savePixels(filename, format, pixels, mW, mH, downsample);
        free_pixels(pixels);
        //restore original viewport:
        glViewport(svp[0], svp[1], svp[2], svp[3]);
    }
    free_pixels(sub);
    return res;
}

//------------------------------------------------------------------------------
#pragma mark - Pixels

int SaveImage::readPixels(GLint X, GLint Y, GLsizei W, GLsizei H, GLvoid *pixels)
{
    //set the alignment to double-words
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

    GLenum glError = glGetError();
    
    if ( glError != GL_NO_ERROR )
    {
        fprintf(err, "Error: could not read pixels (OpenGL error %u)\n", glError);
        return OPENGL_ERROR;
    }
    return 0;
}

#endif

/**
 This will downsample pixelmap `src` and set destination `dst`. The pixel
 array `src` should be of size `3*W*H` with 3 bytes per pixels: R, G, B,
 while `src` will be `bin*bin` times smaller. Pixels are stored in row order
 from the lowest to the highest row, left to right in each row (as in OpenGL).
 The pixels components of `src` are averaged to produce `dst`.
 */
void SaveImage::downsampleRGB(uint8_t dst[], const uint8_t src[],
                              unsigned W, unsigned H, unsigned bin)
{
    const size_t s = bin * bin;
    const size_t sx = W / bin;
    const size_t sy = H / bin;
    
#if ( 0 )
    //reset destination:
    for ( size_t u = 0; u < sx*sy; ++u )
    {
        dst[3*u  ] = 0xFF;
        dst[3*u+1] = 0xFF;
        dst[3*u+2] = 0xFF;
    }
#endif
    
    for ( size_t x = 0; x < sx; ++x )
    for ( size_t y = 0; y < sy; ++y )
    {
        size_t r = 0, g = 0, b = 0;
        for ( size_t dx = 0; dx < bin; ++dx )
        for ( size_t dy = 0; dy < bin; ++dy )
        {
            uint8_t const* p = src + 3 * ( dx + bin*x + W*(dy+bin*y) );
            r += p[0];
            g += p[1];
            b += p[2];
        }
        
        dst[3*(x+sx*y)  ] = (uint8_t)( r / s );
        dst[3*(x+sx*y)+1] = (uint8_t)( g / s );
        dst[3*(x+sx*y)+2] = (uint8_t)( b / s );
    }
}


int SaveImage::savePixels(FILE * file,
                          const char format[],
                          const uint8_t pixels[],
                          int width, int height)
{
    if ( 0 == strcasecmp(format, "ppm") )
        return saveColorPPM(file, pixels, width, height);
    
    if ( 0 == strcasecmp(format, "tga") )
        return saveColorTGA(file, pixels, width, height);

    if ( 0 == strcasecmp(format, "png") )
        return saveColorPNG(file, pixels, width, height);
    
    return UNKNOWN_FORMAT;
}


int SaveImage::savePixels(FILE * file,
                          const char format[],
                          const uint8_t pixels[],
                          int width, int height,
                          int downsample)
{
    if ( downsample > 1 )
    {
        //printf("downsampling %i to : %i %i\n", downsample, mw, mh);
        int W = width / downsample;
        int H = height / downsample;
        uint8_t* img = new_pixels(3*W*H);
        if ( img )
        {
            downsampleRGB(img, pixels, width, height, downsample);
            int res = savePixels(file, format, img, W, H);
            free_pixels(img);
            return res;
        }
    }
    
    return savePixels(file, format, pixels, width, height);
}


FILE * SaveImage::openFile(const char * filename)
{
    if ( filename[0] == '\0' )
        return nullptr;
    FILE * f = fopen(filename, "wb");
    if ( f )
    {
        if ( ferror(f) )
            fclose(f);
        else
            return f;
    }
    return nullptr;
}


int SaveImage::savePixels(const char * filename,
                          const char format[],
                          const uint8_t pixels[],
                          int width, int height,
                          int downsample)
{
    FILE * file = openFile(filename);
    
    if ( file )
    {
        int res = savePixels(file, format, pixels, width, height, downsample);
        
        fclose(file);
        
        if ( res )
        {
            fprintf(err, " error %i while saving %s\n", res, filename);
            // remove file if an error occurred:
            remove(filename);
        }
        
        return res;
    }
    
    return FILE_ERROR;
}

//------------------------------------------------------------------------------
#pragma mark - PPM format

/**
 Write the image in the Portable Pixmap format, also Netpbm format (man -k ppm).
 We use here the 'raw' binary format starting with P6 
 */
int SaveImage::saveColorPPM(FILE* file,
                            const uint8_t pixels[],
                            const int width, const int height)
{
    fprintf(file, "P6\n");
    fprintf(file, "%i %i\n", width, height);
    fprintf(file, "255\n");
    
    //write the pixels binary, line by line:
    for ( int i = height-1; i >= 0; --i )
        fwrite(&pixels[3*i*width], 1, 3*width, file);
    return NO_ERROR;
}


//------------------------------------------------------------------------------
#pragma mark - TGA format

/**
 save RGB Truevision TGA format
 https://en.wikipedia.org/wiki/Truevision_TGA
 */
int SaveImage::saveColorTGA(FILE* file,
                            const uint8_t pixels[],
                            const int width, const int height)
{
    uint8_t header[18] = { 0 };
    // Data code type -- 2 - uncompressed RGB image.
    header[2] = 2;
    // Image width - low byte
    header[12] = width & 0xFF;
    // Image width - high byte
    header[13] = (width >> 8) & 0xFF;
    // Image height - low byte
    header[14] = height & 0xFF;
    // Image height - high byte
    header[15] = (height >> 8) & 0xFF;
    // Color bit depth
    header[16] = 24;

    fwrite(header, sizeof(uint8_t), 18, file);
    fwrite(pixels, sizeof(uint8_t), width * height * 3, file);
    return NO_ERROR;
}


//------------------------------------------------------------------------------
#pragma mark - PNG format

#ifndef HAS_PNG

int SaveImage::savePNG(FILE*, const uint8_t pixels[],
                       const int, const int,
                       const int, const int)
{
    fprintf(err, "PNG format not supported (recompilation needed)\n");
    return UNKNOWN_FORMAT;
}

#else

#include <png.h>

/*
int SaveImage::readPNG(FILE* fp, png_bytep *& row_pointers, int& bit_depth, int& color_mode, int& width, int& height)
{
    if (!fp)
        return 9;
    
    char header[8];    // 8 is the maximum size that can be checked
    fread(header, 1, 8, fp);
    
    if (png_sig_cmp(header, 0, 8))
        return 8;
    
    // initialize stuff 
    png_structp png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    
    if (!png_ptr)
        return 7;
    
    png_infop info_ptr = png_create_info_struct(png_ptr);
    
    if (!info_ptr)
        return 6;
    
    if (setjmp(png_jmpbuf(png_ptr)))
        return 5;
    
    png_init_io(png_ptr, fp);
    png_set_sig_bytes(png_ptr, 8);
    
    png_read_info(png_ptr, info_ptr);
    
    width = png_get_image_width(png_ptr, info_ptr);
    height = png_get_image_height(png_ptr, info_ptr);
    
    color_type = png_get_color_type(png_ptr, info_ptr);
    bit_depth = png_get_bit_depth(png_ptr, info_ptr);
    
    int number_of_passes = png_set_interlace_handling(png_ptr);
    png_read_update_info(png_ptr, info_ptr);
    
    
    // read file
    if (setjmp(png_jmpbuf(png_ptr)))
        return 4;
    
    free(row_pointers);
    
    row_pointers = (png_bytep*) malloc(sizeof(png_bytep) * height);
    
    for (y=0; y<height; y++)
        row_pointers[y] = (png_byte*) malloc(png_get_rowbytes(png_ptr,info_ptr));
    
    png_read_image(png_ptr, row_pointers);
}
*/


static int savePNG(FILE* file,
                   png_bytep row_pointers[],
                   const int bit_depth, const int num_colors,
                   const int width, const int height)
{
    if ( !file )
        return SaveImage::FILE_ERROR;
    
    if ( bit_depth != 8 && bit_depth != 16 )
        return 19;
    
    int color_type = -1;
    
    if ( num_colors == 1 )
        color_type = PNG_COLOR_TYPE_GRAY;
    
    if ( num_colors == 3 )
        color_type = PNG_COLOR_TYPE_RGB;
    
    if ( num_colors == 4 )
        color_type = PNG_COLOR_TYPE_RGBA;
    
    if ( color_type < 0 )
        return 18;
    
    /* initialize stuff */
    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    
    if (!png_ptr)
        return 17;
    
    png_infop info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr)
        return 16;
    
    if (setjmp(png_jmpbuf(png_ptr)))
        return 15;
    
    png_init_io(png_ptr, file);
    
    /* write header */
    if (setjmp(png_jmpbuf(png_ptr)))
        return 14;
    
    png_set_IHDR(png_ptr, info_ptr, width, height,
                 bit_depth, color_type,
                 PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
    
    png_write_info(png_ptr, info_ptr);
    
    /* write bytes */
    if (setjmp(png_jmpbuf(png_ptr)))
        return 13;
    
    if ( num_colors == 1 && bit_depth == 16 )
        png_set_swap(png_ptr);
    
    png_write_image(png_ptr, row_pointers);
    
    
    /* end write */
    if (setjmp(png_jmpbuf(png_ptr)))
        return 12;
    
    png_write_end(png_ptr, NULL);
 
    return 0;
}


int SaveImage::savePNG(FILE* file, const uint8_t pixels[],
                       const int bit_depth, const int num_colors,
                       const int width, const int height)
{
    int res = FAILED_ALLOCATION;
    png_bytep * rows = (png_bytep*)malloc(height*sizeof(png_bytep));
    if ( rows )
    {
        int bytes_per_row = ( bit_depth / 8 ) * num_colors * width;
        
        png_byte * start = (png_byte*)pixels;
        
        for ( int y = 0; y < height; ++y )
            rows[y] = start + bytes_per_row * ( height-y-1 );
        
        res = savePNG(file, rows, bit_depth, num_colors, width, height);
        
        free(rows);
    }
    
    return res;
}

#endif


/**
 Save RGBA image, 4 x 8-bits per pixel
 */
int SaveImage::saveAlphaPNG(FILE* file, const uint8_t pixels[],
                            const int width, const int height)
{
    return savePNG(file, pixels, 8, 4, width, height);
}

/**
 Save RGB image, 3 x 8-bits per pixel
 */
int SaveImage::saveColorPNG(FILE* file, const uint8_t pixels[],
                            const int width, const int height)
{
    return savePNG(file, pixels, 8, 3, width, height);
}

/**
 Save 16-bits gray-level image
 */
int SaveImage::saveGrayPNG(FILE* file, const uint16_t pixels[],
                           const int width, const int height)
{    
    return savePNG(file, (uint8_t*)pixels, 16, 1, width, height);
}

