// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University

#ifndef SAVE_IMAGE_GL_H
#define SAVE_IMAGE_GL_H

#include "save_image.h"

/// Functions to save images obtained from OpenGL
namespace SaveImage
{
    /// save entire viewport in a new file called 'name'. Returns error-code
    int readPixels(int32_t x, int32_t y, uint32_t width, uint32_t height, void * pixels);
    
    /// save entire viewport in a new file called 'name'. Returns error-code
    int readDepthPixels(int32_t x, int32_t y, uint32_t width, uint32_t height, void * pixels);

    /// save entire viewport in a new file called 'name'. Returns error-code
    int saveEntireImage(const char* name, const char format[], int downsample=1);

    /// save a region of the current buffer in a new file called 'name'. Returns error-code
    int saveImage(const char* name, const char format[], const int viewport[4], int downsample=1);
    
    /// save a region of the current depth buffer in a new PNG file called 'name'. Returns error-code
    int saveDepthBuffer(const char* name, const int viewport[4]);

     /// save an image with higher resolution (this is better than saveCompositeImage)
    int saveMagnifiedImage(int mag, const char* name, const char format[], uint32_t width, uint32_t height, void (*display)(int, void *), void* arg, int downsample);

    /// save an image with higher resolution
    int saveCompositeImage(int mag, const char* name, const char format[], uint32_t width, uint32_t height, double pixel_size, void (*display)(int, void *), void* arg, int downsample);
}

#endif
