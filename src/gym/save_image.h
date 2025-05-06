// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University
#ifndef SAVE_IMAGE_H
#define SAVE_IMAGE_H

#include <cstdio>
#include <cstdint>

/// Can save pixel array to files in PNG, TARGA or PPM format
/**
 - PPM files do not require any library, and thus writing is supported on any platform.
   They can be read by various software, in particular the 'pbmplus' toolkit
   However, they are big and usually not supported by most viewers.
   Known viewers on Mac OS X: ImageJ, GraphicConverter, ToyViewer.
 - TARGA is a simple color image format popular on Windows
 - PNG is a modern all-purpose file format that supports RGBA format, and is thus
   very well suited to export OpenGL scenes.
 .
 
 PNG support requires 'libpng', which must be installed separately.
 */
namespace SaveImage
{
    /// error codes
    enum { NO_ERROR=0, FAILED_ALLOCATION=1, FAILED_READ=2, UNKNOWN_FORMAT=3, OPENGL_ERROR=4, FILE_ERROR=5, PNG_ERROR=10 };
    
    /// open a file for binary write (used internally)
    FILE * openFile(const char name[]);
    
    /// true if 'format' is the 3-letter file-entension of a supported image format
    /**
     'ppm' is always supported, and 'png' may be supported depending on compilation
     The 3-letter format string can be lowercase or uppercase.
     */
    bool supported(const char format[]);

    
    /// Netpbm pixel image format, 3 one-byte componentss per pixels (R, G, B)
    int saveColorPPM(FILE*, const uint8_t[], uint32_t width, uint32_t height);
    
    /// uncompressed RGB TGA format (https://en.wikipedia.org/wiki/Truevision_TGA)
    int saveTGA(FILE*, const uint8_t[], bool color, uint32_t width, uint32_t height);
    
    /// save RGB Truevision TGA format
    int saveColorTGA(FILE*, const uint8_t[], uint32_t width, uint32_t height);
    
    /// save Grayscale Truevision TGA format
    int saveGrayTGA(FILE*, const uint8_t[], uint32_t width, uint32_t height);
    
    
    /// write PNG image
    int savePNG(FILE*, const uint8_t[], uint8_t bit_depth, uint8_t num_colors, uint32_t width, uint32_t height);
    
    /// Portable Network Graphic format, 4 one-byte components per pixels (R, G, B, A)
    int saveAlphaPNG(FILE*, const uint8_t[], uint32_t width, uint32_t height);
    
    /// Portable Network Graphic format, 3 one-byte components per pixels (R, G, B)
    int saveColorPNG(FILE*, const uint8_t[], uint32_t width, uint32_t height);
    
    /// save 16-bit grayscale PNG file
    int saveGrayPNG(FILE*, const uint16_t[], uint32_t width, uint32_t height);
    
    
    /// save pixels[] and return error-code
    int savePixels(FILE*, const char format[], const uint8_t[], uint32_t width, uint32_t height);
    
    /// save pixels[] and return error-code
    int savePixels(FILE*, const char format[], uint8_t[], uint32_t width, uint32_t height, int downsample);
    
    /// save pixels[] and return error-code
    int savePixels(const char* name, const char format[], uint8_t[], uint32_t width, uint32_t height, int downsample);
    
}


#endif
