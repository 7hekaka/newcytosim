// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University
#include "save_image.h"
#include <cstdlib>
#include <cstring>
#include <new>

#include "gym_image.h"

/// destination of error messages (set to zero to suppress output)
static FILE * ERF = stderr;


bool SaveImage::supported(const char format[])
{
    if ( 0 == strcasecmp(format, "png") )
        return true;
    if ( 0 == strcasecmp(format, "ppm") )
        return true;
    if ( 0 == strcasecmp(format, "tga") )
        return true;

    return false;
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

//------------------------------------------------------------------------------
#pragma mark - Interface functions

int SaveImage::savePixels(FILE * file,
                          const char format[],
                          const uint8_t pixels[],
                          uint32_t width, uint32_t height)
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
                          uint8_t pixels[],
                          uint32_t width, uint32_t height,
                          int downsample)
{
    if ( downsample > 1 )
    {
        //printf("downsampling %i to : %i %i\n", downsample, mw, mh);
        int W = width / downsample;
        int H = height / downsample;
        gym::downsampleRGB(pixels, W, H, pixels, downsample);
        return savePixels(file, format, pixels, W, H);
    }
    
    return savePixels(file, format, pixels, width, height);
}


int SaveImage::savePixels(const char * filename,
                          const char format[],
                          uint8_t pixels[],
                          uint32_t width, uint32_t height,
                          int downsample)
{
    FILE * file = openFile(filename);
    
    if ( file )
    {
        int res = savePixels(file, format, pixels, width, height, downsample);
        
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

//------------------------------------------------------------------------------
#pragma mark - PPM format

/**
 Write the image in the Portable Pixmap format, also Netpbm format (man -k ppm).
 We use here the 'raw' binary format starting with P6 
 */
int SaveImage::saveColorPPM(FILE* file,
                            const uint8_t pixels[],
                            const uint32_t width, const uint32_t height)
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
 save RGB Truevision TGA format (also known as TARGA)
 https://en.wikipedia.org/wiki/Truevision_TGA
 */
int SaveImage::saveTGA(FILE* file,
                       const uint8_t pixels[], const bool color,
                       const uint32_t width, const uint32_t height)
{
    uint8_t header[18] = { 0 };
    // Data code type -- 2 - uncompressed RGB image.
    header[2] = (color?2:3);
    // Image width - low byte
    header[12] = width & 0xFF;
    // Image width - high byte
    header[13] = (width >> 8) & 0xFF;
    // Image height - low byte
    header[14] = height & 0xFF;
    // Image height - high byte
    header[15] = (height >> 8) & 0xFF;
    // Color bit depth
    header[16] = (color?24:8);

    fwrite(header, 1, 18, file);
    fwrite(pixels, 1, width * height * (color?3:1), file);
    return NO_ERROR;
}


/// save RGB Truevision TGA format
int SaveImage::saveColorTGA(FILE* file, const uint8_t pixels[], uint32_t width, uint32_t height)
{
    return saveTGA(file, pixels, 1, width, height);
}

/// save Grayscale Truevision TGA format
int SaveImage::saveGrayTGA(FILE* file, const uint8_t pixels[], uint32_t width, uint32_t height)
{
    return saveTGA(file, pixels, 0, width, height);
}

//------------------------------------------------------------------------------
#pragma mark - PNG format using libpng

#ifdef HAS_PNG

#include <png.h>

static int savePNG(FILE* file,
                   png_bytep row_pointers[],
                   const uint8_t bit_depth, const uint8_t num_colors,
                   const uint32_t width, const uint32_t height)
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
                       const uint8_t bit_depth, const uint8_t num_colors,
                       const uint32_t width, const uint32_t height)
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

#elif 1

//------------------------------------------------------------------------------
#pragma mark - PNG export using libspng (https://libspng.org)

#include "spng.h"

int SaveImage::savePNG(FILE* file, const uint8_t pixels[],
                       const uint8_t bit_depth, const uint8_t num_colors,
                       const uint32_t width, const uint32_t height)
{
    int res = 0;
    const size_t length = num_colors * width * height * ( bit_depth / 8 );

    uint8_t fmt = SPNG_COLOR_TYPE_TRUECOLOR;
    if ( num_colors == 4 )
        fmt = SPNG_COLOR_TYPE_TRUECOLOR_ALPHA;
    if ( num_colors == 1 )
        fmt = SPNG_COLOR_TYPE_GRAYSCALE;
    
    /* Specify image dimensions, PNG format */
    struct spng_ihdr ihdr =
    {
        .width = width,
        .height = height,
        .bit_depth = bit_depth,
        .color_type = fmt
    };

    /* Creating an encoder context requires a flag */
    spng_ctx *enc = spng_ctx_new(SPNG_CTX_ENCODER);
    spng_set_option(enc, SPNG_IMG_COMPRESSION_LEVEL, 9);

    /* Encode to internal buffer managed by the library */
    //spng_set_option(enc, SPNG_ENCODE_TO_BUFFER, 1);
    
    /* Encode to file directly */
    res = spng_set_png_file(enc, file);
    if (res)
    {
        fprintf(ERF, "spng_set_png_file() error: %s\n", spng_strerror(res));
        goto done;
    }

    /* Image will be encoded according to ihdr.color_type, .bit_depth */
    res = spng_set_ihdr(enc, &ihdr);
    if (res)
    {
        fprintf(ERF, "spng_set_ihdr() error: %s\n", spng_strerror(res));
        goto done;
    }

    /* SPNG_FMT_PNG is a special value that matches the format in ihdr,
       SPNG_ENCODE_FINALIZE will finalize the PNG with the end-of-file marker */
    res = spng_encode_image(enc, pixels, length, SPNG_FMT_PNG, SPNG_ENCODE_FINALIZE|SPNG_ENCODE_FLIP_Y);
    if (res)
        fprintf(ERF, "spng_encode_image() error: %s\n", spng_strerror(res));

done:
    /* Free context memory */
    spng_ctx_free(enc);
    return res;
}

#else

int SaveImage::savePNG(FILE* file, const uint8_t pixels[],
                       const uint8_t bit_depth, const uint8_t num_colors,
                       const uint32_t width, const uint32_t height)
{
    fprintf(ERF, "png export not supported\n");
    return 1;
}

#endif

/**
 Save RGBA image, 4 x 8-bits per pixel
 */
int SaveImage::saveAlphaPNG(FILE* file, const uint8_t pixels[],
                            const uint32_t width, const uint32_t height)
{
    return savePNG(file, pixels, 8, 4, width, height);
}

/**
 Save RGB image, 3 x 8-bits per pixel
 */
int SaveImage::saveColorPNG(FILE* file, const uint8_t pixels[],
                            const uint32_t width, const uint32_t height)
{
    return savePNG(file, pixels, 8, 3, width, height);
}

/**
 Save 16-bits gray-level image
 */
int SaveImage::saveGrayPNG(FILE* file, const uint16_t pixels[],
                           const uint32_t width, const uint32_t height)
{    
    return savePNG(file, (uint8_t*)pixels, 16, 1, width, height);
}

