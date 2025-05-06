/* 
  Save Window's Bitmap version 3 file
  https://www.fileformat.info/format/bmp/egff.htm
  
  The bitmap format BMP allows to write files with 1 or 2 bits per pixels,
  and include a color palette to translate the pixel values to RGB colors.
  This format minimizes the file size, allowing to save very large images.
  
  based on program bmpsuite.c by Jason Summers
  http://entropymine.com/jason/bmpsuite/
  
  based on code by Adam Majewski (fraktal.republika.pl)

  Francois Nedelec, Cambridge University, 27.08.2021
*/        


#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#if 0
/** BMP file header structure **/
typedef struct __attribute__ ((__packed__))
{
    uint16_t  FileType;     /* File type, always 4D42h ("BM") */
    uint32_t  FileSize;     /* Size of the file in bytes */
    uint16_t  Reserved1;    /* Always 0 */
    uint16_t  Reserved2;    /* Always 0 */
    uint32_t  BitmapOffset; /* Starting position of image data in bytes */
} BITMAPFILEHEADER;

#  define BF_TYPE 0x4D42             /* "MB" */

/** BMP file info structure **/
typedef struct __attribute__ ((__packed__))
{
    uint32_t Size;            /* Size of this header in bytes */
    int32_t  Width;           /* Image width in pixels */
    int32_t  Height;          /* Image height in pixels */
    uint16_t Planes;          /* Number of color planes */
    uint16_t BitsPerPixel;    /* Number of bits per pixel */
    uint32_t Compression;     /* Compression methods used */
    uint32_t SizeOfBitmap;    /* Size of bitmap in bytes */
    int32_t  HorzResolution;  /* Horizontal resolution in pixels per meter */
    int32_t  VertResolution;  /* Vertical resolution in pixels per meter */
    uint32_t ColorsUsed;      /* Number of colors in the image */
    uint32_t ColorsImportant; /* Minimum number of important colors */
} BITMAPINFOHEADER;
#endif

static void write_uint16(FILE * f, uint16_t x)
{
    putc(x, f);
    putc(x >> 8, f);
}

static void write_uint32(FILE * f, uint32_t x)
{
    putc(x, f);
    putc(x >> 8, f);
    putc(x >> 16, f);
    putc(x >> 24, f);
}

static void write_sint32(FILE * f, int32_t x)
{
    putc(x, f);
    putc(x >> 8, f);
    putc(x >> 16, f);
    putc(x >> 24, f);
}

static void write_color(FILE * f, uint8_t R, uint8_t G, uint8_t B)
{
    putc(R, f);
    putc(G, f);
    putc(B, f);
    putc(0, f);
}

static uint32_t bytes_per_row(uint32_t w, uint16_t BitsPerPixel)
{
    return 4 * (((w * BitsPerPixel)+31)/32);
}

void save_bitmap(FILE* f, uint8_t bytes[], uint32_t width, uint32_t height, uint16_t BitsPerPixel)
{
    /* in bytes */
    uint32_t FileHeaderSize = 14;
    uint32_t InfoHeaderSize = 40;
    
    uint32_t NumColors = 1 << BitsPerPixel;
    uint32_t PaletteSize = 4 * NumColors;
    uint32_t BytesPerRow = bytes_per_row(width, BitsPerPixel);
    uint32_t BytesSize = BytesPerRow * height;
    uint32_t MapOffset = FileHeaderSize + InfoHeaderSize + PaletteSize;
    uint32_t FileSize = MapOffset + BytesSize;
#if 0
    printf("BytesPerRow = %d\n", BytesPerRow);
    printf("BytesSize = %d\n", BytesSize);
    printf("FileSize = %d\n", FileSize);
    printf("OffBits = %d\n", MapOffset);
#endif
    //------------  BMP file header
    write_uint16(f, 0x4D42);   /* File type, always 4D42h ("BM") */
    write_uint32(f, FileSize); /* Size of the file in bytes */
    write_uint16(f, 0);        /* Always 0 */
    write_uint16(f, 0);        /* Always 0 */
    write_uint32(f, MapOffset);  /* Starting position of image data in bytes */

    //------------  BMP file info structure
    write_uint32(f, InfoHeaderSize); /* Size of this header in bytes */
    write_sint32(f, width);          /* Image width in pixels */
    write_sint32(f, height);         /* Image height in pixels */
    write_uint16(f, 1);              /* Number of color planes */
    write_uint16(f, BitsPerPixel);   /* Number of bits per pixel */
    write_uint32(f, 0);              /* Compression methods used */
    write_uint32(f, 0);              /* Size of bitmap in bytes */
    write_sint32(f, 0);              /* Horizontal resolution in pixels per meter */
    write_sint32(f, 0);              /* Vertical resolution in pixels per meter */
    write_uint32(f, NumColors);      /* Number of colors in the image */
    write_uint32(f, 1);              /* Minimum number of important colors */
    
    //------------  BMP color palette
    
    /*  color table (palette) = RGBX with 1 bytes per component */
    const uint8_t top = ~uint8_t(0) / ( NumColors - 1 );
    for ( uint32_t c = 0; c < NumColors; ++c )
    {
        uint8_t X = c * top;
        write_color(f, X, X, X);
        //printf("color %i : %3i %3i %3i\n", c, X, X, X);
    }
    fwrite(bytes, 1, BytesSize, f);
}


#ifdef __cplusplus


template<uint16_t depth>
class BitMap
{
    uint32_t width;
    uint32_t height;
    uint32_t BPR;
    uint8_t* bytes;

    /**
     set bytes[] at ( colum = x, line = y ) to specified color, assuming corresponding
     bits are zero:
      `x` should be specified in bits: `BitsPerPixel` * pixel_X_coordinate
      `y` should be specified in bytes: `BytesPerRow` * pixel_Y_coordinate
    */
    void set_byte(uint32_t x, uint32_t y, uint8_t color)
    {
        uint32_t i = y + x / 8;
        const uint8_t mask = (uint8_t)(( 1 << depth ) - 1);
        const uint8_t top = depth * ( 8 / depth - 1 );
        int s = top - ( x & 7 );
        bytes[i] = (bytes[i] & ~( mask << s )) | (( color & mask ) << s );
    }

public:
    
    static bool valid_depth(uint16_t d)
    {
        return ( d == 1 || d == 4 || d == 8 || d == 24 );
    }
    
    void resize(uint32_t W, uint32_t H)
    {
        width = W;
        height = H;
        if ( !valid_depth(depth) )
            printf("invalid BitMap depth\n");
        BPR = bytes_per_row(width, depth);
        bytes = (uint8_t*)calloc(height*BPR, 1);
    }
    
    // the matrix is allocated uninitialized
    BitMap(uint32_t W, uint32_t H)
    {
        resize(W, H);
    }
    
    ~BitMap()
    {
        free(bytes);
        bytes = nullptr;
    }
    
    void clear()
    {
        memset(bytes, 0, height*BPR);
    }
    
    // matrix convention: i = line, with 0 at the top of the image
    void set(uint32_t i, uint32_t j, uint8_t color = 1)
    {
        //assert_true(( i < height ) & ( j < width ));
        set_byte(uint32_t(depth)*j, BPR*(height-1-i), color);
    }
    
    // matrix convention: i = line, with 0 at the top of the image
    void set_if(uint32_t i, uint32_t j, uint8_t color = 1)
    {
        if (( i < height ) & ( j < width ))
            set_byte(uint32_t(depth)*j, BPR*(height-1-i), color);
    }

    // x = horizontal, y = vertical, with 0 at the bottom of the image
    void flipset(uint32_t x, uint32_t y, uint8_t color)
    {
        //assert_true(( y < height ) & ( x < width ));
        set_byte(uint32_t(depth)*x, BPR*y, color);
    }

    void save(FILE * f)
    {
        save_bitmap(f, bytes, width, height, depth);
    }
};

#endif

#if 0
int main()
{
    const uint32_t depth = 4;
    uint32_t W=1024, H=512;
    BitMap<depth> bmap(W, H);
    bmap.clear();

    /* paint a disc */
    size_t sup = ( 1 << depth ) - 1;
    size_t radius_square = W * W / 4;
    for ( size_t j = 0; j < W; ++j )
    for ( size_t i = 0; i < H; ++i )
    {
        ssize_t y = j - 400;
        size_t rr = ( i * i + y * y ) * sup / radius_square;
        bmap.set(i, j, rr);
    }
    
    FILE *f = fopen("circle.bmp", "wb");
    bmap.save(f);
    fclose(f);
}
#endif
