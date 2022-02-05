/* 
  Save Window's Bitmap version 3 file
  https://www.fileformat.info/format/bmp/egff.htm
  
  The bitmap format BMP allows to write files with 1 or 2 bits per pixels,
  and include a color palette to translate the pixel values to RGB colors.
  This type of bitmap has a minimal file size.
  
  based on program bmpsuite.c by Jason Summers
  http://entropymine.com/jason/bmpsuite/
  
  based on code by Adam Majewski (fraktal.republika.pl)

  Francois Nedelec, Cambridge University, 27.08.2021
*/        


#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <cstring>


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

static void write_uint16(FILE * fp, uint16_t x)
{
    putc(x, fp);
    putc(x >> 8, fp);
}

static void write_uint32(FILE * fp, uint32_t x)
{
    putc(x, fp);
    putc(x >> 8, fp);
    putc(x >> 16, fp);
    putc(x >> 24, fp);
}

static void write_sint32(FILE * fp, int32_t x)
{
    putc(x, fp);
    putc(x >> 8, fp);
    putc(x >> 16, fp);
    putc(x >> 24, fp);
}

static void write_color(FILE * fp, uint8_t R, uint8_t G, uint8_t B)
{
    putc(R, fp);
    putc(G, fp);
    putc(B, fp);
    putc(255, fp);
}

static void write_color(FILE * fp, uint8_t X)
{
    write_color(fp, X, X, X);
}

static size_t bytes_per_row(unsigned w, uint16_t BitsPerPixel)
{
    return 4 * (((w * BitsPerPixel)+31)/32);
}

void save_bitmap(FILE* fp, uint8_t bytes[], unsigned width, unsigned height, uint16_t BitsPerPixel)
{
    /* in bytes */
    uint32_t FileHeaderSize = 14;
    uint32_t InfoHeaderSize = 40;
    
    uint32_t NumberColors = 1 << BitsPerPixel;
    uint32_t PaletteSize = 4 * NumberColors;
    
    //------------
    uint32_t BytesPerRow = bytes_per_row(width, BitsPerPixel);
    uint32_t BytesSize = BytesPerRow * height;
    uint32_t OffBits = FileHeaderSize + InfoHeaderSize + PaletteSize;
    uint32_t FileSize = OffBits + BytesSize;
#if 0
    printf("BytesPerRow = %d\n", BytesPerRow);
    printf("BytesSize = %d\n", BytesSize);
    printf("FileSize = %d\n", FileSize);
    printf("OffBits = %d\n", OffBits);
#endif
    //------------  BMP file header
    write_uint16(fp, 0x4D42);   /* bfType = "MB" */
    write_uint32(fp, FileSize); /* bfSize */
    write_uint16(fp, 0);        /* bfReserved1 */
    write_uint16(fp, 0);        /* bfReserved2 */
    write_uint32(fp, OffBits);  /* bfOffBits */

    //------------  BMP file info structure
    write_uint32(fp, InfoHeaderSize); /* biSize */
    write_sint32(fp, width);          /* biWidth */
    write_sint32(fp, height);         /* biHeight */
    write_uint16(fp, 1);              /* biPlanes */
    write_uint16(fp, BitsPerPixel);   /* biBitCount */
    write_uint32(fp, 0);              /* biCompression */
    write_uint32(fp, 0);              /* biSizeImage */
    write_sint32(fp, 0);              /* biXPelsPerMeter */
    write_sint32(fp, 0);              /* biYPelsPerMeter */
    write_uint32(fp, NumberColors);   /* biClrUsed */
    write_uint32(fp,0);               /* biClrImportant */
    
    //------------  BMP colormap
    uint8_t colors[256] = { 0, 255 };
    /*  color table (palette) = 2 colors as a RGBA */
    for ( uint32_t c = 0; c < NumberColors; ++c )
        write_color(fp, colors[c]);
    
    fwrite(bytes,1,BytesSize,fp);
}


/**
 set bytes[] at ( colum = x, line = y ) to specified color
 argument `x` should be specified in bits: `BitsPerPixel` * pixel_X_coordinate
 argument `y` should be speicifed in bytes: `BytesPerRow` * pixel_Y_coordinate
*/
static inline void set_bitmap(uint8_t bytes[], uint32_t x, uint32_t y, uint8_t color)
{
    size_t i = y + x / 8;
    bytes[i] |= color << ( 7 - x & 7 );
}


#ifdef __cplusplus


class BitMap
{
    unsigned width;
    unsigned height;
    uint16_t depth;
    uint32_t BPR;
    uint8_t* bytes;

public:
    
    BitMap(unsigned W, unsigned H, uint16_t B)
    {
        width = W;
        height = H;
        depth = B;
        if ( B != 1 && B != 4 && B != 8 && B != 24 )
        {
            printf("invalid BitMap depth\n");
            B = 8;
        }
        BPR = bytes_per_row(width, depth);
        bytes = (uint8_t*)calloc(height*BPR, 1);
    }
    
    ~BitMap()
    {
        free(bytes);
        bytes = nullptr;
    }
    
    void clear()
    {
        std::memset(bytes, 0, height*BPR);
    }
    
    // matrix convention: i = line, with 0 at the top of the image
    void set(uint32_t i, uint32_t j, uint8_t color = 1)
    {
        //assert_true(( i < height ) & ( j < width ));
        set_bitmap(bytes, depth*j, BPR*(height-1-i), color);
    }
    
    // y = line, x = horizontal, x is 0 at the bottom of the image
    void flipset(uint32_t x, uint32_t y, uint8_t color)
    {
        //assert_true(( y < height ) & ( x < width ));
        set_bitmap(bytes, depth*x, BPR*y, color);
    }

    void save(FILE * fp)
    {
        save_bitmap(fp, bytes, width, height, depth);
    }
};

#endif

#if 0
int main()
{
    uint32_t W=1024, H=512;
    size_t imsize = bitmap_size(1, W, H);
    uint8_t * bytes = (uint8_t*)malloc(imsize);
    
    //metset(bytes, 0, imsize);
    for( size_t i = 0; i < imsize; ++i )
        bytes[i] = 0;

    uint32_t BPR = bytes_per_row(1, W);

    /* paint a disc */
    size_t radius_square = W * W / 4;
    for ( size_t j = 0; j < W; ++j )
    for ( size_t i = 0; i < H; ++i )
        if ( i * i + j * j < radius_square )
            set_bitmap(bytes, BPR, i, j, 1);
    
    FILE *fp = fopen("circle.bmp", "wb");
    save_bitmap(fp, bytes, W, H, 1);
    fclose(fp);
    free(bytes);
}
#endif
