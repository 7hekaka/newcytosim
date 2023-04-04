/*
 * fg_font.cc
 *
 * Bitmap and stroke fonts displaying.
 *
 * Copyright (c) 1999-2000 Pawel W. Olszta. All Rights Reserved.
 * Written by Pawel W. Olszta, <olszta@sourceforge.net>
 * Creation date: Thu Dec 16 1999
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
 * PAWEL W. OLSZTA BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
 * IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/*
 Modified March/April 2022 by FJ. Nedelec, to remove dependencies on FreeGLUT
 Keeping only the font-rendering capabilities
 */

#include "fg_font.h"
#include <stdlib.h>
#include <ctype.h>
#include <stdio.h>
#include <string.h>

typedef unsigned char uByte;


/* The bitmap font structure */
struct tagSFG_Font
{
    unsigned      Quantity;     /* Number of chars in font          */
    unsigned      Height;       /* Height of the characters         */
    const uByte** Characters;   /* The characters mapping           */
    float         xorig, yorig; /* Relative origin of the character */
};
typedef struct tagSFG_Font SFG_Font;


#include "fg_font_data.cc"


/// Fonts inherited from GLUT
enum GLUTFontType
{
    BITMAP_9_BY_15 = 2,
    BITMAP_8_BY_13 = 3,
    BITMAP_TIMES_ROMAN_10 = 4,
    BITMAP_TIMES_ROMAN_24 = 5,
    BITMAP_HELVETICA_10 = 6,
    BITMAP_HELVETICA_12 = 7,
    BITMAP_HELVETICA_18 = 8
};


int fgFontHeight(int font)
{
    if ( font == BITMAP_8_BY_13 )        return 15;
    if ( font == BITMAP_9_BY_15 )        return 16;
    if ( font == BITMAP_TIMES_ROMAN_10 ) return 12;
    if ( font == BITMAP_TIMES_ROMAN_24 ) return 27;
    if ( font == BITMAP_HELVETICA_10 )   return 13;
    if ( font == BITMAP_HELVETICA_12 )   return 15;
    if ( font == BITMAP_HELVETICA_18 )   return 22;
    return 13;
}


/*
 * Matches a font ID with a SFG_Font structure pointer.
 * This was changed to match the GLUT header style.
 */
static SFG_Font const* fghFont( int font )
{
    if( font == BITMAP_8_BY_13        ) return &fgFontFixed8x13;
    if( font == BITMAP_9_BY_15        ) return &fgFontFixed9x15;
    if( font == BITMAP_HELVETICA_10   ) return &fgFontHelvetica10;
    if( font == BITMAP_HELVETICA_12   ) return &fgFontHelvetica12;
    if( font == BITMAP_HELVETICA_18   ) return &fgFontHelvetica18;
    if( font == BITMAP_TIMES_ROMAN_10 ) return &fgFontTimesRoman10;
    if( font == BITMAP_TIMES_ROMAN_24 ) return &fgFontTimesRoman24;
    return NULL;
}

/* -- INTERFACE FUNCTIONS -------------------------------------------------- */

void drawBitmap(unsigned W, unsigned H, float X, float Y, float S, const unsigned char* data, const float color[4]);

/// convert binary image into one byte per pixel
void unpackBitmap(unsigned char data[], unsigned W, unsigned H, const unsigned char bits[], unsigned lda);

void drawPixels(unsigned W, unsigned H, float X, float Y, float S, const unsigned char* data, const float color[4]);
void paintBitmap(unsigned W, unsigned H, float X, float Y, float S, const unsigned char* data);


/*
 * Draw a bitmap character
 */
void fgBitmapCharacter(float x, float y, float S, int fontID, const float color[4], int character)
{
    const uByte* face;
    SFG_Font const* font = fghFont( fontID );
    if ( font && character >= 1 && character < 256 )
    {
        /*
         * Find the character we want to draw (???)
         */
        face = font->Characters[ character ];
        drawBitmap(face[0], font->Height, S*x-font->xorig, S*y-font->yorig, S, face+1, color);
    }
}


void fgBitmapString(float X, float Y, float scale, int fontID, const float color[4], const char *string, float vshift)
{
    SFG_Font const* font = fghFont( fontID );
    if ( !font )
        return;

    float gray[4] = { 0.5, 0.5, 0.5, 1 };
    const float* col = color;
    char * str = strdup(string);
    char * token = NULL;
    
    if ( vshift == 0 )
        vshift = font->Height;

    X = scale * ( X - font->xorig );
    Y = scale * ( Y - font->yorig );

    unsigned char * pixels = NULL;
    while ((token = strsep(&str, "\n")) != NULL)
    {
        //printf("%s\n", token);
        if ( token[0] == '%' )
            col = gray;
        else
            col = color;
        //printf("color %3.1f %3.1f %3.1f %3.1f : %c\n", col[0], col[1], col[2], col[3], c);
        const unsigned H = font->Height;
        // calculate total string length in pixels:
        unsigned L = 7;
        for ( char * c = token; *c; ++c )
            L += font->Characters[*c][0];
        pixels = (unsigned char*)realloc(pixels, L*H);
        memset(pixels, 0, L*H);
        unsigned W = 0;
        for ( char * c = token; *c; ++c )
        {
            const uByte* face = font->Characters[*c];
            unpackBitmap(pixels+W, face[0], H, 1+face, L);
            //paintBitmap(face[0], H, X+scale*W, Y, scale, 1+face);
            W += face[0];
        }
        drawPixels(L, H, X, Y, scale, pixels, col);
        // move down one line.
        Y += scale * vshift;
    }
    free(pixels);
    free(str);
}

#ifdef GL_VERSION_2_1
void fgBitmapString0(float x0, float y, int fontID, const float color[4], const char *string, float vshift)
{
    glColor4fv(color);
    glRasterPos2i(0, 0);
    glBitmap(0, 0, 0, 0, x0, y, NULL);
    unsigned char c;
    float x = 0.0f;
    SFG_Font const* font = fghFont( fontID );
    if ( font && string && *string )
    {
        glPushClientAttrib( GL_CLIENT_PIXEL_STORE_BIT );
        glPixelStorei( GL_UNPACK_SWAP_BYTES,  GL_FALSE );
        glPixelStorei( GL_UNPACK_LSB_FIRST,   GL_FALSE );
        glPixelStorei( GL_UNPACK_ROW_LENGTH,  0        );
        glPixelStorei( GL_UNPACK_SKIP_ROWS,   0        );
        glPixelStorei( GL_UNPACK_SKIP_PIXELS, 0        );
        glPixelStorei( GL_UNPACK_ALIGNMENT,   1        );
        
        if ( vshift == 0 )
            vshift = font->Height;
        /*
         * Step through the string, drawing each character.
         * A newline will simply translate the next character's insertion
         * point back to the start of the line and down one line.
         */
        while( ( c = *string++) )
        {
            if( c == '\n' )
            {
                glBitmap ( 0, 0, 0, 0, -x, vshift, NULL );
                x = 0.0f;
            }
            else  /* Not an EOL, draw the bitmap character */
            {
                const uByte* face = font->Characters[ c ];
                
                glBitmap(
                         face[ 0 ], font->Height,     /* Bitmap's width and height    */
                         font->xorig, font->yorig,    /* The origin in the font glyph */
                         ( float )( face[ 0 ] ), 0.0, /* The raster advance; inc. x,y */
                         ( face + 1 )                 /* The packed bitmap data...    */
                         );
                
                x += ( float )( face[ 0 ] );
            }
        }
        glPopClientAttrib( );
    }
}
#endif

/*
 * Returns the width in pixels of a font's character
 */
int fgBitmapWidth( int fontID, unsigned char character )
{
    SFG_Font const* font = fghFont( fontID );
    if ( font && character > 0 && character <= 255 )
        return *( font->Characters[ character ] );
    return 0;
}


/**
 Compute the max width of all the lines in the given text, and the number of lines
 */
int fgTextWidth(int fontID, const char text[], int& lines)
{
    int res = 0;
    int w = 0;
    lines = 0;
    SFG_Font const* font = fghFont(fontID);
    if ( !font )
        return 0;
    for (const char* c = text; *c != '\0' ; ++c)
    {
        if ( *c == '\n' )
        {
            if ( w > res ) res = w;
            ++lines;
            w = 0;
        }
        else if ( isprint(*c) )
        {
            w += font->Characters[(unsigned)*c][0];
        }
    }
    if ( w > res ) res = w;
    if ( res > 0 && lines == 0 )
        lines = 1;
    return res;
}

/*
 * Return the width of a string drawn using a bitmap font
 */
int fgBitmapLength( int fontID, const unsigned char* string )
{
    unsigned char c;
    int length = 0, this_line_length = 0;
    SFG_Font const* font = fghFont( fontID );
    if ( font && string && *string )
    {
        while( ( c = *string++) )
        {
            if( c != '\n' )/* Not an EOL, increment length of line */
                this_line_length += *( font->Characters[ c ]);
            else  /* EOL; reset the length of this line */
            {
                if( length < this_line_length )
                    length = this_line_length;
                this_line_length = 0;
            }
        }
        if ( length < this_line_length )
            length = this_line_length;
    }
    return length;
}

/*
 * Returns the height of a bitmap font
 */
int fgBitmapHeight( int fontID )
{
    SFG_Font const* font = fghFont( fontID );
    if ( font )
        return font->Height;
    return 0;
}

/*** END OF FILE ***/
