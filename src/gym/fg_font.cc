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
    switch ( font )
    {
        case BITMAP_8_BY_13:        return 13;
        case BITMAP_9_BY_15:        return 15;
        case BITMAP_TIMES_ROMAN_10: return 10;
        case BITMAP_TIMES_ROMAN_24: return 24;
        case BITMAP_HELVETICA_10:   return 10;
        case BITMAP_HELVETICA_12:   return 12;
        case BITMAP_HELVETICA_18:   return 18;
    }
    return 13;
}


/*
 * Matches a font ID with a SFG_Font structure pointer.
 * This was changed to match the GLUT header style.
 */
static SFG_Font const* fghFont( int font )
{
    switch ( font )
    {
        case BITMAP_8_BY_13:        return &fgFontFixed8x13;
        case BITMAP_9_BY_15:        return &fgFontFixed9x15;
        case BITMAP_TIMES_ROMAN_10: return &fgFontTimesRoman10;
        case BITMAP_TIMES_ROMAN_24: return &fgFontTimesRoman24;
        case BITMAP_HELVETICA_10:   return &fgFontHelvetica10;
        case BITMAP_HELVETICA_12:   return &fgFontHelvetica12;
        case BITMAP_HELVETICA_18:   return &fgFontHelvetica18;
    }
    return NULL;
}

/* -- INTERFACE FUNCTIONS -------------------------------------------------- */

#include "gym_draw.h"
#include "gym_flute.h"
#include "gym_image.h"
#include "gym_flat.h"



static unsigned countBits(size_t n_bytes, const unsigned char* bytes)
{
    unsigned n = 0;
    for ( unsigned b = 0; b < n_bytes; ++b )
        n += __builtin_popcount(bytes[b]);
    return n;
}


/*
 * Draw one bitmap character
 */
void fgBitmapCharacter(float X, float Y, float S, int fontID, const float col[4], int character)
{
    SFG_Font const* font = fghFont( fontID );
    if ( font && character >= 1 && character < 256 )
    {
        /*
         * Find the character we want to draw (???)
         */
        const uByte* face = font->Characters[character];
        unsigned W = face[0];
        unsigned H = font->Height;
        if ( col ) gym::setColor(col);
#if 0
        unsigned char pixels[W*H+8];
        gym::unpackBitmap(pixels, W, H, face+1, W);
        gym::drawPixels(W, H, S*X-font->xorig, S*Y-font->yorig, S, pixels);
#else
        unsigned n_bits = countBits(H*((W+7)>>3), 1+face);
        flute2* flu = gym::mapBufferV2(6*n_bits);
        unsigned cnt = gym::unpackBitmap(flu, W, H, S*X-font->xorig, S*Y-font->yorig, S, 1+face);
        gym::unmapBufferV2();
        gym::drawTriangleStrip(0, cnt);
#endif
    }
}


/// use gray for text starting with '%'
static inline void setTextColor(char c, char& d, const float color[4])
{
    float gray[4] = { 0.5, 0.5, 0.5, 1 };
    if ( c == '%' )
    {
        if ( d != '%' )
        {
            d = '%';
            gym::setColor(gray);
        }
    }
    else
    {
        if ( d != 'a' )
        {
            d = 'a';
            gym::setColor(color);
        }
    }
}


/** Converts every character into a pixel map (1 byte/pixel) and then render */
void fgBitmapToken0(float X, float Y, float scale, SFG_Font const* font, const char *string)
{
    float H = font->Height;
    // calculate total string length in pixels:
    unsigned L = 7;
    for ( char const* ptr = string; *ptr; ++ptr )
    {
        unsigned char c = *ptr;
        if ( isprint(c) && c < font->Quantity )
            L += font->Characters[c][0];
    }
    unsigned char * pixels = (unsigned char*)malloc(L*H);
    memset(pixels, 0, L*H);
    unsigned W = 0;
    for ( char const* ptr = string; *ptr; ++ptr )
    {
        unsigned char c = *ptr;
        if ( isprint(c) )
        {
            const uByte* face = font->Characters[c];
            gym::unpackBitmap(pixels+W, face[0], H, 1+face, L);
            W += face[0];
        }
    }
    gym::drawPixels(L, H, X, Y, scale, pixels);
    free(pixels);
}



/** Faster: this converts every character into a triangle strip */
void fgBitmapToken(float X, float Y, float scale, SFG_Font const* font, const char *string)
{
    const unsigned H = font->Height;

    X = scale * ( X - font->xorig );
    Y = scale * ( Y - font->yorig );

    unsigned n_bits = 0;
    for ( char const* c = string; *c; ++c )
    {
        const uByte* face = font->Characters[(unsigned char)*c];
        unsigned cW = face[0];
        n_bits += countBits(H*((cW+7)>>3), 1+face);
    }
    const unsigned sup = 4 * n_bits; // empirical upper limit
    flute2* flu = gym::mapBufferV2(sup);
    unsigned cnt = 0;
    unsigned W = 0;
    for ( char const* c = string; *c; ++c )
    {
        const uByte* face = font->Characters[(unsigned char)*c];
        unsigned cW = face[0];
        cnt += gym::unpackBitmap(flu+cnt, cW, H, X+scale*W, Y, scale, 1+face);
        W += cW;
    }
    //if ( cnt > sup )
    //printf("%6u bits %6u vertex: %s\n", n_bits, cnt, string);
    gym::unmapBufferV2();
    gym::drawTriangleStrip(0, cnt);
}


void fgBitmapToken(float X, float Y, float scale, int fontID, const char *string)
{
    SFG_Font const* font = fghFont(fontID);
    if ( font )
        fgBitmapToken(X, Y, scale, font, string);
}


/** Can handle multi-lines strings containing `\n` characters */
void fgBitmapText(float X, float Y, float scale, int fontID, const float color[4], const char *string, float dY)
{
    SFG_Font const* font = fghFont(fontID);
    if ( !font )
        return;

    char * str = strdup(string);
    for ( char * c = str; *c; ++c )
        if ( !isprint(*c) && *c != '\n' ) *c = '*';

    char col = 0;
    char * token = NULL;
    while ((token = strsep(&str, "\n")) != NULL)
    {
        setTextColor(token[0], col, color);
        fgBitmapToken(X, Y, scale, font, token);
        // move down one line:
        Y += dY;
    }
    free(str);
}

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
