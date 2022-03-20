/*
 * fg_font.c
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
 Modified on 20.03.2022 by FJ. Nedelec, to remove dependencies on FreeGLUT
 Keeping only the font-rendering capabilities
 */


/* The bitmap font structure */
typedef struct tagSFG_Font SFG_Font;
struct tagSFG_Font
{
    int             Quantity;     /* Number of chars in font          */
    int             Height;       /* Height of the characters         */
    const GLubyte** Characters;   /* The characters mapping           */
    float           xorig, yorig; /* Relative origin of the character */
};

/* The stroke font structures */
typedef struct tagSFG_StrokeVertex SFG_StrokeVertex;
struct tagSFG_StrokeVertex
{
    GLfloat         X, Y;
};

typedef struct tagSFG_StrokeStrip SFG_StrokeStrip;
struct tagSFG_StrokeStrip
{
    int             Number;
    const SFG_StrokeVertex* Vertices;
};

typedef struct tagSFG_StrokeChar SFG_StrokeChar;
struct tagSFG_StrokeChar
{
    GLfloat         Right;
    int             Number;
    const SFG_StrokeStrip* Strips;
};

typedef struct tagSFG_StrokeFont SFG_StrokeFont;
struct tagSFG_StrokeFont
{
    int             Quantity;                   /* Number of chars in font   */
    GLfloat         Height;                     /* Height of the characters  */
    const SFG_StrokeChar** Characters;          /* The characters mapping    */
};

/*
 * TODO BEFORE THE STABLE RELEASE:
 *
 *  Test things out ...
 */

/* -- IMPORT DECLARATIONS -------------------------------------------------- */

#include "fg_font_data.c"

/*
 * Matches a font ID with a SFG_Font structure pointer.
 * This was changed to match the GLUT header style.
 */
static SFG_Font const* fghFontByID( FontType font )
{
    if( font == BITMAP_8_BY_13        )
        return &fgFontFixed8x13;
    if( font == BITMAP_9_BY_15        )
        return &fgFontFixed9x15;
    if( font == BITMAP_HELVETICA_10   )
        return &fgFontHelvetica10;
    if( font == BITMAP_HELVETICA_12   )
        return &fgFontHelvetica12;
    if( font == BITMAP_HELVETICA_18   )
        return &fgFontHelvetica18;
    if( font == BITMAP_TIMES_ROMAN_10 )
        return &fgFontTimesRoman10;
    if( font == BITMAP_TIMES_ROMAN_24 )
        return &fgFontTimesRoman24;

    return nullptr;
}

/* -- STROKE FONTS ---------------------------------------------------- */

namespace Roman {
#include "fg_stroke_roman.c"
const SFG_StrokeFont fgStrokeFont = {128,152.381f,chars};
}

namespace MonoRoman {
#include "fg_stroke_mono_roman.c"
const SFG_StrokeFont fgStrokeFont = {128,152.381f,chars};
}

/*
 * Matches a font ID with a SFG_StrokeFont structure pointer.
 * This was changed to match the GLUT header style.
 */
static SFG_StrokeFont const* fghStrokeByID( FontType font )
{
    if( font == STROKE_ROMAN      )
        return &Roman::fgStrokeFont;
    if( font == STROKE_MONO_ROMAN )
        return &MonoRoman::fgStrokeFont;
    return nullptr;
}

/* -- INTERFACE FUNCTIONS -------------------------------------------------- */

/*
 * Draw a bitmap character
 */
void fgBitmapCharacter( FontType fontID, int character )
{
    const GLubyte* face;
    SFG_Font const* font = fghFontByID( fontID );
    if ( font && character >= 1 && character < 256 )
    {
        /*
         * Find the character we want to draw (???)
         */
        face = font->Characters[ character ];
        
        glPushClientAttrib( GL_CLIENT_PIXEL_STORE_BIT );
        glPixelStorei( GL_UNPACK_SWAP_BYTES,  GL_FALSE );
        glPixelStorei( GL_UNPACK_LSB_FIRST,   GL_FALSE );
        glPixelStorei( GL_UNPACK_ROW_LENGTH,  0        );
        glPixelStorei( GL_UNPACK_SKIP_ROWS,   0        );
        glPixelStorei( GL_UNPACK_SKIP_PIXELS, 0        );
        glPixelStorei( GL_UNPACK_ALIGNMENT,   1        );
        glBitmap(
                 face[ 0 ], font->Height,      /* The bitmap's width and height  */
                 font->xorig, font->yorig,     /* The origin in the font glyph   */
                 ( float )( face[ 0 ] ), 0.0,  /* The raster advance -- inc. x,y */
                 ( face + 1 )                  /* The packed bitmap data...      */
                 );
        glPopClientAttrib( );
    }
}

void fgBitmapString( FontType fontID, const unsigned char *string, float vshift )
{
    unsigned char c;
    float x = 0.0f ;
    SFG_Font const* font = fghFontByID( fontID );
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
                glBitmap ( 0, 0, 0, 0, -x, -vshift, nullptr );
                x = 0.0f;
            }
            else  /* Not an EOL, draw the bitmap character */
            {
                const GLubyte* face = font->Characters[ c ];
                
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

/*
 * Returns the width in pixels of a font's character
 */
int fgBitmapWidth( FontType fontID, int character )
{
    SFG_Font const* font = fghFontByID( fontID );
    if ( font && character > 0 && character < 256 )
        return *( font->Characters[ character ] );
    return 0;
}

/*
 * Return the width of a string drawn using a bitmap font
 */
int fgBitmapLength( FontType fontID, const unsigned char* string )
{
    unsigned char c;
    int length = 0, this_line_length = 0;
    SFG_Font const* font = fghFontByID( fontID );
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
int fgBitmapHeight( FontType fontID )
{
    SFG_Font const* font = fghFontByID( fontID );
    if ( font )
        return font->Height;
    return 0;
}

/*
 * Draw a stroke character
 */
void fgStrokeCharacter( FontType fontID, int character, bool drawJoinDots )
{
    const SFG_StrokeChar *schar;
    const SFG_StrokeStrip *strip;
    int i, j;
    SFG_StrokeFont const* font = fghStrokeByID( fontID );
    if ( font && character >= 0 && character < font->Quantity )
    {
        schar = font->Characters[ character ];
        if ( schar )
        {
            strip = schar->Strips;
            
            for( i = 0; i < schar->Number; i++, strip++ )
            {
                glBegin( GL_LINE_STRIP );
                for( j = 0; j < strip->Number; j++ )
                    glVertex2f( strip->Vertices[ j ].X, strip->Vertices[ j ].Y );
                glEnd( );
                
                if ( drawJoinDots )
                {
                    glBegin( GL_POINTS );
                    for( j = 0; j < strip->Number; j++ )
                        glVertex2f( strip->Vertices[ j ].X, strip->Vertices[ j ].Y );
                    glEnd( );
                }
            }
            glTranslatef( schar->Right, 0.0, 0.0 );
        }
    }
}

void fgStrokeString( FontType fontID, const unsigned char *string )
{
    unsigned char c;
    int i, j;
    float length = 0.0;
    SFG_StrokeFont const* font = fghStrokeByID( fontID );
    if ( font && string )
    {
        /*
         * Step through the string, drawing each character.
         * A newline will simply translate the next character's insertion
         * point back to the start of the line and down one line.
         */
        while( ( c = *string++) ) if( c < font->Quantity )
        {
            if( c == '\n' )
            {
                glTranslatef ( -length, -( float )( font->Height ), 0.0 );
                length = 0.0;
            }
            else  /* Not an EOL, draw the bitmap character */
            {
                const SFG_StrokeChar *schar = font->Characters[ c ];
                if( schar )
                {
                    const SFG_StrokeStrip *strip = schar->Strips;
                    
                    for( i = 0; i < schar->Number; i++, strip++ )
                    {
                        glBegin( GL_LINE_STRIP );
                        for( j = 0; j < strip->Number; j++ )
                            glVertex2f( strip->Vertices[ j ].X,
                                       strip->Vertices[ j ].Y);
                        
                        glEnd( );
                    }
                    
                    length += schar->Right;
                    glTranslatef( schar->Right, 0.0, 0.0 );
                }
            }
        }
    }
}

/*
 * Return the width in pixels of a stroke character
 */
GLfloat fgStrokeWidthf( FontType fontID, int character )
{
    const SFG_StrokeChar *schar;
    SFG_StrokeFont const* font = fghStrokeByID( fontID );
    if ( font && character >= 0 && character < font->Quantity )
    {
        schar = font->Characters[ character ];
        if ( schar )
            return schar->Right;
    }
    return 0;
}

int fgStrokeWidth(FontType fontID, int character)
{
    return ( int )( fgStrokeWidthf(fontID, character) + 0.5f );
}

/*
 * Return the width of a string drawn using a stroke font
 */
GLfloat fgStrokeLengthf( FontType fontID, const unsigned char* string )
{
    unsigned char c;
    GLfloat length = 0.0;
    GLfloat this_line_length = 0.0;
    SFG_StrokeFont const* font = fghStrokeByID( fontID );
    if ( font && string )
    {
        while( ( c = *string++) )
            if( c < font->Quantity )
            {
                if( c == '\n' ) /* EOL; reset the length of this line */
                {
                    if( length < this_line_length )
                        length = this_line_length;
                    this_line_length = 0.0;
                }
                else  /* Not an EOL, increment the length of this line */
                {
                    const SFG_StrokeChar *schar = font->Characters[ c ];
                    if( schar )
                        this_line_length += schar->Right;
                }
            }
        if( length < this_line_length )
            length = this_line_length;
    }
    return length;
}

int fgStrokeLength( FontType fontID, const unsigned char* string )
{
    return( int )( fgStrokeLengthf(fontID,string) + 0.5f );
}

/*
 * Returns the height of a stroke font
 */
GLfloat fgStrokeHeight( FontType fontID )
{
    SFG_StrokeFont const* font = fghStrokeByID( fontID );
    if ( font )
        return font->Height;
    return 0;
}

/*** END OF FILE ***/
