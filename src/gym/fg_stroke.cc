/*
 * fg_stroke.cc
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

#include "fg_stroke.h"
#include "gym_view.h"
#include "gym_draw.h"
#include "gym_flute.h"

/* -- STROKE FONTS ---------------------------------------------------- */

namespace StrokeRoman
{
    extern "C" {
#include "fg_stroke_roman.c"
    }
}
namespace StrokeMonoRoman
{
    extern "C" {
#include "fg_stroke_mono_roman.c"
    }
}

/*
 * Matches a font ID with a SFG_StrokeFont structure pointer.
 * This was changed to match the GLUT header style.
 */
SFG_StrokeFont const* fghStrokeByID( int mono )
{
    if ( mono )
        return &StrokeMonoRoman::fgStrokeFont;
    return &StrokeRoman::fgStrokeFont;
}

/* -- INTERFACE FUNCTIONS -------------------------------------------------- */

/*
 * Draw a stroke character
 */
void fgStrokeCharacter(float X, float Y, float scale, int mono, unsigned char arg,
                       float line_width, float point_size)
{
    scale *= 0.1;
    const SFG_StrokeChar *schar;
    const SFG_StrokeStrip *strip;
    SFG_StrokeFont const* font = fghStrokeByID( mono );
    if ( font && 0 <= arg && arg < font->Quantity )
    {
        schar = font->Characters[arg];
        if ( schar )
        {
            strip = schar->Strips;
            for( unsigned i = 0; i < schar->Number; i++, strip++ )
            {
                const unsigned num = strip->Number;
                const SFG_StrokeVertex* ptr = strip->Vertices;
                flute2* flu = gym::mapBufferV2(num);
                for ( unsigned j = 0; j < num; ++j )
                    flu[j] = { X+scale*ptr[j].X, Y+scale*ptr[j].Y };
                gym::unmapBufferV2();
                if ( line_width > 0 )
                    gym::drawLineStrip(line_width, 0, num);
                if ( point_size > 0 )
                    gym::drawPoints(point_size, 0, num);
            }
        }
    }
}


void fgStrokeString(float X0, float Y, float scale, int mono, const char *string,
                    float line_width, float point_size, float vshift)
{
    const unsigned ONE = 32;
    const unsigned MAX = 16;
    scale *= 0.1;
    float X = X0;
    SFG_StrokeFont const* font = fghStrokeByID( mono );
    if ( font && string )
    {
        unsigned inx = 0, cnt = 0;
        unsigned start[MAX] = { 0 };
        flute2* flu = gym::mapBufferV2(ONE*MAX);
        /*
         * Step through the string, drawing each character.
         * A newline will simply translate the next character's insertion
         * point back to X0, the alignment X-coordinate and down one line.
         */
        unsigned char c;
        while( ( c = *string++) ) if ( c < font->Quantity )
        {
            if ( c == '\n' )
            {
                X = X0;
                if ( vshift == 0 )
                    vshift = font->Height;
                Y -= scale * vshift;
            }
            else  /* Not an EOL, draw the bitmap character */
            {
                const SFG_StrokeChar *schar = font->Characters[c];
                if ( schar )
                {
                    const SFG_StrokeStrip *strip = schar->Strips;
                    for ( unsigned i = 0; i < schar->Number; i++, strip++ )
                    {
                        const unsigned num = strip->Number;
                        const SFG_StrokeVertex* ptr = strip->Vertices;
                        if ( cnt >= MAX )
                        {
                            gym::unmapBufferV2();
                            if ( line_width > 0 )
                            {
                                gym::drawLineStrip(line_width, 0, start[0]);
                                for ( unsigned s = 1; s < cnt; ++s )
                                    gym::drawLineStrip(start[s-1], start[s]-start[s-1]);
                            }
                            if ( point_size > 0 )
                                gym::drawPoints(point_size, 0, inx);
                            cnt = 0;
                            inx = 0;
                            flu = gym::mapBufferV2(ONE*MAX);
                        }
                        // we translate and scale in CPU:
                        for ( unsigned j = 0; j < num; ++j )
                            flu[inx++] = { X+scale*ptr[j].X, Y+scale*ptr[j].Y };
                        start[cnt++] = inx;
                    }
                    X += scale * schar->Advance;
                }
            }
        }
        gym::unmapBufferV2();
        if ( line_width > 0 )
        {
            gym::drawLineStrip(line_width, 0, start[0]);
            for ( unsigned s = 1; s < cnt; ++s )
                gym::drawLineStrip(start[s-1], start[s]-start[s-1]);
        }
        if ( point_size > 0 )
            gym::drawPoints(point_size, 0, inx);
    }
}

/*
 * Return the width in pixels of a stroke character
 */
float fgStrokeWidthf( int mono, unsigned char character )
{
    const SFG_StrokeChar *schar;
    SFG_StrokeFont const* font = fghStrokeByID( mono );
    if ( font && character >= 0 && character < font->Quantity )
    {
        schar = font->Characters[ character ];
        if ( schar )
            return schar->Advance;
    }
    return 0;
}

int fgStrokeWidth( int mono, unsigned char character )
{
    return ( int )( fgStrokeWidthf(mono, character) + 0.5f );
}

/*
 * Return the width of a string drawn using a stroke font
 */
float fgStrokeLengthf( int mono, const unsigned char* string )
{
    unsigned char c;
    float length = 0.0;
    float this_line_length = 0.0;
    SFG_StrokeFont const* font = fghStrokeByID( mono );
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
                        this_line_length += schar->Advance;
                }
            }
        if( length < this_line_length )
            length = this_line_length;
    }
    return length;
}

int fgStrokeLength( int mono, const unsigned char* string )
{
    return( int )( fgStrokeLengthf(mono,string) + 0.5f );
}

/*
 * Returns the height of a stroke font
 */
float fgStrokeHeight( int mono )
{
    SFG_StrokeFont const* font = fghStrokeByID( mono );
    if ( font )
        return font->Height;
    return 0;
}

/*** END OF FILE ***/
