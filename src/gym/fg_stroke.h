/*
 * fg_stroke.h
 *
 * Copyright (c) 1999-2000 Pawel W. Olszta. All Rights Reserved.
 * Written by Pawel W. Olszta, <olszta@sourceforge.net>
 * Creation date: Thu Dec 16 1999
 */

/*
 Created on 20.03.2022 by FJ. Nedelec, to expose only the font-rendering capabilities
 from the FreeGLUT project, itself derived from GLUT. The two font faces use the same characters.

 NON_MONO FONT:
     A  proportionally  spaced  Roman  Simplex  font  for ASCII characters 32 through 127.
     The maximum top character in the font is 119.05 units; the bottom descends 33.33 units.

 MONO FONT:
     A mono-spaced spaced Roman Simplex font for ASCII characters 32 through 127.
     The maximum top character in the font is 119.05 units; the bottom descends 33.33 units.
     Each character is 104.76 units wide.

 */


/* The stroke font structures */
struct tagSFG_StrokeVertex
{
    float X, Y;
};
typedef struct tagSFG_StrokeVertex SFG_StrokeVertex;

struct tagSFG_StrokeStrip
{
    unsigned Number;
    const SFG_StrokeVertex* Vertices;
};
typedef struct tagSFG_StrokeStrip SFG_StrokeStrip;

struct tagSFG_StrokeChar
{
    float Advance;
    unsigned Number;
    const SFG_StrokeStrip* Strips;
};
typedef struct tagSFG_StrokeChar SFG_StrokeChar;

struct tagSFG_StrokeFont
{
    unsigned Quantity;                /* Number of chars in font   */
    float Height;                     /* Height of the characters  */
    const SFG_StrokeChar** Characters; /* The characters mapping    */
};
typedef struct tagSFG_StrokeFont SFG_StrokeFont;



/*
 * Draw a stroke character. With `scale=pixel_size` the font is ~15 pixels high.
 */
void fgStrokeCharacter(float X, float Y, float scale, int mono, unsigned char character,
                       float stroke_width, float stroke_size);

/*
* Draw a stroke string. With `scale=pixel_size` the font is ~15 pixels high.
*/
void fgStrokeString(float x, float y, float scale, int mono, const char *string,
                    float stroke_width, float stroke_size = 0, float vshift = 0.);

/*
 * Return the width in pixels of a stroke character
 */
int fgStrokeWidth(int mono, unsigned char character);

/*
 * Return the width of a string drawn using a stroke font
 */
int fgStrokeLength(int mono, const unsigned char* string);

/*
 * Returns the height of a stroke font
 */
float fgStrokeHeight(int mono);


