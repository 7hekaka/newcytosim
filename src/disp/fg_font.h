/*
 * fg_font.h
 *
 * Copyright (c) 1999-2000 Pawel W. Olszta. All Rights Reserved.
 * Written by Pawel W. Olszta, <olszta@sourceforge.net>
 * Creation date: Thu Dec 16 1999
 */

/*
 Created on 20.03.2022 by FJ. Nedelec, to expose only the font-rendering capabilities
 from the FreeGLUT project
 */

#include <stdio.h>
#include "opengl.h"


/* The bitmap font structure */
struct tagSFG_Font
{
    int             Quantity;     /* Number of chars in font          */
    int             Height;       /* Height of the characters         */
    const GLubyte** Characters;   /* The characters mapping           */
    float           xorig, yorig; /* Relative origin of the character */
};
typedef struct tagSFG_Font SFG_Font;


/* The stroke font structures */
struct tagSFG_StrokeVertex
{
    GLfloat X, Y;
};
typedef struct tagSFG_StrokeVertex SFG_StrokeVertex;

struct tagSFG_StrokeStrip
{
    int Number;
    const SFG_StrokeVertex* Vertices;
};
typedef struct tagSFG_StrokeStrip SFG_StrokeStrip;

struct tagSFG_StrokeChar
{
    GLfloat Right;
    int     Number;
    const SFG_StrokeStrip* Strips;
};
typedef struct tagSFG_StrokeChar SFG_StrokeChar;

struct tagSFG_StrokeFont
{
    int             Quantity;                   /* Number of chars in font   */
    GLfloat         Height;                     /* Height of the characters  */
    const SFG_StrokeChar** Characters;          /* The characters mapping    */
};
typedef struct tagSFG_StrokeFont SFG_StrokeFont;


/*
* return bitmap font
*/
SFG_Font const* fghFontByID(int font);

/*
* return stroke font
*/
SFG_StrokeFont const* fghStrokeByID(int font);

/*
 * Draw a bitmap character
 */
void fgBitmapCharacter(int font, int character);

/*
 * Draw a bitmap string
 */
void fgBitmapString(int font, const unsigned char *string, float vshift);

/*
 * Returns the width in pixels of a font's character
 */
int fgBitmapWidth(int font, int character);

/*
 * Return the width of a string drawn using a bitmap font
 */
int fgBitmapLength(int font, const unsigned char* string);

/*
 * Returns the height of a bitmap font
 */
int fgBitmapHeight(int font);

/*
 * Draw a stroke character
 */
void fgStrokeCharacter(int font, int character, int drawJoinDots);

/*
* Draw a stroke string
*/
void fgStrokeString(int font, const unsigned char *string);

/*
 * Return the width in pixels of a stroke character
 */
int fgStrokeWidth(int font, int character);

/*
 * Return the width of a string drawn using a stroke font
 */
int fgStrokeLength(int font, const unsigned char* string);

/*
 * Returns the height of a stroke font
 */
GLfloat fgStrokeHeight(int font);


