/*
 * fg_font.h
 *
 * Copyright (c) 1999-2000 Pawel W. Olszta. All Rights Reserved.
 * Written by Pawel W. Olszta, <olszta@sourceforge.net>
 * Creation date: Thu Dec 16 1999
 */

/*
 Modified March--April 2022 by FJ. Nedelec,
 to expose only the font-rendering capabilities from the FreeGLUT project
 */

/*
 * Draw a bitmap character
 */
void fgBitmapCharacter(float x, float y, float S, int font, const float color[4], int character);

/*
 * Draw a bitmap string at position (x, y) with pixel-size S
 */
void fgBitmapToken(float x, float y, float S, int font, const char *string);


/*
 * Draw a bitmap multi-lined string at position (x, y) with pixel-size S
 */
void fgBitmapText(float x, float y, float S, int font, const float color[4], const char *string, float vshift);


/*
 * Returns the width in pixels of a font's character
 */
int fgBitmapWidth(int font, unsigned char character);

/*
* Return the height of a font
*/
int fgFontHeight(int font);

/*
* Return the width of a string drawn using a bitmap font, and the number of lines
*/
int fgTextWidth(int font, const char text[], int& num_lines);

/*
 * Return the width of a string drawn using a bitmap font
 */
int fgBitmapLength(int font, const unsigned char* string);

/*
 * Returns the height of a bitmap font
 */
int fgBitmapHeight(int font);

