// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include <cctype>
#include "point_disp.h"
#include "offscreen.h"
#include "glossary.h"
#include "glut.h"
#include "gle.h"
#include "gle_zoo.h"


/// if this is defined, the pixelmap are stored in graphical memory
#define POINTDISP_USES_PIXEL_BUFFERS 1


void PointDisp::clearPixelmaps()
{
#if POINTDISP_USES_PIXELMAPS
    nPix = 0;
    for ( int i = 0; i < 3; ++i )
    {
        bmp_[i] = nullptr;
        pbo_[i] = 0;
    }
#endif
}


PointDisp::PointDisp(const std::string& k, const std::string& n)
: Property(n)
{
    mKind = k;
    clearPixelmaps();
    clear();
}


PointDisp::PointDisp(PointDisp const& o) : Property(o)
{
    mKind        = o.mKind;
    visible      = o.visible;
    color        = o.color;
    color2       = o.color2;
    coloring     = o.coloring;
    size         = o.size;
    width        = o.width;
    shape        = o.shape;
    style        = o.style;
    symbol       = o.symbol;
    symbol_color = o.symbol_color;

    clearPixelmaps();
}


PointDisp& PointDisp::operator = (PointDisp const& o)
{
    mKind        = o.mKind;
    visible      = o.visible;
    color        = o.color;
    color2       = o.color2;
    coloring     = o.coloring;
    size         = o.size;
    width        = o.width;
    shape        = o.shape;
    style        = o.style;
    symbol       = o.symbol;
    symbol_color = o.symbol_color;
    
    clearPixelmaps();
    
    return *this;
}


PointDisp::~PointDisp()
{
#if POINTDISP_USES_PIXELMAPS
    releasePixelmap();
    
#if POINTDISP_USES_PIXEL_BUFFERS
    if ( pbo_[0] > 0 )
    {
        glDeleteBuffers(3, pbo_);
        pbo_[0] = 0;
        pbo_[1] = 0;
        pbo_[2] = 0;
    }
#endif
#endif
}


void PointDisp::clear()
{
    visible      = 1;
    color        = 0x888888FF;
    color2       = 0x777777FF;
    coloring     = 0;
    size         = 5;
    width        = 2;
    shape        = 'o';
    style        = 7;
    symbol       = 0;
    symbol_color = 0xFFFFFFFF;
}

void PointDisp::strokeA() const
{
    paintShape();
    if ( width > 0.5 )
    {
        //draw a bright rim
        color.darken(2.0).load();
        glLineWidth(3.0);
        strokeShape();
    }
    
    if ( symbol )
    {
        glScalef(1.0f/80, 1.0f/80, 1);
        /*  glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, C)
         character C of width ~104.76 units, and ~150 unit high max
         The translation brings it near the center. */
        if ( islower(symbol) )
            glTranslatef(-52.35f, -35, 0);
        else
            glTranslatef(-52.35f, -50, 0);
        symbol_color.load();
        glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, symbol);
    }
}

void PointDisp::strokeI() const
{
    paintShape();
    
    // radius of the spot that indicate an inactive Hand
    const GLfloat DOT_SIZE = 0.55f;

    // draw a transparent hole in the center:
    GLboolean alpha = glIsEnabled(GL_ALPHA_TEST);
    GLboolean blend = glIsEnabled(GL_BLEND);
    glDisable(GL_ALPHA_TEST);
    glDisable(GL_BLEND);
    glPushMatrix();
    glScalef(DOT_SIZE, DOT_SIZE, DOT_SIZE);
    glColor4f(0, 0, 0, 0);
    gle::disc(); //strokeShape();
    glPopMatrix();
    if ( alpha ) glEnable(GL_ALPHA_TEST);
    if ( blend ) glEnable(GL_BLEND);
}


#pragma mark - Bitmaps

#if POINTDISP_USES_PIXELMAPS


/**
 Allocate memory for pixelmaps of size `pixSize x pixSize`
 3 bitmaps x 4 colors x pixSize x pixSize pixels
 \todo use Graphic card memory for the bitmaps
 */
void PointDisp::allocatePixelmap()
{
    mOffs = -0.5f * pixSize;
    
    // allocate only if needed
    if ( !bmp_[0] || pixSize != nPix )
    {
        if ( bmp_[0] )
            delete(bmp_[0]);
    
        size_t dd = pixSize * pixSize;
        uint8_t * mem = new uint8_t[12*dd];
        bmp_[0] = mem;
        bmp_[1] = mem + 4*dd;
        bmp_[2] = mem + 8*dd;
        
        for ( size_t y = 0; y < 12*dd; ++y )
            mem[y] = 0;
        
        nPix = pixSize;
    }
}


void PointDisp::releasePixelmap()
{
    if ( bmp_[0] )
    {
        delete(bmp_[0]);
        bmp_[0] = nullptr;
    }
    nPix = 0;
}


/// print pixel map in ASCII
void printPixels(uint8_t const* pix, unsigned sx, unsigned sy)
{
    for ( size_t y = 0; y < sy; ++y )
    {
        for ( size_t x = 0; x < sx; ++x )
            printf("%02X", pix[4*(x+sx*y)+3]);
        printf("\n");
    }
}

/**
 This will downsample pixelmap `src` and set destination `dst`. The pixel
 array `dst` should be of size `4*sx*sy` with 4 bytes per pixels: R, G, B and A,
 while `src` should be `bin*bin` times larger. Pixels are stored in row order
 from the lowest to the highest row, left to right in each row (as in OpenGL).
 The pixels components of `src` are averaged to produce `dst`.
 */
void PointDisp::downsampleRGBA(uint8_t dst[], unsigned sx, unsigned sy,
                               uint8_t const src[], unsigned bin)
{
    const size_t s = bin * bin;

#if ( 0 )
    //reset destination:
    for ( size_t u = 0; u < sx*sy; ++u )
    {
        dst[4*u  ] = 0xFF;
        dst[4*u+1] = 0xFF;
        dst[4*u+2] = 0xFF;
        dst[4*u+3] = 0xFF;
    }
#endif
    
    for ( unsigned x = 0; x < sx; ++x )
    for ( unsigned y = 0; y < sy; ++y )
    {
        size_t r = 0, g = 0, b = 0, a = 0;
        for ( unsigned dx = 0; dx < bin; ++dx )
        for ( unsigned dy = 0; dy < bin; ++dy )
        {
            uint8_t const* p = src + 4 * ( dx+bin*(x+sx*(dy+bin*y)) );
            r += p[0];
            g += p[1];
            b += p[2];
            a += p[3];
        }
            
        dst[4*(x+sx*y)  ] = (uint8_t)( r / s );
        dst[4*(x+sx*y)+1] = (uint8_t)( g / s );
        dst[4*(x+sx*y)+2] = (uint8_t)( b / s );
        dst[4*(x+sx*y)+3] = (uint8_t)( a / s );
    }
}


void PointDisp::storePixelmap(uint8_t* bitmap, unsigned dim, GLuint pbi) const
{
#if POINTDISP_USES_PIXEL_BUFFERS
    assert_true(pbi);
    glBindBuffer(GL_PIXEL_PACK_BUFFER_ARB, pbi);
    glBufferData(GL_PIXEL_PACK_BUFFER_ARB, 4*dim*dim, bitmap, GL_STATIC_DRAW);
    glBindBuffer(GL_PIXEL_PACK_BUFFER_ARB, 0);
    assert_true(glIsBuffer(pbi));
    CHECK_GL_ERROR("PointDisp::storePixelmap");
#endif
}


#if ( 0 )

#include "saveimage.h"

// Export bitmap to file in PNG format
void PointDisp::savePixelmap(uint8_t* bitmap, unsigned dim, unsigned id) const
{
    if ( SaveImage::supported("png") )return given name of property
    {
       char str[32];
        snprintf(str, sizeof(str), "bitmap_%s_%02u.png", name_str(), id);
        FILE * f = fopen(str, "w");
        if ( f )
        {
            if ( !ferror(f) )
            {
                SaveImage::saveAlphaPNG(f, bitmap, dim, dim);
                fclose(f);
                std::clog << "PointDisp saved " << str << '\n';
            }
            fclose(f);
        }
    }
}
#endif


void PointDisp::drawPixelmap(size_t inx) const
{
    //translate to center the bitmap:
    glBitmap(0,0,0,0,mOffs,mOffs,nullptr);
#if POINTDISP_USES_PIXEL_BUFFERS
    glBindBuffer(GL_PIXEL_UNPACK_BUFFER, pbo_[inx]);
    glDrawPixels(nPix, nPix, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
    //glBindBuffer(GL_PIXEL_UNPACK_BUFFER, 0);
#else
    glDrawPixels(nPix, nPix, GL_RGBA, GL_UNSIGNED_BYTE, bmp_[inx]);
#endif
    CHECK_GL_ERROR("PointDisp::drawPixelmap");
}


/**
 `sampling` defines the level of oversampling used to improve the quality of bitmaps
 */
void PointDisp::makePixelmaps(GLfloat unit_value, unsigned sampling)
{
    assert_true(pixSize==nPix);
    CHECK_GL_ERROR("1 PointDisp::makePixelmaps");
    
    unsigned dim = sampling * nPix;
    GLfloat s = 0.5f * sampling * size * unit_value;
    GLfloat t = 0.5f * dim;
    GLfloat w = width * sampling * unit_value;
    
    GLboolean alpha = glIsEnabled(GL_ALPHA_TEST);
    GLboolean light = glIsEnabled(GL_LIGHTING);
    GLboolean blend = glIsEnabled(GL_BLEND);
    glDisable(GL_ALPHA_TEST);
    glDisable(GL_LIGHTING);
    glDisable(GL_BLEND);

    GLint svp[4];
    glGetIntegerv(GL_VIEWPORT, svp);

    //match projection to viewport:
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();

#ifdef __APPLE__
    OffScreen::createBuffer(dim, dim, 0);
    glOrtho(0, dim, 0, dim, 0, 1);
#else
    glPushAttrib(GL_PIXEL_MODE_BIT|GL_VIEWPORT_BIT|GL_ENABLE_BIT|GL_COLOR_BUFFER_BIT);
    glOrtho(0, svp[2], 0, svp[3], 0, 1);
#endif
    
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glTranslatef(t, t, 0);
    glScalef(s,s,s);
    if ( w > 0 ) glLineWidth(w);
    // we use a transparent background, because points will overlap
    glClearColor(0,0,0,0);

    for ( int i = 0; i < 3; ++i )
    {
        glClear(GL_COLOR_BUFFER_BIT);
        switch ( i )
        {
            case 0:
                color2.load();
                strokeI();
                break;
            case 1:
                color2.load();
                strokeA();
                break;
            case 2:
                color.load();
                strokeA();
                break;
        }
        if ( sampling > 1 )
        {
            uint8_t * tmp = new uint8_t[4*dim*dim];
            glReadPixels(0, 0, dim, dim, GL_RGBA, GL_UNSIGNED_BYTE, tmp);
            downsampleRGBA(bmp_[i], nPix, nPix, tmp, sampling);
            delete[] tmp;
#if ( 0 )
            //savePixelmap(tmp, dim, i);
            std::clog << name() << i << "\n";
            printPixels(bmp_[i], nPix, nPix);
#endif
        }
        else
        {
            glReadPixels(0, 0, nPix, nPix, GL_RGBA, GL_UNSIGNED_BYTE, bmp_[i]);
            //savePixelmap(bmp_[i], nPix, i+10);
        }
        CHECK_GL_ERROR("5 PointDisp::makePixelmaps");
        storePixelmap(bmp_[i], nPix, pbo_[i]);
    }
 
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
    
    if ( alpha ) glEnable(GL_ALPHA_TEST);
    if ( light ) glEnable(GL_LIGHTING);
    if ( blend ) glEnable(GL_BLEND);

#ifdef __APPLE__
    OffScreen::releaseBuffer();
    glViewport(svp[0], svp[1], svp[2], svp[3]);
#else
    glPopAttrib();
#endif
}

#endif


void PointDisp::createPixelmaps(GLfloat uf)
{
#if POINTDISP_USES_PIXEL_BUFFERS
    if ( pbo_[0] == 0 )
        glGenBuffers(3, pbo_);
#endif
    
    if ( pixSize != nPix )
    {
        CHECK_GL_ERROR("1 PointDisp::prepare");
        allocatePixelmap();
        //fprintf(stderr, " new %i bitmap for %s\n", pixSize, name_str());
        CHECK_GL_ERROR("2 PointDisp::prepare");
    }
    
    makePixelmaps(uf, 3);
}


void PointDisp::prepare(GLfloat uf, GLfloat sf, bool pixelmaps)
{
    realSize    = size * sf;
    unsigned sz = static_cast<unsigned>(std::ceil(uf*(size+width)));
    // make it a multiple of 4:
    pixSize     = ( sz + 4U ) & ~3U;
    perceptible = visible && ( uf*(size+width) > 0.25 );
    
#if POINTDISP_USES_PIXELMAPS
    if ( pixelmaps )
        createPixelmaps(uf);
#endif
}


#pragma mark - I/O


void PointDisp::read(Glossary& glos)
{
    glos.set(visible,      "visible");
    
    // set 'color2' as a darker tone of 'color':
    if ( glos.set(color,   "color") )
        color2 = color.alpha_scaled(DIM==2?0.25f:0.5f);
    glos.set(color2,       "color", 1, "back_color", 0);
    glos.set(coloring,     "coloring");
    
    // if 'size' is specified, width is set accordingly:
    if ( glos.set(size,    "size") )
        width = 2 * size / 3;
    else
        glos.set(size,     "point_size");
#if BACKWARD_COMPATIBILITY < 100
    glos.set(size,         "points");
    glos.set(shape,        "points", 1);
#endif

    glos.set(width,        "width");
    glos.set(style,        "style");
    glos.set(shape,        "shape");
    glos.set(symbol,       "symbol");
    glos.set(symbol_color, "symbol", 1);
    
    if ( ! isprint(symbol) )
        symbol = 0;
    shape = tolower(shape);
    
#if POINTDISP_USES_PIXELMAPS
    releasePixelmap();
#endif
}


void PointDisp::write_values(std::ostream& os) const
{
    write_value(os, "visible",     visible);
    if ( color2 != color.alpha_scaled(DIM==2?0.25f:0.5f) )
        write_value(os, "color",   color, color2);
    else
        write_value(os, "color",   color);
    write_value(os, "coloring",    coloring);
    write_value(os, "size",        size);
    write_value(os, "width",       width);
    write_value(os, "shape",       shape);
    write_value(os, "style",       style);
    write_value(os, "symbol",      symbol, symbol_color);
}

