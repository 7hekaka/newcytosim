// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include <cctype>
#include "point_disp.h"
#include "offscreen.h"
#include "glossary.h"
#include "gle.h"
#include "gym_text.h"
#include "gym_check.h"
#include "gym_view.h"
#include "gym_draw.h"
#include "gym_image.h"


/// if this is defined, the pixelmap are stored in graphical memory
#define POINTDISP_USES_PIXEL_BUFFERS 0


void PointDisp::clearPixelmaps()
{
#if POINTDISP_USES_PIXELMAPS
    pixSize = 0;
    pixAlloc_ = 0;
    for ( int i = 0; i < 3; ++i )
    {
        pixels_ = nullptr;
        pbo_ = 0;
    }
#endif
}


PointDisp::PointDisp(const std::string& k, const std::string& n)
: Property(n)
{
    mKind = k;
    clearPixelmaps();
    clear();
    sizeR = 0;
    sizeX = 0;
    widthX = 0;
}


PointDisp::PointDisp(PointDisp const& o) : Property(o)
{
    mKind    = o.mKind;
    visible  = o.visible;
    color    = o.color;
    color2   = o.color2;
    coloring = o.coloring;
    size     = o.size;
    width    = o.width;
    shape    = o.shape;
    style    = o.style;
    symbol   = o.symbol;
    colorS   = o.colorS;

    clearPixelmaps();
}


PointDisp& PointDisp::operator = (PointDisp const& o)
{
    mKind    = o.mKind;
    visible  = o.visible;
    color    = o.color;
    color2   = o.color2;
    coloring = o.coloring;
    size     = o.size;
    width    = o.width;
    shape    = o.shape;
    style    = o.style;
    symbol   = o.symbol;
    colorS   = o.colorS;
    
    clearPixelmaps();
    return *this;
}


PointDisp::~PointDisp()
{
#if POINTDISP_USES_PIXELMAPS
    releasePixelmap();
#endif
}


void PointDisp::clear()
{
    visible  = 1;
    color    = 0x888888FF;
    color2   = 0x777777FF;
    coloring = 0;
    size   = 5;
    width  = 2;
    shape  = 'o';
    style  = 7;
    symbol = 0;
    colorS = 0xFFFFFFFF;
}

void PointDisp::strokeA(float w) const
{
    paintShape();
    if ( w > 0 )
    {
        gym::color(color.darken(2.0));
        strokeShape(w);
        if ( symbol )
        {
            /* Character C of width ~104.76 units, and ~150 unit high max
             The translation brings it near the center. */
            const float G = 0.0125;
            const float x = -52.35 * G;
            const float y = ( islower(symbol) ? -35 : -50 ) * G;
            gym::strokeCharacter(x, y, G, 1, colorS, symbol, w);
        }
    }
}

void PointDisp::strokeI() const
{
    paintShape();
    
    // draw a transparent hole in the center:
    gym::scale(0.5f);
    gym::color(0,0,0,0);
    gle::disc();
    gym::scale(2.0f);
}


#pragma mark - Bitmaps


#if ( 0 )
#include "save_image.h"

void savePixelmap(uint8_t* bitmap, unsigned dim, unsigned id, char const* name)
{
    if ( SaveImage::supported("png") )
    {
        char str[64];
        snprintf(str, sizeof(str), "bitmap_%s_%02u.png", name, id);
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


#if POINTDISP_USES_PIXELMAPS


/**
 Allocate memory for pixelmaps of size `pixSize x pixSize`
 3 bitmaps x 4 colors x pixSize x pixSize pixels
 */
void PointDisp::allocatePixelmap(unsigned dim)
{
    pixAlloc_ = pixSize;
    pixStride = ( 4 * pixSize * pixSize + 7U ) & ~7U;
    delete(pixels_);
    pixels_ = new uint8_t[3*pixStride]();
}


void PointDisp::releasePixelmap()
{
    delete(pixels_);
    pixels_ = nullptr;
    pixAlloc_ = 0;
#if POINTDISP_USES_PIXEL_BUFFERS
    if ( pbo_ > 0 )
    {
        glDeleteBuffers(1, &pbo_);
        pbo_ = 0;
    }
#endif
}


void PointDisp::drawPixelmap(size_t inx) const
{
    //translate to center the bitmap:
    glBitmap(0,0,0,0,pixOffset_,pixOffset_,nullptr);
#if POINTDISP_USES_PIXEL_BUFFERS
    glBindBuffer(GL_PIXEL_UNPACK_BUFFER, pbo_);
    glDrawPixels(pixSize, pixSize, GL_RGBA, GL_UNSIGNED_BYTE, (void*)(inx*pixStride));
    //glBindBuffer(GL_PIXEL_UNPACK_BUFFER, 0);
#else
    glDrawPixels(pixSize, pixSize, GL_RGBA, GL_UNSIGNED_BYTE, pixels_+inx*pixStride);
#endif
    CHECK_GL_ERROR("PointDisp::drawPixelmap");
}


/**
 `sampling` defines the level of oversampling used to improve the quality of bitmaps
 */
void PointDisp::makePixelmaps(float unit_value, unsigned sampling, unsigned dim)
{
    float s = size * unit_value * 0.5f;
    float w = width * unit_value * 0.5f;
    gym::translate(dim*0.5, dim*0.5, 0);
    gym::scale(s);

    if ( width > 0 ) glLineWidth(w);

    for ( int i = 0; i < 3; ++i )
    {
        uint8_t * pix = pixels_ + i * pixStride;
        // we use a transparent background, because points will overlap
        gym::clearPixels(0,0,0,0);
        switch ( i )
        {
            case 0:
                gym::color(color2);
                strokeI();
                break;
            case 1:
                gym::color(color2);
                strokeA(w);
                break;
            case 2:
                gym::color(color);
                strokeA(w);
                break;
        }
        if ( sampling > 1 )
        {
            uint8_t * tmp = new uint8_t[4*dim*dim];
            glReadPixels(0, 0, dim, dim, GL_RGBA, GL_UNSIGNED_BYTE, tmp);
            gym::downsampleRGBA(pix, pixSize, pixSize, tmp, sampling);
            //gym::printPixels(stdout, tmp, dim, dim);
            delete[] tmp;
        }
        else
        {
            glReadPixels(0, 0, pixSize, pixSize, GL_RGBA, GL_UNSIGNED_BYTE, pix);
        }
        //for ( size_t u = 0; u < pixStride; ++u ) pix[u] = 255;
#if ( 0 )
        //savePixelmap(pix, pixSize, i, name_str());
        std::clog << name() << i << " " << size << " " << width << "\n";
        gym::printPixels(stdout, pix, pixSize, pixSize);
#endif
        CHECK_GL_ERROR("5 PointDisp::makePixelmaps");
    }
}

/**
 `sampling` defines the level of oversampling used to improve the quality of bitmaps
 */
void PointDisp::makePixelmaps(float unit_value, unsigned sampling)
{
    CHECK_GL_ERROR("1 PointDisp::makePixelmaps");
    GLint svp[4];
    glGetIntegerv(GL_VIEWPORT, svp);
    
    unsigned dim = sampling * pixSize;
    unsigned buf = OffScreen::openBuffer(dim, dim, 0);
    if ( buf )
        gym::one_view(dim, dim);
    else
        gym::one_view(svp[2], svp[3]);

    gym::disableLighting();
    gym::disableBlending();
    gym::disableAlphaTest();
    //gym::printCaps();
    makePixelmaps(unit_value, sampling, dim);
    gym::restoreAlphaTest();
    gym::restoreBlending();
    gym::restoreLighting();

    if ( buf )
        OffScreen::releaseBuffer();
    glViewport(svp[0], svp[1], svp[2], svp[3]);
}


void PointDisp::createPixelmaps(float uv)
{
    if ( pixSize > pixAlloc_ )
    {
        CHECK_GL_ERROR("1 PointDisp::prepare");
        allocatePixelmap(pixSize);
        //fprintf(stderr, " new %i bitmap for %s\n", pixSize, name_str());
        CHECK_GL_ERROR("2 PointDisp::prepare");
    }
    
    pixOffset_ = -0.5f * pixSize;
    makePixelmaps(3*uv, 3);

#if POINTDISP_USES_PIXEL_BUFFERS
    if ( pbo_ == 0 )
        glGenBuffers(1, &pbo_);
    glBindBuffer(GL_PIXEL_PACK_BUFFER_ARB, pbo_);
    glBufferData(GL_PIXEL_PACK_BUFFER_ARB, 3*pixStride, pixels_, GL_STATIC_DRAW);
    glBindBuffer(GL_PIXEL_PACK_BUFFER_ARB, 0);
    CHECK_GL_ERROR("PointDisp::storePixelmap");
#endif
}

#endif


void PointDisp::setPixels(float ps, float uv, bool make_maps)
{
    float sw = size + width;
    perceptible = visible && ( uv*sw > 0.5 );
    
    sizeR = size * uv * ps * 0.5f;
    sizeX = std::max(size * uv, 0.25f);;
    widthX = std::max(width * uv, 0.25f);;
    //printf("widthX %6.3f sizeX %6.3f\n", widthX, sizeX);

#if POINTDISP_USES_PIXELMAPS
    unsigned S = static_cast<unsigned>(std::ceil(uv*sw));
    // make it a multiple of 2:
    pixSize = ( S + 2U ) & ~1U;
    //pixSize = S;

    if ( make_maps )
        createPixelmaps(uv);
#endif
}


#pragma mark - I/O


void PointDisp::read(Glossary& glos)
{
    glos.set(visible, "visible");
    
    // set 'color2' as a darker tone of 'color':
    if ( glos.set(color, "color") )
    {
        color2 = color.alpha_scaled(0.5f);
        colorS = color.inverted();
    }
    glos.set(color2, "color", 1, "back_color", 0);
    glos.set(coloring, "coloring");
    
    // if 'size' is specified, width is set accordingly:
    if ( glos.set(size, "size") )
        width = size / 4;
    else
        glos.set(size, "point_size");
    // harmless backward compatibility
    glos.set(size, "points");
    glos.set(shape, "points", 1);

    glos.set(width, "width");
    glos.set(style, "style");
    glos.set(shape, "shape");
    glos.set(symbol, "symbol");
    glos.set(colorS, "symbol", 1);
    
    if ( ! isprint(symbol) )
        symbol = 0;
    shape = tolower(shape);
    
#if POINTDISP_USES_PIXELMAPS
    releasePixelmap();
#endif
}


void PointDisp::write_values(std::ostream& os) const
{
    write_value(os, "visible", visible);
    if ( color2 != color.alpha_scaled(0.5f) )
        write_value(os, "color", color, color2);
    else
        write_value(os, "color", color);
    write_value(os, "coloring", coloring);
    write_value(os, "size", size);
    write_value(os, "width", width);
    write_value(os, "shape", shape);
    write_value(os, "style", style);
    write_value(os, "symbol", symbol, colorS);
}

