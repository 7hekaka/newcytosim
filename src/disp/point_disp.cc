// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include <cctype>
#include "point_disp.h"
#include "offscreen.h"
#include "glossary.h"
#include "gle.h"
#include "fg_stroke.h"
#include "gym_check.h"
#include "gym_view.h"
#include "gym_draw.h"
#include "gym_image.h"
#include "gym_flute.h"


void PointDisp::clearPixelmaps()
{
#if POINTDISP_USES_PIXELMAPS
    pixSize = 0;
    pixAlloc_ = 0;
    pixels_ = nullptr;
    texture_ = 0;
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
    scale    = o.scale;
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
    scale    = o.scale;
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
    scale  = 1;
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
            const float X = -52.35 * G;
            const float Y = ( islower(symbol) ? -35 : -50 ) * G;
            gym::color(colorS);
            fgStrokeCharacter(X, Y, G, 1, symbol, 3, 3);
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
    delete(pixels_);
    pixels_ = new uint8_t[16*pixSize*pixSize]();
}


void PointDisp::releasePixelmap()
{
    delete(pixels_);
    pixels_ = nullptr;
    pixAlloc_ = 0;
    if ( texture_ > 0 )
    {
        glDeleteTextures(1, &texture_);
        texture_ = 0;
    }
}


void PointDisp::drawPixelmap(float X, float Y, float Z, size_t inx) const
{
    CHECK_GL_ERROR("drawPixelmap0");
    float S = sizeR;
    float T = 0.25 * inx;
    float U = 0.25 + T;
    gym::ref_view();
    gym::color(0,1,0,1);
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, texture_);
    flute6* flu = gym::mapBufferV4T2(4);
    flu[0] = { X-S, Y+S, Z, 0, 0, T };
    flu[1] = { X-S, Y-S, Z, 0, 0, U };
    flu[2] = { X+S, Y+S, Z, 0, 1, T };
    flu[3] = { X+S, Y-S, Z, 0, 1, U };
    gym::unmapBufferV4T2();
    gym::drawTriangleStrip(0, 4);
    gym::clenupTexture();
    CHECK_GL_ERROR("drawPixelmap1");
}


/**
 `sampling` defines the level of oversampling used to improve the quality of bitmaps
 */
void PointDisp::makePixelmaps(float unit_value, unsigned sampling, unsigned dim)
{
    float S = size * unit_value * 0.5f;
    float W = width * unit_value * 0.5f;
    gym::transScale(dim*0.5, dim*0.5, 0, S);

    if ( width > 0 ) glLineWidth(W);

    for ( int i = 0; i < 3; ++i )
    {
        uint8_t * pix = pixels_ + i * pixSize * pixSize;
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
                strokeA(W);
                break;
            case 2:
                gym::color(color);
                strokeA(W);
                break;
        }
        if ( sampling > 1 )
        {
            uint8_t * tmp = new uint8_t[4*dim*dim];
            glReadPixels(0, 0, dim, dim, GL_RGBA, GL_UNSIGNED_BYTE, tmp);
            gym::downsampleRGBA(pix, pixSize, pixSize, tmp, sampling);
            gym::printPixels(stdout, tmp, dim, dim);
            delete[] tmp;
        }
        else
        {
            glReadPixels(0, 0, pixSize, pixSize, GL_RGBA, GL_UNSIGNED_BYTE, pix);
        }
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
    //gym::printCaps("P");
    makePixelmaps(unit_value, sampling, dim);
    gym::restoreAlphaTest();
    gym::restoreBlending();
    gym::restoreLighting();

    if ( buf )
        OffScreen::releaseBuffer();
    glViewport(svp[0], svp[1], svp[2], svp[3]);
}

static void setRGBA(size_t cnt, uint8_t ptr[], uint8_t R, uint8_t G, uint8_t B, uint8_t A)
{
    for ( size_t i = 0; i < cnt; ++i )
    {
        ptr[4*i+3] = R;
        ptr[4*i+2] = G;
        ptr[4*i+1] = B;
        ptr[4*i+0] = A;
    }
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
    
    makePixelmaps(3*uv, 3);

    if ( ! texture_ )
        glGenTextures(1, &texture_);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glBindTexture(GL_TEXTURE_2D, texture_);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    if ( 1 )
    {
        size_t S = pixSize*pixSize;
        setRGBA(S, pixels_, 100, 90, 128, 255);
        setRGBA(S, pixels_+4*S, 100, 200, 0, 255);
        setRGBA(S, pixels_+8*S, 0, 0, 255, 255);
        setRGBA(S, pixels_+12*S, 255, 0, 0, 255);
    }
    glTexImage2D(GL_TEXTURE_2D, 0, 4, pixSize, pixSize, 0, GL_RGBA, GL_UNSIGNED_INT_8_8_8_8, pixels_);
    //printf("%10s:texture %i x %i\n", name().c_str(), pixSize, pixSize);
    CHECK_GL_ERROR("PointDisp::createPixelmaps");
}

#endif


/// return smallest power of 2 that is greater or equal to `x`
static unsigned next_power(unsigned x)
{
    if ( x > 0 )
    {
        --x;
        x |= x >> 1;
        x |= x >> 2;
        x |= x >> 4;
        x |= x >> 8;
        x |= x >> 16;
    }
    return x+1;
}

/// arguments are ps = pixel_size; uv = unit_value
void PointDisp::setPixels(float ps, float uv, bool make_maps)
{
    float sw = size + width;
    // object is 'perceptible' if it covers more than half a pixel:
    perceptible = visible && ( uv*sw > 0.5 );
    
    sizeR = size * uv * ps * 0.5f;
    sizeX = std::max(size * uv, 0.25f);
    widthX = std::max(width * uv, 0.25f);
    //printf("widthX %6.3f sizeX %6.3f\n", widthX, sizeX);

#if POINTDISP_USES_PIXELMAPS
    // make it a power of 2:
    pixSize = next_power(std::ceil(uv*sw));
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

    glos.set(width, "width") || glos.set(width, "size", 1);
    glos.set(scale, "scale");
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
    write_value(os, "scale", scale);
    write_value(os, "shape", shape);
    write_value(os, "style", style);
    if ( isprint(symbol) )
        write_value(os, "symbol", symbol, colorS);
}

