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
#include "gym_flat.h"


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
    ulna_ = 0;
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


#pragma mark - I/O


void PointDisp::read(Glossary& glos)
{
    glos.set(visible, "visible");
    
    // set 'color2' as a darker tone of 'color':
    if ( glos.set(color, "color") )
    {
        color2 = color.alpha_scaled(0.1875f);
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


#pragma mark - Stroke

void PointDisp::strokeA(float w) const
{
    gle::disc();
    if ( w > 0 )
    {
        gym::color(color.darken(2.0));
        gle::thinRing();  /*gym::zoo_stroke(shape);*/ 
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
    gle::disc();
    
    // punch a transparent hole:
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
    float S = ulna_;
    float T = 0.25 * inx;
    float U = 0.25 + T;
    gym::ref_view();
    gym::color(1,1,1,1);
    gym::enableTexture(texture_);
    flute6* flu = gym::mapBufferV4T2(4);
    flu[0] = { X-S, Y+S, Z, 0, 0, T };
    flu[1] = { X-S, Y-S, Z, 0, 0, U };
    flu[2] = { X+S, Y+S, Z, 0, 1, T };
    flu[3] = { X+S, Y-S, Z, 0, 1, U };
    gym::unmapBufferV4T2();
    gym::drawTriangleStrip(0, 4);
    gym::cleanupTexture();
    CHECK_GL_ERROR("drawPixelmap1");
}


/**
 `sampling` defines the level of oversampling used to improve the quality of bitmaps
 */
void PointDisp::makePixelmaps(unsigned sampling, unsigned dim)
{
    //float col[4] = {1,1,1,1};
    for ( int i = 0; i < 3; ++i )
    {
        uint8_t * pix = pixels_ + i * pixSize * pixSize * 4;
        // we use a transparent background, because points will overlap
        gym::clearPixels(0,0,0,0);
        //gym::fillRectangle(0, 0, dim, dim, 0, col);
        switch ( i )
        {
            case 0:
                //gle::disc(); //gym::fillRectangle(0, 0, dim, dim, 0, col);
                gym::color(color2);
                strokeI();
                break;
            case 1:
                gym::color(color2);
                strokeA(widthX*sampling);
                break;
            case 2:
                gym::color(color);
                strokeA(widthX*sampling);
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
            /*
            glBindFramebuffer(GL_READ_FRAMEBUFFER, 0);
            glBindFramebuffer(GL_DRAW_FRAMEBUFFER, texture_);
            glBlitFramebuffer(0, 0, dim, dim, 0, 0, dim, dim, GL_COLOR_BUFFER_BIT, GL_NEAREST);
            glBindFramebuffer(GL_READ_FRAMEBUFFER, 0);
            glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 1);
             */
            glReadPixels(0, 0, pixSize, pixSize, GL_RGBA, GL_UNSIGNED_BYTE, pix);
        }
        if ( 0 )
        {
            //savePixelmap(pix, pixSize, i, name_str());
            printf("%s-%i %3.1f %3.1f %u\n", name().c_str(), i, size, width, pixSize);
            gym::printPixels(stdout, pix, pixSize, pixSize);
        }
        CHECK_GL_ERROR("5 PointDisp::makePixelmaps");
    }
}

/**
 `sampling` defines the level of oversampling used to improve the quality of bitmaps
 */
void PointDisp::makePixelmaps(unsigned sampling)
{
    CHECK_GL_ERROR("1 PointDisp::makePixelmaps");
    GLint svp[4];
    glGetIntegerv(GL_VIEWPORT, svp);
    
    unsigned dim = sampling * pixSize;
    unsigned buf = 0; //OffScreen::openBuffer(dim, dim, 0);
    if ( buf )
        gym::one_view(dim, dim);
    else
        gym::one_view(svp[2], svp[3]);
    gym::translate_scale(0.5*dim, 0.5*dim, 0, 0.5*dim*sizeX/pixSize);

    gym::disableLighting();
    gym::disableAlphaTest();
    gym::disableDepthTest();
    gym::disableBlending();
    //gym::printCaps("P");
    makePixelmaps(sampling, dim);
    gym::restoreBlending();
    gym::restoreDepthTest();
    gym::restoreAlphaTest();
    gym::restoreLighting();

    if ( buf )
        OffScreen::releaseBuffer();
    glViewport(svp[0], svp[1], svp[2], svp[3]);
}

static void setRGBA(size_t cnt, uint8_t ptr[], uint8_t R, uint8_t G, uint8_t B, uint8_t A)
{
    for ( size_t i = 0; i < cnt; ++i )
    {
        ptr[4*i+0] = R;
        ptr[4*i+1] = G;
        ptr[4*i+2] = B;
        ptr[4*i+3] = A;
    }
}
                           
void PointDisp::createPixelmaps()
{
    if ( pixSize > pixAlloc_ )
    {
        CHECK_GL_ERROR("1 PointDisp::createPixelmaps");
        allocatePixelmap(pixSize);
        //fprintf(stderr, " new %i bitmap for %s\n", pixSize, name_str());
        CHECK_GL_ERROR("2 PointDisp::createPixelmaps");
    }
    
    makePixelmaps(3);

    if ( 0 )
    {
        size_t S = pixSize * pixSize;
        setRGBA(S/2, pixels_, 0, 90, 255, 200);
        setRGBA(S, pixels_+4*S, 0, 200, 0, 255);
        setRGBA(S, pixels_+8*S, 0, 0, 255, 255);
        setRGBA(S, pixels_+12*S, 0, 255, 0, 255);
    }
    if ( ! texture_ )
        glGenTextures(1, &texture_);
    glBindTexture(GL_TEXTURE_2D, texture_);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, pixSize, 4*pixSize, 0, GL_RGBA, GL_UNSIGNED_BYTE, pixels_);
    //glGenerateMipmap(GL_TEXTURE_2D);
    //printf("%10s:texture %i x %i\n", name().c_str(), pixSize, pixSize);
    CHECK_GL_ERROR("3 PointDisp::createPixelmaps");
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
    perceptible = visible && ( uv*sw > 0.25 );
    
    sizeX = std::max(size * uv, 0.25f);
    widthX = std::max(width * uv, 0.25f);
    //printf("%s: sizeX %6.3f widthX %6.3f\n", name.c_str(), sizeX, widthX);

#if POINTDISP_USES_PIXELMAPS
    // make it a power of 2:
    if ( make_maps )
    {
        pixSize = next_power(std::ceil(uv*size*M_SQRT2));
        ulna_ = pixSize * ps * 0.5f;
        createPixelmaps();
    }
#endif
}
