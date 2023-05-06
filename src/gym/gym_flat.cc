// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

#include "gym_flat.h"
#include "gym_view.h"
#include "gym_draw.h"
#include "gym_flute.h"
#include "gym_check.h"
#include "gym_draw.h"
#include "gym_image.h"

/// name of OpenGL texture buffer
GLuint gym_font_texture_ = 0;

/// global function accessible from C
void drawPixels(unsigned W, unsigned H, float X, float Y, float S, const unsigned char* bytes, const float color[4])
{
    gym::color(color);
    gym::drawPixels(W, H, X, Y, S, bytes);
}

/// global function accessible from C
///\todo: we should unpack the whole font data only once!
void drawBitmap(unsigned W, unsigned H, float X, float Y, float S, const unsigned char* bits, const float color[4])
{
    gym::color(color);
    unsigned char pixels[W*H+8];
    unpackBitmap(pixels, W, H, bits, W);
    gym::drawPixels(W, H, X, Y, S, pixels);
}


void gym::printPixels(unsigned W, unsigned H, const unsigned char* pixels, unsigned lda)
{
    static const char LOOKUP[17] = " 123456789ABCDEF";

    for ( unsigned i = 0; i < H; ++i )
    {
        const unsigned char * row = pixels + i * lda;
        for ( unsigned j = 0; j < W; ++j )
            putchar(LOOKUP[row[j]&15]);
        putchar('\n');
    }
}


/** This is drawing `bits` by using a texture over [X, X+S*W]x[Y, Y+S*H].
 S is the pixel dimension */
void gym::drawPixels(unsigned W, unsigned H, float X, float Y, float S, const unsigned char* pixels)
{
    //printPixels(W, H, pixels, W);
    CHECK_GL_ERROR("drawPixels0");
    glEnable(GL_TEXTURE_2D);
    if ( ! gym_font_texture_ )
        glGenTextures(1, &gym_font_texture_);

    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glBindTexture(GL_TEXTURE_2D, gym_font_texture_);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_ALPHA, W, H, 0, GL_ALPHA, GL_UNSIGNED_BYTE, pixels);

    flute4* flu = gym::mapBufferV2T2(4);
    flu[0] = { X,     Y+S*H, 0, 0 };
    flu[1] = { X,     Y,     0, 1 };
    flu[2] = { X+S*W, Y+S*H, 1, 0 };
    flu[3] = { X+S*W, Y,     1, 1 };
    gym::unmapBufferV2T2();
    CHECK_GL_ERROR("drawPixels1");

    gym::drawTriangleStrip(0, 4);
    glDisableClientState(GL_TEXTURE_COORD_ARRAY);
    glDisable(GL_TEXTURE_2D);
    CHECK_GL_ERROR("drawPixels2");
}

void paintBitmap(unsigned W, unsigned H, float X0, float Y0, float S, const unsigned char* bits, const float col[4])
{
    gym::color(col);
    gym::paintBitmap(W, H, X0, Y0, S, bits);
}

/** This is drawing pixels of `bits`, line by line, using triangle strips */
void gym::paintBitmap(unsigned W, unsigned H, float X0, float Y0, float S, const unsigned char* bits)
{
    const unsigned Wb = ( W + 7 ) >> 3;
    unsigned n_bits = 0;
    for ( unsigned b = 0; b < H*Wb; ++b )
        n_bits += __builtin_popcount(bits[b]);
    if ( !n_bits )
        return;
    flute2* flu = gym::mapBufferV2(6*n_bits);
    flute2* ptr = flu;
    for ( unsigned i = 0; i < H; ++i )
    {
        float X = X0;
        float Y = Y0 + S * i, T = Y + S;
        const unsigned char* row = bits + i * Wb;
        unsigned char old = 0;
        for ( unsigned j = 0; j < Wb; ++j )
        {
            for ( int k = 7; k >= 0; --k )
            {
                unsigned char bit = ( row[j] >> k ) & 1;
                if ( bit != old )
                {
                    old = bit;
                    ptr[0] = { X, Y };
                    ptr[1] = { X, (bit?Y:T) };
                    ptr[2] = { X, T };
                    ptr += 3;
                }
                X += S;
            }
        }
        if ( old )
        {
            ptr[0] = { X, Y };
            ptr[1] = { X, T };
            ptr[2] = { X, T };
            ptr += 3;
        }
    }
    //printf("%4lu %4u\n", ptr-flu, 6*n_bits);
    gym::unmapBufferV2();
    gym::drawTriangleStrip(0, ptr-flu);
    CHECK_GL_ERROR("paintBitmap");
}


/** This is drawing squares of dimension WxH for every '1' in str[] */
void gym::paintSequence(float X, float Y, float W, float H, const char str[])
{
    size_t n = strlen(str);
    X -= n * W * 0.5;
    Y -= H * 0.5;
    float T = Y + H;
    flute2* flu = gym::mapBufferV2(3*n+2);
    flute2* ptr = flu;
    char d = '0';
    for ( char const* c = str; *c; ++c )
    {
        X += W;
        if ( *c != d )
        {
            d = *c;
            ptr[0] = { X, Y };
            ptr[1] = { X, (d=='0'?T:Y) };
            ptr[2] = { X, T };
            ptr += 3;
        }
    }
    if ( d != '0' )
    {
        ptr[0] = { X, Y };
        ptr[1] = { X, T };
        ptr += 2;
    }
    gym::unmapBufferV2();
    //assert_true( ptr-flu <= 3*n+2 );
    gym::drawTriangleStrip(0, ptr-flu);
    CHECK_GL_ERROR("paintSequence");
}

/**
 rectangle should be specified as [ left, bottom, right, top ]
 The rectangle will be drawn counter-clockwise
 */
void gym::drawRectangle(const int rec[4], float width)
{
    float L(rec[0]), B(rec[1]), R(rec[2]), T(rec[3]);
    flute2 * flu = gym::mapBufferV2(4);
    flu[0] = {L, B};
    flu[1] = {R, B};
    flu[2] = {R, T};
    flu[3] = {L, T};
    flu[4] = {L, B};
    gym::unmapBufferV2();
    gym::drawLineStrip(width, 0, 5);
}

void gym::drawRectangle(float L, float B, float R, float T, float Z, float width)
{
    flute3 * flu = gym::mapBufferV3(5);
    flu[0] = {L, B, Z};
    flu[1] = {R, B, Z};
    flu[2] = {R, T, Z};
    flu[3] = {L, T, Z};
    flu[4] = {L, B, Z};
    gym::unmapBufferV3();
    gym::drawLineStrip(width, 0, 5);
}

void gym::fillRectangle(float L, float B, float R, float T, float Z, const float col[4])
{
    gym::color(col);
    flute3 * flu = gym::mapBufferV3(4);
    flu[0] = {L, B, Z};
    flu[1] = {R, B, Z};
    flu[2] = {L, T, Z};
    flu[3] = {R, T, Z};
    gym::unmapBufferV3();
    gym::drawTriangleStrip(0, 4);
}


//-----------------------------------------------------------------------

int copy_parity(const int a, const int b)
{
    return a + (( std::abs(a) + b ) & 1 );
}

static void setRectangle(flute3*& flu, float L, float B, float R, float T, float Z)
{
    flu[0] = {L, B, Z};
    flu[1] = {L, B, Z};
    flu[2] = {L, T, Z};
    flu[3] = {R, B, Z};
    flu[4] = {R, T, Z};
    flu[5] = {R, T, Z};
    flu += 6;
}

void gym::drawTiledFloor(int R, float T, float Z)
{
    float H = T * 0.5;
    int Q = std::floor( double(R) * M_SQRT1_2 );
    
    int x = R;
    int RX = 2 * x - 3;
    int RY = 0;
    flute3 * flu = gym::mapBufferV3(24*Q*Q);
    flute3 * ptr = flu;
    
    for ( int y = 0; y <= x; ++y )
    {
        /*
         using the Midpoint circle algorithm
         https://en.wikipedia.org/wiki/Midpoint_circle_algorithm
        */
        if ( RY > RX )
        {
            RX += 4 * ( x - 1 );
            --x;
        }
        RY += 4 * y + 8;
        for ( int i = copy_parity(-x,y); i <= x; i+=2 )
        {
            float X = i * T;
            float Y = y * T;
            setRectangle(ptr,  X-H, Y-H, X+H, Y+H, Z);
            setRectangle(ptr, -X+H,-Y+H,-X-H,-Y-H, Z);
        }
        for ( int i = copy_parity(Q,y); i <= x; i+=2 )
        {
            float X = y * T;
            float Y = i * T;
            setRectangle(ptr, X-H, Y-H, X+H, Y+H, Z);
            setRectangle(ptr,-X+H,-Y+H,-X-H,-Y-H, Z);
        }
        for ( int i = copy_parity(Q,y); i <= x; i+=2 )
        {
            float X = y * T;
            float Y = i * T;
            setRectangle(ptr,-X-H, Y-H,-X+H, Y+H, Z);
            setRectangle(ptr, X+H,-Y+H, X-H,-Y-H, Z);
        }
    }
    gym::unmapBufferV3();
    gym::drawTriangleStrip(0, ptr-flu);
}


void gym::paintOctagon(float L, float B, float R, float T, const float col[4], float D)
{
    gym::color(col);
    flute2 * flu = gym::mapBufferV2(8);
    flu[0] = {R, B+D};
    flu[1] = {R, T-D};
    flu[2] = {R-D, B};
    flu[3] = {R-D, T};
    flu[4] = {L+D, B};
    flu[5] = {L+D, T};
    flu[6] = {L, B+D};
    flu[7] = {L, T-D};
    gym::unmapBufferV2();
    glDrawArrays(GL_TRIANGLE_STRIP, 0, 8);
}

void gym::drawOctagon(float L, float B, float R, float T, const float col[4], float D, float W)
{
    gym::color(col);
    flute2 * flu = gym::mapBufferV2(10);
    flu[0] = {L, B+D};
    flu[1] = {L+D, B};
    flu[2] = {R-D, B};
    flu[3] = {R, B+D};
    flu[4] = {R, T-D};
    flu[5] = {R-D, T};
    flu[6] = {L+D, T};
    flu[7] = {L, T-D};
    flu[8] = {L, B+D};
    gym::unmapBufferV2();
    gym::drawLineStrip(W, 0, 9);
}

