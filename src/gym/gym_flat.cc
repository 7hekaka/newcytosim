// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

#include "gym_flat.h"
#include "gym_view.h"
#include "gym_draw.h"
#include "gym_flute.h"
#include "gym_check.h"
#include "gym_draw.h"
#include "gym_image.h"


/// global function accessible from C
void drawBitmap(unsigned W, unsigned H, float X, float Y, float S, const unsigned char* bits, const float color[4])
{
    gym::color(color);
    gym::drawBitmap(W, H, X, Y, S, bits);
}

///\todo: we should unpack the whole font data only once!
/** This is drawing `bits` by using a texture over a square */
void gym::drawPixels(unsigned W, unsigned H, float X, float Y, float S, const unsigned char* pixels)
{
    static GLuint tex = 0;
    if ( ! tex ) glGenTextures(1, &tex);

    glEnable(GL_TEXTURE_2D);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glBindTexture(GL_TEXTURE_2D, tex);

    glTexImage2D(GL_TEXTURE_2D, 0, GL_ALPHA, W, H, 0, GL_ALPHA, GL_UNSIGNED_BYTE, pixels);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

    flute4* flu = gym::mapBufferV2T2(4);
    flu[0] = { X,     Y+S*H, 1, 1 };
    flu[1] = { X,     Y,     1, 0 };
    flu[2] = { X+S*W, Y+S*H, 0, 1 };
    flu[3] = { X+S*W, Y,     0, 0 };
    gym::unmapBufferV2T2();
    
    gym::drawTriangleStrip(0, 4);
    glDisableClientState(GL_TEXTURE_COORD_ARRAY);
    //glDeleteTextures(1, &tex);
    glDisable(GL_TEXTURE_2D);
    CHECK_GL_ERROR("drawBitmap");
}


///\todo: we should unpack the whole font data only once!
/** This is drawing `bits` by using a texture over a square */
void gym::drawBitmap(unsigned W, unsigned H, float X, float Y, float S, const unsigned char* bits)
{
    static GLuint tex = 0;
    if ( ! tex ) glGenTextures(1, &tex);
    
    unsigned char pixels[W*H+8];
    gym::unpackBitmap(pixels, W, H, bits);
    drawPixels(W, H, X, Y, S, pixels);
}


/** This is drawing every pixels of `bits` with triangle strips */
void gym::paintBitmap(unsigned W, unsigned H, float X0, float Y0, float S, const unsigned char* bytes)
{
    const size_t Wb = ( ( W + 7 ) & ~7 ) / 8;
    for ( size_t i = 0; i < H; ++i )
    {
        float X = X0;
        float Y = Y0 + S * i, T = Y + S;
        const unsigned char* row = bytes + i * Wb;
        flute2* flu = gym::mapBufferV2(24*Wb+2);
        flute2* ptr = flu;
        unsigned char old = 0;
        for ( size_t j = 0; j < Wb; ++j )
        {
            for ( size_t k = 0; k < 8; ++k )
            {
                X += S;
                unsigned char bit = !( ( row[j] >> (7-k) ) & 1 );
                if ( bit != old )
                {
                    old = bit;
                    ptr[0] = { X, T };
                    ptr[1] = { X, (bit?Y:T) };
                    ptr[2] = { X, Y };
                    ptr += 3;
                }
            }
        }
        if ( old )
        {
            ptr[0] = { X, T };
            ptr[1] = { X, Y };
            ptr += 2;
        }
        gym::unmapBufferV2();
        gym::drawTriangleStrip(0, ptr-flu);
    }
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
            ptr[0] = { X, T };
            ptr[1] = { X, (d=='0'?Y:T) };
            ptr[2] = { X, Y };
            ptr += 3;
        }
    }
    if ( d != '0' )
    {
        ptr[0] = { X, T };
        ptr[1] = { X, Y };
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
    flu[3] = {L, T, Z};
    flu[2] = {R, T, Z};
    gym::unmapBufferV3();
    gym::drawTriangleStrip(0, 4);
}


//-----------------------------------------------------------------------

int copy_parity(const int a, const int b)
{
    return a + (( std::abs(a) + b ) & 1 );
}

void gym::drawTiledFloor(int R, float T, float Z, const float col[4], const float back[4])
{
    float H = T * 0.5;
    int Q = std::floor( double(R) * M_SQRT1_2 );
    
    if ( back && back[3] > 0 )
    {
        float U = R * T;
        fillRectangle(-U, -U, U, U, Z, back);
    }
    
    int x = R;
    int RX = 2 * x - 3;
    int RY = 0;
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
            fillRectangle( X-H, Y-H, X+H, Y+H, Z, col);
            fillRectangle(-X+H,-Y+H,-X-H,-Y-H, Z, col);
        }
        for ( int i = copy_parity(Q,y); i <= x; i+=2 )
        {
            float X = y * T;
            float Y = i * T;
            fillRectangle( X-H, Y-H, X+H, Y+H, Z, col);
            fillRectangle(-X+H,-Y+H,-X-H,-Y-H, Z, col);
        }
        for ( int i = copy_parity(Q,y); i <= x; i+=2 )
        {
            float X = y * T;
            float Y = i * T;
            fillRectangle(-X-H, Y-H,-X+H, Y+H, Z, col);
            fillRectangle( X+H,-Y+H, X-H,-Y-H, Z, col);
        }
    }
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

