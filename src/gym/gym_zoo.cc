// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University

#include <cmath>
#include "opengl.h"
#include "gym_zoo.h"
#include "vector.h"
#include "gym_color.h"
#include "gym_flute.h"
#include "gym_draw.h"
#include "flute.h"

/*
 This is a set of 2D shapes, which is unfinished
 */

static GLuint zoo_buffer_ = 0;

static GLsizei zoo_[16] = { 0 };

static inline unsigned zoo_index(char c)
{
    switch ( c )
    {
        case 't': return 0;
        case 'v': return 1;
        case 'q': return 2;
        case 'r': return 3;
        case '+': return 4;
        case 'p': return 5;
        case 'h': return 6;
        case 's': return 7;
        default:
        case 'c': return 8;
    }
}

void gym::zoo_stroke(char c, float w)
{
    unsigned i = zoo_index(c);
    gym::bindBufferV2(zoo_buffer_);
    gym::drawLineStrip(w, 1+zoo_[i], zoo_[i+1]-zoo_[i]-1);
}


void gym::zoo_paint(char c)
{
    unsigned i = zoo_index(c);
    gym::bindBufferV2(zoo_buffer_);
    gym::drawArrays(GL_TRIANGLE_FAN, zoo_[i], zoo_[i+1]-zoo_[i]);
}


size_t gym::zoo_init(flute2* flt, flute2* const ori)
{
    size_t j = 0;
    // triangle
    {
        zoo_[j++] = flt - ori;
        const float B(-0.5f); //std::sqrt(3)/2;
        const float T(0.8660254037844386f); //std::sqrt(3)/2;
        *flt++ = { 0, 0 };
        *flt++ = { 0, 1 };
        *flt++ = {-T, B };
        *flt++ = { T, B };
        *flt++ = { 0, 1 };
    }
    // inverted triangle
    {
        zoo_[j++] = flt - ori;
        const float T(0.8660254037844386f); //std::sqrt(3)/2;
        *flt++ = { 0, 0 };
        *flt++ = { 0, -1.f };
        *flt++ = { T,  0.5 };
        *flt++ = {-T,  0.5 };
        *flt++ = { 0, -1.f };
    }
    // square
    {
        zoo_[j++] = flt - ori;
        *flt++ = { 0, 0 };
        *flt++ = { 1,  1};
        *flt++ = {-1,  1};
        *flt++ = {-1, -1};
        *flt++ = { 1, -1};
        *flt++ = { 1,  1};
    }
    // rectangle
    {
        zoo_[j++] = flt - ori;
        *flt++ = { 0, 0 };
        *flt++ = { 1,  0.5};
        *flt++ = {-1,  0.5};
        *flt++ = {-1, -0.5};
        *flt++ = { 1, -0.5};
        *flt++ = { 1,  0.5};
    }
    // plus
    {
        zoo_[j++] = flt - ori;
        const float R = 1.1f;
        const float C = 0.4f;
        *flt++ = { 0, 0 };
        *flt++ = { R,  C};
        *flt++ = { C,  C};
        *flt++ = { C,  R};
        *flt++ = {-C,  R};
        *flt++ = {-C,  C};
        *flt++ = {-R,  C};
        *flt++ = {-R, -C};
        *flt++ = {-C, -C};
        *flt++ = {-C, -R};
        *flt++ = { C, -R};
        *flt++ = { C, -C};
        *flt++ = { R, -C};
        *flt++ = { R,  C};
    }
    // pentagon
    {
        zoo_[j++] = flt - ori;
        const float A(M_PI * 0.1);
        const float B(M_PI * 0.3);
        const float R = 1.3512958724134987f; //sqrt(4*M_PI/sqrt(25+10*std::sqrt(5)));
        const float C1 = R * cosf(A), S1 = R * sinf(A);
        const float C3 = R * cosf(B), S3 = R * sinf(B);
        *flt++ = {  0,  0};
        *flt++ = {  0,  R};
        *flt++ = {-C1,  S1};
        *flt++ = {-C3, -S3};
        *flt++ = { C3, -S3};
        *flt++ = { C1,  S1};
        *flt++ = {  0,  R};
    }
    // hexagon
    {
        zoo_[j++] = flt - ori;
        const float R = 1.0996361107912678f; //sqrt( 2 * M_PI / ( 3 * sqrt(3) ));
        const float H = R * 0.8660254037844386f; // sqrt(3)/2;
        const float X = R * 0.5f;
        *flt++ = { 0,  0};
        *flt++ = { R,  0};
        *flt++ = { X,  H};
        *flt++ = {-X,  H};
        *flt++ = {-R,  0};
        *flt++ = {-X, -H};
        *flt++ = { X, -H};
        *flt++ = { R,  0};
    }
    // star
    {
        zoo_[j++] = flt - ori;
        const float A(M_PI * 0.1);
        const float B(M_PI * 0.3);
        const float R  = 1.2f, H = -0.6f;
        const float C1 = R * cosf(A), S1 = R * sinf(A);
        const float C3 = R * cosf(B), S3 = R * sinf(B);
        *flt++ = {    0,     0};
        *flt++ = {    0,     R};
        *flt++ = { H*C3, -H*S3};
        *flt++ = {  -C1,    S1};
        *flt++ = { H*C1,  H*S1};
        *flt++ = {  -C3,   -S3};
        *flt++ = {    0,   H*R};
        *flt++ = {   C3,   -S3};
        *flt++ = {-H*C1,  H*S1};
        *flt++ = {   C1,    S1};
        *flt++ = {-H*C3, -H*S3};
        *flt++ = {    0,     R};
    }
    // nearly a circle
    {
        zoo_[j++] = flt - ori;
        float a(M_PI/6.0);
        *flt++ = { 1, 0 };
        for ( int u = 1; u < 12; ++u )
            *flt++ = { cosf(u*a),  sinf(u*a) };
        *flt++ = { 1, 0 };
    }
    zoo_[j] = flt - ori;
    return flt - ori;
}

