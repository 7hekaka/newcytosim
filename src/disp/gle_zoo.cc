// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University

#include <cmath>
#include "opengl.h"
#include "gle_zoo.h"
#include "vector.h"
#include "gle_color.h"
#include "flute.h"

/*
 This is a set of 2D shape, which is unfinished
 */

static GLuint buffer_ = 0;

static GLsizei zoo_[16] = { 0 };


static void zoo_draw(int i)
{
    glBindBuffer(GL_ARRAY_BUFFER, buffer_);
    glVertexPointer(2, GL_FLOAT, 0, nullptr);
    glDrawArrays(GL_LINE_STRIP, 1+zoo_[i], zoo_[i+1]-zoo_[i]-1);
}

static void zoo_fill(int i)
{
    glBindBuffer(GL_ARRAY_BUFFER, buffer_);
    glVertexPointer(2, GL_FLOAT, 0, nullptr);
    glDrawArrays(GL_TRIANGLE_FAN, zoo_[i], zoo_[i+1]-zoo_[i]);
}

void gle::zoo_stroke(char c)
{
    return;
    switch ( c )
    {
        case 't': zoo_draw(0); break;
        case 'v': zoo_draw(1); break;
        case 'q': zoo_draw(2); break;
        case 'r': zoo_draw(3); break;
        case '+': zoo_draw(4); break;
        case 'p': zoo_draw(5); break;
        case 'h': zoo_draw(6); break;
        case 's': zoo_draw(7); break;
        default:
        case 'c': zoo_draw(8); break;
    }
}


void gle::zoo_paint(char c)
{
    return;
    switch ( c )
    {
        case 't': zoo_fill(0); break;
        case 'v': zoo_fill(1); break;
        case 'q': zoo_fill(2); break;
        case 'r': zoo_fill(3); break;
        case '+': zoo_fill(4); break;
        case 'p': zoo_fill(5); break;
        case 'h': zoo_fill(6); break;
        case 's': zoo_fill(7); break;
        default:
        case 'c': zoo_fill(8); break;
    }
}


void gle::zoo_init(flute2* flt, flute2* const ori)
{
    size_t j = 0;
    // triangle
    {
        zoo_[j++] = flt - ori;
        constexpr GLfloat B(-0.5f); //std::sqrt(3)/2;
        constexpr GLfloat T(0.8660254037844386f); //std::sqrt(3)/2;
        *flt++ = { 0, 0 };
        *flt++ = { 0, 1 };
        *flt++ = {-T, B };
        *flt++ = { T, B };
        *flt++ = { 0, 1 };
    }
    // inverted triangle
    {
        zoo_[j++] = flt - ori;
        constexpr GLfloat T(0.8660254037844386f); //std::sqrt(3)/2;
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
        const GLfloat R = 1.1f;
        const GLfloat C = 0.4f;
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
        const GLfloat A(M_PI * 0.1);
        const GLfloat B(M_PI * 0.3);
        const GLfloat R  = 1.3512958724134987f; //std::sqrt( 4 * M_PI / std::sqrt( 25 + 10 * std::sqrt(5)) );
        const GLfloat C1 = R * cosf(A), S1 = R * sinf(A);
        const GLfloat C3 = R * cosf(B), S3 = R * sinf(B);
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
        const GLfloat R = 1.0996361107912678f; //std::sqrt( 2 * M_PI / ( 3 * std::sqrt(3) ));
        const GLfloat H = R * 0.8660254037844386f; // sqrtf(3)/2;
        const GLfloat X = R * 0.5f;
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
        const GLfloat A(M_PI * 0.1);
        const GLfloat B(M_PI * 0.3);
        const GLfloat R  = 1.2f, H = -0.6f;
        const GLfloat C1 = R * cosf(A), S1 = R * sinf(A);
        const GLfloat C3 = R * cosf(B), S3 = R * sinf(B);
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
        GLfloat a(M_PI/6.0);
        *flt++ = { 1, 0 };
        for ( int u = 1; u < 12; ++u )
            *flt++ = { cosf(u*a),  sinf(u*a) };
        *flt++ = { 1, 0 };
    }
    zoo_[j] = flt - ori;
}

