// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University

#include <cmath>
#include "opengl.h"
#include "gle_zoo.h"

void gle::triangleS()
{
    constexpr GLfloat H = 0.8660254037844386f; //std::sqrt(3)/2;
    constexpr GLfloat pts[] = { 0, 1, 0, -H, -.5f, 0, H, -.5f, 0 };
    constexpr GLfloat dir[] = { 0, 0, 1, 0, 0, 1, 0, 0, 1 };
    
    glEnableClientState(GL_NORMAL_ARRAY);
    glVertexPointer(3, GL_FLOAT, 0, pts);
    glNormalPointer(GL_FLOAT, 0, dir);
    glDrawArrays(GL_TRIANGLES, 0, 1);
    glEnableClientState(GL_NORMAL_ARRAY);
}

void gle::triangleL()
{
    glBegin(GL_LINE_LOOP);
    glNormal3f(0, 0, 1);
    const GLfloat H = 0.8660254037844386f; //std::sqrt(3)/2;
    glVertex2f( 0,  1.0);
    glVertex2f(-H, -0.5);
    glVertex2f( H, -0.5);
    glEnd();
}

//-----------------------------------------------------------------------

void gle::nablaS()
{
    glBegin(GL_TRIANGLES);
    glNormal3f(0, 0, 1);
    const GLfloat H = 0.8660254037844386f; //std::sqrt(3)/2;
    glVertex2f( 0, -1.0);
    glVertex2f( H,  0.5);
    glVertex2f(-H,  0.5);
    glEnd();
}

void gle::nablaL()
{
    glBegin(GL_LINE_LOOP);
    glNormal3f(0, 0, 1);
    const GLfloat H = 0.8660254037844386f; //std::sqrt(3)/2;
    glVertex2f( 0, -1.0);
    glVertex2f( H,  0.5);
    glVertex2f(-H,  0.5);
    glEnd();
}

//-----------------------------------------------------------------------
void gle::squareS()
{
    glBegin(GL_TRIANGLE_FAN);
    glNormal3f(0, 0, 1);
    glVertex2f( 1,  1);
    glVertex2f(-1,  1);
    glVertex2f(-1, -1);
    glVertex2f( 1, -1);
    glEnd();
}

void gle::squareL()
{
    glBegin(GL_LINE_LOOP);
    glNormal3f(0, 0, 1);
    glVertex2f( 1,  1);
    glVertex2f(-1,  1);
    glVertex2f(-1, -1);
    glVertex2f( 1, -1);
    glEnd();
}

//-----------------------------------------------------------------------
void gle::rectangleS()
{
    glBegin(GL_TRIANGLE_FAN);
    glNormal3f(0, 0, 1);
    glVertex2f( 1,  0.5);
    glVertex2f(-1,  0.5);
    glVertex2f(-1, -0.5);
    glVertex2f( 1, -0.5);
    glEnd();
}

void gle::rectangleL()
{
    glBegin(GL_LINE_LOOP);
    glNormal3f(0, 0, 1);
    glVertex2f( 1,  0.5);
    glVertex2f(-1,  0.5);
    glVertex2f(-1, -0.5);
    glVertex2f( 1, -0.5);
    glEnd();
}

//-----------------------------------------------------------------------

void gle::plusS()
{
    const GLfloat R = 1.1f;
    const GLfloat C = 0.4f;
    
    glBegin(GL_TRIANGLE_FAN);
    glNormal3f(0, 0, 1);
    glVertex2f( R,  C);
    glVertex2f(-R,  C);
    glVertex2f(-R, -C);
    glVertex2f( R, -C);
    glEnd();
    
    glBegin(GL_TRIANGLE_FAN);
    glNormal3f(0, 0, 1);
    glVertex2f( C,  R);
    glVertex2f(-C,  R);
    glVertex2f(-C,  C);
    glVertex2f( C,  C);
    glEnd();
    
    glBegin(GL_TRIANGLE_FAN);
    glNormal3f(0, 0, 1);
    glVertex2f( C, -C);
    glVertex2f(-C, -C);
    glVertex2f(-C, -R);
    glVertex2f( C, -R);
    glEnd();
}

void gle::plusL()
{
    const GLfloat R = 1.2f;
    const GLfloat C = 0.6f;
    
    glBegin(GL_LINE_LOOP);
    glNormal3f(0, 0, 1);
    glVertex2f( C,  R);
    glVertex2f(-C,  R);
    glVertex2f(-C,  C);
    glVertex2f(-R,  C);
    glVertex2f(-R, -C);
    glVertex2f(-C, -C);
    glVertex2f(-C, -R);
    glVertex2f( C, -R);
    glVertex2f( C, -C);
    glVertex2f( R, -C);
    glVertex2f( R,  C);
    glVertex2f( C,  C);
    glEnd();
}
//-----------------------------------------------------------------------
/// draw pentagon that has the same surface as a disc of radius 1.
void gle::pentagonS()
{
    const GLfloat A = GLfloat(M_PI * 0.1);
    const GLfloat B = GLfloat(M_PI * 0.3);
    glBegin(GL_TRIANGLE_FAN);
    glNormal3f(0, 0, 1);
    glVertex2f(0, 0);
    const GLfloat R  = 1.3512958724134987f; //std::sqrt( 4 * M_PI / std::sqrt( 25 + 10 * std::sqrt(5)) );
    const GLfloat C1 = R * cosf(A), S1 = R * sinf(A);
    const GLfloat C3 = R * cosf(B), S3 = R * sinf(B);
    
    glVertex2f(  0,  R);
    glVertex2f(-C1,  S1);
    glVertex2f(-C3, -S3);
    glVertex2f( C3, -S3);
    glVertex2f( C1,  S1);
    glVertex2f(  0,  R);
    glEnd();
}

void gle::pentagonL()
{
    const GLfloat A = GLfloat(M_PI * 0.1);
    const GLfloat B = GLfloat(M_PI * 0.3);
    glBegin(GL_LINE_LOOP);
    glNormal3f(0, 0, 1);
    const GLfloat R  = 1.3512958724134987f; //std::sqrt( 4 * M_PI / std::sqrt( 25 + 10 * std::sqrt(5)) );
    const GLfloat C1 = R * cosf(A), S1 = R * sinf(A);
    const GLfloat C3 = R * cosf(B), S3 = R * sinf(B);
    
    glVertex2f(  0,  R);
    glVertex2f(-C1,  S1);
    glVertex2f(-C3, -S3);
    glVertex2f( C3, -S3);
    glVertex2f( C1,  S1);
    glEnd();
}

//-----------------------------------------------------------------------
/// draw hexagon that has the same surface as a disc of radius 1.
void gle::hexagonS()
{
    glBegin(GL_TRIANGLE_FAN);
    glNormal3f(0, 0, 1);
    glVertex2f(0, 0);
    const GLfloat R = 1.0996361107912678f; //std::sqrt( 2 * M_PI / ( 3 * std::sqrt(3) ));
    const GLfloat H = R * 0.8660254037844386f; // sqrtf(3)/2;
    const GLfloat X = R * 0.5f;
    glVertex2f( R,  0);
    glVertex2f( X,  H);
    glVertex2f(-X,  H);
    glVertex2f(-R,  0);
    glVertex2f(-X, -H);
    glVertex2f( X, -H);
    glVertex2f( R,  0);
    glEnd();
}

void gle::hexagonL()
{
    glBegin(GL_LINE_LOOP);
    glNormal3f(0, 0, 1);
    const GLfloat R = 1.0996361107912678f; //std::sqrt( 2 * M_PI / ( 3 * std::sqrt(3) ));
    const GLfloat H = R * 0.8660254037844386f; // sqrtf(3)/2;
    const GLfloat X = R * 0.5f;
    glVertex2f( R,  0);
    glVertex2f( X,  H);
    glVertex2f(-X,  H);
    glVertex2f(-R,  0);
    glVertex2f(-X, -H);
    glVertex2f( X, -H);
    glEnd();
}


//-----------------------------------------------------------------------

void gle::starS()
{
    const GLfloat A = GLfloat(M_PI * 0.1);
    const GLfloat B = GLfloat(M_PI * 0.3);
    const GLfloat R  = 1.2f, H = -0.6f;
    const GLfloat C1 = R * cosf(A), S1 = R * sinf(A);
    const GLfloat C3 = R * cosf(B), S3 = R * sinf(B);

    glBegin(GL_TRIANGLE_FAN);
    glNormal3f(0, 0, 1);
    glVertex2f(0, 0);
    glVertex2f(    0,     R);
    glVertex2f( H*C3, -H*S3);
    glVertex2f(  -C1,    S1);
    glVertex2f( H*C1,  H*S1);
    glVertex2f(  -C3,   -S3);
    glVertex2f(    0,   H*R);
    glVertex2f(   C3,   -S3);
    glVertex2f(-H*C1,  H*S1);
    glVertex2f(   C1,    S1);
    glVertex2f(-H*C3, -H*S3);
    glVertex2f(    0,     R);
    glEnd();
}

void gle::starL()
{
    const GLfloat A = GLfloat(M_PI * 0.1);
    const GLfloat B = GLfloat(M_PI * 0.3);
            const GLfloat R  = 1.2f, H = -0.6f;
    const GLfloat C1 = R * cosf(A), S1 = R * sinf(A);
    const GLfloat C3 = R * cosf(B), S3 = R * sinf(B);
    glBegin(GL_LINE_LOOP);
    glNormal3f(0, 0, 1);
    glVertex2f(    0,     R);
    glVertex2f( H*C3, -H*S3);
    glVertex2f(  -C1,    S1);
    glVertex2f( H*C1,  H*S1);
    glVertex2f(  -C3,   -S3);
    glVertex2f(    0,   H*R);
    glVertex2f(   C3,   -S3);
    glVertex2f(-H*C1,  H*S1);
    glVertex2f(   C1,    S1);
    glVertex2f(-H*C3, -H*S3);
    glEnd();
}

