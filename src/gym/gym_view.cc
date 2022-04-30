// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

#include "gym_view.h"
#include <cmath>

/// current modelview matrix
GLfloat gym::mvp_[16] = { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1 };

/// modelview representing the current view
GLfloat gym::ref_[16] = { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1 };


GLboolean gym::light_ = 0;
GLboolean gym::alpha_ = 0;


void gym::set_projection(GLfloat mat[16])
{
    glMatrixMode(GL_PROJECTION);
    glLoadMatrixf(mat);
    glMatrixMode(GL_MODELVIEW);
}

void gym::one_view(int W, int H)
{
    glMatrixMode(GL_PROJECTION);
    gym::mat_ortho(mvp_, 0, W, 0, H, 0, 1);
    load();
    
    glMatrixMode(GL_MODELVIEW);
    gym::mat_diagonal(mvp_, 1);
    load();
}


// rotate to align Z with X and translate to center 'P'
void gym::transAlignZX(float P, float R, float X)
{
    float Y = std::copysign(R, X);
    float mat[16] = {
        0, Y, 0, 0,
        0, 0, R, 0,
        Y, 0, 0, 0,
        P, 0, 0, 1};
    apply(mat);
}

// rotate to align Z with X and translate to center 'A'
void gym::stretchAlignZX(float A, float B, float R)
{
    float X = B - A;
    float Y = std::copysign(R, X);
    //warning! this matrix appears here transposed
    float mat[16] = {
        0, Y, 0, 0,
        0, 0, R, 0,
        X, 0, 0, 0,
        A, 0, 0, 1 };
    apply(mat);
}

// rotate to align Z with Y and translate to center 'A'
void gym::stretchAlignZY(float A, float B, float R)
{
    float X = B - A;
    float Y = std::copysign(R, X);
    //warning! this matrix appears here transposed
    float mat[16] = {
        0, 0, Y, 0,
        R, 0, 0, 0,
        0, X, 0, 0,
        0, A, 0, 1 };
    apply(mat);
}

//-----------------------------------------------------------------------
#pragma mark - Clip Planes

void gym::enableLineStipple(short pattern)
{
#ifdef GL_VERSION_2_1
    glLineStipple(1, pattern);
    glEnable(GL_LINE_STIPPLE);
#endif
}

void gym::disableLineStipple()
{
#ifdef GL_VERSION_2_1
    glDisable(GL_LINE_STIPPLE);
#endif
}

void gym::setClipPlane(unsigned glp, double X, double Y, double Z, double S)
{
    GLdouble eq[4] = { X, Y, Z, S };
    glClipPlane(GL_CLIP_PLANE0+glp, eq);
}

