// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

#ifndef GYM_VIEW_H
#define GYM_VIEW_H

#include "opengl.h"
#include "gym_matrix.h"

namespace gym
{
    /// current modelview matrix
    extern GLfloat mvp_[16];
    
    /// modelview matrix of the current View
    extern GLfloat ref_[16];

    /// flag for Lighting effects
    extern GLboolean light_, alpha_;
    
#pragma mark -

    /// load current matrix to OpenGL
    inline void load() { glLoadMatrixf(mvp_); }
    
    /// load reference matrix
    inline void load_ref() { glLoadMatrixf(ref_); }

    /// copy reference to current
    inline void pull_ref() { gym::mat_copy(mvp_, ref_); }

#pragma mark - Modifying the reference view

    /// get current reference modelview matrix
    inline void get_view(float mat[16]) { gym::mat_copy(mat, ref_); }
    
    /// replace reference modelview matrix
    inline void set_view(const float mat[16]) { gym::mat_copy(ref_, mat); load_ref(); }
    
    /// make reference matrix current and load
    inline void ref_view() { pull_ref(); load_ref(); }
    
    /// change projection matrix
    void set_projection(GLfloat mat[16]);

#pragma mark - Set the current view

    /// set Identity transformation and load (reference view is not changed)
    inline void eye_view() { gym::mat_diagonal(mvp_, 1); load(); }
    
    /// center view on (X, Y, Z) and scale by S (reference view is not changed)
    inline void abs_view(float X, float Y, float Z, float S) { gym::mat_diagonal(mvp_, 1); gym::mat_transscale(mvp_, X, Y, Z, S); load(); }
    
    /// make one-to-one correspondance between pixel and model coordinates
    void one_view(int W, int H);

#pragma mark - Modifying the current view

    /// multiply current matrix by 'mat'
    inline void apply(const float mat[16]) { gym::mat_mul(mvp_, ref_, mat); load(); }
    
    /// translate current view
    inline void translate(float x, float y, float z) { gym::mat_translate(mvp_, x, y, z); load(); }

    /// scale current view
    inline void scale(float S) { gym::mat_scale(mvp_, S, S, S); load(); }
    
    /// scale current view in each direction separately
    inline void scale(float X, float Y, float Z) { gym::mat_scale(mvp_, X, Y, Z); load(); }
    
    /// rotate current view around axis (X, Y, Z) by angle defined by (C, S)
    inline void rotate(float X, float Y, float Z, float C, float S) { GLfloat T[16];  gym::mat_rotation(T, X, Y, Z, C, S); apply(T); }

    /// rotate current view around axis X by angle defined by (C, S)
    inline void rotateX(float C, float S) { gym::mat_rotateX(mvp_, C, S); load(); }
    
    /// rotate current view around axis Y by angle defined by (C, S)
    inline void rotateY(float C, float S) { gym::mat_rotateY(mvp_, C, S); load(); }
    
    /// rotate current view around axis Z by angle defined by (C, S)
    inline void rotateZ(float C, float S) { gym::mat_rotateZ(mvp_, C, S); load(); }

#pragma mark -
    
    /// translate by (X, Y, Z) and then scale by S
    inline void transScale(float X, float Y, float Z, float S) { gym::mat_copy(mvp_, ref_); gym::mat_transscale(mvp_, X, Y, Z, S); load(); }

    /// translate by pos; rotate to align X to Z, scale by rad
    void transAlignZX(float pos, float rad, float dir);
    /// translate by A; rotate to align X to Z, scale XY by rad and Z by B-A
    void stretchAlignZX(float A, float B, float rad);
    /// translate by A; rotate to align X to Z, scale XY by rad and Z by B-A
    void stretchAlignZY(float A, float B, float rad);

#pragma makr - deprecated features
    
    /// enable Lighting effects
    inline void enableLighting() { light_ = glIsEnabled(GL_LIGHTING); glEnable(GL_LIGHTING); }
    /// disable Lighting
    inline void disableLighting() { light_ = glIsEnabled(GL_LIGHTING); glDisable(GL_LIGHTING); }
    /// restore previous Lighting state
    inline void restoreLighting() { if ( light_ ) glEnable(GL_LIGHTING); else glDisable(GL_LIGHTING);  }
               
    /// enable Lighting effects
    inline void enableAlphaTest() { alpha_ = glIsEnabled(GL_ALPHA_TEST); glEnable(GL_ALPHA_TEST); }
    /// disable Lighting
    inline void disableAlphaTest() { alpha_ = glIsEnabled(GL_ALPHA_TEST); glDisable(GL_ALPHA_TEST); }
    /// restore previous Lighting state
    inline void restoreAlphaTest() { if ( alpha_ ) glEnable(GL_ALPHA_TEST); else glDisable(GL_ALPHA_TEST); }

    /// enable Line Stipple
    void enableLineStipple(short);
    
    /// disable Line Stipple
    void disableLineStipple();
    
#pragma mark - Clip Planes
    
    /// enable
    inline void enableClipPlane(unsigned i) { glEnable(GL_CLIP_PLANE0+i); }
    
    /// disable clip plane
    inline void disableClipPlane(unsigned i) { glDisable(GL_CLIP_PLANE0+i); }

    /// define cliping equation
    void setClipPlane(unsigned, double, double, double, double);

    inline void enableClipPlane(unsigned i, double X, double Y, double Z, double S)
    {
        setClipPlane(i, X, Y, Z, S);
        enableClipPlane(i);
    }

}

#endif
