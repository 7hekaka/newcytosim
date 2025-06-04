// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

#ifndef GYM_VIEW_H
#define GYM_VIEW_H

#include "opengl.h"
#include "gym_matrix.h"

namespace gym
{
    /// current modelview matrix
    extern GLfloat mvp_[16];
    
    /// modelview matrix of the reference View
    extern GLfloat ref_[16];
    
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
    inline void set_view(const float mat[16]) { gym::mat_copy(ref_, mat); gym::mat_copy(mvp_, mat); load_ref(); }
    
    /// replace reference modelview matrix by `mat` translated by (X, Y, Z)
    inline void set_view(const float mat[16], float X, float Y, float Z) { gym::mat_copytrans(ref_, mat, X, Y, Z); gym::mat_copy(mvp_, ref_); load_ref(); }
    
    /// replace reference modelview matrix by `mat` scaled by (S, S, S)
    inline void set_view(const float mat[16], float S) { gym::mat_copy(ref_, mat); gym::mat_scale(ref_, S, S, S); load_ref(); }

    /// make reference matrix current and loaded
    inline void ref_view() { pull_ref(); load_ref(); }
    
    /// change projection matrix
    void set_projection(GLfloat mat[16]);

#pragma mark - Set the current view

    /// set Identity transformation and load (reference view is not changed)
    inline void eye_view(float Z, float S) { gym::mat_diagonal(mvp_, S, -Z); load(); }
    
    /// center view on (X, Y, Z) and scale by S (reference view is not changed)
    inline void eye_view(float X, float Y, float Z, float S) { gym::mat_diagonal(mvp_, S); gym::mat_translate(mvp_, X/S, Y/S, Z/S); load(); }
    
    /// make one-to-one correspondance between pixel and model coordinates
    void one_view(int W, int H);

#pragma mark - Modifying the current view
    
    /// keep scale and translation, but remove rotation component
    inline void cancelRotation() { gym::mat_unrotate(mvp_, mvp_); }
    
    /// translate by (X, Y, Z), relative to reference
    inline void translate_ref(float X, float Y, float Z) { gym::mat_translate(mvp_, ref_, X, Y, Z); load(); }

    /// multiply current matrix by 'mat'
    inline void apply(const float mat[16]) { gym::mat_multiply(mvp_, ref_, mat); load(); }
    
    /// translate current view
    inline void translate(float x, float y, float z) { gym::mat_translate(mvp_, x, y, z); load(); }
    
    /// translate current view, but do not load
    inline void shift(float x, float y, float z) { gym::mat_translate(mvp_, x, y, z); }

    /// scale current view uniformly
    inline void scale(float S) { gym::mat_scale(mvp_, S); load(); }
    
    /// scale current view in each direction separately
    inline void scale(float X, float Y, float Z) { gym::mat_scale(mvp_, X, Y, Z); load(); }
    
    /// translate by (X, Y, Z) and then scale by S
    inline void translate_scale(float X, float Y, float Z, float S) { gym::mat_transscale(mvp_, X, Y, Z, S); load(); }

    /// rotate current view around axis (X, Y, Z) by angle defined by (C, S)
    inline void rotate(float X, float Y, float Z, float C, float S) { GLfloat T[16]; gym::mat_rotation(T, X, Y, Z, C, S); apply(T); }

    /// rotate current view around axis X by angle defined by (C, S)
    inline void rotateX(float C, float S) { gym::mat_rotateX(mvp_, C, S); load(); }
    
    /// rotate current view around axis Y by angle defined by (C, S)
    inline void rotateY(float C, float S) { gym::mat_rotateY(mvp_, C, S); load(); }
    
    /// rotate current view around axis Z by angle defined by (C, S)
    inline void rotateZ(float C, float S) { gym::mat_rotateZ(mvp_, C, S); load(); }

#pragma mark -
    
    /// translate by (X, Y, Z) and then scale by S
    inline void transScale(float X, float Y, float Z, float S) { gym::mat_transscale(mvp_, ref_, X, Y, Z, S); load(); }

    /// translate by pos; rotate to align X to Z, scale by rad in Z and trans in the other directions
    void transAlignZX(float pos, float rad, float trans);
    
    /// translate by A; rotate to align X to Z, scale XY by rad and Z by B-A
    void stretchAlignZX(float A, float B, float rad);
    
    /// translate by A; rotate to align X to Z, scale XY by rad and Z by B-A
    void stretchAlignZY(float A, float B, float rad);

    
    void print_view();
}

#endif
