// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

#ifndef GYM_DRAW_H
#define GYM_DRAW_H

#include <string.h>
#include <algorithm>
#include "opengl.h"

namespace gym
{
    inline void drawPoints(float size, int off, int cnt)
    {
        if ( size > 0 )
        {
            glPointSize(size);
            glDrawArrays(GL_POINTS, off, cnt);
        }
    }
    
    inline void drawSquarePoints(float size, int off, int cnt)
    {
        if ( size > 0 )
        {
            glPointSize(size);
            glDisable(GL_POINT_SMOOTH);
            glDrawArrays(GL_POINTS, off, cnt);
            glEnable(GL_POINT_SMOOTH);
        }
    }
    
    inline void drawLines(float width, int off, int cnt)
    {
        if ( width > 0 )
        {
            glLineWidth(width);
            glDrawArrays(GL_LINES, off, cnt);
        }
    }
    
    inline void drawLineStrip(float width, int off, int cnt)
    {
        if ( width > 0 )
        {
            glLineWidth(width);
            glDrawArrays(GL_LINE_STRIP, off, cnt);
        }
    }
    
    inline void drawLineStrip(int off, int cnt)
    {
        glDrawArrays(GL_LINE_STRIP, off, cnt);
    }

    inline void drawTriangles(int off, int cnt)
    {
        glDrawArrays(GL_TRIANGLES, off, cnt);
    }

    inline void drawTriangleStrip(int off, int cnt)
    {
        glDrawArrays(GL_TRIANGLE_STRIP, off, cnt);
    }
    
    inline void drawArrays(GLenum mode, int first, int cnt)
    {
        glDrawArrays(mode, first, cnt);
    }
    
    
#pragma mark - Functions to set the Current Color
    
    /// make RGBA color current
    void setColor(const float col[]);

    /// make RGBA color current
    inline void color(const float col[]) { glColor4fv(col); }

    /// make RGB color current
    inline void color(float R, float G, float B) { glColor3f(R,G,B); }

    /// make RGBA color current
    inline void color(float R, float G, float B, float A) { glColor4f(R,G,B,A); }

    /// return value clamped to [0, 1]
    static inline float clamp(float s) { return std::max(float(0), std::min(s, float(1))); }

    inline static void no_emission(GLenum face)
    {
        /*
        float blk[4] = { 0, 0, 0, 1 };
        glMaterialfv(face, GL_EMISSION, blk);
         */
    }
    
    /// set FRONT and BACK material property for lighting
    inline void set_emission(const float col[])
    {
        float blk[4] = { 0, 0, 0, 1 };
        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, blk);
        glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, col);
    }
    
    /// set FRONT material property for lighting
    inline void color_front(const float col[])
    {
        glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, col);
        no_emission(GL_FRONT);
    }
    
    /// set FRONT material property for lighting, and current color
    inline void color_load(const float col[])
    {
        glColor4fv(col);
        glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, col);
        no_emission(GL_FRONT);
    }
    
    /// set front OpenGL color, with `a` as alpha component
    inline void color_front(const float col[], float a)
    {
        float mat[4] = { col[0], col[1], col[2], clamp(a) };
        glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, mat);
        no_emission(GL_FRONT);
    }
    
    inline void color_front(float R, float G, float B, float A=1.0)
    {
        float col[4] = { clamp(R), clamp(G), clamp(B), clamp(A) };
        glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, col);
        no_emission(GL_FRONT);
    }
    
    /// set BACK material property for lighting
    inline void color_back(const float col[])
    {
        //glMaterialfv(GL_BACK, GL_AMBIENT_AND_DIFFUSE, col_);
        //float blk[4] = { 0, 0, 0, 1 };
        glMaterialfv(GL_BACK, GL_AMBIENT_AND_DIFFUSE, col);
        no_emission(GL_BACK);
    }
    
    /// set BACK material property for lighting
    inline void color_back(const float col[], float a)
    {
        float mat[4] = { col[0], col[1], col[2], clamp(a) };
        glMaterialfv(GL_BACK, GL_AMBIENT_AND_DIFFUSE, mat);
        no_emission(GL_BACK);
    }

    inline void color_back(float R, float G, float B, float A)
    {
        float col[4] = { clamp(R), clamp(G), clamp(B), clamp(A) };
        glMaterialfv(GL_BACK, GL_AMBIENT_AND_DIFFUSE, col);
        no_emission(GL_BACK);
    }

    /// set FRONT and BACK material property for lighting
    inline void color_both(const float col[])
    {
#if 0
        float blk[4] = { 0, 0, 0, 0 };
        glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, col_);
        glMaterialfv(GL_BACK, GL_AMBIENT, col_);
        glMaterialfv(GL_BACK, GL_DIFFUSE, blk);
#endif
        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, col);
        no_emission(GL_FRONT_AND_BACK);
    }
    
    /// set FRONT and BACK material property for lighting
    inline void color_both(const float col[], float a)
    {
        float blk[4] = { 0, 0, 0, 1 };
        float mat[4] = { col[0], col[1], col[2], clamp(a) };
        glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, mat);
        glMaterialfv(GL_BACK, GL_AMBIENT, mat);
        glMaterialfv(GL_BACK, GL_DIFFUSE, blk);
        no_emission(GL_FRONT_AND_BACK);
    }
    
    /// set FRONT and BACK material property for lighting
    inline void color_both(float R, float G, float B, float A)
    {
        float mat[4] = { R, G, B, A };
        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, mat);
        no_emission(GL_FRONT_AND_BACK);
    }

#pragma mark -

    inline void clearPixels() { glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT); }
    inline void clearPixels(const float c[4]) { glClearColor(c[0],c[1],c[2],c[3]); clearPixels(); }
    inline void clearPixels(float R, float G, float B, float A) { glClearColor(R,G,B,A); clearPixels(); }

    inline void clearStencil(GLint x) { glClearStencil(x); glClear(GL_STENCIL_BUFFER_BIT); }

    /// display back facing triangles followed by front facing ones
    void dualPass(void primitive());
}

#endif
