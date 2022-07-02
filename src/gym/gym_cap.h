// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

#ifndef GYM_CAP_H
#define GYM_CAP_H

#include <string.h>
#include <algorithm>
#include "opengl.h"

namespace gym
{
    /// flag for OpenGL effects
    extern GLboolean depth_;
    extern GLboolean cull_;
    extern GLboolean blend_;
    
    /// flag for Lighting effects
    extern GLboolean light_;
    extern GLboolean alpha_;
    
#pragma mark - Features supported in OpenGL ES 2.0

    /// enable Blending effects
    inline void enableBlending() { blend_ = glIsEnabled(GL_BLEND); glEnable(GL_BLEND); }
    /// disable Blending
    inline void disableBlending() { blend_ = glIsEnabled(GL_BLEND); glDisable(GL_BLEND); }
    /// restore previous Blending state
    inline void restoreBlending() { if ( blend_ ) glEnable(GL_BLEND); else glDisable(GL_BLEND); }

    /// enable depth test
    inline void enableDepthTest() { depth_ = glIsEnabled(GL_DEPTH_TEST); glEnable(GL_DEPTH_TEST); }
    /// disable depth test
    inline void disableDepthTest() { depth_ = glIsEnabled(GL_DEPTH_TEST); glDisable(GL_DEPTH_TEST); }
    /// restore previous depth test state
    inline void restoreDepthTest() { if ( depth_ ) glEnable(GL_DEPTH_TEST); else glDisable(GL_DEPTH_TEST); }
    
    inline void openDepthMask() { glDepthMask(GL_TRUE); }
    inline void closeDepthMask() { glDepthMask(GL_FALSE); }
    
    /// enable  Cull Face
    inline void enableCullFace(GLenum face) { cull_ = glIsEnabled(GL_CULL_FACE); glEnable(GL_CULL_FACE); glCullFace(face); }
    /// disable  Cull Face
    inline void disableCullFace() { cull_ = glIsEnabled(GL_CULL_FACE); glDisable(GL_CULL_FACE); }
    /// restore previous Cull Face state
    inline void restoreCullFace() { if ( cull_ ) glEnable(GL_CULL_FACE); else glDisable(GL_CULL_FACE); }
    /// change Cull Face parameter
    inline void switchCullFace(GLenum face) { glCullFace(face); }
    
#pragma mark - deprecated features
    
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
    
    /// define clipping equation
    void setClipPlane(unsigned, double, double, double, double);
    
    inline void enableClipPlane(unsigned i, double X, double Y, double Z, double S)
    {
        setClipPlane(i, X, Y, Z, S);
        enableClipPlane(i);
    }

}

#endif
