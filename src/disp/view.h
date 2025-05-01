// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University
#ifndef VIEW_H
#define VIEW_H

#include "opengl.h"
#include "vector1.h"
#include "vector2.h"
#include "view_prop.h"


/// Handles the viewing angle, projection and other aspects of an OpenGL display
/**
 ViewProp does not depend on the window-system (GLUT),
 but only on the rendering engine (OpenGL)
 */
class View : public ViewProp
{
private:
    
    /// modelview matrix
    mutable float modelview_[16];
    
    /// projection matrix
    mutable float projection_[16];
    
    /// dimension of camera window if `zoom=1`
    mutable float viewScale[2];
    
    /// translation between origin and camera plane
    mutable float eyeDistance;
    
    /// OpenGL viewport
    GLint viewport_[4];

    /// window number in GLUT
    int window_;
    
    /// flag for Region of Interest
    bool hasROI;
    
    /// Region of interest
    Vector3 mROI[2];
    
    /// text displayed near top right corner of window
    std::string top_message;
    
    /// text displayed near bottom left corner of window
    std::string full_label;
    
    /// text display in upper-right corner to provide user feedback
    mutable std::string flash;
    
    /// time beyond which `flash` is not displayed
    double flash_end;

    /// display callback
    int (*drawCallback)(View&);
    
    /// display callback
    void (*drawMagFunc)(View&);

    /// set OpenGL Fog, with mode (GL_EXP or GL_LINEAR), intensity and color
    void setFog(GLint mode, float density, gym_color) const;

public:
    
    /// constructor
    explicit View(const std::string& n, int depth);
    
    /// destructor
    ~View();
    
    /// return window-id
    int window() const { return window_; }
    
    /// set window-id
    void window(int w) { window_ = w; }
    
    /// adjust viewport size to window_size[]
    void resize(int);

    /// handle window resize events
    void reshape(int, int);
    
    /// reset viewport
    void loadViewport() const;
    
    /// reset viewport
    void setViewport(int, int, unsigned, unsigned) const;
    
    /// height of display area in pixels
    void copyViewport(GLint vp[4]) const { for(int i=0; i<4; ++i) vp[i] = viewport_[i]; }

    /// adjust parameters of projections, given window size
    void adjust(int, int) const;

    /// set OpenGL Projection and ModelView matrices
    void adjust() const { adjust(viewport_[2], viewport_[3]); }

    /// upload OpenGL Projection matrix
    void setProjection() const;
    
    /// set OpenGL Projection matrix
    void setPickProjection(float X, float Y, float W, float H) const;

    /// upload Model-View matrix
    void setModelView() const;
    
    /// upload Projection and ModelView matrices to OpenGL
    void loadView() const;

    /// adjust view to only show a slice of the world
    void sliceView(int) const;
    
    /// reset the view (no-rotation, zoom=1), and enable auto_scale
    void reset();
    
    //---------------------------------------------------------------------------

    /// width of display area in pixels
    int width() const { return viewport_[2]; }
    
    /// height of display area in pixels
    int height() const { return viewport_[3]; }
    
    /// distance between origin and camera plate
    float depth() const { return eyeDistance; }

    /// size of pixel in drawing units
    float pixelSize() const;
    
    /// half visible range in X or Y
    float range(int i) const { return 0.5f * viewScale[i] / zoom; }
    
    /// return first direction of view parallel to display screen
    Vector3 firstAxis() const;

    /// return direction of view that is orthogonal to display screen
    Vector3 depthAxis() const;
    
    /// calculate projection of {X, Y, Z} onto screen
    void project(float& H, float& V, const real XYZ[3]) const;

    /// transform window coordinates to 3D world-coordinates
    Vector3 unproject(float x, float y, float z) const;

    //---------------------------------------------------------------------------
    
    /// set display callback
    void setDisplayFunc(int (*f)(View&)) { if (f) drawCallback = f; }

    /// set second display callback
    void setDrawMagFunc(void (*f)(View&)) { if (f) drawMagFunc = f; }

    /// clear pixels and set clipping planes and fog parameters
    void openDisplay() const;
    
    /// unset clipping planes and fog parameters, display axes and scale bar
    void closeDisplay() const;

    /// display frame-per-seconds
    void drawFPS() const;

    /// draw bottom-left label
    void drawTextAndScale(std::string const&) const;
    
    /// display scale bar, info text, etc.
    void drawInteractiveFeatures() const;
    
    /// call drawCallback
    int callDraw() { return drawCallback(*this); }

    /// clear color and depth buffer
    void clearPixels();

    //---------------------------------------------------------------------------
    
    /// init OpenGL parameters
    void initGL() const;

    /// toggle depth clamp GL capability
    void toggleDepthClamp();
    
    /// set OpenGL Lights for lighting effects
    void setLights() const;
    
    /// set OpenGL Lights for lighting effects
    void setLightsEye() const;

    /// set text displayed near bottom-left corner of window
    void setLabel(std::string const& arg) { full_label = label + " " + arg; }
    
    /// set message displayed in center of window
    void setMemo(std::string const& arg) { memo = arg; };

    /// set text displayed near top left corner of window
    void setMessage(std::string const& arg) { top_message = arg; }

    /// set text displayed near top left corner of window
    void setSubtitle(std::string const& arg) { subtitle = arg; }

    /// set OpenGL Fog, with mode (GL_EXP or GL_LINEAR), intensity and color
    void enableFog(GLint mode, float param, gym_color);
    void enableFog(GLint mode, float param);
    
    /// establish all enabled clipping planes
    void setClipping() const;
    
    /// disable clipping planes
    void endClipping() const;
    
    /// set equations for a clipping plane, and enable it in View
    void enableClipPlane(int, real X, real Y, real Z, real S, bool absolute=true);
    
    /// disable clipping plane in View
    void disableClipPlane(int);
    
    /// return enable/disable state
    int hasClipPlane(int) const;
    
    //---------------------------------------------------------------------------
    
    /// position 'pos' in the center of the display
    void move_to(Vector3 const&);
    
    /// set additional translation of focal point
    void move_shift(Vector3 const& vec) { focus_shift = vec;}

    /// translate view
    void move_by(Vector3 const& trans) { move_to( focus - trans ); }

    //---------------------------------------------------------------------------
    
    /// set rotation to given Quaternion
    void rotate_to(const Quaternion<real>&);
    
    /// rotate to have `dir` aligned with the X-axis
    void align_with(const real dir[3], bool can_flip, real angle_scale);

    /// rotate view
    void rotate_by(const Quaternion<real>& Q) { rotate_to( rotation * Q ); }
    
    //---------------------------------------------------------------------------

    /// adjust the size of the visible area with zoom==1.
    void set_scale(float);
    
    /// set absolute zoom (zoom <- z)
    void zoom_to(float z);
    
    /// increase zoom (zoom <- zoom * z)
    void zoom_in(float z) { zoom_to( zoom * z ); }
    
    /// decrease zoom (zoom <- zoom / z)
    void zoom_out(float z) { zoom_to( zoom / z ); }
    
    //---------------------------------------------------------------------------
    
    /// return ROI i-th point, i in { 0, 1 }
    Vector3 roi(size_t i) const { return mROI[i]; }
    
    /// adjust zoom and focus to match the ROI specificed by two corner points
    void matchROI();
    
    /// set ROI to match the current view
    void adjustROI(float Z);
    
    /// define ROI
    void setROI(Vector3, Vector3);
    
    /// return 'true' if point is inside ROI
    bool insideROI(Vector3) const;
    
    /// display zoomed in regions around position (mX, mY)
    void drawMagnifier(float factor, Vector3 foc, Vector3 cen, int mX, int mY, int R) const;
    
    /// draw text for a brief moment near top-right corner
    bool flashText(std::string const& str);

    //---------------------------------------------------------------------------

    /// Fonts inherited from GLUT
    enum FontType
    {
        BITMAP_9_BY_15 = 2,
        BITMAP_8_BY_13 = 3,
        BITMAP_TIMES_ROMAN_10 = 4,
        BITMAP_TIMES_ROMAN_24 = 5,
        BITMAP_HELVETICA_10 = 6,
        BITMAP_HELVETICA_12 = 7,
        BITMAP_HELVETICA_18 = 8
    };
    
    /// display text in top-left corner
    void strokeString(const char[], float width) const;

    /// stroke text at given position
    void strokeString(float X, float Y, float Z, const char str[], float scale=1.f) const;
    
    /// draw `text` at position `pos`
    void drawText(float X, float Y, float Z, const float color[4], const char text[], float offset=0.5, FontType = BITMAP_9_BY_15) const;

    /// draw `text` at position `pos`
    void drawText(Vector3 const& pos, const float color[4], const char text[], float offset=0.5, FontType font = BITMAP_9_BY_15) const
    {
        drawText(pos.XX, pos.YY, pos.ZZ, color, text, offset, font);
    }

    /// draw `text` at position `pos`
    void drawText(Vector2 const& pos, const float color[4], const char text[], float offset=0.5, FontType font = BITMAP_9_BY_15) const
    {
        drawText(pos.XX, pos.YY, 0, color, text, offset, font);
    }
    
    /// draw `text` at position `pos`
    void drawText(Vector1 const& pos, const float color[4], const char text[], float offset=0.5, FontType font = BITMAP_9_BY_15) const
    {
        drawText(pos.XX, 0, 0, color, text, offset, font);
    }

    /// display text on a rectangle of color `bcol`, in a corner of the center of the display window
    void frameText(int position, FontType, const float color[4], const char text[], const float back[4], int width, int height) const;

    /// display a scale bar vertical or horizontal
    void drawScaleHV(float, float, float, void (*func)(float*, int cnt, float, float, float), float=0.f) const;
    
    /// display crossed scale bars
    void drawScaleX(float, float=0.f) const;

    /// display a scale bar (mode is vertical, horizontal, centered)
    void drawScaleBar(int mode, float, const float[4]) const;
    
};

#endif
