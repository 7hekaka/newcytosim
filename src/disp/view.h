// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University
#ifndef VIEW_H
#define VIEW_H

#include "opengl.h"
#include "view_prop.h"
#include "gym_text.h"

/// Handles the viewing angle, projection and other aspects of an OpenGL display
/**
 ViewProp does not depend on the window-system (GLUT),
 but only on the rendering engine (OpenGL)
 */
class View : public ViewProp
{
private:
    
    /// viewport obtained by getMatrices()
    GLint mViewport[4];
    
    /// modelview obtained by getMatrices()
    GLfloat mModelview[16];
    
    /// projection obtained by getMatrices()
    GLfloat mProjection[16];
    
    /// half-size of the OpenGL visible region in OpenGL units
    GLfloat visRegion[4];
    
    /// translation between center of volume and camera
    GLfloat eyePosition[4];
    
    /// window number in GLUT
    int mWindowId;
    
    /// flag for Region of Interest
    bool hasROI;
    
    /// Region of interest
    Vector3 mROI[2];

    /// used to check that getMatrices() was called
    bool hasMatrices;
    
    /// text displayed near top right corner of window
    std::string top_message;
    
    /// text displayed near bottom left corner of window
    std::string full_label;

    /// display callback
    void (*displayCallback)(View&, int);
    
    /// set OpenGL Fog, with mode (GL_EXP or GL_LINEAR), intensity and color
    void setFog(GLint mode, GLfloat density, gle_color) const;

public:
    
    /// constructor
    explicit View(const std::string& n);
    
    /// destructor
    ~View();
    
    /// return window-id
    int window() const { return mWindowId; }
    
    /// set window-id
    void window(int w) { mWindowId = w; }
    
    /// handle window resize events
    void reshape(int, int);
    
    /// adjust parameters of projections
    void adjust();
    
    /// set OpenGL Projection matrix
    void setProjection();
    
    /// set OpenGL Model-View matrix
    void setModelView() const;
    
    /// set OpenGL Projection and ModelView matrices
    void load();

    /// adjust view to only show a slice of the world
    void sliceView(int) const;
    
    /// reset the view (no-rotation, zoom=1), and enable auto_scale
    void reset();
    
    /// width of display area in pixels
    int width() const { return window_size[0]; }
    
    /// height of display area in pixels
    int height() const { return window_size[1]; }
    
    /// size of pixel in drawing units
    float pixelSize() const { return view_size / ( zoom * std::max(width(), height()) ); }
    
    /// return direction of view that is orthogonal to display screen
    Vector3 depthAxis() const;

    //---------------------------------------------------------------------------
    
    /// set display callback
    void setDisplayFunc(void (*f)(View&, int)) { displayCallback = f; }

    /// clear pixels and set clipping planes and fog parameters
    void openDisplay();
    
    /// unset clipping planes and fog parameters, display axes and scale bar
    void closeDisplay() const;
    
    /// display scale bar, info text, etc.
    void drawInteractiveFeatures() const;
    
    /// call displayCallback
    void display() { displayCallback(*this, 1); }
    
    //---------------------------------------------------------------------------
    
    /// init OpenGL parameters
    void initGL();

    /// toggle depth clamp GL capability
    void toggleDepthClamp();
    
    /// set OpenGL Lights for lighting effects
    void setLights(bool local = false) const;
    
    /// set text displayed in center of window
    void setLabel(std::string const& arg) { full_label = label + " " + arg; }
    
    /// set message displayed on current window
    void setMemo(std::string const& arg) { memo = arg; };

    /// set text displayed near top of window
    void setMessage(std::string const& arg) { top_message = arg; }

    /// set OpenGL Fog, with mode (GL_EXP or GL_LINEAR), intensity and color
    void enableFog(GLint mode, GLfloat param, gle_color);
    void enableFog(GLint mode, GLfloat param);
    
    /// enable cliping plane in OpenGL
    void setClipPlane(GLenum glp, Vector3 dir, real sca) const;

    /// enable cliping plane in OpenGL
    void setClipPlaneEye(GLenum glp, Vector3 dir, real sca) const;
    
    /// call setClipPlane(int) for all enabled clipping planes
    void setClipping() const;
    
    /// disable cliping planes in OpenGL
    void endClipping() const;
    
    /// set equations for a clipping plane, and enable it in View
    void enableClipPlane(int, Vector3 dir, real scal, bool absolute=true);
    
    /// disable cliping plane in View
    void disableClipPlane(int);
    
    /// return enable/disable state
    int hasClipPlane(int) const;
    
    //---------------------------------------------------------------------------

    /// store the matrices defining the current OpenGL Model-View and Projection
    void getMatrices();
    
    /// transform window coordinates to 3D world-coordinates
    Vector3 unproject(GLfloat x, GLfloat y, GLfloat z);
    
    //---------------------------------------------------------------------------
    
    /// position 'pos' in the center of the display
    void move_to(const Vector3 & pos);
    
    /// set additional translation of focal point
    void move_shift(const Vector3 & pos);
    
    /// translate view
    void move_by(const Vector3 & trans) { move_to( focus - trans ); }

    //---------------------------------------------------------------------------
    
    /// set rotation to given Quaternion
    void rotate_to(const Quaternion<real>&);
    
    /// rotate to have `dir` aligned with the X-axis
    void align_with(const Vector3& dir);

    /// rotate view
    void rotate_by(const Quaternion<real> &q) { rotate_to( rotation * q ); }
    
    //---------------------------------------------------------------------------

    /// set absolute zoom
    void zoom_to(GLfloat z);
    
    /// increase zoom (multiplicative)
    void zoom_in(GLfloat z) { zoom_to( zoom * z ); }
    
    /// decrease zoom (multiplicative)
    void zoom_out(GLfloat z) { zoom_to( zoom / z ); }
    
    //---------------------------------------------------------------------------
    
    /// return ROI i-th point, i in { 0, 1 }
    Vector3 roi(size_t i) const { return mROI[i]; }
    
    /// adjust zoom and focus to match the ROI specificed by two corner points
    void matchROI();
    
    /// set ROI to match the current view
    void adjustROI(real Z);
    
    /// define ROI
    void setROI(Vector3, Vector3);
    
    /// return 'true' if point is inside ROI
    bool insideROI(Vector3) const;
    
    /// draw cuboid
    static void drawCuboid(Vector3 const&, Vector3 const&);

    /// draw ROI
    void drawROI() const;
    
    /// display zoomed in regions around position (mX, mY)
    void displayMagnifier(GLint factor, Vector3 foc, GLint mX, GLint mY) const;

    //---------------------------------------------------------------------------
    
    /// draw text
    void drawText(FontType, std::string const& str, gle_color, int pos) const;

    /// display a scale bar vertical or horizontal
    void drawScaleH(float, float, float) const;
    
    /// display a scale bar vertical or horizontal
    void drawScaleV(float, float, float) const;
    
    /// display a scale bar vertical or horizontal
    void drawScaleX(float) const;

    /// display a scale bar (mode is vertical, horizontal, centered)
    void drawScaleBar(int mode, float) const;
    
};

#endif
