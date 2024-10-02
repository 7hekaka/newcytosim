// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University
#ifndef VIEW_PROP_H
#define VIEW_PROP_H

#include "real.h"
#include "vector3.h"
#include "vector4.h"
#include "quaternion.h"
#include "property.h"
#include "gym_color.h"

///properties needed to define a view
class ViewProp : public Property
{
public:
    
    /// number of user-accessible OpenGL clipping planes
    static constexpr int NB_CLIP_PLANES = 4U;
    
    /**
     @defgroup ViewPar Display Parameters: View
     @ingroup DisplayParameters
     @{
     */
    
    /// zoom factor = ratio between visible area and `view_scale`
    float zoom;
    
    /// pixel magnification
    unsigned magnify;

    /// range that is visible if `zoom=1`, in sim-units (default=10)
    float view_scale;
    
    /// enables the display area to be set from the size of the simulation space
    /**
     If ( `auto_scale` > 0 ), `view_scale` is set automatically to match the simulation space.
     This is on by default.
     */
    unsigned auto_scale;
    
    /// the point that is in the center of the window in real-world coordinates
    Vector3 focus;
    
    /// additional translation used by autoFocus()
    Vector3 focus_shift;
    
    /// orientation of display
    Quaternion<real> rotation;
    
    /// flag to enable perspective view in 3D
    /**
     By default, cytosim uses a orthographic projection to view the 3D space,
     but it will use a 3D perspective if 'perspective==true'.
     This is only meaningful in 3D mode.
     */
    int perspective;
    
    /// modifies the display to show only the front, the back or a slice of the world
    /**
     possible values are:
     - `off`    (0)
     - `front`  (1)
     - `back`   (2)
     - `slice`  (3)
     .
     */
    unsigned int slice;

    /// color of background
    gym_color back_color;
    
    /// color used to highlight objects
    gym_color front_color;

    /// flag to use a double buffer for smoother rendering (default=1)
    /**
     http://en.wikipedia.org/wiki/Multiple_buffering#Double_buffering_in_computer_graphics
     */
    bool buffered;

    /// flag to enable OpenGL depth buffer (default=1)
    /**
     This is useful for 3D rendering.
     http://en.wikipedia.org/wiki/Z-buffering
     */
    int depth_test;
    
    /// flag to perform depth-clamp (default=false)
    /** http://www.opengl.org/registry/specs/NV/depth_clamp.txt */
    int depth_clamp;

    /// flag to enable native device resolution on mac osx
    /**
     This works only if you use Renaud Blanch's modified GLUT
     http://iihm.imag.fr/blanch/software/glut-macosx
     */
    int retina;
    
    /// flag to enable OpenGL stencil buffer (default=0)
    int stencil;
    
    /// if > 0, enables OpenGL full scene anti-aliasing (default=0)
    /**
     This defines the number of samples used to build an image.
     Higher values result in nicer (but slower) display.
     http://en.wikipedia.org/wiki/Multisample_anti-aliasing
     Many graphic cards only support 8 samples max, so try 4 or 8.
     */
    int multisample;
    
    
    /// string at start of `message` (if `none` is specified, no message is shown)
    std::string label;
    
    /// text shown near the bottom of window
    std::string subtitle;

    /// automatically adjust view to keep fibers in window
    /**
     Possible values:
     - 0 : off
     - 1 : translate to track the center of gravity of the cloud of fiber-points
     - 2 : rotate to align the principal direction of the fiber
     - 3 : translate and rotate ( 1 and 2 are combined )
     - 4 : rotate to align two principal directions
     - 5 : translate and rotate ( 1 and 4 are combined )
     .
     The translation defined by focus is applied after this adjustment.
     */
    unsigned track_fibers;
    
    /// position of window on screen (top-left corner, in pixels)
    int window_position[2];
    
    /// desired size of window in pixels (also known as `size`)
    int window_size[2];
    
    /// display flag for a scalebar (default=0)
    unsigned scalebar;

    /// length of scale-bar in sim-world units (set as `scalebar[1]`)
    real scalebar_length;
    
    /// color of scale-bar (set as `scalebar[2]`)
    gym_color scalebar_color;

    /// display flag for displaying X-Y-Z axes
    unsigned axes;
    
    /// length of axes (set a `axes[1]`, default=1)
    float axes_size;

    /// on/off flags for clipping (defined as `clip_plane?`)
    /**
     Up to 4 clipping planes can be defined: clip_plane0 to clip_plane3
     
     Syntax:
     
         clip_plane? = BOOL, VECTOR, REAL
     
     The Boolean enables the clipping plane.
     The plane is specified by a normal vector `n` (VECTOR) and a scalar `a` (REAL).
     The visible half-space corresponds to <em> n.pos + a > 0 </em>
     
     Example:
     To define a slice perpendicular to the X-axis of width 2: 
     
         set system display
         {
            clip_plane1 = 1,  1 0 0, 1
            clip_plane2 = 1, -1 0 0, 1
         }
     
     */
    int clip_plane_mode[NB_CLIP_PLANES];

    /// direction perpendicular to clipping plane (defined as `clip_plane?[1]`)
    Vector4 clip_plane[NB_CLIP_PLANES];

    /// characteristics of OpenGL fog (also known as `fog[0]`)
    int fog_type;
    
    /// density of fog (also known as `fog[1]`)
    float fog_param;
    
    /// color of fog (also known as `fog[2]`)
    gym_color fog_color;

    /// draw 'tiled' floor
    int floor_radius;
    
    /// parameters for 'tiled' floor
    float floor_tile, floor_height;
    
    /// color of floor
    gym_color floor_color;

    /// @}
    
    /// text displayed in center of window
    std::string memo;
    
    /// flag to display information on screen
    int draw_memo;
    
public:
   
    /// constructor
    ViewProp(const std::string& n) : Property(n)  { clear(); }
    
    /// destructor
    ~ViewProp()  { }
    
    /// identifies the property
    std::string category() const { return "view"; }
     
    /// set default values
    void clear();
    
    /// set from a Glossary
    void read(Glossary&);
    
    /// return a carbon copy of object
    Property* clone() const { return new ViewProp(*this); }

    /// write all values
    void write_values(std::ostream&) const;

    /// invert front and back colors
    void invertColors();
    
    /// use back and white for front/back colors
    void blackAndWhite();

};

#endif
