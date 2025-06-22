// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef FIBER_DISP_H
#define FIBER_DISP_H

#include "real.h"
#include "assert_macro.h"
#include "gym_color.h"
#include "property.h"
#include "vector3.h"


/// Display parameters for a class of Fiber
/**
 Holds the display attributes for a certain class of Fiber.
 
 There is one FiberDisp for each FiberProp.
 */
class FiberDisp : public Property
{
public:

    /// possible values for fiber:coloring
    enum ColoringModes {
        COLORING_OFF,
        COLORING_RANDOM,
        COLORING_DIRECTION,
        COLORING_MARK,
        COLORING_FLAG,
        COLORING_FAMILY,
        COLORING_CLUSTER,
        COLORING_AGE,
        COLORING_PSTATE
    };
    
public:
    
    /**
     @defgroup FiberDispPar Display parameters: Fibers
     @ingroup DisplayParameters
     @{
     */
    
    /// general rendering style
    /**
     Possible values of `style`:
     - 0 or 'line'        : line or cylinders for style=3 (this is the default),
     - 1 or 'backbone'    : thin broken line
     - 2 or 'striped'     : tube with colored stripes every 24nm.
     - 3 or 'filament'    : a stretch of beads rendering a proto-filament.
     - 4 or 'actin'       : actin-like rendering using beads for monomers,
     - 5 or 'microtubule' : microtubule-like rendering using beads for monomers.
     .
     */
    int style;
    
    /// visibility flag
    int visible;
    
    /// color of fiber
    gym_color color;
    
    /// color of inner surfaces of cylinder in 3D display (set as color[1])
    gym_color back_color;
    
    /// color for unselected objects, default=invisible (set as color[2])
    gym_color hide_color;

    /// if true, vary the colors used to display the fibers
    /**
     This option is used to attribute a different color to each fiber,
     painting all the segments of the fiber with the same color.
     
     Effects of `coloring`:
     - 0 : no coloring,
     - 1 : color fibers randomly,
     - 2 : color fibers depending on direction, relative to `right_direction`,
     - 3 : color fibers depending on the `mark`,
     - 4 : color clusters defined by couple-connectivity,
     - 5 : color fibers according to age.
     .
     */
    int coloring;

    /// width of lines (also known as `line[0]` or `width`)
    float line_width;

    /// style for lines (also known as `line[1]`)
    /**
     Possible values of `line_style`:
     - 0 : hide,
     - 1 : plain lines,
     - 2 : color rendering of longitudinal tensions (2)
     - 3 : color rendering of longitudinal tensions, using a jet color scheme (2)
     - 4 : color rendering of local curvature (1)
     - 5 : color rendering of the angular orientation relative to the X-axis
     - 6 : color rendering based on distance from minus end (1)
     - 7 : color rendering based on distance from plus end (1)
     - 8 : color rendering based on distance to edge of Space (1)
     .
     The color value is calculated using (1): 'length_scale' or (2): 'tension_scale'.
     */
    int line_style;
    
    /// if true, close the end of the fiber (valid only for style==3)
    /**
     Possible values of `line_caps`:
     - 0: leave fibers open (unfinished),
     - 1: use a disc to make a flat end,
     - 2: use a hemisphere to make a round end.
     This is only valid for style==3
     */
    int line_caps;
    
    /// width of lines used to draw backbone
    float bone_width;

    /// if set, draw fiber's outline with the background color (supported only for DIM==3)
    float outline_width;

    /// diameter of points (also known as `point[0]` or `size`)
    /**
     `point_size` and `line_width` are normally set in pixels, 
     but if `display`:point_value is set, their value is understood 
     in multiples of `point_value`, which itself is a distance.
     
     For example, if you set line_width=2.5 and point_value=0.01,
     the fibers will be displayed with a diameter of 0.025.
     */
    float point_size;
    
    /// style for display of points (also known as `point[1]`)
    /**
     Possible values for `point_style`:
     - 0 : show nothing,
     - 1 : show vertices,
     - 2 : show arrowheads separated by `point_gap`,
     - 3 : draw chevrons separated by `point_gap`,
     - 4 : show middle point of each fiber
     .
     */
    int point_style;
    
    /// distance between arrows for `point_style=2` (also known as `point[2]`)
    real point_gap;
    

    /// style of fiber tips for { plus end, minus end }
    /**
     `end_style[0]` determines the style of the plus end,
     and `end_style[1]` the style of the minus end.
     
     Possible end_style:
     - 0 : hide,
     - 1 : disc/sphere,
     - 2 : cone,
     - 3 : cylinder centered on the fiber end,
     - 4 : arrow fins,
     - 5 : arrow fins in the reversed direction
     - 6 : cube
     - 7 : cylinder aligned with the fiber end
     - 8 : hemisphere
     .
     */
    int end_style[2];
    
    /// size of fiber tips for { plus end, minus end }
    /**
     You can also specify:
     
         plus_end  = SIZE, STYLE
         minus_end = SIZE, STYLE
         
     */
    float end_size[2];
    
    /// colors of the different FiberTip states
    /**
     This determines the set of color that are used to display the fiber tips,
     according to their assembly state, Fiber::endState():
     - end_colors[0]: static ends,
     - end_colors[1]: growing ends,
     - end_colors[2] and [3]: intermediate states,
     - end_colors[4]: shrinking ends
     .
     The default colors are: white, green, yellow, orange, red
     */
    gym_color end_colors[6];
    
    /// draw fresh polymer assembly
    int growth_style;

    /// if true, specify the style for displaying lattice content (also known as `lattice[0]`)
    int lattice_style;
    
    /// defines the range of colors when displaying the lattice (also known as `lattice[1]`)
    real lattice_scale;
    
    /// rescale concentration for the cells at the edge of reduced length
    bool lattice_rescale;
    
    /// style of labels
    /**
     Possible `label_style`:
     - 0 : hide,
     - 1 or 2 : name of fiber and index of vertices
     - 4 : abscissa along fiber
     .
     */
    int label_style;
    
    
    /// size for speckle display (also know as `speckles`)
    float speckle_size;

    /// style for speckle display (also know as `speckles[1]`)
    /**
     Possible values for `speckle_style`:
     - 0 : hide,
     - 1 : random speckles, separated on average by `interval`,
     - 2 : regular speckes, separated by `interval`.
     .
     */
    int speckle_style;

    /// average distance between speckles (also known as `speckles[2]`)
    real speckle_gap;

    /// a bit-field to hide certain categories of fibers
    /**
     Possible values for `exclude`:
     - 0 : all fibers are displayed,
     - 1 : show only right-pointing fibers,
     - 2 : show only left-pointing fibers,
     - 4 : show only counter-clockwise fibers,
     - 8 : show only clockwise fibers.
     .
     
     You may also address each bit directly, knowning that:
     - 1st bit true: hide left-pointing fibers
     - 2nd bit true: hide right-pointing fibers
     - 3rd bit true: hide clockwise fibers
     - 4th bit true: hide counter-clockwise fibers
     .
     */
    int hide;
    
    /// the direction used for hiding left- or right-pointing fibers, etc. (known as `exclude[1]`)
    Vector3 hide_axis;
    
    /// hide filament with specified state
    unsigned hide_state;
    
    /// hide filament which do not have specified mark (default: disabled)
    unsigned show_marked;

    /// number of bits equal to `1` in the mask bitfield
    /**
     This parameter can be used to hide a fraction of the fibers.
     Each fiber will be visible with a probability `1/2^mask`.
     `mask_bitfield' is set randomly with `mask` bits set to 1, 
     When the parameter is read.
     */
    unsigned mask;
    
    /// selection bitfield used to hide some fibers (known as `mask[1]`)
    /**
     A 32-bitfield is compared to the signature of the object,
     itself a random bitfield. The Object is hidden if the result is non-zero.
     So if the mask bitfield has many 1s, fewer filaments will be visible.
     Note that `mask_bitfield' is set automatically if `mask` is given.
     */
    unsigned mask_bitfield;
    
    
    /// conversion coefficient from length to color, for some line styles
    real length_scale;

    /// conversion coefficient from tension to color, for `line_style==3` (tension)
    /**
     The values of `tension_scale` determines how longitudinal tensions are displayed:
     - tension_scale < 0 : compressive forces are highlighted,
     - tension_scale > 0 : tensile forces are highlighted.
     .
     A longitudinal tension equal to `tension_scale` will be displayed with a blue tint,
     while a value three times higher will be displayed red.
     Lower tension_scale values will yield brighter colors for the same force in the fiber.
     */
    real tension_scale;

    /// ( if != 0 ) display the net forces acting on vertices (known as `force`)
    int force_style;

    /**
     if `force_style != 0), a force `F` acting on a vertex is displayed as a
     segment of length `force_scale * F`.
     ( known as force[1], default = 1 )
     */
    real force_scale;
    
    /// this color is specified as forces[2]
    gym_color force_color;   
    
    /// if true, display the average fiber
    /**
     The 'average fiber' is calculated from the centroid of the fiber tips,
     and the centroid of the polymer mass.
     It is useful to evaluate the amount of order in the network.
     */
    int draw_average;

    /// @}
    
public:
    
    /// constructor
    FiberDisp(const std::string& n) : Property(n)  { clear(); }
    
    /// destructor
    ~FiberDisp() { }
    
    /// identifies the property
    std::string category() const { return "fiber:display"; }

    /// clear to default values
    void clear();
    
    /// set from glossary
    void read(Glossary&);
    
    /// return a carbon copy of object
    Property* clone() const { return new FiberDisp(*this); }

    /// write all values
    void write_values(std::ostream&) const;

};


#endif

