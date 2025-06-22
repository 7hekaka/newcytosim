// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef DISPLAY_PROP_H
#define DISPLAY_PROP_H

#include "property.h"


/// Parameters for Play
class DisplayProp : public Property
{

public:
    
    /**
     @defgroup DisplayPar Display parameters: World
     @ingroup DisplayParameters
     @{
     */
    
    /// style of display { 1, 2, 3 }
    /**
     3 styles are implemented:
     - style 1 used OpenGL lines and points. It is the best style for 2D work.
     - style 2 uses pixelmap to draw Hands, and is slower than style 1, also made for 2D.
     - style 3 is for 3D. It draws 3D tubes and uses OpenGL lighting for rendering.
     .
     */
    unsigned style;

    /// if true, repeat the display along periodic boundary directions
    /**
    This is only useful is the main Space is periodic, strip or cylinderP
    In this case, the whole system will be displayed multiple times,
    shifted appropriately in the directions that are periodic.
    */
    int tile;
    
    /// if true, translate objects to place them in the root cell for periodic boundary conditions
    int fold;
    
    /// default diameter of points
    float point_size;

    /// default width of lines
    float line_width;
    
    /// default width of links
    float link_width;

    /// if set > 0, this defines the unit scale used for `point_size` and `line_width`
    /**
     Set this parameter to specify the fiber radius and point size in length units.

     `point_size` and `line_width` are usually specified in pixels, but if `point_value` is set,
     they are interpreted as multiples of `point_value`, which itself is given in length unit.
     
     For example, by setting `line_width=25` and `point_value=0.001`,
     the fibers will be displayed with a diameter of 0.025 (i.e. 25 nanometers).
     
     <em> default = 0 </em>
     */
    float point_value;
    
    /// if `true`, unattached Couples are display randomly with one or the other Hand (default=false)
    unsigned couple_flip;
    
    /// selection bitfield for Couples
    /**
     This is a bitfield:
     
         bit 1 (value 1) means `show free couple`
         bit 2 (value 2) means `show bound couple`
         bit 3 (value 4) means `show bridge couple`
     
     Hence value 7 (the default) will show all couples, while 3 will show free and bound couples but not bridge couple, etc.
     */
    unsigned couple_select;
    
    /// selection bitfield for Singles
    /**
     This is a bitfield:
     
         bit 1 (value 1) means `show free single`
         bit 2 (value 2) means `show bound single`
     
     Hence value 3 (the default) will show all singles.
     */
    unsigned single_select;
    
    /// flag to display Meca's links
    bool draw_links;
    
    
    /// the 'explosion' effect shift the fibers in space
    /**
     This can be useful to visualize dense regions,
     but is only implemented for style=2
     */
    int explode_style;
    
    /// amount of lateral shift to separate fibers when display is exploded (known as `explode[1]`)
    float explode_range;
 
    /// @}

public:
    
    /// constructor
    DisplayProp(const std::string& n) : Property(n) { clear(); }
    
    /// destructor
    ~DisplayProp() { }
    
    /// identifies the property
    std::string category() const { return "simul:display"; }
        
    /// set default values
    void clear();
    
    /// set from a Glossary
    void read(Glossary&);
    
    /// return a carbon copy of object
    Property* clone() const { return new DisplayProp(*this); }
    
    /// write all values
    void write_values(std::ostream&) const;

};


#endif


