// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.
#ifndef POINT_DISP_H
#define POINT_DISP_H

#include "gym_color.h"
#include "property.h"
#include "gle.h"
#include "gym_zoo.h"
#include "gym_view.h"
#include "gym_vect.h"
#include "gym_draw.h"


/** \todo: replace PIXELMAPS by a texture, and draw each particle as a textured Square */
#define POINTDISP_USES_PIXELMAPS 1


/// the parameters necessary to display a point-like object
/**
 This is used by Hands, Sphere, Solid, Beads
 Note that these parameters may be interpreted differently when displaying different classes. 
 For example `coloring` is only implemented for Sphere, Beads but not for Hands,
 while `shape` and `symbol` are only implemented for Hands.
 */
class PointDisp : public Property
{
private:

#if POINTDISP_USES_PIXELMAPS
    
    /// pointer to 3 square bitmaps with RGBA components
    uint8_t * pixels_;
    
    /// texture ID
    GLuint texture_;
    
    /// size in pixels, must be a power of 2
    unsigned pixSize;

    /// allocated size for square bitmaps, in pixels
    unsigned pixAlloc_;

    /// allocate pixelmap memory
    void allocatePixelmap(unsigned);
    
    /// release pixelmap memory
    void releasePixelmap();
    
    /// create the pixelmaps
    void makePixelmaps(unsigned supersampling, unsigned dim);
    
    /// create the pixelmaps
    void makePixelmaps(unsigned supersampling);

    /// create the pixelmaps
    void createPixelmaps();

    /// draw pixel map
    void drawPixelmap(float X, float Y, float Z, size_t) const;

#endif

    /// used to differentiate between different uses of the class
    std::string mKind;

    /// draw active state with OpenGL vector primitives
    void strokeA(float) const;
    
    /// draw inactive state with OpenGL vector primitives
    void strokeI() const;
    
    /// clear pointers
    void clearPixelmaps();

public:
    
    /**
     @defgroup PointDispPar Display parameters: Points
     @ingroup DisplayParameters
     @{
     */
    
    
    /// visibility flag : 0=hidden, 1=opaque
    /**
     For a Space, bits are used to enable the display of front/back sides,
     such that one can use these values:
     - 1 : display only faces that are facing outside,
     - 2 : display only faces that are facing inside,
     - 3 : display both sides.
     .
     */
    int visible;
    
    /// color of object (in 3D display, the color of outer surfaces)
    gym_color color;
    
    /// second color (set as color[1])
    /**
     This is used to display unattached Single and unbridging Couple, 
     and the inner surfaces of objects such as Sphere, Solid, Bead and Space.
     If it is not defined, `color2` is set to be a darker tone of `color`.
     */
    gym_color color2;

    /// if true, attribute random colors to individual objects
    int coloring;
    
    /// display diameter of points in pixel units
    float size;
    
    /// display width of lines
    float width;
    
    /// miscellaneous display parameter
    float scale;

    /// 'c' for circle, 'h' for hexagon, 's' for star, etc.
    char shape;
    
    /// a bitfield to set different display options
    int style;
    
    /// character displayed (do not set, or set as 0 to disable this feature)
    char symbol;
    
    /// color of symbol (set as symbol[1])
    gym_color colorS;
    
    /// @}
    
    /// visible and big enough to be seen
    bool perceptible;
    
    /// radius of feature in simulation units
    float ulna_;
    
    /// rescaled point size in pixels
    float sizeX;
    
    /// rescaled line width in pixels
    float widthX;

public:
    
    /// constructor
    PointDisp(const std::string& k, const std::string& n);
    
    /// copy constructor
    PointDisp(PointDisp const&);
    
    /// copy assignment
    PointDisp& operator = (PointDisp const&);

    /// destructor
    ~PointDisp();
    
    /// identifies the property
    std::string category() const { return mKind; }

    /// clear to default values
    void clear();
    
    /// set from glossary
    void read(Glossary&);
    
    /// return a carbon copy of object
    Property* clone() const { return new PointDisp(*this); }

    /// write all values
    void write_values(std::ostream&) const;
    
    /// recalculate bitmaps
    void setPixels(float ps, float uv, bool make_maps);
    
    /// draw inactive state
    template < typename VECTOR >
    void drawI(VECTOR const& vec) const
    {
        if ( perceptible )
        {
    #if POINTDISP_USES_PIXELMAPS
            drawPixelmap(vec.XX, vec.y(), vec.z(), 0);
    #else
            gym::transScale(vec, radius);
            gym::color(color2);
            gle::disc();
    #endif
        }
    }

    /// draw active state, unattached
    template < typename VECTOR >
    void drawF(VECTOR const& vec) const
    {
        if ( perceptible )
        {
    #if POINTDISP_USES_PIXELMAPS
            drawPixelmap(vec.XX, vec.y(), vec.z(), 1);
    #else
            gym::transScale(vec, radius);
            gym::color(color2);
            strokeA(widthX);
    #endif
        }
    }

    /// draw active state, attached
    template < typename VECTOR >
    void drawA(VECTOR const& vec) const
    {
        if ( perceptible )
        {
    #if POINTDISP_USES_PIXELMAPS
            drawPixelmap(vec.XX, vec.y(), vec.z(), 2);
    #else
            gym::transScale(vec, radius);
            gym::color(color);
            strokeA(widthX);
    #endif
        }
    }

};


#endif
