// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University

#ifndef PLAYER_PROP_H
#define PLAYER_PROP_H

#include "real.h"
#include "property.h"
#include "quaternion.h"


/// Parameters for the Player
class PlayerProp : public Property
{
    
public:
    
    /// number of programmable keys
    static constexpr int NB_MAGIC_KEYS = 4;

public:
    
    /**
     @defgroup PlayPar Parameters of Play
     @ingroup DisplayParameters
     @{
     */
    
    /// direction of replay: 1: forward and -1: reverse
    int replay;
    
    /// if true, jump to first frame after last frame
    unsigned loop;
    
    /// number of simulation steps done between two drawings
    /**
     if period==2, only half of the frames will be displayed
     */
    unsigned period;
    
    /// number of milli-seconds between refresh
    unsigned delay;
    
    /// number of images to export
    unsigned save_images;
    
    /// if > 1, downsample images before writing them out
    /**
     This can be used to reduce pixelation artifacts:
     Specify a larger image size than desired, and the equivalent downsampling,
     For example, to produce a 512x256 image:
     
     window_size = 1024, 512
     downsample = 2
     
     */
    unsigned downsample;

    /// specifies information displayed near the bottom left corner of window
    std::string report;

    /// format of exported images [png, ppm]
    std::string image_format;
    
    /// name of image to be saved
    mutable std::string image_name;

    /// directory where images are exported
    std::string image_dir;

    /// associate a piece of custom code to a key
    /**
     Example:
     
     % define a magic key to delete fibers:
     set system display
     {
     magic_key1 = m, ( delete 10 microtubule )
     magic_key2 = C, ( cut microtubule { plane = 1 0 0, 0 } )
     label = (Press 'm' to delete fibers!)
     }
     
     up to 4 keys (magic_key, magic_key1 ... 3) can be defined.
     */
    char magic_key[NB_MAGIC_KEYS];

    /// if true, program will quit when end-of-file is reached
    int auto_pilot;

    /// if true, program will quit when end-of-file is reached
    unsigned auto_exit;
    
    /// zoom applied after reading
    float auto_zoom;
    
    /// rotation applied after reading
    Quaternion<real> auto_rotate;
    
    /** @} */

    /// index of report which is displayed
    unsigned report_index;

    /// index used to build the name of the exported image
    unsigned image_index;
    
    /// index used to build the name of the exported poster
    unsigned poster_index;
    
    /// time of last image exported
    double saved_image_time;
    
    /// the piece of cytosim code executed when `magic_key` is pressed (set as `magic_key[1]`)
    std::string magic_code[NB_MAGIC_KEYS];
    
public:

    /// constructor
    PlayerProp(const std::string& n) : Property(n)  { clear(); }
    
    /// destructor
    ~PlayerProp() { }
    
    /// identifies the property
    std::string category() const { return "simul:display"; }

    /// set default values
    void clear();
    
    /// set from a Glossary
    void read(Glossary&);
    
    /// return a carbon copy of object
    Property* clone() const { return new PlayerProp(*this); }

    /// write all values
    void write_values(std::ostream&) const;
    
    /// change `report` that is displayed onscreen
    void toggleReport(int alt);
};


#endif


