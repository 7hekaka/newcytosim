// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
// Created by Francois Nedelec on 08/08/2010.


#include "gle_color.h"


namespace gle
{
    /// Holds the name and the RGBA components of a color
    class named_color
    {
    public:
        /// name of color
        char const * name;
        /// components of color (RGBA)
        uint32_t hex;
    };
    
    /// a small set of contrasted colors (indx is wrapped to the number of colors)
    gle_color nice_color(size_t indx);
  
    /// a set of standard colors (indx is wrapped to the number of colors)
    gle_color std_color(size_t indx);

    /// a set of standard named html colors
    gle_color std_color(const std::string& name);

    /// a large set of colors from Crayola crayons
    gle_color alt_color(size_t indx);
    
    /// extract colors having a brightness between `minb` and `maxb`
    int       select_colors(gle_color list[], size_t list_size, GLfloat minb, GLfloat maxb);
    
    /// one of the crayola color, with a brightness() > threshold
    gle_color bright_color(size_t indx, GLfloat minb = 0.6f, GLfloat maxb = 3.0f);

    /// print a list of colors with names and values
    void      print_colors(std::ostream&, named_color* list, size_t list_size);
    
    /// print the list of standard colors
    void      print_std_colors(std::ostream&);

}

