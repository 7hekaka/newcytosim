// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
// Created by Francois Nedelec on 08/08/2010.

#ifndef GLE_COLOR_LIST_H
#define GLE_COLOR_LIST_H

#include "gle_color.h"


namespace gle
{
    class named_color;

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
    
    /// one of the crayola color, with a brightness() in [minb, maxb]
    gle_color bright_color(size_t indx, GLfloat minb, GLfloat maxb);

    /// one of the crayola color, not too dark
    inline gle_color bright_color(size_t indx) { return bright_color(indx, 0.6, 3.0); }

    /// print a list of colors with names and values
    void      print_colors(std::ostream&, named_color* list, size_t list_size);
    
    /// print the list of standard colors
    void      print_std_colors(std::ostream&);

}

#endif

