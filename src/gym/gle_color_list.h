// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.
// Created by Francois Nedelec on 08/08/2010.

#ifndef GLE_COLOR_LIST_H
#define GLE_COLOR_LIST_H

#include "gle_color.h"


namespace gle
{
    class named_color;

    /// a small set of contrasted colors (indx is wrapped to the number of colors)
    gle_color nice_color(size_t indx);
  
    /// a set of ~150 standard colors (indx is wrapped to the number of colors)
    gle_color std_color(size_t indx);

    /// a set of ~150 standard named html colors
    gle_color std_color(const std::string& name);

    /// number of colors from the XKCD project
    size_t nb_alt_color();

    /// one of People's 949 popular colors, from the XKCD project
    gle_color alt_color(size_t indx);
    
    /// set a list of colors that are different from `back`
    size_t filter_colors(gle_color list[], size_t list_size, const gle_color back);
    
    /// color used to select a 'bright_color()'
    extern gle_color background_color;
    
    /// one of the ~260 crayola color, that differs significantly from `back`
    gle_color bright_color(size_t indx, gle_color back=background_color);

    /// print a list of colors with names and values
    void print_colors(std::ostream&, named_color* list, size_t list_size);
    
    /// print the list of standard colors
    void print_std_colors(std::ostream&);

}

#endif

