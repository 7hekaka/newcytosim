// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.
// Created by Francois Nedelec on 08/08/2010.

#ifndef GYM_COLOR_LIST_H
#define GYM_COLOR_LIST_H

#include "gym_color.h"


namespace gym
{
    class named_color;
    
    /// color used to select a 'bright_color()'
    extern gym_color background_color;

    /// return one color amongs a small set of contrasted colors
    gym_color get_color(unsigned indx);
  
    /// change indx-th color
    void set_color(unsigned indx, gym_color const&);

    /// a set of ~150 standard colors (indx is wrapped to the number of colors)
    gym_color std_color(unsigned indx);

    /// a set of ~150 standard named html colors
    gym_color std_color(const std::string& name);

    /// number of colors from the XKCD project
    unsigned nb_alt_color();

    /// one of People's 949 popular colors, from the XKCD project
    gym_color alt_color(unsigned indx);
    
    /// set a list of colors that are different from `back`
    unsigned filter_colors(gym_color list[], unsigned list_size, const gym_color back);

    /// one of the ~260 crayola color, that differs significantly from `back`
    gym_color bright_color(unsigned indx, gym_color back=background_color);

    /// print a list of colors with names and values
    void print_colors(std::ostream&, named_color* list, unsigned list_size);
    
    /// print the list of standard colors
    void print_std_colors(std::ostream&);

}

#endif

