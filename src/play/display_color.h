// Cytosim was created by Francois Nedelec.
// Copyright 2020 Sainsbury Laboratory, Cambridge University

#include "point_disp.h"
#include "gle_color_list.h"



/**
 returns the front color used to display Object
 if `coloring` is enabled, this loads the N-th bright color,
 otherwise load the object's display color
 */
template < typename T >
gle_color bodyColorF(T const& obj)
{
    const PointDisp * disp = obj.prop->disp;
    if ( disp->coloring )
    {
        size_t i = ( disp->coloring == 2 ? obj.mark() : obj.signature());
        return gle::bright_color(i).match_a(disp->color);
    }
    return disp->color;
}

/**
 returns the front color used to display Object
 if `coloring` is enabled, this loads the N-th bright color,
 otherwise load the object's display color
 */
inline gle_color bodyColorF(PointDisp const* disp, size_t s)
{
    if ( disp->coloring )
        return gle::bright_color(s).match_a(disp->color);
    return disp->color;
}

/**
Sets color material for lighting mode
if `coloring` is enabled, this loads a bright color,
otherwise load the object's display color
*/
template < typename T >
void bodyColor(T const& obj)
{
    const PointDisp * disp = obj.prop->disp;
    if ( disp->coloring )
    {
        size_t i = ( disp->coloring == 2 ? obj.mark() : obj.signature());
        gle_color col = gle::bright_color(i);
        col.load_front();
        col.darken(0.5).load_back();
    }
    else
    {
        disp->color.load_front(1.0);
        disp->color2.load_back();
    }
}

/**
 Sets color material for lighting mode
 if `coloring` is enabled, this loads the N-th bright color,
 otherwise load the object's display color
 */
inline void bodyColor(PointDisp const* disp, size_t s)
{
    if ( disp->coloring )
    {
        gle_color col = gle::bright_color(s);
        col.load_front();
        col.darken(0.5).load_back();
    }
    else
    {
        disp->color.load_front(1.0);
        disp->color2.load_back();
    }
}

