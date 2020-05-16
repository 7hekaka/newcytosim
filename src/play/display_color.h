// Cytosim was created by Francois Nedelec.
// Copyright 2020 Sainsbury Laboratory, Cambridge University

#include "point_disp.h"
#include "gle_color_list.h"


/**
if `coloring` is enabled, this loads a bright color,
otherwise load the object's display color
*/
template < typename T >
void bodyColor(T const& obj)
{
    const PointDisp * disp = obj.prop->disp;
    if ( disp->coloring )
    {
        unsigned long i = ( disp->coloring == 2 ? obj.mark() : obj.signature());
        gle_color col = gle::bright_color(i);
        col.load();
        col.load_front();
        col.darken(0.5).load_back();
    }
    else
    {
        disp->color.load_load(1.0);
        disp->color2.load_back();
    }
}


/**
 if `coloring` is enabled, this loads the N-th bright color,
 otherwise load the object's display color
 */
template < typename T >
void bodyColor2(T const& obj)
{
    const PointDisp * disp = obj.prop->disp;
    if ( disp->coloring )
    {
        unsigned long i = ( disp->coloring == 2 ? obj.mark() : obj.signature());
        gle::bright_color(i).match_a(disp->color).load();
    }
    else
        disp->color.load();
}

/**
 This is used for transparent objects.
 if `coloring` is enabled, this loads the N-th bright color,
 with an alpha value matched to the one of the object's display color.
 */
template < typename T >
void bodyColorT(T const& obj)
{
    //assert_true(disp->color.transparent());
    const PointDisp * disp = obj.prop->disp;

    if ( disp->coloring )
    {
        unsigned long i = ( disp->coloring == 2 ? obj.mark() : obj.signature());
        gle_color col = gle::bright_color(i).match_a(disp->color);
        col.load();
        col.load_both();
    }
    else
    {
        disp->color.load();
        disp->color.load_both();
    }
}

/**
 if `coloring` is enabled, this loads the N-th bright color,
 otherwise load the object's display color
 */
inline void bodyColor(PointDisp const* disp, size_t s)
{
    if ( disp->coloring )
    {
        gle_color col = gle::bright_color(s);
        col.load();
        col.load_front();
        col.darken(0.5).load_back();
    }
    else
    {
        disp->color.load_load(1.0);
        disp->color2.load_back();
    }
}


/**
 if `coloring` is enabled, this loads the N-th bright color,
 otherwise load the object's display color
 */
inline void bodyColor2(PointDisp const* disp, size_t s)
{
    if ( disp->coloring )
        gle::bright_color(s).match_a(disp->color).load();
    else
        disp->color.load();
}


/**
 This is used for transparent objects.
 if `coloring` is enabled, this loads the N-th bright color,
 with an alpha value matched to the one of the object's display color.
 */
inline void bodyColorT(PointDisp const* disp, size_t s)
{
    //assert_true(disp->color.transparent());
    
    if ( disp->coloring )
    {
        gle_color col = gle::bright_color(s).match_a(disp->color);
        col.load();
        col.load_both();
    }
    else
    {
        disp->color.load();
        disp->color.load_both();
    }
}

