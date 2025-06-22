// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include "fiber_disp.h"
#include "glossary.h"

#include "random_pcg.h"
using namespace PCG32;

// this controls compatibility in display parameter, which is not critical
#define BACKWARD_COMPATIBLE 1


void FiberDisp::clear()
{
    style      = 0;
    visible    = 1;
    color      = 0xFFFFFFFF;
    back_color = 0x777777FF;
    hide_color = 0xFFFFFF00;
    coloring   = 0;
    
    line_style = 1;
    line_width = 2;
    line_caps  = 1;
    bone_width = 1;
    outline_width = 0;
    
    point_style = 0;
    point_size  = 5;
    point_gap   = 1;

    end_style[0] = 0;
    end_style[1] = 0;
    end_size[0] = 6;
    end_size[1] = 6;
    growth_style = 0;
    
    end_colors[0] = 0xFFFFFFFF;  // white
    end_colors[1] = 0x00FF00FF;  // green
    end_colors[2] = 0xFFFF00FF;  // yellow
    end_colors[3] = 0xFF7538FF;  // orange
    end_colors[4] = 0xFF0000FF;  // red
    end_colors[5] = 0x777777FF;  // gray

    lattice_style   = 0;
    lattice_scale   = 1;
    lattice_rescale = 0;
    
    label_style   = 0;
    speckle_size  = 3;
    speckle_style = 0;
    speckle_gap   = 1;
    
    hide = 0;
    hide_axis.set(1,0,0);
    
    mask          = 0;
    mask_bitfield = 0;
    hide_state    = 7;  // not a possible state
    show_marked   = ~0U;
    
    length_scale  = 1;
    tension_scale = 1;
    
    force_style = 0;
    force_scale = 1;
    force_color = 0xFF00FFFF;
    
    draw_average  = 0;
}


void FiberDisp::read(Glossary& glos)
{
    glos.set(style, "style", {{"line", 0}, {"backbone", 1}, {"stripe", 2},
        {"striped", 2}, {"filament", 3}, {"actin", 4}, {"microtubule", 5}, {"wide", 8}});
    glos.set(visible, "visible");
    if ( glos.set(color, "color") )
        back_color = color.darken(0.625);
    glos.set(back_color, "color", 1);
    glos.set(hide_color, "color", 2);

    glos.set(back_color, "back_color");
    glos.set(hide_color, "hide_color");
    glos.set(coloring, "coloring");
    
    glos.set(line_width, "line_width", 0, "width", 0);
    std::string key = glos.has_key("line") ? "line" : "lines";
    glos.set(line_width, key);
    glos.set(line_style, "line_style", 0, key, 1, {{"off", 0}, {"line", 1}, {"tension", 2},
                                        {"rainbow_tension", 3}, {"curvature", 4}, {"orientation", 5},
                                        {"minus_end", 6}, {"plus_end", 7}, {"height", 8}, {"grid", 9}});
    glos.set(line_caps, "line_caps", 0, key, 2);
    glos.set(bone_width, "bone_width");
    glos.set(outline_width, "outlines");
    
    glos.set(point_size, "point_size", 0, "size", 0);
    key = glos.has_key("point") ? "point" : "points";
    if ( glos.set(point_size, key, 0) && !point_style )
        point_style = 1;
    glos.set(point_style, "point_style", 0, key, 1, {{"off", 0}, {"point", 1}, {"arrow", 2}, {"chevron", 3}, {"center", 4}});
    glos.set(point_gap, "point_gap", 0, key, 2);

    if ( point_gap <= 0 )
        point_gap = 1;

    if ( glos.set(end_size[0], "plus_end") )
        end_style[0] = 2;
    glos.set(end_style[0], "plus_end", 1, {{"off", 0}, {"sphere", 1}, {"cone", 2},
        {"cylinder", 3}, {"arrow", 4}, {"fins", 5}, {"cube", 6}, {"tip", 7}, {"hemisphere", 8}});
    
    if ( glos.set(end_size[1], "minus_end") )
        end_style[1] = 1;
    glos.set(end_style[1], "minus_end", 1, {{"off", 0}, {"sphere", 1}, {"cone", 2},
        {"cylinder", 3}, {"arrow", 4}, {"fins", 5}, {"cube", 6}, {"tip", 7}, {"hemisphere", 8}});
    
    glos.set(end_size,  2, "end_size");
    glos.set(end_style, 2, "end_style");
    glos.set(end_colors, 5, "end_color");
    
    glos.set(growth_style, "growth");
    
#if BACKWARD_COMPATIBLE
    glos.set(lattice_style, "draw_lattice");
    glos.set(lattice_scale, "lattice_max");
    glos.set(tension_scale, "tension");
#endif
    
    glos.set(lattice_style,   "lattice_style", 0, "lattice", 0);
    glos.set(lattice_scale,   "lattice_scale", 0, "lattice", 1);
    glos.set(lattice_rescale, "lattice", 2);

    glos.set(label_style, "label_style", 0, "labels", 0);

    key = glos.has_key("speckle") ? "speckle" : "speckles";
#if BACKWARD_COMPATIBLE
    if ( glos.num_values(key) == 2 )
    {
        speckle_size = line_width * 2;
        glos.set(speckle_style, key);
        glos.set(speckle_gap, key, 1);
    }
    else
#endif
    {
        glos.set(speckle_size,  "speckle_size", 0, key, 0);
        glos.set(speckle_style, "speckle_style", 0, key, 1, {{"off", 0}, {"random", 1}, {"regular", 2}});
        glos.set(speckle_gap,   "speckle_gap", 0, key, 2) || glos.set(speckle_gap, "interval");
    }
    if ( speckle_gap <= 0 )
        speckle_gap = 1;

    glos.set(hide, "hide", "exclude");
    glos.set(hide_axis, "hide_axis") || glos.set(hide_axis, "hide", 1, "exclude", 1);
    glos.set(hide_state, "hide_state");
    glos.set(show_marked, "show_marked");

    glos.set(mask, "mask");
    glos.set(mask_bitfield, "mask", 1);
    if ( mask && !mask_bitfield )
        mask_bitfield = distribute_bits(mask, pcg32_state);

    glos.set(length_scale,  "length_scale");
    glos.set(tension_scale, "tension_scale");

    glos.set(force_style, "forces", "force");
    glos.set(force_scale, "forces", 1, "force", 1);
    glos.set(force_color, "forces", 2, "force", 2);
    
    glos.set(draw_average, "draw_average");
}


void FiberDisp::write_values(std::ostream& os) const
{
    write_value(os, "style",        style);
    write_value(os, "visible",      visible);
    write_value(os, "color",        color, back_color, hide_color);
    write_value(os, "coloring",     coloring);
    
    write_value(os, "points",       point_size, point_style, point_gap);
    write_value(os, "lines",        line_width, line_style, line_caps);
    write_value(os, "bone_width",   bone_width);
    write_value(os, "outlines",     outline_width);
    write_value(os, "plus_end",     end_size[0], end_style[0]);
    write_value(os, "minus_end",    end_size[1], end_style[1]);
    write_value(os, "end_color",    end_colors, 6);
    write_value(os, "growth",       growth_style);
    
    write_value(os, "lattice",      lattice_style, lattice_scale, lattice_rescale);
    write_value(os, "labels",       label_style);
    write_value(os, "speckles",     speckle_size, speckle_style, speckle_gap);
    write_value(os, "hide",         hide, hide_axis);
    write_value(os, "hide_state",   hide_state);
    write_value(os, "show_marked",  show_marked);
    write_value(os, "mask",         mask, mask_bitfield);
    
    write_value(os, "length_scale", length_scale);
    write_value(os, "tension_scale",tension_scale);
    write_value(os, "forces",       force_style, force_scale, force_color);
    write_value(os, "draw_average", draw_average);
}

