// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include "display_prop.h"
#include "gym_color_list.h"
#include "glossary.h"

//------------------------------------------------------------------------------
void DisplayProp::clear()
{
    style          = 1;
    tile           = 8;
    fold           = 1;
    draw_links     = true;

    couple_flip    = 0;
    couple_select  = 7;
    single_select  = 3;
    
    point_value    = 0;
    point_size     = 5;
    link_width     = 4;
    line_width     = 2;

    explode_style = 0;
    explode_range = 0;
}

//------------------------------------------------------------------------------
void DisplayProp::read(Glossary& glos)
{
    gym_color col;
    for ( unsigned i = 0; glos.set(col, "colors", i); ++i )
        gym::set_color(i, col);

    glos.set(style, "style", {{"fast", 1}, {"detailed", 2}, {"nice", 3}});
    glos.set(tile, "tile");
    glos.set(fold, "fold");
    glos.set(fold, "tile", 1);

    glos.set(tile, "tiled", "periodic");
    glos.set(fold, "tiled", 1, "periodic", 1);

    glos.set(draw_links, "draw_links");

    glos.set(couple_flip, "couple_flip");
    glos.set(couple_select, "couple_select");
    glos.set(single_select, "single_select");
    
    glos.set(point_value, "point_value");
    glos.set(point_size, "point_size");
    // unless specified, `link_width` will be equal to `line_width`:
    if ( glos.set(line_width, "line_width") )
        link_width = line_width;
    glos.set(link_width, "link_width", "link_size");
    
    glos.set(explode_style, "explode");
    glos.set(explode_range, "explode", 1);
    
#if BACKWARD_COMPATIBLE
    if ( glos.set(explode_range, "display_shift") )
        explode_style = 1;
#endif
}


//------------------------------------------------------------------------------

void DisplayProp::write_values(std::ostream& os) const
{
    write_value(os, "style",         style);
    write_value(os, "tile",          tile, fold);
    write_value(os, "draw_links",    draw_links);
    write_value(os, "couple_flip",   couple_flip);
    write_value(os, "couple_select", couple_select);
    write_value(os, "single_select", single_select);
    write_value(os, "point_value",   point_value);
    write_value(os, "point_size",    point_size);
    write_value(os, "link_width",    link_width);
    write_value(os, "line_width",    line_width);
    write_value(os, "explode",       explode_style, explode_range);
}


