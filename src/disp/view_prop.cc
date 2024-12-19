// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University
#include "view_prop.h"
#include "glossary.h"

//------------------------------------------------------------------------------
void ViewProp::clear()
{
    zoom         = 1;
    magnify      = 1;
    view_scale   = 10;
    auto_scale   = 1;
    focus.reset();
    focus_shift.reset();
    rotation.set(1,0,0,0);
    perspective  = 0;
    slice        = 0;
    
    back_color   = 0x000000FF;
    front_color  = 0xFFFFFFFF;
    
    buffered     = 1;
    depth_test   = 1;
    depth_clamp  = 0;
    retina       = 0;
    stencil      = 0;
    multisample  = 0;
    
    label        = "Cytosim";
    subtitle     = "";
    memo         = "Please, visit www.cytosim.org";
    draw_memo    = 0;

    track_fibers = 0;
    
    window_size[0] = 800;
    window_size[1] = 800;
    window_position[0] = 0;
    window_position[1] = 32;
    
    scalebar        = 0;
    scalebar_length = 10;
    scalebar_color  = 0xFFFF88AA;
    
    axes      = 0;
    axes_size = 1;

    for ( int k = 0; k < NB_CLIP_PLANES; ++k )
    {
        clip_plane_mode[k] = 0;
        clip_plane[k].set(1,0,0,0);
    }
    
    fog_type     = 0;
    fog_param    = 1;
    fog_color    = 0x000000FF;
    
    floor_radius = 0;
    floor_tile   = 1;
    floor_height = 0;
    floor_color  = 0x242424FF;
}


void ViewProp::invertColors()
{
    back_color = back_color.inverted();
    front_color = front_color.inverted();
    fog_color = fog_color.inverted();
    floor_color = floor_color.inverted();
}

void ViewProp::blackAndWhite()
{
    if ( back_color.brightness() < 0.5 )
    {
        back_color.set(0,0,0);
        front_color.set(1,1,1);
    }
    else
    {
        back_color.set(1,1,1);
        front_color.set(0,0,0);
    }
}


void ViewProp::read(Glossary& glos)
{
    glos.set(zoom,        "zoom");
    glos.set(magnify,     "magnify");
    if ( glos.set(view_scale, "view_scale") )
        auto_scale = 0;
    glos.set(auto_scale,  "auto_scale");
    glos.set(auto_scale,  "autoscale");
    glos.set(focus,       "focus");
    glos.set(rotation,    "rotation");
    rotation.normalize();

    glos.set(perspective, "perspective");
    glos.set(slice, "slice", {{"off", 0},{"front", 1},{"back", 2},{"slice", 3}});

    if ( glos.set(back_color, "back_color", "background_color") )
    {
        fog_color = back_color;
        front_color = back_color.inverted();
    }
    glos.set(front_color, "front_color");
    glos.set(buffered,    "buffered", "buffer");
    glos.set(depth_test,  "depth_test");
    glos.set(depth_clamp, "depth_clamp");
    glos.set(retina,      "retina");
    if ( glos.use_key("+") )
        retina = 1;
    glos.set(stencil,     "stencil");
    
    glos.set(multisample, "multisample");
    glos.set(multisample, "samples");
    // harmless backward compatibility
    glos.set(multisample, "gl_samples");
    glos.set(label,       "label");
    glos.set(draw_memo,   "draw_memo");
    glos.set(subtitle,    "subtitle");
    
    glos.set(track_fibers,        "track_fibers");
    glos.set(window_position, 2,  "window_position");
#if 1
    // A square window is made if only one value is given.
    if ( glos.set(window_size, 2, "window_size") == 1 )
        window_size[1] = window_size[0];

    // A square window is made if only one value is given.
    if ( glos.set(window_size, 2, "image_size") == 1 )
        window_size[1] = window_size[0];
#endif
    // 'size' is an alias to set the size of the window.
    if ( glos.set(window_size, 2, "size") == 1 )
        window_size[1] = window_size[0];
    /*
     window_size can be changed here, but we cannot resize
     the window, since we do not have access to GLUT 
     */
    //std::clog << this << " window " << window_size[0] << "x" << window_size[1] << '\n';
    
    glos.set(scalebar,        "scalebar");
    glos.set(scalebar_length, "scale_bar", 1);
    glos.set(scalebar_color,  "scale_bar", 2);

    glos.set(axes,      "axes");
    glos.set(axes_size, "axes", 1);

    for ( int k = 0; k < NB_CLIP_PLANES; ++k )
    {
        std::string var = "clip_plane" + std::to_string(k);
        glos.set(clip_plane_mode[k], var);
        Vector3 V(0,0,1);
        real S = 0;
        glos.set(V, var, 1);
        glos.set(S, var, 2);
        clip_plane[k].set(V.XX, V.YY, V.ZZ, S);
    }
    
    Glossary::dict_type<int> keys{{"off", 0}, {"linear", 1}, {"exponential", 2}, {"exponential2", 3}};

    glos.set(fog_type,  "fog_type", keys);
    glos.set(fog_param, "fog_param");
    glos.set(fog_color, "fog_color");
 
    glos.set(fog_type,  "fog", keys);
    glos.set(fog_param, "fog", 1);
    glos.set(fog_color, "fog", 2);
    
    glos.set(floor_radius, "floor");
    glos.set(floor_tile, "floor", 1);
    glos.set(floor_height, "floor", 2);
    glos.set(floor_color, "floor", 3);
}

//------------------------------------------------------------------------------

void ViewProp::write_values(std::ostream& os) const
{
    write_value(os, "zoom",          zoom);
    write_value(os, "magnify",       magnify);
    write_value(os, "view_scale",    view_scale);
    write_value(os, "auto_scale",    auto_scale);
    write_value(os, "focus",         focus+focus_shift);
    write_value(os, "rotation",      rotation);
    write_value(os, "perspective",   perspective);
    write_value(os, "slice",         slice);
    write_value(os, "back_color",    back_color);
    write_value(os, "front_color",   front_color);
    write_value(os, "buffered",      buffered);
    write_value(os, "depth_test",    depth_test);
    write_value(os, "depth_clamp",   depth_clamp);
    write_value(os, "stencil",       stencil);
    write_value(os, "multisample",   multisample);
    write_value(os, "label",         "("+label+")");
    write_value(os, "track_fibers",  track_fibers);
    //write_value(os, "window_position", window_position, 2);
    write_value(os, "window_size",   window_size, 2);
    write_value(os, "scalebar",      scalebar, scalebar_length, scalebar_color);
    write_value(os, "axes",          axes, axes_size);
    write_value(os, "fog",           fog_type, fog_param, fog_color);
    write_value(os, "floor",         floor_radius, floor_tile, floor_height, floor_color);
    
    for ( int k = 0; k < NB_CLIP_PLANES; ++k )
    {
        std::string var = "clip_plane" + std::to_string(k);
        write_value(os, var, clip_plane_mode[k], Vector3(clip_plane[k]), clip_plane[k].ZZ);
    }
}

