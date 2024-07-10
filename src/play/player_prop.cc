// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University

#include "player_prop.h"
#include "glossary.h"
#include "save_image.h"


void PlayerProp::clear()
{
    replay = 0;
    loop = 0;
    period = 1;
    delay = 32;
    auto_exit = 0;
    auto_zoom = 1;
    auto_rotate.set(1,0,0,0);
    auto_pilot = 0;
    
    for ( int k = 0; k < NB_MAGIC_KEYS; ++k )
    {
        magic_key[k] = 0;
        magic_code[k] = "";
    }
    
    save_images = 0;
    downsample = 1;
    if ( SaveImage::supported("png") )
        image_format = "png";
    else
        image_format = "ppm";

    image_dir    = ".";
    image_name   = "image%04i";
    report_index = 0;
    image_index  = 0;
    poster_index = 0;
    saved_image_time = -1;
}


void PlayerProp::read(Glossary& glos)
{
    glos.set(replay, "play");
    glos.set(loop, "loop");
    if ( glos.set(period, "period") )
        period = std::max(1u, period);
    if ( glos.set(delay, "delay") )
        delay = std::max(2u, delay);
    glos.set(save_images,  "save_images", "save_image");
    glos.set(image_format, "image_format");
    glos.set(image_name,   "image_name");
    glos.set(image_dir,    "image_dir");
    glos.set(downsample,   "downsample");
    glos.set(report,       "report");
    glos.set(auto_exit,    "auto_exit");
    glos.set(auto_zoom,    "auto_zoom");
    glos.set(auto_rotate,  "auto_rotate");
    auto_pilot = ( auto_zoom != 1 || auto_rotate != Quaternion<real>(1,0,0,0) );
    
    if ( ! SaveImage::supported(image_format.c_str()) )
        throw InvalidParameter("unsupported image format");
    
    std::string var = "magic_key";
    for ( int k = 0; k < NB_MAGIC_KEYS; ++k )
    {
        glos.set(magic_key[k], var);
        glos.set(magic_code[k], var, 1);
        var = "magic_key" + std::to_string(k+1);
    }
}


void PlayerProp::write_values(std::ostream& os) const
{
    write_value(os, "play", replay);
    write_value(os, "loop",   loop);
    write_value(os, "period", period);
    write_value(os, "delay",  delay);
    write_value(os, "report", report);
    write_value(os, "save_images", save_images);
    write_value(os, "image_format", image_format);
    write_value(os, "image_dir", image_dir);
    write_value(os, "downsample", downsample);
    write_value(os, "auto_zoom", auto_zoom);
    write_value(os, "auto_rotate", auto_rotate);

    for ( int k = 0; k < NB_MAGIC_KEYS; ++k )
    {
        std::string var = "magic_key" + std::to_string(k+1);
        write_value(os, var, magic_key[k], "("+magic_code[k]+")");
    }
}

//------------------------------------------------------------------------------

static std::string standardReport(unsigned inx)
{
    switch( inx )
    {
        case 0: return "";
        case 1: return "fiber:lengths,fiber:mark";
        case 2: return "fiber:dynamics{split=1}";
        case 3: return "single,couple";
        case 4: return "single:force,couple:force";
        case 5: return "couple:configuration{split=1}";
        case 6: return "simul";
        case 7: return "simul:inventory";
        case 8: return "fiber:energy";
        case 9: return "space";
        case 10: return "platelet";
        case 11: return "fiber:mesh,field";
        case 12: return "fiber:lattice";
        case 13: return "fiber:cluster{couple=1}";
        case 14: return "fiber:age";
        case 15: return "fiber:distribution";
        case 16: return "fiber:distribution{interval=0.25}";
    }
    return "";
}


void PlayerProp::toggleReport(int alt)
{
    if ( alt == 0 )
    {
        report_index = 0;
        report = "";
    }
    else
    {
        report_index = ( report_index + alt + 18 ) % 18;
        report = standardReport(report_index);
    }
}

