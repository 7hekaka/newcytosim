// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University

#include "player_prop.h"
#include "glossary.h"
#include "saveimage.h"


void PlayerProp::clear()
{
    play      = 0;
    loop      = 0;
    auto_exit = 0;
    period    = 1;
    delay     = 32;

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

    image_dir    = "";
    image_name   = "image";
    report_index = 0;
    image_index  = 0;
    poster_index = 0;
    saved_image  = -1;
    goLive = 0;
}


void PlayerProp::read(Glossary& glos)
{
    glos.set(play,         "play");
    glos.set(loop,         "loop");
    if ( glos.set(period,  "period") )
        period = std::max(1u, period);
    if ( glos.set(delay,   "delay") )
        delay = std::max(2u, delay);
    glos.set(save_images,  "save_images") || glos.set(save_images,  "save_image");
    glos.set(image_format, "image_format");
    glos.set(image_name,   "image_name");
    glos.set(image_dir,    "image_dir");
    glos.set(downsample,   "downsample");
    glos.set(report,       "report");
    glos.set(auto_exit,    "auto_exit");

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
    write_value(os, "play",   play);
    write_value(os, "loop",   loop);
    write_value(os, "period", period);
    write_value(os, "delay",  delay);
    write_value(os, "report", report);
    write_value(os, "save_images", save_images);
    write_value(os, "image_format", image_format);
    write_value(os, "image_dir", image_dir);
    write_value(os, "downsample", downsample);

    for ( int k = 0; k < NB_MAGIC_KEYS; ++k )
    {
        std::string var = "magic_key" + std::to_string(k);
        write_value(os, var, magic_key[k], "("+magic_code[k]+")");
    }
}

//------------------------------------------------------------------------------

std::string standardReport(unsigned inx)
{
    switch( inx )
    {
        case 0: return "";
        case 1: return "fiber:lengths";
        case 2: return "fiber:dynamics";
        case 3: return "single";
        case 4: return "single:force";
        case 5: return "couple";
        case 6: return "couple:force";
        case 7: return "couple:configuration";
        case 8: return "platelet";
        case 9: return "space";
    }
    return "";
}

std::string alternativeReport(unsigned inx)
{
    switch( inx )
    {
        case 0: return "";
        case 1: return "system:inventory";
        case 2: return "system";
        case 3: return "fiber:energy";
        case 4: return "fiber:mesh,field";
        case 5: return "fiber:lattice";
        case 6: return "fiber:segment";
        case 7: return "fiber:cluster{couple=1}";
        case 8: return "fiber:age";
        case 9: return "fiber:distribution";
    }
    return "";
}


void PlayerProp::toggleReport(bool alt)
{
    report_index = ( report_index + 1 ) % 10;
    
    if ( alt )
        report = alternativeReport(report_index);
    else
        report = standardReport(report_index);
}

