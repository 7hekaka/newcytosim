// Cytosim 3.0 -  Copyright Francois Nedelec et al.  EMBL 2015

#include "mechoui_param.h"
#include "glossary.h"
#include <cmath>
#include <iomanip>


void MechouiParam::clear()
{
    back_color  = 0xDDDDDDFF;
    point_color = 0x000000FF;
    face_color  = 0xFFFFFF55;
    point_size  = 1;
    point_style = 0;
    file        = "";
    dir         = ".";
    config      = "";
    delay       = 250;
    selected    = 0;
}



void MechouiParam::read(Glossary& glos)
{
    glos.set(point_size,    "point_size");
    glos.set(back_color,    "back_color");
    glos.set(face_color,    "face_color");
    glos.set(point_color,   "point_color");
    glos.set(point_style,   "point_style");

    glos.set(file,          "file");
    glos.set(dir,           "dir");
    glos.set(config,        "config");
    glos.set(delay,         "delay");
}


/// formatted output of one line
template<typename T>
static  void write_param(std::ostream& os, std::string const& name, T const& c)
{
    os << " " << std::left << std::setw(20) << name << " = " << c << ";" << std::endl;
}



void MechouiParam::write(std::ostream& os) const
{
    write_param(os, "back_color",  back_color);
    write_param(os, "face_color",  face_color);
    write_param(os, "point_size",  point_size);
    write_param(os, "point_color", point_color);
    write_param(os, "point_style", point_style);
    write_param(os, "file",        file);
    write_param(os, "dir",         dir);
    write_param(os, "config",      config);
    write_param(os, "delay",       delay);
}


