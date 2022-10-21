// Cytosim 3.0 -  Copyright Francois Nedelec et al.  EMBL 2015

#include "mechoui_param.h"
#include "glossary.h"
#include <cmath>
#include <iomanip>


void MechouiParam::clear()
{
    face_color  = 0xFFFFFF55;
    point_color = 0x0000FFFF;
    point_size  = 3;
    point_style = 0;
    face_style = 1;
    path     = ".";
    config   = "";
    delay    = 250;
    selected = 0;
}


void MechouiParam::read(Glossary& glos)
{
    glos.set(point_size,  "point_size");
    glos.set(face_color,  "face_color");
    glos.set(point_color, "point_color");
    glos.set(point_style, "point_style");
    glos.set(face_style, "face_style");
    
    glos.set(path,   "path");
    glos.set(config, "config");
    glos.set(delay,  "delay");
}


/// formatted output of one line
template<typename T>
static  void write_param(std::ostream& os, std::string const& name, T const& c)
{
    os << " " << std::left << std::setw(20) << name << " = " << c << ";\n";
}


void MechouiParam::write(std::ostream& os) const
{
    write_param(os, "face_color",  face_color);
    write_param(os, "point_size",  point_size);
    write_param(os, "point_color", point_color);
    write_param(os, "point_style", point_style);
    write_param(os, "face_style", face_style);
    write_param(os, "config",      config);
    write_param(os, "delay",       delay);
    os.flush();
}


