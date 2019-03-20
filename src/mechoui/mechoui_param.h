// Cytosim 3.0 -  Copyright Francois Nedelec et al.  EMBL 2015

#ifndef MECHOUI_PARAM_H
#define MECHOUI_PARAM_H

#include "real.h"
#include <string>
#include "gle_color.h"

class Glossary;

/// set of parameters for Mechoui
class MechouiParam
{
public:
    
    /// color of background
    gle_color    back_color;
    
    /// color of faces
    gle_color    face_color;
    
    /// color of points
    gle_color    point_color;
    
    /// size of points
    float        point_size;
    
    /// flag to display points
    int          point_style;
    
    /// selected cell
    int          selected;
    
    /// input file containing mesh
    std::string  file;
    
    /// patter for files
    std::string  dir;
    
    /// parameter
    std::string  config;
    
    /// refresh speed in milliseconds
    unsigned int delay;

    //----------------------------------------------------------
    
    /// constructor
    MechouiParam()  { clear(); }
    
    /// destructor
    ~MechouiParam() {}
    
    /// set default values
    void clear();
    
    /// read
    void read(Glossary& glos);

    /// write all values
    void write(std::ostream &) const;
};

#endif
