// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University.

#ifndef PRIMITIVES_H
#define PRIMITIVES_H

#include <iostream>
#include "vector.h"
#include "isometry.h"

class Space;

namespace Cytosim
{
    /// read a position primitives, such as 'circle 5', etc.
    Vector readPositionPrimitive(std::string const&, size_t&, Space const*);
    
    /// read a direction primitives, such as 'horizontal', etc.
    Vector readDirectionPrimitive(std::string const&, size_t&, Vector const&, Space const*);

    /// read a position in space
    Vector readPosition(std::string const&, size_t&, Space const*);
    
    /// modify a position in space
    void modifyPosition(std::string const&, size_t&, Space const*, Vector&);

    /// convert string to a position
    Vector readPosition(std::string const&, Space const*);
    
    /// attempt to read a position in string, multiple times
    Vector findPosition(std::string const&, Space const*);
    
    /// read an orientation, and return a normalized vector
    Vector readDirection(std::string const&, size_t&, Vector const&, Space const*);
    
    /// modify a position in space
    void modifyDirection(std::string const&, size_t&, Space const*, Vector&);

    /// read an orientation, and return a normalized vector
    Vector readDirection(std::string const&, Vector const&, Space const*);

    /// read a rotation specified in stream
    Rotation readRotation(std::string const&, size_t&);
    
    /// read a rotation specified in string
    Rotation readRotation(std::string const&);

    /// read a rotation specified in stream, at position `pos`
    Rotation readOrientation(std::string const&, size_t&,  Vector const&, Space const*);
    
    /// skip space and new-line if `eat_line==true`, returning the next character
    int skip_space(std::string const&, size_t&, bool eat_line);

    /// true if string has significant material after 'pos'
    size_t has_trail(std::string const&, size_t pos);

}
#endif

