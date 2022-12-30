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
    Vector readPositionPrimitive(std::istream&, Space const*);
    
    /// read a direction primitives, such as 'horizontal', etc.
    Vector readDirectionPrimitive(std::istream&, Vector const&, Space const*);
    
    /// read a position in space
    Vector readPosition(std::istream&, Space const*);
    
    /// convert string to a position
    Vector readPosition(std::string const&, Space const*);
    
    /// read a position in space
    Vector readPosition(std::string const&);
    
    /// read an orientation, and return a normalized vector
    Vector readDirection(std::istream&, Vector const&, Space const*);
    
    /// read a rotation specified in stream
    Rotation readRotation(std::istream&);
    
    /// read a rotation specified in stream, at position `pos`
    Rotation readOrientation(std::istream&, Vector const&, Space const*);
    
}
#endif

