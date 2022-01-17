// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.
/**
 @file
 @brief Global Compile Switches and common definitions: FiberEnd, etc.
*/

#ifndef CYMDEF_H
#define CYMDEF_H

/**
 Enable code to be able to read old trajectory files.
 The value indicates the earliest file format that should still be readable.
 
 default = 50
 */
#define BACKWARD_COMPATIBILITY 50


// enable compact trajectory format created on 23/06/2021
/** This makes the trajectory files smaller, at the expense of precision
 This option is normally not desired, expect to compress old simulations files */
#define NEW_COMPACT_STORAGE 0

/// Enables support for some advanced Space class
#define NEW_SPACES 1


/// Option to not use constraints on the fiber's segment lengths
#define NEW_UNCONSTRAINED_LENGTH 0

/// Option to enable rigidity terms to link the two ends of a fiber
#define NEW_FIBER_LOOP 0


/// Designates the tip of a Fiber, but also the origin and center points
enum FiberEnd
{
    NO_END      = 0,   ///< not an end
    PLUS_END    = 1,   ///< highest abscissa = last vertex
    MINUS_END   = 2,   ///< lowest abscissa = fist vertex at index 0
    BOTH_ENDS   = 3,   ///< used to designate any of the two ends
    ORIGIN      = 7,   ///< refers to the origin of abscissa
    CENTER      = 8    ///< the mid-point between the two ends
};


/// Possible dynamic states for the tip of a Fiber [dynamic instability]
/**
 The naming is intentionally vague and does not refer to the nature of the states,
 since their actual interpretation may be be different in different types of Fiber.
 */
enum AssemblyState
{
    STATE_WHITE  = 0,   ///<  Used to indicate a non-dynamic end
    STATE_GREEN  = 1,   ///<  First dynamic state: usually growing
    STATE_YELLOW = 2,   ///<  Intermediate dynamic state
    STATE_ORANGE = 3,   ///<  Intermediate dynamic state
    STATE_RED    = 4,   ///<  Third dynamic state: usually shrinking
    STATE_BLACK  = 7    ///<  used by cutter to indicate deletion
};


/// used as function argument to define an AssemblyState
/** This is needed as ENUM are treated by default as signed int */
typedef unsigned state_t;


/// Possible modes of confinements
enum Confinement
{
    CONFINE_OFF           = 0,   ///< not confined
    CONFINE_INSIDE        = 1,   ///< confine vertices inside the Space
    CONFINE_OUTSIDE       = 2,   ///< confine vertices outside the Space
    CONFINE_ON            = 3,   ///< confine all vertices to the surface of the Space (always apply force)
    CONFINE_ALL_INSIDE    = 4,   ///< confine the entire sphere inside the Space
    CONFINE_ALL_OUTSIDE   = 5,   ///< confine the entire sphere outside the Space
    CONFINE_PLUS_END      = 10,  ///< confine PLUS_END of fibers to the surface of the Space
    CONFINE_MINUS_END     = 11,  ///< confine MINUS_END of fibers to the surface of the Space
    CONFINE_BOTH_ENDS     = 12,  ///< confine both ends of fibers to the surface of the Space
    CONFINE_PLUS_OUT      = 13,  ///< confine the PLUS_END outside the Space
    CONFINE_POINT         = 14,  ///< confine first vertex of a Solid to the surface of the Space
    CONFINE_POINT_INSIDE  = 15,  ///< confine first vertex of a Solid, inside the Space
    CONFINE_RANGE         = 17   ///< confine vertices within a range of abscissa to the surface
};


#endif
