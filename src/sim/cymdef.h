// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.
/**
 @file
 @brief Global Compile Switches and common definitions: FiberEnd, etc.
*/

#ifndef CYMDEF_H
#define CYMDEF_H

/**
 Enable code to read trajectory files writen in older format.
 The value indicates the earliest file format that this aimed to be supported
 
 default = 50
 */
#define BACKWARD_COMPATIBILITY 50

/// Enables support for some advanced Space class
#define NEW_SPACES 1

/// Option to enable bending elasticity terms to link the two ends of a fiber
#define NEW_FIBER_LOOP 0


/// Designates the tip of a Fiber, but also the origin and center points
enum FiberEnd
{
    NO_END      = 0,   ///< not an end
    PLUS_END    = 1,   ///< highest abscissa = last vertex
    MINUS_END   = 2,   ///< lowest abscissa = fist vertex at index 0
    BOTH_ENDS   = 3,   ///< used to designate any of the two ends
    ORIGIN      = 4,   ///< refers to the origin of abscissa
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
    STATE_RED    = 4    ///<  Third dynamic state: usually shrinking
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
    CONFINE_IN_OUT        = 5,   ///< confine vertices with different stiffness inside/outside
    CONFINE_ALL_INSIDE    = 8,   ///< confine entire sphere inside the Space
    CONFINE_ALL_OUTSIDE   = 9,   ///< confine entire sphere outside the Space
    CONFINE_PLUS_END      = 16,  ///< confine plus end of fibers to the surface of the Space
    CONFINE_MINUS_END     = 17,  ///< confine minus end of fibers to the surface of the Space
    CONFINE_BOTH_ENDS     = 18,  ///< confine both ends of fibers to the surface of the Space
    CONFINE_MINUS_OUT     = 19,  ///< confine the plus end outside the Space
    CONFINE_PLUS_OUT      = 20,  ///< confine the plus end outside the Space
    CONFINE_POINT         = 24,  ///< confine first vertex of a Solid to the surface of the Space
    CONFINE_POINT_INSIDE  = 25,  ///< confine first vertex of a Solid, inside the Space
    CONFINE_RANGE         = 26   ///< confine vertices within a range of abscissa to the surface
};


#endif
