// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
/**
@file vector.h
 `Vector` is defined as an alias to Vector1, 2 or 3, according to DIM
 The types Vector1, Vector2 and Vector3 are always available
*/

#ifndef VECTOR_H
#define VECTOR_H

#include "dim.h"
#include "vector1.h"
#include "vector2.h"
#include "vector3.h"
#include "vector4.h"

#if ( DIM == 1 )

typedef Vector1 Vector;
typedef real    Torque;
constexpr Torque nullTorque(0.0);

/// helper function to normalize a 'Torque'
inline Torque normalize(Torque x) { return sign_real(x); }

#elif ( DIM == 2 )

typedef Vector2 Vector;
typedef real    Torque;
constexpr Torque nullTorque(0.0);

/// helper function to normalize a 'Torque'
inline Torque normalize(Torque x) { return sign_real(x); }

#elif ( DIM == 3 )

typedef Vector3 Vector;
typedef Vector3 Torque;
constexpr Torque nullTorque(0.0,0.0,0.0);

#endif

#endif

