// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University
#ifndef SPHERICAL_CODE_H
#define SPHERICAL_CODE_H

#include "real.h"
#include <cstdio>

///\todo we could replace here the Coulomb repulsive interaction by a linear force
/* 
 Idea: A linear forces would allow us to solve an associated linear system on the
 coordinates of the points, using an iterative solver, that might be quite fast.
 The repulsive interaction only needs to take the first neighbours into account,
 so we could have local forces only, which would scale better than having all points
 interact, as is the case with the Coulomb energy.
*/

/// Distribute points on the unit sphere, minimizing the 'electrostatic' energy
/**  The number of points is arbitrary, see
http://mathworld.wolfram.com/SphericalCode.html

Algorithm:
 -# The points are distributed randomly on the sphere
 -# A 1/r^3 repulsive force is assumed for all points, 
 to which corresponds a certain potential energy in 1/r^2
 -# New positions are calculated form the current one, the forces
 and an adjustable scaling factor:  dx = scale * ( force at x )
 -# The potential energy of the new configuration is calculated, and:
     - If the total energy is lower, the move is accepted and the scaling factor is increased,
     - If the energy is higher, the move is rejected and the scaling factor is reduced.
     .
 .
 
The procedure (steps 2-4) is continues until convergence.\n

The main method is the class constructor, or equivalently distributePoints(),
which take the number of points as argument and performs the calculation.
The coordinates of points can then be retrieved using either:

    - copyPositionsForAllPoints()
    - copyCoordinatesOfPoint()

\author FJN, created August 2002, last modified Avril 2021
*/
class SphericalCode
{
public:
    
    /// a-priori expected distance between neighboring points, as a function of number of points
    static real expectedDistance(size_t);
    
    /// set coordinates point P randomly on the sphere
    static void randomize(real[3]);

    /// distribute point randomly
    static void randomize(size_t, real*);
    
    /// distribute point regularly
    static void distribute(size_t, real*);
    
    /// project W on the sphere
    static bool project(real P[3], const real W[3]);
    
    /// Calculate distance between point given their coordinates P and Q (3-dim)
    static real distance3(const real P[], const real Q[]);
    
    /// Calculate distance between point given their coordinates P and Q (3-dim)
    static real distance3Sqr(const real P[], const real Q[]);
    
    /// calculate Coulomb energy
    static real coulombEnergy(size_t, const real P[]);
    
    /// write points coordinates
    static void printPoints(size_t, real vec[], FILE* file = stdout);

    /// default constructor, does nothing
    SphericalCode();
    
    /// constructor that also calls distributePoints(),
    SphericalCode(size_t nbp);
    
    /// constructor that also calls distributePoints(), 
    SphericalCode(size_t nbp, real precision, size_t mx_nb_iterations);
    
    /// move point from old to new coordinates
    size_t refinePoints(real precision, size_t mx_nb_iterations);

    /// distribute the nbp points on the sphere and store their coordinates
    size_t distributePoints(size_t nbp, real precision, size_t mx_iterations);

    /// default destructor
    virtual ~SphericalCode();
    
    /// number of points in the configuration
    size_t nbPoints()  const { return num_points_;  }
    
    /// the 'virtual' total energy of the configuration
    real finalEnergy() const { return energy_; }
    
    /// minimum distance in the actual configuration, in 3D space
    real minimumDistance();
    
    /// multiply all coordinates by `factor`
    void scale(real factor);
    
    /// address where the coordinates for point `inx` are stored
    const real* addr(const size_t inx) const { return &coord_[3 * inx]; }
    
    /// copy the coordinates from point `inx` onto the given 3-dim array ptr
    void putPoint(real ptr[3], size_t inx) const;
    
    /// copy the coordinates from point `inx` onto x,y,z
    void putPoint(double* x, double* y, double* z, size_t inx) const;
    
    /// copy the coordinates from point `inx` onto x,y,z
    void putPoint(float* x, float* y, float* z, size_t inx) const;
    
    /// copy the points coordinates onto `x[]`, allocated to hold `sup` elements
    void putPoints(real ptr[], size_t sup) const;
    
private:

    /// This number affects convergence speed but not the result
    static constexpr unsigned SEVEN = 7U;
    
    /// number of point on the sphere
    size_t num_points_;
    
    /// coordinates of the points in a array
    /** in the array all the coordinates are together (x,y,z) point 1, (x,y,z) point 2, etc.
     so the coordinates for the first point are:
     x = coord_[0], y = coord_[1], z = coord_[3]
     the coordinates of point `i` are:
     x = coord_[3*i+0], y = coord_[3*i+1], z = coord_[3*i+2]
     */
    real* coord_;
    
    /// Coulomb energy of current configuration
    real energy_;
    
    /// calculate Coulomb forces
    void setForces(real forces[], real threshold);
    
    /// move point from old to new coordinates
    void movePoints(real Pnew[], const real Pold[], real forces[], real S);

    /// allocated memory
    void allocate(size_t);

};

#endif
