// Cytosim 3.0 - F. Nedelec and Laboratory, Copyright EMBL 2007

#ifndef SPACE_DYNAMIC_ELLIPSE_H
#define SPACE_DYNAMIC_ELLIPSE_H

#include "dim.h"
#include "space.h"
#include "space_ellipse.h"
#include "meca.h"

/// ellipse in 2D, ellipsoid or spheroid in 3D 
/**
 Space `ellipse` is aligned with the principal axes X, Y and Z.
 In 2D, there are two principal axis, and this is called an ellipse.
 In 3D, it is called an ellipsoid.
 
    ellipse sizeX sizeY sizeZ

 With:
 - sizeX = half length of X axis
 - sizeY = half length of Y axis
 - sizeZ = half length of Z axis
 .

 This Space has no corners: setInteraction() relies on project()
 */

class SpaceDynamicEllipse : public SpaceEllipse
{
private:
    
    /// Orientation matrix
    MatrixD  mat;
    
    /// Inverse of mat
    MatrixD  inv;
    
    /// pressure : lagrange parameter for volume conservation
    real     pressure;
    
    /// Forces from interactions
    Vector   inter_forces;
    
    /// Radial and circular forces
    mutable Vector Rforces;
    mutable Torque Torques;
    
    /// Reset forces
    void        reset_forces() const;
    
    /// Decompose forces between radial and circular components
    void        decompose_force(Vector const &, Vector const &, Vector const &) const;
    
    /// Add radial component
    void        add_radial_force(Vector const &, Vector const &) const;
    
    /// Find the optimal value of the pressure
    real        compute_pressure(Vector const &, Vector const&) const;
    
    /// Return forces corresponding to given pressure
    Vector      pressure_forces(real pressure) const;
    
    /// return normalized forces caused by surface tension
    Vector      tension_forces() const;
    
    /// Gives the i-th eigenvector of the ellipsoid
    Vector      director(unsigned i) const;
    /// Surface area of an ellipse of given axis length
    static real surfaceEllipse(const Vector &);
    
    /// Volume of an ellipse of given axis length
    static real volumeEllipse(const Vector &);
    
public:
     
    /// constructor
    SpaceDynamicEllipse(const SpaceProp*);
    
    /// add interactions to a Meca
    void    setInteractions(Meca &, FiberSet const&) const;

    /// setInteraction and changes the forces the ellipse undergoes
    void    setInteraction(Vector const &pos, Mecapoint const& pe, Meca & meca, real stiff) const;

    ///    ContractEllipse has a step function to adjust shape
    void    step();
    
    /// project point on the closest edge of the Space
    Vector  project(Vector const& pos) const
    {
        Vector p = inv * pos;
        return mat * SpaceEllipse::project(p);
    }
    
    /// true if the point is inside the Space 
    bool    inside(Vector const& pos) const
    {
        Vector p = inv * pos;
        return SpaceEllipse::inside(p);
    }

    ///
    void    resize();
    
    /// read from file
    void    read(Inputter& , Simul&, ObjectTag);
    
    /// write to file
    void    write(Outputter&) const;
    
    /// return forces
    void    report(std::ostream&) const;

    /// OpenGL display function; returns true if successful
    bool    draw() const;
};

#endif

