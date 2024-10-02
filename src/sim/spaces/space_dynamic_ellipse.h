// Cytosim 3.0 - F. Nedelec and Laboratory, Copyright 2020 Cambridge University.

#ifndef SPACE_DYNAMIC_ELLIPSE_H
#define SPACE_DYNAMIC_ELLIPSE_H

#include "dim.h"
#include "space.h"
#include "space_ellipse.h"
#include "space_dynamic_prop.h"
#include "meca.h"


/// deformable ellipse in 2D, ellipsoid or spheroid in 3D
/**
 Space `ellipse` change its axis length and rotate
 
 Parameters:
     - length = total length in X, Y and Z directions
     .

 This Space has no corners: setConfinement() relies on project()
 */

class SpaceDynamicEllipse : public SpaceEllipse
{
private:
    
    /// Orientation matrix
    MatrixD mat;
    
    /// Inverse of mat
    MatrixD inv;
    
    /// pressure : lagrange parameter for volume conservation
    real pressure;
    
    /// external radial forces
    Vector rad_forces;
    
    /// Radial and circular forces
    mutable Vector Rforces;
    mutable Torque Torques;
    
    /// Reset forces
    void reset_forces() const;
    
    /// Decompose forces between radial and circular components
    void decompose_force(Vector const&, Vector const&, Vector const&) const;
    
    /// Add radial component
    void add_radial_force(Vector const&, Vector const&) const;
    
    /// Find the optimal value of the pressure
    real compute_pressure(Vector const&, Vector const&) const;
    
    /// Return forces corresponding to given pressure
    Vector pressure_forces(real pressure) const;
    
    /// return normalized forces caused by surface tension
    Vector tension_forces() const;
    
    /// Gives the i-th eigenvector of the ellipsoid
    Vector director(index_t i) const;
    
    /// Surface area of an ellipse of given axis length
    static real surfaceEllipse(Vector const&);
    
    /// Volume of an ellipse of given axis length
    static real volumeEllipse(Vector const&);
    
public:
     
    /// constructor
    SpaceDynamicEllipse(SpaceDynamicProp const*);
    
    /// Property
    SpaceDynamicProp const* prop() const { return static_cast<SpaceDynamicProp const*>(Space::prop); }

    /// change dimensions
    void resize(Glossary& opt);

    /// add interactions to a Meca
    void setInteractions(Meca&, Simul const&) const;

    /// this also register the forces the ellipse undergoes
    void setConfinement(Vector const&pos, Mecapoint const& mp, Meca& meca, real stiff) const;

    ///    ContractEllipse has a step function to adjust shape
    void step();
    
    /// project point on the closest edge of the Space
    Vector project(Vector const& pos) const
    {
        Vector p = inv * pos;
        return mat * SpaceEllipse::project(p);
    }
    
    /// true if the point is inside the Space 
    bool inside(Vector const& pos) const
    {
        Vector p = inv * pos;
        return SpaceEllipse::inside(p);
    }
    
    /// debug printout
    void dump(std::ostream&) const;
    
    /// descriptive printout
    void report(std::ostream&) const;
    
    /// write to file
    void write(Outputter&) const;

    /// read from file
    void read(Inputter&, Simul&, ObjectTag);

    
    /// OpenGL display function
    void draw3D() const;
    
    /// OpenGL display function
    void draw2D(float) const { draw3D(); }

};

#endif

