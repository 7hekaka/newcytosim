// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef SPACE_DYNAMIC_SPHERE_H
#define SPACE_DYNAMIC_SPHERE_H

#include "space_sphere.h"
#include "space_dynamic_prop.h"

/// A disc centered at the origin, with variable radius.
/**
 Space `dynamic_sphere` is a disc or a sphere centered around the origin
 Forces registered with 'setInteractions' are added, and used to update the
 radius of the Space. How fast the radius changes is set by the value 'mobility'
 in SpaceProp.
 
 Parameters:
     - radius

 @ingroup SpaceGroup
 
 FJN, Strasbourg 29.01.2017
 */

class SpaceDynamicSphere : public SpaceSphere
{
private:
    
    /// radial force
    mutable real force_;
    
public:
    
    /// constructor
    SpaceDynamicSphere(SpaceDynamicProp const*);

    /// Property
    SpaceDynamicProp const* prop() const { return static_cast<SpaceDynamicProp const*>(Space::prop); }

    /// add interactions to a Meca
    void setInteractions(Meca&, Simul const&) const;

    /// apply a force directed towards the edge of the Space
    void setConfinement(Vector const& pos, Mecapoint const&, Meca&, real stiff) const;

    /// apply a force directed towards the edge of the Space
    void setConfinement(Vector const& pos, Mecapoint const&, real rad, Meca&, real stiff) const;
    
    /// the step function can change the radius
    void step();

};

#endif

