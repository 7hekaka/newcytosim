// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University
#ifndef SINGLE_PROP_H
#define SINGLE_PROP_H

#include "real.h"
#include "vector.h"
#include "property.h"
#include "hand_prop.h"

class Mecable;
class Single;
class Wrist;
class Space;

#define NEW_MOBILE_SINGLE 0
#define NEW_ANCHOR_STIFFNESS 1

/// Property for Single
/**
 @ingroup Properties
*/
class SingleProp : public Property
{
    friend class Single;
    friend class Wrist;
    friend class WristLong;

public:
    
    /**
     @defgroup SinglePar Parameters of Single
     @ingroup Parameters
     @{
     */
    
    
    /// name of Hand
    std::string hand;
    
    /// stiffness of link (pN/um)
    real stiffness;
    
#if NEW_ANCHOR_STIFFNESS
    /// A stiffness to anchor fibers on Solids with angular constraints
    real anchor_stiff;
#endif
    
    /// resting length of link (um)
    real length;
    
    /// diffusion coefficient
    real diffusion;

    /// if set > 0, assumes uniform concentration of diffusing Single
    /**
     The possible values for `fast_diffusion` are:
     - 0: disabled.
     Every Single is explicitly represented by a diffusing point. The
     unattached Single move randomly, with the specified diffusion constant,
     and they may accumulate at certain regions of the simulation. That
     happens in particular with asters when the motors move inward. In fact,
     the aster can really act like black hole for the motors, as seen in
     experiments!
     
     - 1: attach Single along the length of fibers
     It is assumed that free Singles are uniformly distributed, and they are
     thus not explicitly represented. In this mode, Cytosim will only keep a
     count of the number of free motors, and directly attach these Singles by
     one of their Hand, at random positions of the fibers, irrespective of the
     manner in which the filaments are distributed in the simulation.
     The diffusion constant is not relevant.
     
     - 2: attach Single only near growing PLUS_ENDS
     Similar to mode 1, but Hand1 of Single is directly attached
     at random positions on the growing ends of the fibers.
     
     - 3: attach Single only near the edge of the Space.
     .

     */
    int fast_diffusion;
    
    /// if > 0, the number of candidates for binding considered for `fast_diffusion`
    size_t fast_reservoir;

#if NEW_MOBILE_SINGLE
    /// constant drift
    Vector speed;
#endif
    
    /// Confinement can be `none`, `inside` (default) or `surface`
    Confinement confine;
    
    /// Unused Parameter: confinement stiffness (also known as `confine[1]`)
    real confine_stiff;
    
    /// name of space used for confinement (also known as `confine[2]`)
    std::string confine_label;
    
    /// specialization
    /**
     @copydetails SingleGroup
     */
    std::string activity;
    
    /// @}
    
    /// derived variable: Property of associated Hand
    HandProp * hand_prop;

protected:
    
    /// pointer to actual confinement Space, derived from `confine_label`
    Space const* confine_space;

    /// displacement in one time step
    real diffusion_dt;
    
#if NEW_MOBILE_SINGLE
    /// movement in one time step
    Vector speed_dt;
#endif

public:
    
    /// constructor
    SingleProp(const std::string& n) : Property(n)  { clear(); }

    /// destructor
    ~SingleProp() { }
    
    /// create a Single with this property
    virtual Single * newSingle() const;
    
    /// create a Single at given position
    Single * newSingle(Vector const& pos) const;

    /// create a Wrist attached to a Mecable vertex
    Wrist * newWrist(Mecable const*, size_t inx) const;
    
    /// create a Wrist attached to a point interpolated from Mecable's vertices
    Wrist * newWrist(Mecable const*, size_t ref, Vector const&) const;

    /// identifies the property
    std::string category() const { return "single"; }
    
    /// set default values
    void clear();
    
    /// set from a Glossary
    void read(Glossary&);
    
    /// compute derived parameter values
    void complete(Simul const&);

    /// return a carbon copy of object
    Property* clone() const { return new SingleProp(*this); }
    
    /// write all values
    void write_values(std::ostream&) const;

    /// return volume of confine_label
    real spaceVolume() const;

};

#endif

