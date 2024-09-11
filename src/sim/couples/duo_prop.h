// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

#ifndef DUO_PROP_H
#define DUO_PROP_H

#include "vector2.h"
#include "couple_prop.h"
class Space;

// if enabled, a Duo can exert torque
#define NEW_DUO_HAS_TORQUE 1

/// Additional Property for Duo and DuoLong
/**
 @ingroup Properties
*/
class DuoProp : public CoupleProp
{
    
    friend class Duo;
    
public:
    
    /**
     @defgroup DuoPar Parameters of Duo
     @ingroup Parameters
     Inherits @ref CouplePar.
     @{
     */
    
    /// name of the Space inside which the Duo is activated
    std::string activation;

    /// rate of deactivation
    real deactivation_rate;
    
    /// type of deactivation: can lead to object deletion
    int deactivation_mode;
    
    /// if true, the deactivation clock runs at all time
    bool vulnerable;
    
#if NEW_DUO_HAS_TORQUE
    /// Resting angle in radian (set as `torque[0]`)
    real rest_angle;
    
    /// Stiffness of the angular link, in Torque per radians (pN.um/radian) (set as `torque[1]`)
    real angular_stiffness;
    
    /// Allow the angle to flip in 2D
    bool flip;
#endif
    
    /// @}

    /// last message from splash()
    mutable std::string splashed;

    /// deactivation_rate * time_step
    float deactivation_rate_dt;
    
    // Space inside which the Duo is activated
    Space const* activation_space;
    
#if NEW_DUO_HAS_TORQUE
    /// derived variable: [cos(angle), sin(angle)]
    Vector2 rest_dir;
#endif

public:
    
    /// constructor
    DuoProp(const std::string& n) : CoupleProp(n)  { clear(); }
    
    /// destructor
    ~DuoProp() { }
    
    /// return a Duo or a DuoLong with this property
    Couple * newCouple() const;
    
    /// set default values
    void clear();
    
    /// set from a Glossary
    void read(Glossary&);
    
    /// print some info
    void splash(std::ostream&) const;

    /// compute values derived from the parameters
    void complete(Simul const&);
    
    /// return a carbon copy of object
    Property* clone() const { return new DuoProp(*this); }

    /// write all values
    void write_values(std::ostream&) const;

};

#endif

