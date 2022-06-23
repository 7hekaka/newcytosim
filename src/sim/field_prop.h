// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

#ifndef FIELD_PROP_H
#define FIELD_PROP_H

#include "real.h"
#include "property.h"
#include "vector.h"

class Space;

/// Property for Field
/**
 @ingroup Properties
*/
class FieldProp : public Property
{
    friend class FieldSet;
    
public:
    
    /**
     @defgroup FieldPar Parameters of Field
     @ingroup Parameters
     @{
     */
    
    /// name of confining Space
    std::string field_space;

    /// size of square unit cell
    real step;
    
    /// use periodic boundary conditions for diffusion
    bool field_periodic;
    
    /// diffusion constant
    real slow_diffusion;
    
    /// diffusion constant
    real full_diffusion;
    
    /// type of boundary condition
    /*
     can be:
     - `flux`   (default) zero-flux at edges
     - `edge`   fixed values on all edges
     - `edgeX`  fixed values on X-edges
     - `edgeY`  fixed values on Y-edges
     - `edgeZ`  fixed values on Z-edges
     - `edgeXY` fixed values on XY-edges
     .
     */
    int boundary_condition;
    
    /// value of Dirichlet boundary conditions:
    real boundary_value;
    
    /// decay rate per unit time
    real decay_rate;
    
    /// coefficient for transport along Fiber (experimental)
    real transport_strength;
    
    /// length of average transport event (experimental)
    real transport_length;
    
    /// flag for cutting of Fiber
    int cut_fibers;
    
    /// flag for depolymerization at Fiber tips (experimental)
    int chew_fibers;
    
    /// if > 0, the simulation is stopped in any cell becomes negative
    /** The test is done just after writing the field to file */
    int positive;
    
    /// if false, the field is not recorded in the trajectory file
    bool save;
    
    
    /// a scale used to display the field
    real display_scale;
    
    /// a display flag
    int visible;
    /// @}
    
public:
    
    /// pointer to actual confinement Space, derived from `field_space`
    Space const*  field_space_ptr;

public:
    
    /// constructor
    FieldProp(const std::string& n) : Property(n)  { clear(); }

    /// destructor
    ~FieldProp() { }
    
    /// identifies the property
    std::string category() const { return "field"; }
        
    /// set default values
    void clear();
    
    /// set from a Glossary
    void read(Glossary&);
    
    /// derive parameters
    void complete(Simul const&);
    
    /// return a carbon copy of object
    Property* clone() const { return new FieldProp(*this); }

    /// write all values
    void write_values(std::ostream&) const;
    
};

#endif

