// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef SPACE_DYNAMIC_ELLIPSE_PROP_H
#define SPACE_DYNAMIC_ELLIPSE_PROP_H

#include "dynamic_space_prop.h"
#include "property.h"

//class SpaceDynamicEllipse;
//class DynamicSpaceProp;
//class Space;


class SpaceDynamicEllipseProp : public DynamicSpaceProp 
{
    
    friend class SpaceDynamicEllipse;
    
public:
    
    // tension of the ellipse
    real    tension;
    
    // volume of the ellipse (mutable because changed by const method)
    mutable real    volume;
	
public:
    
	
    /// constructor
    SpaceDynamicEllipseProp(const std::string& n) : DynamicSpaceProp(n)  { clear(); }
    
    /// destructor
    ~SpaceDynamicEllipseProp() { }
    
    /// create a new, uninitialized, Space
    Space * newSpace() const;

	/// create a new Space according to specifications
	Space * newSpace(Glossary&) const;
        
    /// set default values
    void clear();
    
    /// set from a Glossary
    void read(Glossary&);
	
    /// check and derive more parameters
    void complete(Simul const&);
	
    /// write all values
    void write_values(std::ostream&) const;
    
	
};

#endif

