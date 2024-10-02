// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef FAKE_H
#define FAKE_H

#include "object.h"
#include "organizer.h"
#include "fake_prop.h"
#include "mecapoint.h"
#include "solid.h"

//------------------------------------------------------------------------------

///a set of two asters held together by a Solid 
/**
 This object cannot handle the destruction of the Asters
 
 The Fake should just link two Solid, without reference to the Asters
 
 @ingroup OrganizerGroup
 */
class Fake : public Organizer
{

private:

    /// Solid on which Fake is build
    Solid * fkSolid;
    
    /// Property
    FakeProp const* prop;
    
    /// connections
    std::vector<Mecapoint> asterPoints, solidPoints;

public:
    
    /// constructor
    Fake(FakeProp const* p) : fkSolid(nullptr), prop(p) { }
 
    /// destructor
    ~Fake();

    /// construct all the dependent Objects of the Organizer
    ObjectList build(Glossary&, Simul&);

    /// perform one Monte-Carlo step
    void step();
    
    /// add interactions to a Meca
    void setInteractions(Meca&) const;
    
    /// return pointer to central Solid
    Solid * solid() const { return fkSolid; }

    
    /// number of links to be displayed using getLink()
    index_t nbLinks() const { return (index_t)asterPoints.size(); }

    /// retrieve link between Solid and Aster's core
    bool getLink(index_t, Vector&, Vector&) const;
    
    /// return unique character identifying the class
    ObjectTag tag() const { return FAKE_TAG; }
    
    /// return associated Property
    Property const* property() const { return prop; }
    
    /// read from IO
    void read(Inputter&, Simul&, ObjectTag);
    
    /// write to IO
    void write(Outputter&) const;

 };


#endif

