// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef ORGANIZER_H
#define ORGANIZER_H

#include <vector>

#include "assert_macro.h"
#include "object.h"
#include "buddy.h"

class Meca;
class Simul;
class Mecable;
class Glossary;
class PointDisp;
class Display;
class Solid;
class Sphere;

/// An assemblage of Mecable
/** 
 The Organizer contains an Array of pointers of type `Mecable*`, allowing
 `setInteractions()` to create specific mechanical links and typically to arrange
 the Mecables relative to each other in Space following a certain pattern.
 The functions are implemented in derived classes: Bundle, Aster, Nucleus, etc.
*/
class Organizer: public Object, protected Buddy
{

protected:

    typedef std::vector<Mecable*> MecableList;
    
    /// list of Objects that are `organized`
    MecableList mObjects;
    
public:

    /// default constructor
    Organizer() { }
    
    /// destructor
    virtual ~Organizer();
    
    /// create all the Objects of the Organizer, and return in list
    virtual ObjectList build(Glossary&, Simul&) = 0;

    //--------------------------------------------------------------------------

    /// number of objects currently organized
    size_t nbOrganized() const { return mObjects.size(); }
    
    /// set number of objects
    void nbOrganized(size_t n) { mObjects.resize(n, nullptr); }
    
    /// return Mecable at index `n`
    Mecable * organized(size_t n) const { assert_true(n<mObjects.size()); return mObjects[n]; }
    
    /// add Mecable at end of list
    void grasp(Mecable *);

    /// add Mecable at given index
    void grasp(Mecable *, size_t);
    
    /// returns true if Mecables is one of the organized object
    bool check(Mecable const*) const;

    /// handles the disapearance of one of the organized object
    void goodbye(Buddy const*);
    
    /// add objects to Simul if they are not already linked
    void checkOrganized() const;

    //--------------------------------------------------------------------------

    /// Organizer cannot be moved
    int mobile() const { return 0; }
    
    /// return the center of gravity
    virtual Vector position() const;

    /// return the average of all vertices
    virtual Vector positionP(index_t) const;
/*
    /// move all associated objects
    void translate(Vector const& T);
    
    /// rotate all associated objects
    void rotate(Rotation const& T);
*/
    /// Stochastic simulation
    virtual void step() {}
    
    /// add interactions to a Meca
    virtual void setInteractions(Meca&) const {}
    
    /// sum the drag coefficient of all objects
    real sumDragCoefficient() const;
    
    /// number of links to be displayed using getLink()
    virtual index_t nbLinks() const { return 0; }

    /// retrieve end positions of link number `inx`, or returns zero if this link does not exist
    virtual bool getLink(index_t inx, Vector&, Vector&) const { return false; }
    
    /// object from that gives its display parameters
    virtual Solid * solid() const { return nullptr; }
    
    /// object from that gives its display parameters
    virtual Sphere * sphere() const { return nullptr; }

    //--------------------------------------------------------------------------
    
    /// character identifying each class
    static const ObjectTag ASTER_TAG = 'a';
    static const ObjectTag BUNDLE_TAG = 'u';
    static const ObjectTag NUCLEUS_TAG = 'n';
    static const ObjectTag FAKE_TAG = 'k';
    
    //--------------------------------------------------------------------------

    /// a static_cast<> of Object::next()
    Organizer * next() const { return static_cast<Organizer*>(next_); }
    
    /// a static_cast<> of Object::prev()
    Organizer * prev() const { return static_cast<Organizer*>(prev_); }
    
    //--------------------------------------------------------------------------
    
    /// write list of Mecables
    void writeOrganized(Outputter&) const;

    /// read list of Mecables
    void readOrganized(Inputter&, Simul&, size_t);

    /// read
    void read(Inputter&, Simul&, ObjectTag);
};


#endif
