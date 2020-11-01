// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef SPACE_BEADS_H
#define SPACE_BEADS_H

#include "space.h"
#include "array.h"
#include "bead_prop.h"

class Bead;


/// a volume defined as the union of Beads
/** 
 Space `bead` is a volume defined as the union of Beads
 
 Parameters:
 
    bead = BEAD_NAME
 
 This code in inherited from an old (unfinished) mitotic spindle project.
 The implementation is quite incomplete, as only inside() works.
 */

class SpaceBeads : public Space
{    
    typedef Array<Bead*> BeadList;
    
    BeadProp *  mBeadProp;
    
    BeadList    mBeads;
    
    real        bbMin[3];
    real        bbMax[3];

    std::string mBeadName;
    
    /// set bounding box
    void        setBoundaries();
    /// set mBeads
    void        setBeads();
    
public:

    /// constructor
    SpaceBeads(SpaceProp const*);

    /// change dimensions
    void        resize(Glossary& opt);
 
    /// return bounding box in `inf` and `sup`
    void        boundaries(Vector& inf, Vector& sup) const;
    
    /// the volume inside
    real        volume() const;
    
    /// true if the point is inside the Space
    bool        inside(Vector const&) const;
    
    /// set `proj` as the point on the edge that is closest to `point`
    Vector      project(Vector const& pos) const;
    
    /// apply a force directed towards the edge of the Space
    void        setInteraction(Vector const& pos, Mecapoint const&, Meca&, real stiff) const;

    /// apply a force directed towards the edge of the Space
    void        setInteraction(Vector const& pos, Mecapoint const&, real rad, Meca&, real stiff) const;

    /// find the beads
    void        step();

};

#endif

