// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SOLID_H
#define SOLID_H

#include "dim.h"
#include "array.h"
#include "object.h"
#include "mecable.h"
#include "matrix33.h"
#include "solid_prop.h"

class Meca;
class Fiber;
class SingleProp;
class Wrist;



/// Undeformable set of points
/**
 This is a Mecable behaving like a undeformable cloud of points.
 Each point can have its own radius and together they define the viscous drag
 of the Solid in the medium.
 
 \par Geometry:
 
 The ensemble can rotate and translate like a rigid body under external forces,
 but the relative configuration of the points in space is fixed: 
 the distance between any two points is constant.  

 A snapshot of the current geometry is saved in soShape[] by fixShape().
 This configuration is reapplied to the current points by reshape(),
 finding the best isometric transformation of soShape[] into the current configuration,
 thus maintaining the current position and orientation of the solid.
 
 \par Viscous Drag:
 
 The distance between the points, and their radii define a total drag
 coefficient according to Stokes' law applied to individual spheres.
 Points that have a radius = 0 do not induce viscous drag.
 The hydrodynamic interactions between the beads in the ensemble,
 and more advanced hydrodynamic effects are neglected.
 The drag coefficent for translation is simply the sum of Stokes' law,
 for all points that have a radius > 0.
 The rotational drag coefficient involves the second momentum of the configuration.
 
 \par Related classes:
 
 Solid is an extension of Bead. 
 A Solid with only one point is equivalent to a Bead, but slower to simulate.
*/
class Solid : public Mecable
{
private:
    
#if ( DIM > 2 )
    /// matrix containing the reduced momentum of inertia for 3D
    Matrix33 soMomentum;
#endif
    
    /// the dimensions used in Stokes' law to calculate overall mobility
    real * soRadius;
    
    /// array of coordinates for the reference shape of the solid
    real * soShape;
    
    /// the mean of the vertices weighted by their drag coefficients
    Vector soCenter;

    /// the reduced total (all points summed) drag coefficient for translation
    real soDrag;
    
    /// the reduced total drag coefficient for rotation
    real soDragRot;

    /// second moment of the reference shape
    real soVariance;
    
    /// the number of points when fixShape() was last called, used for verifications.
    SIZE_T soAmount;
    
#if NEW_SOLID_CLAMP
    /// position of clamp
    Vector clamp_place;
    
    /// stiffness of clamping force (known as `clamp[1]`)
    real clamp_stiff;
#endif

    /// reset private variables
    void reset();

    void reshape1D(real const*);
    void reshape2D(real const*);
    void reshape3D(real const*);
    
    /// part of build()
    index_t makePoint(ObjectList&, Glossary&, std::string const&, Simul&);
    
    /// part of build()
    index_t makeSphere(ObjectList&, Glossary&, std::string const&, Simul&);
    
    /// part of build()
    index_t makeBall(ObjectList&, Glossary&, std::string const&, Simul&);

    /// part of build()
    void makeWrist(ObjectList&, Glossary&, std::string const&, Simul&);
    
    /// part of build()
    Fiber* makeFiber(ObjectList&, Glossary&, std::string const&, Simul&);

    /// part of build()
    void addWrists(ObjectList&, index_t num, SingleProp const*, index_t ref, std::string const&);

public:
    
    /// Property
    SolidProp const* prop;
    
    /// default constructor
    Solid();

    /// constructor
    Solid(SolidProp const*);
    
    /// Copy constructor
    Solid(const Solid&);
    
    /// Assignement operator
    Solid& operator = (const Solid&);

    /// destructor
    virtual ~Solid();
    
    /// initialize according to options given in Glossary
    ObjectList build(Glossary&, Simul&);

    //------------------------------- Mecable ----------------------------------
    
    /// allocate memory
    void allocateMecable(index_t);
    
    /// free allocated memory
    void release();

    /// prepare for Meca
    void prepareMecable();

    /// sets the mobility
    void setDragCoefficient();

    /// total translation drag-coefficient (force = drag * speed)
    real dragCoefficient() const { return ( 6 * M_PI ) * prop->viscosity * sumRadius(); }
    
    /// The mobility of a model vertex ( speed = mobility * point_force )
    real pointMobility() const { return nbPoints() / dragCoefficient(); }
    
    /// Number of distance constraints applied to the movements of vertices
    index_t nbConstraints() const { return DIM * nPoints - ( DIM + (DIM-1) * (nPoints>1) ); }

    /// add the interactions due to confinement
    void setInteractions(Meca&) const;

    /// prepare for constrained projection
    void makeProjection();
    
    /// calculates the speed of points in Y, for the forces given in X
    void projectForces(const real* X, real* Y) const;
    
    /// calculates the speed of points in Y, for the forces given in X
    void projectForces0(const real* X, real* Y) const;

    /// add contribution of Brownian forces
    real addBrownianForces(real const* fce, real, real* rhs) const;
    
    /// Stochastic simulation
    void step();
    
    //------------------------------- Shaping ----------------------------------

    /// set the reference shape as a copy of the current one
    void fixShape();
    
    /// scale the reference shape
    void scaleShape(const real[DIM]);
    
    /// scale current shape to match the reference set in fixShape()
    void rescale();
    
    /// change coordinate values
    void getPoints(real const*);

    /// add a new point with a sphere (extends Mecable::addPoint)
    index_t addSphere(Vector const&, real radius);
    
    /// change radius of the sphere around point `i`
    void setRadius(index_t i, real val);

    /// add DIM points separated by `len`, to make a coordinate system around the last point
    index_t addTriad(real len);

    /// rotate the 3 points to align the diagonal with the X-axis
    void rotateTriad(index_t, Rotation const&);

    /// check if a coordinate system already exist at given index
    real hasTriad(index_t, real epsilon = 0.001) const;
    
    //--------------------------------------------------------------------------
    
    /// radius of the sphere around point `i`
    real radius(const index_t i) const { return abs_real(soRadius[i]); }
    
    /// sum of all sphere's radius
    real sumRadius() const;

    /// set index of Spheres that are nearest neighbors to `inx`; return number of values set
    index_t closestSpheres(index_t inx, index_t&, index_t&, index_t&) const;

    /// mean of all spheres weighted with their drag coefficients (or equivalently radius)
    Vector centroid() const;
    
    /// direction calculated from all the points
    Vector orientation() const;

    /// Position of center of gravity
    Vector position() const { return centroid(); }

#if NEW_SOLID_CLAMP
    /// returns clamping stiffness
    real clampStiffness() const { return clamp_stiff; }

    /// returns clamp position
    Vector clampPosition() const { return clamp_place; }
    
    /// set clamp stiffness and position
    void setClamp(Vector const& pos, real val) { clamp_place=pos; clamp_stiff=val; }

    /// force due to clamping
    Vector clampForce() const { return clamp_stiff * ( clamp_place - posPoint(0) ); }
#endif

    //--------------------------------------------------------------------------

    /// a static_cast<> of Object::next()
    Solid * next() const { return static_cast<Solid*>(next_); }
    
    /// a static_cast<> of Object::prev()
    Solid * prev() const { return static_cast<Solid*>(prev_); }
    
    //--------------------------------------------------------------------------

    /// a unique character identifying the class
    static const ObjectTag TAG = 'd';
    
    /// to store info relating the Solid
    static const ObjectTag SOLID_TAG = 'D';
    
    /// to store the clamp stiffness and position
    static const ObjectTag CLAMP_TAG = 'C';

    /// return unique character identifying the class
    ObjectTag tag() const { return TAG; }
    
    /// return associated Property
    Property const* property() const { return prop; }
    
    /// convert pointer to Solid* if the conversion seems valid; returns 0 otherwise
    static Solid* toSolid(Object * obj)
    {
        if ( obj  &&  obj->tag() == TAG )
            return static_cast<Solid*>(obj);
        return nullptr;
    }
    
    /// convert pointer to Solid* if the conversion seems valid; returns 0 otherwise
    static Solid const* toSolid(Object const* obj)
    {
        if ( obj  &&  obj->tag() == TAG )
            return static_cast<Solid const*>(obj);
        return nullptr;
    }

    //--------------------------------------------------------------------------

    /// read from file
    void read(Inputter&, Simul&, ObjectTag);
    
    /// write to file
    void write(Outputter&) const;

    /// Human friendly ouput
    void print(std::ostream&, bool write_shape = false) const;
};

/// output operator:
std::ostream& operator << (std::ostream& os, Solid const&);

#endif
