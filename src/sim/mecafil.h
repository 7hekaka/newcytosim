// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef MECAFIL_H
#define MECAFIL_H

#include "chain.h"
#include "fiber_prop.h"  // needed for NEW_FIBER_LOOP

/**
 If the keyword below is defined, the viscous drag of the fibers
 will be different in the transverse and parallel directions, such that
 it will be 2x easier to move a fiber along it longitudinal direction.
 
 This is unpublished development, and you should set to zero
 */
#define NEW_ANISOTROPIC_FIBER_DRAG 0

/**
 Enable this option to build the projection matrix explicitly.
 Alternatively, the projection is calculated directly using vectors only.
 Having two methods is useful for cross-validation, but the matrix version is SLOWER
 
 Conclusion : do not enable this normally
*/
#define PROJECT_WITH_MATRIX 0

class Matrix;

/// incompressible Filament with bending elasticity
/**
 Implements the methods of a Mecable for the Chain:
 
 -# projectForces() includes longitudinal incompressibility,
 which means keeping successive points equidistants:
 norm( point(p+1) - point(p) ) = segmentation()
 
 -# addRigidity() implements bending elasticity.
 .
 
 \todo Rename Mecafil -> Filament
*/
class Mecafil : public Chain
{
private:
    
    /// normalized differences of successive vertices
    real   *    iDir;

    /// Lagrange multipliers associated with longitudinal imcompressibility
    real   *    iLag;
    
    /// work array allocated to hold `DIM*nbPoints` scalar values
    real   *    iLLG;
    
#if NEW_ANISOTROPIC_FIBER_DRAG
    /// local filament direction vectors used to calculate anisotropic drag
    real   *    iAni;
#endif

#if PROJECT_WITH_MATRIX
    
    /* variables used for projecting with a matrix ( mecafil_projectmat.cc ) */
    
    /// projection matrix
    real   *    iProj;
    
    /// differential of projection matrix
    real   *    iDProj;
    
    /// part of the projection matrix
    real   *    iJJtiJ;
    
    /// intermediate of calculus
    real   *    iTMP;

#else
    
    /// J*J' is a tridiagonal symmetric matrix of size (nbPoints-1).
    /** iJJt[] holds the diagonal elements and iJJtU[] the off-diagonal ones. */
    real   *    iJJt, * iJJtU;

#endif
    
    /// vector for the projection correction of size nbSegments
    real   *    iJJtiJforce;
    
    /// true if all elements of iJJtiJforce[] are null
    bool        useProjectionDiff;
    
protected:
    
    /// mobility of the points (all points have the same drag coefficient)
    real        iPointMobility;
    
    /// rigidity scaling factor used in addRigidity()
    real        iRigidity;
    
#if NEW_FIBER_LOOP
    /// link filament into a loop
    bool        iRigidityLoop;
#endif
    
    /// calculate the normalized difference of successive vertices in iDir[]
    void        storeDirections();

private:
    
    /// reset the memory pointers for the projection
    void        buildProjection();
    
    /// allocate memory for the projection
    void        allocateProjection(size_t);
    
    /// free the memory for the projection
    void        destroyProjection();

public:
    
    /// Constructor
    Mecafil();
    
    /// copy constructor
    Mecafil(Mecafil const&);
    
    /// copy assignment
    Mecafil& operator=(Mecafil const&);
    
    /// Destructor
    virtual    ~Mecafil();
    
    /// allocate memory
    size_t      allocateMecable(size_t);
    
    /// free allocated memory
    void        release();

    /// compute longitudinal tensions in the segments
    void        computeTensions(const real* force);
    
    /// copy Lagrange multipliers computed in projectForces()
    void        storeTensions(const real* force);
    
    /// debug output
    void        printTensions(std::ostream&) const;

    /// longitudinal force dipole between vertices `p` and `p+1`
    /**
     Tensions are calculated as the Lagrange multipliers associated with the
     distance between neigboring vertices, i.e. the fiber segment's lengths.
     This tension is:
     - positive when the segment is being pulled
     - negative when the segment is under compression
     .
     It is given in units of force (pico-Newton, if all quantitites use our units).
     */
    real        tension(size_t p) const { assert_true(p+1<nPoints); return iLag[p]; }
    
    /// total drag-coefficient of object (force = drag * speed)
    real        dragCoefficient() const { return nPoints / iPointMobility; }
    
    /// drag coefficient of one point
    real        leftoverMobility() const { return iPointMobility; }
    
    //--------------------- Projection  / Dynamics
    
    /// prepare for projection
    void        makeProjection();

    /// prepare the correction to the projection
    void        makeProjectionDiff(const real* );
    
    /// add the contribution from the projection correction
    void        addProjectionDiff(const real*, real*) const;
    
    /// true if addProjectionDiff() does something
    bool        hasProjectionDiff() const { return useProjectionDiff; }

    /// add displacements due to the Brownian motion to rhs[]
    real        addBrownianForces(real const* rnd, real alpha, real* rhs) const;

    /// calculate the speeds from the forces, including projection
    void        projectForces(const real* X, real* Y) const;
    
    /// print projection matrix
    void        printProjection(std::ostream&) const;

    //--------------------- Rigidity
    
    /// return fiber rigidity
    int         hasRigidity() const { return ( nPoints > 2 ) & ( iRigidity > 0 ); }

    /// return fiber rigidity
    real        fiberRigidity() const { return iRigidity; }

    /// add the rigidity force corresponding to configuration X into vector Y
    void        addRigidity(const real* X, real* Y) const;

};


#endif
