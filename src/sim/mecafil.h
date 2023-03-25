// Cytosim was created by Francois Nedelec. Copyright Cambridge University 2020
#ifndef MECAFIL_H
#define MECAFIL_H

#include "chain.h"
#include "cymdef.h"  // needed for NEW_FIBER_LOOP

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
#define OLD_PROJECT_WITH_MATRIX 0

class Matrix;

/// incompressible Filament with bending elasticity
/**
 Implements the methods of a Mecable for the Chain:
 
 -# projectForces() includes longitudinal incompressibility,
 which means keeping successive points equidistants:
 norm( point(p+1) - point(p) ) = segmentation()
 
 -# addRigidity() implements bending elasticity.
 .
 
*/
class Mecafil : public Chain
{
private:
    
    /// normalized differences of successive vertices
    real * iDir;

    /// Lagrange multipliers associated with longitudinal imcompressibility
    real * iLag;
    
    /// work array allocated to hold `DIM*nbPoints` scalar values
    real * iLLG;
    
#if NEW_ANISOTROPIC_FIBER_DRAG
    /// local filament direction vectors used to calculate anisotropic drag
    real * iAni;
#endif

#if OLD_PROJECT_WITH_MATRIX
    
    /* variables used for projecting with a matrix ( mecafil_projectmat.cc ) */
    
    /// projection matrix
    real * iProj;
    
    /// differential of projection matrix
    real * iDProj;
    
    /// part of the projection matrix
    real * iJJtiJ;
    
    /// intermediate of calculus
    real * iTMP;

#else
    
    /// J*J' is a tridiagonal symmetric matrix of size (nbPoints-1).
    /** iJJt[] holds the diagonal elements and iJJtU[] the off-diagonal ones. */
    real * iJJt, * iJJtU;

#endif
#if ADD_PROJECTION_DIFF
    /// vector for the correction to the projection of size (nbPoints-1)
    real * iJJtJF;
    
    /// true if not all elements of iJJtJF[] are null
    bool useProjectionDiff;
#endif
protected:
    
    /// mobility of the points (all points have the same drag coefficient)
    real iPointMobility;
    
    /// rigidity scaling factor used in addRigidity()
    real iRigidity;
    
#if NEW_FIBER_LOOP
    /// link filament into a loop
    bool iRigidityLoop;
#endif
    
    /// calculate the normalized difference of successive vertices in iDir[]
    void storeDirections();

private:
    
    /// reset the memory pointers for the projection
    void initProjection();
    
    /// allocate memory for the projection
    void allocateProjection(size_t, real*);
    
    /// free the memory for the projection
    void destroyProjection();

public:
    
    /// Constructor
    Mecafil();
    
    /// copy constructor
    Mecafil(Mecafil const&);
    
    /// copy assignment
    Mecafil& operator = (Mecafil const&);
    
    /// Destructor
    virtual ~Mecafil();
    
    //------------------------------- Mecable ----------------------------------

    /// allocate memory
    void allocateMecable(size_t);
    
    /// free allocated memory
    void release();

    /// compute Lagrange multipliers associated with length constraints, given the force
    void computeTensions(const real* force);
    
    /// set Lagrange multipliers
    void setTensions(const real* ptr) { copy_real(nbSegments(), ptr, iLag); }
    
    /// debug output
    void printTensions(FILE *, char = ' ') const;
    
    /// replace current forces by the ones provided as argument, and compute tensions
    void getForces(const real* ptr);

    /// longitudinal force dipole between vertices `p` and `p+1`
    /**
     Tensions are calculated as the Lagrange multipliers associated with the constrains
     of distance between adjacent vertices, i.e. the segment length.
     This tension is:
     - positive when the segment is being pulled
     - negative when the segment is under compression
     .
     It is given in units of force (picoNewton, in Cytosim's standard units).
     */
    real tension(size_t p) const { assert_true(p+1<nPoints); return iLag[p]; }
    
    /// total drag-coefficient of object (force = drag * speed)
    real dragCoefficient() const { return nPoints / iPointMobility; }
    
    /// The mobility of a model vertex ( speed = mobility * point_force )
    real pointMobility() const { return iPointMobility; }

    /// drag coefficient of one point
    real leftoverMobility() const { return iPointMobility; }

    //--------------------- Projection  / Dynamics
    
    /// prepare for projection by building projection matrix
    void makeProjection();
    
#if NEW_ANISOTROPIC_FIBER_DRAG
    /// build projection matrix corresponding to current configuration
    void makeProjectionAnisotropic();
#endif
#if ADD_PROJECTION_DIFF
    /// select corrections to the projection, from the Lagrange multipliers
    void setProjectionDiff(real threshold);
    
    /// prepare the correction to the projection
    void makeProjectionDiff(const real* );
    
    /// add the contribution from the projection correction
    void addProjectionDiff(const real*, real*) const;
    
    /// add projection correction matrix
    void addProjectionDiff(real*) const;

    /// true if addProjectionDiff() does something
    bool hasProjectionDiff() const { return useProjectionDiff; }
#endif
    
    /// add displacements due to the Brownian motion to rhs[]
    real addBrownianForces(real const* rnd, real, real* rhs) const;

    /// calculate the speeds from the forces, including projection
    void projectForces(const real* X, real* Y) const;
    
    /// print projection matrix
    void printProjection(std::ostream&) const;

    //--------------------- Rigidity
    
    /// return fiber rigidity
    int hasRigidity() const { return ( nPoints > 2 ) & ( iRigidity > 0 ); }

    /// return fiber rigidity
    real fiberRigidity() const { return iRigidity; }

    /// add the rigidity force corresponding to configuration X into vector Y
    void addRigidity(const real* X, real* Y) const;

};


#endif
