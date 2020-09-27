// Cytosim was created by Francois Nedelec.  Copyright 2020 Cambridge University.

#ifndef MECA1D_H
#define MECA1D_H

#include "array.h"
#include "mecable.h"
#include "sparmatsym1.h"
#include "monitor.h"
#include "allocator.h"
#include "bicgstab.h"
#include "simul.h"


/// Solves the motion of Objects along the X axis
/**
 This class is used to solve the motion of Mecables effectively in 1D, along the X axis.
 It works in 2D and 3D, but each Mecable is represented here by only one coordinate X.
 The Mecables are then translated in the X direction by the solution of the system.
 
 If you are curious about knowning how cytosim works, this is a good place to start!
 This is a bare-bone solver, which should be easy to understand.
 Meca does essentially the same work in higher dimension, with bending elasticity.
 */
class Meca1D
{
    size_t allocated_;           ///< allocated size of vectors
    
    int    ready_;               ///< true if the solution is contained in 'vSOL'
    
public:
   
    Array<Mecable *> mecables;   ///< list of mobile objects

    real * vSOL;                 ///< position of the points
    real * vBAS;                 ///< base points of forces and intermediate of calculus
    real * vMOB;                 ///< the mobility coefficients of the objects
    real * vRHS;                 ///< right-hand side term of the equation

    /// matrix containing the elasticity coefficients
    SparMatSym1   mA;
    
    /// working memory allocator
    LinearSolvers::Allocator allocator;

    Meca1D()
    {
        allocated_ = 0;
        ready_ = -1;
        vSOL = nullptr;
        vBAS = nullptr;
        vMOB = nullptr;
        vRHS = nullptr;
    }
    
    ~Meca1D()
    {
        //std::cerr << "~Meca1D\n";
        deallocate();
    }

    void deallocate()
    {
        free_real(vBAS);
        free_real(vSOL);
        free_real(vMOB);
        free_real(vRHS);
        allocated_ = 0;
        vBAS = nullptr;
        vSOL = nullptr;
        vMOB = nullptr;
        vRHS = nullptr;
    }

    void prepare(Simul const* sim, real time_step, real kT)
    {
        ready_ = 0;
        mecables.clear();
        
        for(Fiber * fib = sim->fibers.first(); fib; fib=fib->next())
            mecables.push_back(fib);

        size_t dim = mecables.size();
        if ( dim > allocated_ )
        {
            // make a multiple of chunk to align memory:
            allocated_ = chunk_real(dim);
            
            free_real(vBAS);
            free_real(vSOL);
            free_real(vMOB);
            free_real(vRHS);
            
            vBAS = new_real(allocated_);
            vSOL = new_real(allocated_);
            vMOB = new_real(allocated_);
            vRHS = new_real(allocated_);
        }
        
        mA.resize(dim);
        mA.reset();

        zero_real(dim, vBAS);
        zero_real(dim, vRHS);

        size_t ii = 0;
        for ( Mecable * mec : mecables )
        {
            mec->setIndex(ii);
            mec->setDragCoefficient();
            // Put the x coordinate of the origin in vSOL[ii]
            vSOL[ii] = mec->posPoint(0).XX;
            vMOB[ii] = time_step / mec->dragCoefficient();
            ++ii;
        }
    }
    
    /// add a clamp between point at index 'ii' to position 'dx'
    void addClamp(size_t ii, real w, real dx)
    {
        mA(ii, ii) -= w;
        vBAS[ii]   += w * dx;
    }
    
    /// add a link between points 'ii' and 'jj' with a position shift 'dx'
    void addLink(size_t ii, size_t jj, real w, real dx)
    {
        mA(ii, ii) -= w;
        mA(ii, jj) += w;
        mA(jj, jj) -= w;
        vBAS[ii] += w * dx;
        vBAS[jj] -= w * dx;
    }
    
    /// add Brownian motion and return smallest magnitude
    real setRightHandSide(real kT)
    {
        real res = INFINITY;
        for ( size_t ii = 0; ii < mecables.size(); ++ii )
        {
            real b = std::sqrt( 2 * kT * vMOB[ii] );
            vRHS[ii] = vMOB[ii] * vBAS[ii] + b * RNG.gauss();
            if ( b < res )
                res = b;
        }
        return res;
    }
    
    /// solve the system into 'vSOL', given 'vRHS' and matrix m
    /**
     Given that the force can be expressed as:
     
         F(POS) = mA * ( POS + vBAS )
     
     Where mA is a matrix and vBAS is a vector.
     This solves the discretized equation:

         ( newPOS - oldPOS ) / time_step = MOB * mA * ( newPOS + vBAS ) + Noise
  
     where MOB is a diagonal matrix of mobility coefficients (ie just a vector).
     Importantly, we used 'newPOS' on the right hand side, for implicit integration!
     We define
     
         vMOB = time_step * MOB
         vRHS = time_step * MOB * mA * ( oldPOS + vBAS ) + time_step * Noise
     
     leading to:
     
         ( I - vMOB * mA ) ( newPOS - oldPOS ) = vRHS
     
     We solve the linear system to determine:
     
         vSOL = newPOS - oldPOS

     */
    size_t solve(real precision)
    {
        assert_true(ready_==0);
        mA.prepareForMultiply(1);
        LinearSolvers::Monitor monitor(dimension(), precision);
        // solve without preconditionner
        LinearSolvers::BCGS(*this, vRHS, vSOL, monitor, allocator);
        ready_ = monitor.converged();
        return monitor.count();
    }
    
    /// apply translation to registered Mecables based on 'vSOL'
    void apply()
    {
        if ( ready_ )
        {
            size_t ii = 0;
            for ( Mecable * mec : mecables )
            {
                // Move the Mecable along the X direction as calculated
                mec->translate(Vector(vSOL[ii], 0, 0));
                ++ii;
            }
        }
    }
    
    /// Implements the LinearOperator
    size_t dimension() const { return mecables.size(); }
    
    /// Implements the LinearOperator:  Y <- X - vMOB * mA * X
    void multiply(const real * X, real * Y) const
    {
        assert_true( X != Y  &&  X != vBAS  &&  Y != vBAS );
        
        zero_real(mecables.size(), vBAS);
        mA.vecMulAdd(X, vBAS);
        
        for( size_t ii = 0; ii < mecables.size(); ++ii )
            Y[ii] = X[ii] - vMOB[ii] * vBAS[ii];
    }
};

#endif

