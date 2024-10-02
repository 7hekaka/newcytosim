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
    index_t allocated_;          ///< allocated size of vectors

    real * vSOL;                 ///< position of the points
    real * vBAS;                 ///< base points of forces and intermediate of calculus
    real * vMOB;                 ///< the mobility coefficients of the objects
    real * vRHS;                 ///< right-hand side term of the equation
    
    Array<Mecable *> mecables;   ///< list of mobile objects

    SparMatSym1 matX;            ///< matrix containing the elasticity coefficients

    LinearSolvers::Allocator allocator_;  ///< working memory allocator
    
    int verbose_;
    int ready_;                  ///< true if the solution is contained in 'vSOL'

public:

    Meca1D()
    {
        allocated_ = 0;
        verbose_ = 1;
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

    void allocate(size_t dim)
    {
        if ( dim > allocated_ )
        {
            free_real(vBAS);
            free_real(vSOL);
            free_real(vMOB);
            free_real(vRHS);
            allocated_ = chunk_real(dim);
            vBAS = new_real(allocated_);
            vSOL = new_real(allocated_);
            vMOB = new_real(allocated_);
            vRHS = new_real(allocated_);
        }
    }
    
    void getReady(Simul const& sim)
    {
        ready_ = 0;
        mecables.clear();
        
        for(Fiber * fib = sim.fibers.first(); fib; fib=fib->next())
            mecables.push_back(fib);
        
        index_t dim = mecables.size();
        
        allocate(dim);
        matX.resize(dim);
        matX.reset();

        zero_real(dim, vBAS);
        zero_real(dim, vRHS);

        index_t inx = 0;
        for ( Mecable * mec : mecables )
        {
            mec->setIndex(inx);
            mec->setDragCoefficient();
            // Put the x coordinate of the origin in vSOL[inx]
            vSOL[inx] = mec->posPoint(0).XX;
            vMOB[inx] = sim.prop.time_step / mec->dragCoefficient();
            ++inx;
        }
    }
    
    /// clamp point at index 'i' to position 'pos' with weight w
    void addClamp(index_t i, real w, real pos)
    {
        matX(i, i) -= w;
        vBAS[i] += w * pos;
    }
    
    /// link points 'i' and 'j' with a position offset 'delta' and weight w
    void addLink(index_t i, index_t j, real w, real delta)
    {
        matX(i, i) -= w;
        matX(i, j) += w;
        matX(j, j) -= w;
        vBAS[i] += w * delta;
        vBAS[j] -= w * delta;
    }
    
    /// add Brownian motion and return smallest magnitude
    real setRightHandSide(real kT)
    {
        real res = INFINITY;
        for ( size_t i = 0; i < mecables.size(); ++i )
        {
            real b = std::sqrt( 2 * kT * vMOB[i] );
            vRHS[i] = vMOB[i] * vBAS[i] + b * RNG.gauss();
            res = std::min(res, b);
        }
        return res;
    }
    
    /// solve the system into 'vSOL', given 'vRHS' and matrix m
    /**
     Given that the force can be expressed as:
     
         F(POS) = MAT * ( POS + vBAS )
     
     Where MAT is a matrix and vBAS is a vector.
     This solves the discretized equation:

         ( newPOS - oldPOS ) / time_step = MOB * MAT * ( newPOS + vBAS ) + Noise
  
     where MOB is a diagonal matrix of mobility coefficients (ie just a vector).
     Importantly, we used 'newPOS' on the right hand side, for implicit integration!
     We define
     
         vMOB = time_step * MOB
         vRHS = time_step * MOB * MAT * ( oldPOS + vBAS ) + time_step * Noise
     
     leading to:
     
         ( I - vMOB * MAT ) ( newPOS - oldPOS ) = vRHS
     
     We solve the linear system to determine:
     
         vSOL = newPOS - oldPOS

     */
    index_t solve(real tol)
    {
        assert_true(ready_==0);
        matX.prepareForMultiply(1);
        LinearSolvers::Monitor monitor(2*dimension(), tol);
        // Solve linear system using Bi-Conjugate Gradient Stabilized Method:
        LinearSolvers::BCGS(*this, vRHS, vSOL, monitor, allocator_);
        ready_ = monitor.converged();
        if ( verbose_ )
        {
            std::stringstream oss;
            oss << "\tuniaxial size " << dimension();
            oss << " " << matX.what();
            oss << " count " << std::setw(4) << monitor.count();
            oss << " residual " << std::setw(11) << std::left << monitor.residual();
            if ( !ready_ ) oss << " FAILED!";
            Cytosim::out << oss.str() << std::endl;
            //std::clog << oss.str() << std::endl;
        }
        return monitor.count();
    }
    
    /// apply translation to registered Mecables based on 'vSOL'
    void apply()
    {
        if ( ready_ )
        {
            size_t i = 0;
            for ( Mecable * mec : mecables )
            {
                // Move the Mecable along the X direction as calculated
                mec->translate(Vector(vSOL[i], 0, 0));
                ++i;
            }
        }
    }
    
    /// Implements the LinearOperator
    size_t dimension() const { return mecables.size(); }
    
    /// Implements the LinearOperator:  Y <- X - vMOB * MAT * X
    void multiply(const real * X, real * Y) const
    {
        assert_true( X != Y  &&  X != vBAS  &&  Y != vBAS );
        
        zero_real(mecables.size(), vBAS);
        matX.vecMulAdd(X, vBAS);
        
        for( size_t i = 0; i < mecables.size(); ++i )
            Y[i] = X[i] - vMOB[i] * vBAS[i];
    }
};

#endif

