// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef LINEAR_SOLVERS_H
#define LINEAR_SOLVERS_H

#include "real.h"
#include "linear_operator.h"
#include "allocator.h"
#include "monitor.h"


/// Iterative methods to solve a system of linear equations
/**
 The linear system is defined by class LinearOperator.
 
 The iterative solver is monitored by class Monitor,
 where the desired convergence criteria can be specified,
 and that will also keep track of iteration counts.
 
 F. Nedelec, 22.03.2012
*/
namespace LinearSolvers
{
    
    /// Conjugate Gradient
    void CG(const LinearOperator&, const real* rhs, real* sol, Monitor&, Allocator&);
    
    /// Conjugate Gradient, with Preconditionning
    void CGP(const LinearOperator&, const real* rhs, real* sol, Monitor&, Allocator&);
    
    /// Bi-Conjugate Gradient
    void BCG(const LinearOperator&, const real* rhs, real* sol, Monitor&, Allocator&);
    
    /// Bi-Conjugate Gradient Stabilized
    void BCGS(const LinearOperator&, const real* rhs, real* sol, Monitor&, Allocator&);
    
    /// Bi-Conjugate Gradient Stabilized with Preconditionning
    void BCGSP(const LinearOperator&, const real* rhs, real* sol, Monitor&, Allocator&);
        
};

#endif

