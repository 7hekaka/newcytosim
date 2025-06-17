// Cytosim was created by Francois Nedelec.  Copyright 2023 Cambridge University.

#ifndef BICGSTAB_H
#define BICGSTAB_H

#include "real.h"
#include "blas.h"
#include "cytoblas.h"
#include "allocator.h"
#include "monitor.h"


#ifndef __FAST_MATH__
[[ maybe_unused ]]
static void check_numbers(real const* vec, size_t len, size_t cnt)
{
    size_t nan = 0;
    for ( int i = 0; i < len; ++i )
        nan += std::isnan(vec[i]);
    if ( nan )
        fprintf(stderr, "BCGS %lu : %lu NaNs\n", cnt, nan);
}
#endif


/// using BLAS requires functions calls, and more loops are required.
#define BICGSTAB_USES_BLAS 0

/// assumes that vector 'sol' is zero initially
#define BICGSTAB_ZERO_INITIALIZED 0

/// Bi-Conjugate Gradient Stabilized method to solve a system of linear equations
/**
 F. Nedelec, 27.03.2012 - 13.03.2017
*/
namespace LinearSolvers
{
    /// Bi-Conjugate Gradient Stabilized without Preconditionning
    /*
     This solves `mat * S = rhs` with a tolerance specified in 'monitor'
     */
    template < typename LinearOperator, typename Monitor, typename Allocator >
    void BCGS(const LinearOperator& mat, const real* rhs, real* S,
              Monitor& monitor, Allocator& allocator)
    {
        double rho = 1.0, rho_old = 1.0, alpha = 0.0, beta = 0.0, omega = 1.0;
        double upsilon = 0.0;
        
        const int dim = mat.dimension();
        allocator.allocate(dim, 5);
        real * R  = allocator.bind(0);
        real * R0 = allocator.bind(1);
        real * P  = allocator.bind(2);
        real * T  = allocator.bind(3);
        real * V  = allocator.bind(4);

#if BICGSTAB_ZERO_INITIALIZED
        #pragma ivdep
        for ( int i = 0; i < dim; ++i )
        {
            S[i] = 0;
            real x = rhs[i];
            R0[i] = x;
            R[i] = x;
            P[i] = x;
        }
#else
        mat.multiply(S, R0);                  // r0 = MAT * s
#if BICGSTAB_USES_BLAS
        blas::xcopy(dim, rhs, 1, R, 1);       // r = rhs
        blas::xaxpy(dim, -1.0, R0, 1, R, 1);  // r = rhs - MAT * s
        if ( monitor.finished(dim, R) )
            return;
        blas::xcopy(dim, R, 1, R0, 1);        // r0 = r
        blas::xcopy(dim, R, 1, P, 1);
#else
        #pragma ivdep
        for ( int i = 0; i < dim; ++i )
        {
            real x = rhs[i] - R0[i];
            R0[i] = x;
            R[i] = x;
            P[i] = x;
        }
        if ( monitor.finished(dim, R) )
            return;
#endif
#endif
        rho = blas::dot(dim, R, R);
#if 0
        if ( rho != rho )
        {
            fprintf(stderr, "BCGSP called with invalid RHS argument\n");
            return;
        }
#endif
        goto start;
        
        while ( ! monitor.finished(dim, R) )
        {
            rho_old = rho;
            rho = blas::dot(dim, R0, R);
            
            if ( rho == 0.0 )
            {
#if ( 1 )
                /* The residual vector became nearly orthogonal to the
                 arbitrarily chosen direction r0, and we restart with a new r0 */
                blas::xcopy(dim, rhs, 1, R, 1);       // r = rhs
                mat.multiply(S, R0);                  // r0 = A*x
                blas::xaxpy(dim, -1.0, R0, 1, R, 1);  // r = rhs - MAT * x
                blas::xcopy(dim, R, 1, R0, 1);        // r0 = r
                rho = blas::dot(dim, R0, R0);
#else
                monitor.finish(2, dim, R);
                break;
#endif
            }
            
            beta = ( rho / rho_old ) * ( alpha / omega );
            // p = r + beta * ( p - omega * v )
#if BICGSTAB_USES_BLAS
            blas::xaxpy(dim, -omega, V, 1, P, 1);  // p = p - omega * v
            blas::xpay(dim, R, beta, P);           // p = r + beta * p
#else
            upsilon = beta * omega;
            for ( int i = 0; i < dim; ++i )
                P[i] = R[i] + beta * P[i] - upsilon * V[i];
#endif

        start:
            
            mat.multiply(P, V);                   // v = MAT * p;
            alpha = rho / blas::dot(dim, R0, V);

#if BICGSTAB_USES_BLAS
            blas::xaxpy(dim, -alpha, V, 1, R, 1); // r = r - alpha * v;
            blas::xaxpy(dim,  alpha, P, 1, S, 1); // s = s + alpha * p;
#else
            #pragma ivdep
            for ( int i = 0; i < dim; ++i )
            {
                R[i] -= alpha * V[i];
                S[i] += alpha * P[i];
            }
#endif

            //if ( monitor.finished(dim, R) )
            //    break;
            
            mat.multiply(R, T);                   // t = MAT * r;
            monitor += 2;

            double tdt = blas::dot(dim, T, T);
            
            if ( tdt > 0.0 )
            {
                omega = blas::dot(dim, T, R) / tdt;
                
                if ( omega == 0.0 )
                {
                    monitor.finish(3, dim, R);
                    break;
                }
#if BICGSTAB_USES_BLAS
                blas::xaxpy(dim,  omega, R, 1, S, 1); // s = s + omega * r
                blas::xaxpy(dim, -omega, T, 1, R, 1); // r = r - omega * t
#else
                #pragma ivdep
                for ( int i = 0; i < dim; ++i )
                {
                    S[i] += omega * R[i];
                    R[i] -= omega * T[i];
                }
#endif
            }
            else
                omega = 0.0;
        }
        allocator.release();
    }
    
    
    /// Bi-Conjugate Gradient Stabilized with right-sided Preconditionning
    /*
     This solves `mat * S = rhs` with a tolerance specified in 'monitor'
     */
    template < typename LinearOperator, typename Monitor, typename Allocator >
    void BCGSP(const LinearOperator& mat, const real* rhs, real* S,
               Monitor& monitor, Allocator& allocator)
    {
        double rho = 1.0, rho_old = 1.0, alpha = 0.0, beta = 0.0, omega = 1.0;
        double delta = 0.0, upsilon = 0.0;
        
        const int dim = mat.dimension();
        allocator.allocate(dim, 7);
        real * R    = allocator.bind(0);
        real * R0   = allocator.bind(1);
        real * P    = allocator.bind(2);
        real * T    = allocator.bind(3);
        real * V    = allocator.bind(4);
        real * Phat = allocator.bind(5);
        real * Shat = allocator.bind(6);

#if BICGSTAB_ZERO_INITIALIZED
        for ( int i = 0; i < dim; ++i )
        {
            S[i] = 0;
            real x = rhs[i];
            R0[i] = x;
            R[i] = x;
            P[i] = x;
        }
#else
        mat.multiply(S, R0);                  // r0 = MAT * s
#if BICGSTAB_USES_BLAS
        blas::xcopy(dim, rhs, 1, R, 1);       // r = rhs
        blas::xaxpy(dim, -1.0, R0, 1, R, 1);  // r = rhs - r0 = rhs - MAT * s
        if ( monitor.finished(dim, R) )
            return;
        blas::xcopy(dim, R, 1, R0, 1);        // r0 = r
        blas::xcopy(dim, R, 1, P, 1);
#else
        #pragma ivdep
        for ( int i = 0; i < dim; ++i )
        {
            real x = rhs[i] - R0[i];
            R0[i] = x;
            R[i] = x;
            P[i] = x;
        }
        if ( monitor.finished(dim, R) )
            return;
#endif
#endif
        rho = blas::dot(dim, R, R);
        goto start;

        while ( ! monitor.finished(dim, R) )
        {
            //fprintf(stderr, "BCGSP %4i res %14.9f\n", monitor.count(), monitor.residual());
            rho_old = rho;
            rho = blas::dot(dim, R0, R);
            
            if ( rho == 0.0 )
            {
#if ( 1 )
                /* The residual vector became nearly orthogonal to the
                 arbitrarily chosen direction r0, and we restart with a new r0 */
                mat.multiply(S, R);
                blas::xaxpy(dim, -1.0, rhs, 1, R, 1); // r = rhs - MAT * x
                blas::xcopy(dim, R, 1, R0, 1);        // r0 = r
                rho = blas::dot(dim, R0, R0);
#else
                monitor.finish(2, dim, R);
                break;
#endif
            }
            
            beta = ( rho * alpha ) / ( rho_old * omega );
            // p = r + beta * ( p - omega * v )
#if BICGSTAB_USES_BLAS
            blas::xaxpy(dim, -omega, V, 1, P, 1);  // p = p - omega * v
            blas::xpay(dim, R, beta, P);           // p = r + beta * p
#else
            upsilon = beta * omega;
            #pragma ivdep
            for ( int i = 0; i < dim; ++i )
                P[i] = R[i] + beta * P[i] - upsilon * V[i];
#endif
        start:
            
            mat.precondition(P, Phat);             // phat = PC * p;
            mat.multiply(Phat, V);                 // v = MAT * PC * p;
            
            delta = blas::dot(dim, R0, V);
            if ( delta == 0.0 )
            {
                ++monitor;
                monitor.finish(4, dim, R);
                break;
            }
            
            alpha = rho / delta;
#if BICGSTAB_USES_BLAS
            blas::xaxpy(dim, -alpha,    V, 1, R, 1); // r = r - alpha * v;
            blas::xaxpy(dim,  alpha, Phat, 1, S, 1); // s = s + alpha * phat;
#else
            #pragma ivdep
            for ( int i = 0; i < dim; ++i )
            {
                R[i] -= alpha * V[i];
                S[i] += alpha * Phat[i];
            }
#endif
            mat.precondition(R, Shat);             // shat = PC * r
            mat.multiply(Shat, T);                 // t = MAT * PC * r

            monitor += 2;

            double tdt = blas::dot(dim, T, T);
            
            if ( tdt > 0.0 )
            {
                omega = blas::dot(dim, T, R) / tdt;
            
                if ( omega == 0.0 )
                {
                    monitor.finish(3, dim, R);
                    break;
                }
#if BICGSTAB_USES_BLAS
                blas::xaxpy(dim,  omega, Shat, 1, S, 1); // s = s + omega * shat
                blas::xaxpy(dim, -omega,    T, 1, R, 1); // r = r - omega * t
#else
                #pragma ivdep
                for ( int i = 0; i < dim; ++i )
                {
                    S[i] += omega * Shat[i];
                    R[i] -= omega * T[i];
                }
#endif
            }
            else
                omega = 0.0;
        }
#if 0
        /* recalculate the true residual r0 = rhs - MAT * s */
        real est = monitor.residual();
        mat.multiply(S, R);
        blas::xaxpy(dim, -1.0, rhs, 1, R, 1);
        monitor.finished(dim, R);
        fprintf(stderr, "[BCGSP %4i res %14.9f %14.9f]", monitor.count(), est, monitor.residual());
#endif
        //fprintf(stderr, "BCGSP %4i res %14.9f - end\n", monitor.count(), monitor.residual());
        allocator.release();
    }
    
    
    /**
     This is an alternative implementation adapted from the CUSP project
     https://cusplibrary.github.io/index.html
     */
    template < typename LinearOperator, typename Monitor, typename Allocator >
    void bicgstab(const LinearOperator& mat, const real* rhs, real* sol,
                  Monitor& monitor, Allocator& allocator)
    {
        const size_t dim = mat.dimension();
        
        allocator.allocate(dim, 8);
        real * p   = allocator.bind(0);
        real * r   = allocator.bind(1);
        real * r0  = allocator.bind(2);
        real * s   = allocator.bind(3);
        real * Mp  = allocator.bind(4);
        real * AMp = allocator.bind(5);
        real * Ms  = allocator.bind(6);
        real * AMs = allocator.bind(7);
        
        // r <- A*x
        mat.multiply(sol, p);
        
        // r <- b - A*x
        blas::xcopy(dim, rhs, 1, r, 1);
        blas::xaxpy(dim, -1.0, p, 1, r, 1);
        
        // p <- r
        blas::xcopy(dim, r, 1, p, 1);
        
        // r_star <- r
        blas::xcopy(dim, r, 1, r0, 1);
        
        double r0_old = blas::dot(dim, r0, r);
        
        while ( !monitor.finished(dim, r) )
        {
            // Mp = M*p
            mat.precondition(p, Mp);
            
            // AMp = A*Mp
            mat.multiply(Mp, AMp);

            // alpha = (r_j, r_star) / (A*M*p, r_star)
            double alpha = r0_old / blas::dot(dim, r0, AMp);
            
            // s_j = r_j - alpha * AMp
            blas::xcopy(dim, r, 1, s, 1);
            blas::xaxpy(dim, -alpha, AMp, 1, s, 1);
            
            if (monitor.finished(dim, s))
            {
                // x += alpha*M*p_j
                blas::xaxpy(dim, alpha, Mp, 1, sol, 1);
                ++monitor;
                break;
            }
            
            // Ms = M*s_j
            mat.precondition(s, Ms);
            
            // AMs = A*Ms
            mat.multiply(Ms, AMs);
            monitor += 2;

            // omega = (AMs, s) / (AMs, AMs)
            double omega = blas::dot(dim, AMs, s) / blas::dot(dim, AMs, AMs);
            
            // x_{j+1} = x_j + alpha*M*p_j + omega*M*s_j
            blas::xaxpy(dim, alpha, Mp, 1, sol, 1);
            blas::xaxpy(dim, omega, Ms, 1, sol, 1);
            
            // r_{j+1} = s_j - omega*A*M*s
            blas::xcopy(dim, s, 1, r, 1);
            blas::xaxpy(dim, -omega, AMs, 1, r, 1);

            // beta_j = (r_{j+1}, r_star) / (r_j, r_star) * (alpha/omega)
            double r0_new = blas::dot(dim, r0, r);
            double beta = (r0_new / r0_old) * (alpha / omega);
            r0_old = r0_new;
            
            // p_{j+1} = r_{j+1} + beta*(p_j - omega*A*M*p)
            blas::xaxpy(dim, -omega, AMp, 1, p, 1);
            blas::xscal(dim, beta, p, 1);
            blas::xaxpy(dim, 1.0, r, 1, p, 1);
        }
        allocator.release();
    }

}

#endif

