// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef QMR_H
#define QMR_H

#include "real.h"
#include "cblas.h"
#include "allocator.h"
#include "monitor.h"

/// QMR method to solve a system of linear equations
/**
 F. Nedelec, 20.03.2018 ----  UNFINISHED
*/
namespace LinearSolvers
{
    
    /// Quasi-Minimal Residual Method
    /*
     This solves `mat * sol = rhs` with a tolerance specified in 'monitor'
     */
    template < typename LinearOperator, typename Monitor, typename Allocator >
    void QMR(const LinearOperator& mat, const real* rhs, real* sol, Monitor& monitor, Allocator& allocator)
    {
        const int dim = mat.dimension();
        real resid, eta, delta, ep, beta;
        real rho, rho_1, xi, gamma, gamma_1, theta, theta_1;
        
        //allocate workspace
        allocator.allocate(dim, 14);
        
        real * r     = allocator.bind(0);
        real * v_tld = allocator.bind(1);
        real * y     = allocator.bind(2);
        real * w_tld = allocator.bind(3);
        real * z     = allocator.bind(4);
        real * v     = allocator.bind(5);
        real * w     = allocator.bind(6);
        real * y_tld = allocator.bind(7);
        real * z_tld = allocator.bind(8);
        real * p     = allocator.bind(9);
        real * q     = allocator.bind(10);
        real * p_tld = allocator.bind(11);
        real * d     = allocator.bind(12);
        real * s     = allocator.bind(13);

        //std::clog << "norm_rhs = " << blas::nrm2(dim, rhs) << '\n';
        // compute initial residual:
        blas::xcopy(dim, rhs, 1, r, 1);
        mat.multiply(sol, s);
        ++monitor;
        
        blas::xaxpy(dim, -1.0, s, 1, r, 1); // r = b - A*x

        blas::xcopy(dim, r, 1, v_tld, 1);
        mat.precondition(v_tld, y);
        rho = blas::nrm2(dim, y);
        
        blas::xcopy(dim, r, 1, w_tld, 1);
        blas::xcopy(dim, w_tld, 1, z, 1);
        xi = blas::nrm2(dim, z);

        gamma = 1.0;
        eta = -1.0;
        theta = 0.0;
        
        int it = 1;
        do {
            
            if (rho == 0.0)
            {
                monitor.flag(2);
                return;
            }
            
            if (xi == 0.0)
            {
                monitor.flag(7);
                return;
            }
            
            // v = (1.0 / rho) * v_tld;
            blas::xcopy(dim, v_tld, 1, v, 1);
            blas::xscal(dim, 1.0/rho, v_tld, 1);
            //y = (1.0 / rho) * y;
            blas::xscal(dim, 1.0/rho, y, 1);
        
            //w = (1.0 / xi) * w_tld;
            blas::xcopy(dim, w_tld, 1, w, 1);
            blas::xscal(dim, 1.0/xi, w_tld, 1);
            //z = (1.0 / xi) * z;
            blas::xscal(dim, 1.0/xi, z, 1);
            
            delta = blas::dot(dim, z, y);
            if (delta == 0.0)
            {
                monitor.flag(5);
                return;
            }

            //y_tld = M2.solve(y);
            //z_tld = M1.trans_solve(z);
            
            mat.precondition(y, y_tld);
            blas::xcopy(dim, z, 1, z_tld, 1);

            if (it > 1)
            {
                //p = y_tld - (xi * delta / ep) * p;
                blas::xscal(dim, -xi * delta / ep, p, 1);
                blas::xaxpy(dim, 1.0, y_tld, 1, p, 1);

                //q = z_tld - (rho * delta / ep) * q;
                blas::xscal(dim, -rho * delta / ep, q, 1);
                blas::xaxpy(dim, 1.0, z_tld, 1, q, 1);
            }
            else
            {
                //p = y_tld;
                blas::xcopy(dim, y_tld, 1, p, 1);
                //q = z_tld;
                blas::xcopy(dim, z_tld, 1, q, 1);
            }
            
            //p_tld = A * p;
            mat.multiply(p, p_tld);
            ++monitor;
            
            ep = blas::dot(dim, q, p_tld);
            if (ep == 0.0)
            {
                monitor.flag(6);
                return;
            }
            
            beta = ep / delta;
            if (beta == 0.0)
            {
                monitor.flag(3);
                return;
            }
            
            //v_tld = p_tld - beta * v;
            blas::xcopy(dim, p_tld, 1, v_tld, 1);
            blas::xaxpy(dim, -beta, v, 1, v_tld, 1);

            //y = M1.solve(v_tld);
            mat.precondition(v_tld, y);
            
            rho_1 = rho;
            rho = blas::nrm2(dim, y);
            
            //w_tld = A.trans_mult(q) - beta * w;
            mat.trans_multiply(q, w_tld);
            ++monitor;
            
            blas::xaxpy(dim, -beta, w, 1, w_tld, 1);
            
            //z = M2.trans_solve(w_tld);
            blas::xcopy(dim, w_tld, 1, z, 1);

            xi = blas::nrm2(dim, z);
            
            gamma_1 = gamma;
            theta_1 = theta;
            
            theta = rho / ( gamma_1 * beta );
            gamma = 1.0 / sqrt( 1.0 + theta * theta );
            
            if (gamma == 0.0)
            {
                monitor.flag(4);
                return;
            }
           
            eta = -eta * rho_1 * gamma * gamma / (beta * gamma_1 * gamma_1);
            
            if (it > 1)
            {
                //d = eta * p     + (theta_1 * theta_1 * gamma * gamma) * d;
                blas::xscal(dim, theta_1 * theta_1 * gamma * gamma, d, 1);
                blas::xaxpy(dim, eta, p, 1, d, 1);

                //s = eta * p_tld + (theta_1 * theta_1 * gamma * gamma) * s;
                blas::xscal(dim, theta_1 * theta_1 * gamma * gamma, s, 1);
                blas::xaxpy(dim, eta, p_tld, 1, s, 1);
            }
            else
            {
                //d = eta * p;
                blas::xcopy(dim, p, 1, d, 1);
                blas::xscal(dim, eta, d, 1);

                //s = eta * p_tld;
                blas::xcopy(dim, p_tld, 1, s, 1);
                blas::xscal(dim, eta, s, 1);
            }
            
            //x += d;
            blas::xaxpy(dim, 1.0, d, 1, sol, 1);

            //r -= s;
            blas::xaxpy(dim, -1.0, s, 1, r, 1);
           
        } while ( ! monitor.finished(dim, r) );
        
#if ( 1 )
        // calculate true residual = rhs - A * x
        mat.multiply(sol, r);
        blas::xaxpy(dim, -1.0, rhs, 1, r, 1);
        resid = blas::nrm2(dim, r);
        fprintf(stderr, "QMR %4i iteration %4i residual %10.6f\n", dim, monitor.count(), resid);
#endif

        allocator.release();
    }
};

#endif

