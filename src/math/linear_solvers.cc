// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

/*
 Conjugate gradient and related iterative methods
 to solve linear systems: http://www.netlib.org/templates
*/
#include "linear_solvers.h"
#include "cblas.h"


#pragma mark - Iterative Methods

/**
 Conjugate Gradient, no Preconditionning
 */

void LinearSolvers::CG(const LinearOperator& mat, const real* rhs, real* x, Monitor & monitor, Allocator & allocator)
{
    const int dim = mat.dimension();
    allocator.allocate(dim, 4);
    real * d = allocator.bind(0);
    real * s = allocator.bind(1);
    real * r = allocator.bind(2);
    real * q = allocator.bind(3);
    
    double alpha, beta, dold, dnew;
    
    blas::xcopy(dim, rhs, 1, r, 1 );
    mat.multiply( x, s );
    blas::xaxpy(dim, -1.0, s, 1, r, 1);          //   r <- rhs - A * x
    blas::xcopy(dim, r, 1, d, 1 );               //   d <- r
    dnew = blas::dot(dim, r, r );
        
    while ( ! monitor.finished(dim, r) )
    {
        mat.multiply( d, q );                     //   q = A * d
        ++monitor;

        alpha = dnew / blas::dot(dim, d, q);
        blas::xaxpy(dim,  alpha, d, 1, x, 1 );   //   x += alpha * d
        blas::xaxpy(dim, -alpha, q, 1, r, 1 );   //   r -= alpha * q
        
        dold = dnew;
        dnew = blas::dot(dim, r, r);
        beta = dnew / dold;
        blas::xscal(dim, beta, d, 1 );
        blas::xaxpy(dim, 1.0, r, 1, d, 1 );      //   d = beta * d + r
    }
    
    allocator.release();
}


/**
Conjugate Gradient, with Preconditioning
*/

void LinearSolvers::CGP(const LinearOperator& mat, const real* rhs, real* x, Monitor & monitor, Allocator & allocator)
{    
    const int dim = mat.dimension();
    allocator.allocate(dim, 4);
    real * d  = allocator.bind(0);
    real * s  = allocator.bind(1);
    real * r  = allocator.bind(2);
    real * q  = allocator.bind(3);

    double alpha, beta, dold, dnew;
    
    blas::xcopy(dim, rhs, 1, r, 1 );
    mat.multiply( x, s );
    blas::xaxpy(dim, -1.0, s, 1, r, 1);          //   r = rhs - M * x
    
    mat.precondition( r, d );                      //   d <- inv(M) * r
    
    dnew = blas::dot(dim, r, d);
    
    while ( ! monitor.finished(dim, r) )
    {
        mat.multiply( d, q );                      //   q = M * d
        ++monitor;

        alpha = dnew / blas::dot(dim, d, q);
        blas::xaxpy(dim,  alpha, d, 1, x, 1 );    //   x += alpha * d
        blas::xaxpy(dim, -alpha, q, 1, r, 1 );    //   r -= alpha * q
        
        mat.precondition( r, s );                  //   s = inv(M) * r;
        
        dold = dnew;
        dnew = blas::dot(dim, r, s);
        beta = dnew / dold;
        blas::xscal(dim, beta, d, 1 );
        blas::xaxpy(dim, 1.0, s, 1, d, 1 );      //   d = beta * d + s
    }
    
    allocator.release();
}


/**
 Bi-Conjugate Gradient
*/

void LinearSolvers::BCG(const LinearOperator& mat, const real* rhs, real* x, Monitor & monitor, Allocator & allocator)
{
    const int dim = mat.dimension();
    allocator.allocate(dim, 6);
    real * r  = allocator.bind(0);
    real * rb = allocator.bind(1);
    real * p  = allocator.bind(2);
    real * pb = allocator.bind(3);
    real * q  = allocator.bind(4);
    real * qb = allocator.bind(5);
    
    double alpha, beta, dold, dnew;
    
    blas::xcopy(dim, rhs, 1, r, 1 );
    mat.multiply( x, rb );
    blas::xaxpy(dim, -1.0, rb, 1, r, 1);          //   r = rhs - A * x
    
    blas::xcopy(dim, r, 1, p, 1 );
    blas::xcopy(dim, r, 1, rb, 1 );
    blas::xcopy(dim, r, 1, pb, 1 );
    
    dnew = blas::dot(dim, rb, r);
    
    while ( ! monitor.finished(dim, r) )
    {
        mat.multiply( p, q );                       //   q = A * p
        mat.trans_multiply( pb, qb );               //   qb = A' * pb
        monitor += 2;
       
        alpha = dnew / blas::dot(dim, pb, q);
        blas::xaxpy(dim,  alpha, p, 1, x, 1 );    //   x  += alpha * p
        blas::xaxpy(dim, -alpha, q, 1, r, 1 );    //   r  -= alpha * q
        blas::xaxpy(dim, -alpha, qb, 1, rb, 1 );  //   rb -= alpha * qb
        
        dold = dnew;
        dnew = blas::dot(dim, r, rb);
        beta = dnew / dold;
        blas::xscal(dim, beta, p, 1 );
        blas::xaxpy(dim, 1.0, r, 1, p, 1 );       //   p  = beta * p  + r
        blas::xscal(dim, beta, pb, 1 );
        blas::xaxpy(dim, 1.0, rb, 1, pb, 1 );     //   pb = beta * pb + rb
    
    }
    
    allocator.release();
}


/**
 Bi-Conjugate Gradient Stabilized
*/

void LinearSolvers::BCGS(const LinearOperator& mat, const real* rhs, real* x, Monitor & monitor, Allocator & allocator)
{
    const int dim = mat.dimension();
    allocator.allocate(dim, 5);
    real * r      = allocator.bind(0);
    real * rtilde = allocator.bind(1);
    real * p      = allocator.bind(2);
    real * t      = allocator.bind(3);
    real * v      = allocator.bind(4);

    double rho_1 = 1, rho_2, alpha = 0, beta = 0, omega = 1;
    
    blas::xcopy(dim, rhs, 1, r, 1 );
    mat.multiply( x, rtilde );
    blas::xaxpy(dim, -1.0, rtilde, 1, r, 1);       // r = rhs - A * x
    blas::xcopy(dim, r, 1, rtilde, 1 );
    
    while ( !monitor.finished(dim, r) )
    {
        rho_2 = rho_1;
        rho_1 = blas::dot(dim, rtilde, r);
        
        if ( rho_1 == 0.0 )
        {
            monitor.finish(2, dim, r);
            break;
        }
        
        beta = ( rho_1 / rho_2 ) * ( alpha / omega );
        if ( beta == 0 )
        {
            // p = r;
            blas::xcopy(dim, r, 1, p, 1 );
        }
        else {
            // p = r + beta * ( p - omega * v )
            blas::xaxpy(dim, -omega, v, 1, p, 1);
            blas::xpay(dim, r, beta, p);
        }
        
        mat.multiply( p, v );                      // v = A * p;
        alpha = rho_1 / blas::dot(dim, rtilde, v);
        
        blas::xaxpy(dim, -alpha, v, 1, r, 1);     // r = r - alpha * v;
        blas::xaxpy(dim,  alpha, p, 1, x, 1);     // x = x + alpha * p;
        
        mat.multiply( r, t );                      // t = A * r;
        ++monitor;
        omega = blas::dot(dim, t, r) / blas::dot(dim, t, t);
        
        if ( omega == 0.0 )
        {
            monitor.finish(3, dim, r);
            break;
        }
        
        blas::xaxpy(dim,  omega, r, 1, x, 1);     // x = x + omega * r;
        blas::xaxpy(dim, -omega, t, 1, r, 1);     // r = r - omega * t;
    }

    allocator.release();
}


/**
 Bi-Conjugate Gradient Stabilized with Preconditionning
*/

void LinearSolvers::BCGSP(const LinearOperator& mat, const real* rhs, real* x, Monitor & monitor, Allocator & allocator)
{
    const int dim = mat.dimension();
    allocator.allocate(dim, 7);
    real * r      = allocator.bind(0);
    real * rtilde = allocator.bind(1);
    real * p      = allocator.bind(2);
    real * t      = allocator.bind(3);
    real * v      = allocator.bind(4);
    real * phat   = allocator.bind(5);
    real * shat   = allocator.bind(6);

    double rho_1 = 1, rho_2, alpha = 0, beta = 0, omega = 1.0, delta;
    
    blas::xcopy(dim, rhs, 1, r, 1 );
    mat.multiply( x, rtilde );
    blas::xaxpy(dim, -1.0, rtilde, 1, r, 1);         // r = rhs - A * x
    blas::xcopy(dim, r, 1, rtilde, 1 );              // r_tilde = r
    
    while ( ! monitor.finished(dim, r) )
    {
        rho_2 = rho_1;
        rho_1 = blas::dot(dim, rtilde, r);
        
        if ( rho_1 == 0.0 )
        {
            monitor.finish(2, dim, r);
            break;
        }
        
        beta = ( rho_1 / rho_2 ) * ( alpha / omega );
        if ( beta == 0.0 )
        {
            // p = r;
            blas::xcopy(dim, r, 1, p, 1 );
        }
        else {
            // p = r + beta * ( p - omega * v )
            blas::xaxpy(dim, -omega, v, 1, p, 1);
            blas::xpay(dim, r, beta, p);
        }
        
        mat.precondition( p, phat );                // phat = PC * p;
        mat.multiply( phat, v );                    // v = M * phat;
        ++monitor;

        delta = blas::dot(dim, rtilde, v);
        if ( delta == 0.0 )
        {
            monitor.finish(4, dim, r);
            break;
        }
        
        alpha = rho_1 / delta;
        blas::xaxpy(dim, -alpha, v, 1, r, 1);       // r = r - alpha * v;
        blas::xaxpy(dim,  alpha, phat, 1, x, 1);    // x = x + alpha * phat;
        
        if ( monitor.finished(dim, r) )
            break;

        mat.precondition( r, shat );                // shat = PC * r
        mat.multiply( shat, t );                    // t = M * shat
        ++monitor;
        
        omega = blas::dot(dim, t, r) / blas::dot(dim, t, t);
        
        if ( omega == 0.0 )
        {
            monitor.finish(3, dim, r);
            break;
        }
        
        blas::xaxpy(dim,  omega, shat, 1, x, 1);    // x = x + omega * shat
        blas::xaxpy(dim, -omega, t, 1, r, 1);       // r = r - omega * t
        
    }
    
    allocator.release();
}

