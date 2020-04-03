// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
/*
 Conjugate gradient and related iterative methods
 to solve linear systems: http://www.netlib.org/templates
*/

#include "conjgradient.h"
#include "blas.h"

/*
 The linear system is defined by functions provided as arguments:
 int dim
 matVect(const real* X, real* Y)
 precond(const real* X, real* Y)
 */

//-------all purpose allocation function:

void SolverC::allocate(size_t dim,
                       real** vec1, real** vec2, real** vec3, real** vec4,
                       real** vec5, real** vec6, real** vec7, real** vec8)
{
    real ** vec[] = { vec1, vec2, vec3, vec4, vec5, vec6, vec7, vec8 };
    
    for ( size_t ii = 0; ii < 8; ++ii )
    {
        if ( vec[ii] )
        {
            free_real(*vec[ii]);
            *vec[ii] = new_real(dim);
            //zero_real(dim, *vec[ii]);
        }
    }
}

void SolverC::release(real** vec1, real** vec2, real** vec3, real** vec4,
                      real** vec5, real** vec6, real** vec7, real** vec8)
{
    real* * vec[] = { vec1, vec2, vec3, vec4, vec5, vec6, vec7, vec8 };
    
    for ( size_t ii = 0; ii < 8; ++ii )
    {
        if ( vec[ii] )
            free_real(*vec[ii]);
    }
}


//=============================================================================
//              Conjugate Gradient, no Preconditionning
//=============================================================================


void SolverC::CG(int dim, const real* rhs, real* x,
                 void (*matVect)( const real*, real* ),
                 int& max_itr, real& max_res)
{
    real* d = nullptr, *s=nullptr, *r=nullptr, *q=nullptr;
    
    allocate(dim, &d, &s, &r, &q );
    
    double alpha, beta, dold, dnew;
    
    blas::xcopy(dim, rhs, 1, r, 1 );
    matVect( x, s );
    blas::xaxpy(dim, -1.0, s, 1, r, 1);          //   r <- rhs - A * x
    blas::xcopy(dim, r, 1, d, 1 );               //   d <- r
    dnew = blas::dot(dim, r, r );
    
    real res = 0;
    int  itr = 0;
    
    for ( ; itr <= max_itr; ++itr )
    {
        matVect( d, q );                          //   q = A * d
        
        alpha = dnew / blas::dot(dim, d, q);
        blas::xaxpy(dim,  alpha, d, 1, x, 1 );   //   x += alpha * d
        blas::xaxpy(dim, -alpha, q, 1, r, 1 );   //   r -= alpha * q
        
        dold = dnew;
        dnew = blas::dot(dim, r, r);
        beta = dnew / dold;
        blas::xscal(dim, beta, d, 1 );
        blas::xaxpy(dim, 1.0, r, 1, d, 1 );      //   d = beta * d + r
        
        res = blas::nrm2(dim, r );
        if ( res < max_res )
            break;
    }
    
    max_res = res;
    max_itr = itr;
    release(&d, &s, &r, &q);
}


//=============================================================================
//              Conjugate Gradient, with Preconditioning
//=============================================================================


void SolverC::CGP(int dim, const real* rhs, real* x,
                  void (*matVect)( const real*, real* ),
                  void (*precond)( const real*, real* ),
                  int& max_itr, real& max_res)
{
    real* d = nullptr, *s = nullptr, *r = nullptr, *q = nullptr;
    allocate(dim, &d, &s, &r, &q );

    double alpha, beta, dold, dnew;
    
    blas::xcopy(dim, rhs, 1, r, 1 );
    matVect( x, s );
    blas::xaxpy(dim, -1.0, s, 1, r, 1);          //   r = rhs - M * x
    
    precond( r, d );                               //   d <- inv(M) * r
    
    dnew = blas::dot(dim, r, d);
    
    real res = 0;
    int  itr = 0;
    
    for ( ; itr <= max_itr; ++itr )
    {
        matVect( d, q );                           //   q = M * d
        
        alpha = dnew / blas::dot(dim, d, q);
        blas::xaxpy(dim,  alpha, d, 1, x, 1 );   //   x += alpha * d
        blas::xaxpy(dim, -alpha, q, 1, r, 1 );   //   r -= alpha * q
        
        precond( r, s );                           //   s = inv(M) * r;
        
        dold = dnew;
        dnew = blas::dot(dim, r, s);
        beta = dnew / dold;
        blas::xscal(dim, beta, d, 1 );
        blas::xaxpy(dim, 1.0, s, 1, d, 1 );      //   d = beta * d + s
        
        res = blas::nrm2(dim, r );
        if ( res < max_res )
            break;
    }

    max_res = res;
    max_itr = itr;
    release(&d, &s, &r, &q);
}


//=============================================================================
//                      Bi-Conjugate Gradient
//=============================================================================


void SolverC::BCG(int dim, const real* rhs, real* x,
                   void (*matVect)( const real*, real* ),
                   void (*matVectTrans)( const real*, real* ),
                   int& max_itr, real& max_res)
{
    real *r = nullptr, *rb = nullptr, *p = nullptr;
    real *pb = nullptr, *q = nullptr, *qb = nullptr;
    allocate(dim, &r, &rb, &p, &pb, &q, &qb );
    
    real alpha, beta, dold, dnew;
    
    blas::xcopy(dim, rhs, 1, r, 1 );
    matVect( x, rb );
    blas::xaxpy(dim, -1.0, rb, 1, r, 1);          //   r = rhs - A * x
    
    blas::xcopy(dim, r, 1, p, 1 );
    blas::xcopy(dim, r, 1, rb, 1 );
    blas::xcopy(dim, r, 1, pb, 1 );
    
    dnew = blas::dot(dim, rb, r);
    
    real res = 0;
    int  itr = 0;
    
    for ( ; itr <= max_itr; ++itr )
    {
        matVect( p, q );                           //   q = A * p
        matVectTrans( pb, qb );                    //   qb = A' * pb
        
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
        
        res = blas::nrm2(dim, r );
        if ( res < max_res )
            break;
    }

    max_res = res;
    max_itr = itr;
    release(&r, &rb, &p, &pb, &q, &qb );
}


//=============================================================================
//                 Bi-Conjugate Gradient Stabilized
//=============================================================================


int SolverC::BCGS(int dim, const real* rhs, real* x,
                  void (*matVect)( const real*, real* ),
                  int& max_itr, real& max_res)
{
    double rho_1, rho_2 = 1, alpha = 0, beta, omega = 1.0;
    real* r = nullptr, *rtilde = nullptr, *p = nullptr, *t = nullptr, *v = nullptr;
    
    allocate(dim, &r, &rtilde, &p, &t, &v );
    
    blas::xcopy(dim, rhs, 1, r, 1 );
    matVect( x, rtilde );
    blas::xaxpy(dim, -1.0, rtilde, 1, r, 1);       // r = rhs - A * x
    blas::xcopy(dim, r, 1, rtilde, 1 );
    
    real res = 0;
    int  ret = 1;
    int  itr = 0;
    
    for ( ; itr <= max_itr; ++itr )
    {
        rho_1 = blas::dot(dim, rtilde, r);
        
        if ( rho_1 == 0 )
        {
            res = blas::nrm2(dim, r);
            ret = 2;
            break;
        }
        
        if ( itr == 1 )
            blas::xcopy(dim, r, 1, p, 1 );        // p = r;
        else {
            beta = (rho_1/rho_2) * (alpha/omega);
            blas::xaxpy(dim, -omega, v, 1, p, 1);
            blas::xscal(dim, beta, p, 1);
            blas::xaxpy(dim, 1.0, r, 1, p, 1);    // p = r + beta*(p-omega*v);
        }
        
        matVect( p, v );                          // v = A * p;
        alpha = rho_1 / blas::dot(dim, rtilde, v);
        blas::xaxpy(dim, -alpha, v, 1, r, 1);     // r = r - alpha * v;
        blas::xaxpy(dim,  alpha, p, 1, x, 1);     // x = x + alpha * p;
        
        res = blas::nrm2(dim, r);
        if ( res < max_res )
        {
            ret = 0;
            break;
        }
        
        matVect( r, t );                          // t = A * s;
        
        omega = blas::dot(dim, t, r) / blas::dot(dim, t, t);
        blas::xaxpy(dim,  omega, r, 1, x, 1);     // x = x + omega * r;
        blas::xaxpy(dim, -omega, t, 1, r, 1);     // r = r - omega * t;
        
        res = blas::nrm2(dim, r);
        
        if ( res < max_res )
        {
            ret = 0;
            break;
        }
        
        if ( omega == 0 )
        {
            ret = 3;
            break;
        }
        
        rho_2 = rho_1;
    }

    max_res = res;
    max_itr = itr;
    release(&r, &rtilde, &p, &t, &v);
    return ret;
}


//=============================================================================
//        Bi-Conjugate Gradient Stabilized with Preconditionning
//=============================================================================


int SolverC::BCGSP(int dim, const real* rhs, real* x,
                   void (*matVect)( const real*, real* ),
                   void (*precond)( const real*, real* ),
                   int& max_itr, real& max_res)
{
    double rho_1, rho_2 = 1, alpha = 0, beta, omega = 1.0, delta;
    real *r = nullptr, *rtilde = nullptr, *p = nullptr, *phat = nullptr;
    real *shat = nullptr, *t = nullptr, *v = nullptr;
    
    allocate(dim, &r, &rtilde, &p, &t, &v, &phat, &shat );
    
    blas::xcopy(dim, rhs, 1, r, 1 );
    matVect( x, rtilde );
    blas::xaxpy(dim, -1.0, rtilde, 1, r, 1);         // r = rhs - A * x
    
    blas::xcopy(dim, r, 1, rtilde, 1 );              // r_tilde = r
     
    real res = 0;
    int  ret = 1;
    int  itr = 0;
    
    for ( ; itr <= max_itr; ++itr )
    {
        rho_1 = blas::dot(dim, rtilde, r);
        
        if ( rho_1 == 0 )
        {
            res = blas::nrm2(dim, r);
            ret = 2;
            break;
        }
        
        if ( itr == 1 )
            blas::xcopy(dim, r, 1, p, 1 );          // p = r;
        else {
            beta = (rho_1/rho_2) * (alpha/omega);
            //we should test here the value of beta, which is scalar
            blas::xaxpy(dim, -omega, v, 1, p, 1);
            blas::xscal(dim, beta, p, 1);
            blas::xaxpy(dim, 1.0, r, 1, p, 1);      // p = r + beta*(p-omega*v);
        }
        
        precond( p, phat );                         // phat = inv(M) * p;
        matVect( phat, v );                         // v = M * phat;
        
        //added test for failure detected by D. Foethke, Jan 2005
        delta = blas::dot(dim, rtilde, v);
        if ( delta == 0 )
        {
            res = blas::nrm2(dim, r);
            ret = 4;
            break;
        }
        alpha = rho_1 / delta;
        blas::xaxpy(dim, -alpha, v, 1, r, 1);       // r = r - alpha * v;
        blas::xaxpy(dim,  alpha, phat, 1, x, 1);    // x = x + alpha * phat;
        
        res = blas::nrm2(dim, r);
        if ( res < max_res )
        {
            ret = 0;
            break;
        }
        
        precond( r, shat );                         // shat = inv(M) * r
        matVect( shat, t );                         // t = M * shat
        
        omega = blas::dot(dim, t, r) / blas::dot(dim, t, t);
        blas::xaxpy(dim,  omega, shat, 1, x, 1);    // x = x + omega * shat
        blas::xaxpy(dim, -omega, t, 1, r, 1);       // r = r - omega * t
        
        res = blas::nrm2(dim, r);
        if ( res < max_res )
        {
            ret = 0;
            break;
        }
        if ( omega == 0 )
        {
            ret = 3;
            break;
        }
        rho_2 = rho_1;
    }
    
    max_res = res;
    max_itr = itr;
    release(&r, &rtilde, &p, &t, &v, &phat, &shat);
    return ret;
}

