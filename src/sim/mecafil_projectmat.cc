// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
/**
 projections performed with explicit matrices
 This is a slow method, but it can be useful to compare with other methods
*/
//#include "vecprint.h"


void Mecafil::buildProjection()
{
    //reset all variables for the projections:
    mtJJAlloc    = 0;
    mtP          = 0;
    mtDiffP      = 0;
    mtJJtiJ      = 0;
    mtJJtiJforce = 0;
}


void Mecafil::allocateProjection(const size_t nbp)
{
    if ( mtJJAlloc < nbp )
    {
        //std::clog << reference() << "allocateProjection(" << nbp << ")\n";
        if ( mtP )     free_real(mtP);
        if ( mtDiffP ) free_real(mtDiffP);
        if ( mtJJtiJ ) free_real(mtJJtiJ);
        mtJJAlloc = nbp;
        mtP          = new_real(DIM*nbp*DIM*nbp);
        mtDiffP      = new_real(DIM*nbp*DIM*nbp);
        mtJJtiJforce = new_real(nbp);
        mtJJtiJ      = new_real(DIM*nbp*nbp);
    }
}


void Mecafil::destroyProjection()
{
    //std::clog << reference() << "destroyProjection\n";
    if ( mtP )          free_real(mtP);
    if ( mtDiffP )      free_real(mtDiffP);
    if ( mtJJtiJforce ) free_real(mtJJtiJforce);
    if ( mtJJtiJ )      free_real(mtJJtiJ);
    mtP          = 0;
    mtDiffP      = 0;
    mtJJtiJforce = 0;
    mtJJtiJ      = 0;
}


/*
 Computes the projection matrix

     P = I - J' ( J J' )^-1 J

 that is associated with the length constraints:

     | point(p+1) - point(p) |^2 = lambda^2

 */
void Mecafil::makeProjection()
{
    const unsigned int nbc = nbSegments();             //number of constraints
    const unsigned int nbv = DIM * nbPoints();         //number of variables
    assert_true( nbc > 0 );
    
    assert_true( mtJJAlloc >= nbPoints() );
    
    //----- we allocate the arrays needed:
    real* J   = new_real(nbv*nbc);
    real* JJt = new_real(2*nbc);
    
    //------------compute the projection matrix
    real x;
    Vector v, w, dv, dw;
    int ofs, jj, kk;
    int info=0;
    
    zero_real(nbv*nbc, J);
    
    //set up the Jacobian matrix J and the diagonals of J * Jt
    w  = posP(0);
    for ( jj = 0; jj < nbc ; ++jj )
    {
        //set J:
        for ( ofs = 0; ofs < DIM ; ++ofs )
        {
            kk = DIM * jj + ofs;
            x = pPos[kk+DIM] - pPos[kk];
            J[jj+nbc*kk      ] = -x;
            J[jj+nbc*(kk+DIM)] =  x;
        }
        
        //set the diagonal and off-diagonal term of JJt:
        v  = w;
        w  = posP(jj+1);
        dv = dw;
        dw = w - v;
        JJt[jj] = 2 * dw.normSqr();     //diagonal term
        JJt[jj+nbc] = - dv * dw;        //off-diagonal term
    }
    
    // JJtiJ <- J
    //blas::xcopy( nbc * nbv, J, 1, mtJJtiJ, 1 );
    copy_real(nbc*nbv, J, mtJJtiJ);
    
    // JJtiJ <- inv( JJt ) * J
    lapack::xptsv(nbc, nbv, JJt, JJt+nbc+1, mtJJtiJ, nbc, &info);
    if ( info ) ABORT_NOW("lapack::ptsv() failed");
    
    // mtP <-  - Jt * JJtiJ
    blas::xgemm('T', 'N', nbv, nbv, nbc, -1., J, nbc,
               mtJJtiJ, nbc, 0., mtP, nbv );
    
    // mtP <- mtP + I
    for ( jj = 0; jj < nbv*nbv; jj += nbv+1 )
        mtP[jj] += 1.0;
    
    //printf(" m%lu\n", name); VecPrint::print( nbv, nbv, mtP );
    free_real(J);
    free_real(JJt);
}


void Mecafil::setSpeedsFromForces(const real* X, real alpha, real* Y) const
{
    const unsigned nbv = DIM * nbPoints();
    blas::xsymv('U', nbv, alpha, mtP, nbv, X, 1, 0.0, Y, 1);
}


void Mecafil::computeTensions(const real* force)
{
    const unsigned nbs = nbSegments();
    const unsigned nbv = DIM * nbPoints();
    
    // calculate the lagrangian coefficients:
    blas::xgemv('N', nbs, nbv, 1., mtJJtiJ, nbs, force, 1, 0., rfLag, 1);
}


void Mecafil::storeTensions(const real* force)
{
    computeTensions(force);
}

//------------------------------------------------------------------------------


void Mecafil::makeProjectionDiff(const real* force)
{
    const unsigned nbs = nbSegments();             //number of constraints
    const unsigned nbv = DIM * nbPoints();         //number of variables
    assert_true( nbs > 0 );
    
    // calculate the lagrangian coefficients:
    blas::xgemv('N', nbs, nbv, 1., mtJJtiJ, nbs, force, 1, 0., rfLag, 1);
    
    //printf("Lagrange: "); VecPrint::print(std::clog, nbc, rfLag);
    
    // select expensive forces ( lagrangian > 0 )
    for ( int ii = 0; ii < nbs; ++ii )
    {
        if ( rfLag[ii] > 0 )
            mtJJtiJforce[ii] = rfLag[ii];
        else
            mtJJtiJforce[ii] = 0.0;
    }
    
    //printf("diffP ");VecPrint::print(std::clog, nbs, mtJJtiJforce);
    
    //set up the first term in the derivative of J with respect to variable x[ii]
    //set up term  P * (DJ)t (JJti) J force:
    for ( int jj = 0; jj < nbv; ++jj )
    {
        real* coljj = mtDiffP + nbv * jj;
        zero_real(nbv, coljj);
        int lin = jj / DIM;
        if ( lin > 0 ) {
            coljj[jj-DIM] = +mtJJtiJforce[lin-1];
            coljj[jj    ] = -mtJJtiJforce[lin-1];
        }
        if ( lin < nbs ) {
            coljj[jj    ] += -mtJJtiJforce[lin];
            coljj[jj+DIM]  = +mtJJtiJforce[lin];
        }
    }

    /*
     The final matrix is symmetric, for any force,
     as can be seen from the above relations to set its columns
     */
    //printf("projectionDiff\n");
    //VecPrint::print(std::clog, nbv, nbv, mtDiffP);
}


void Mecafil::addProjectionDiff( const real* X, real* Y ) const
{
    unsigned int nbv = DIM * nbPoints();
    blas::xsymv('U', nbv, 1.0, mtDiffP, nbv, X, 1, 1.0, Y, 1);
}

