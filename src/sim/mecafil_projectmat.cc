// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
/**
 projections performed with explicit matrices
 This is a slow method, but it can be useful to compare with other methods
*/
#include "vecprint.h"


void Mecafil::buildProjection()
{
    //reset all variables for the projections:
    iProj       = nullptr;
    iDProj      = nullptr;
    iJJtiJ      = nullptr;
    iTMP        = nullptr;
    iJJtiJforce = nullptr;
}


void Mecafil::allocateProjection(const size_t ms)
{
    //std::clog << reference() << "allocateProjection(" << nbp << ")\n";
    free_real(iProj);
    free_real(iDProj);
    free_real(iJJtiJ);
    free_real(iJJtiJforce);
    free_real(iTMP);
    const size_t N = DIM * ms;
    iProj       = new_real(N*N);
    iDProj      = new_real(N*N);
    iJJtiJ      = new_real(N*ms);
    iJJtiJforce = new_real(ms);
    iTMP        = new_real(N);
}


void Mecafil::destroyProjection()
{
    //std::clog << reference() << "destroyProjection\n";
    free_real(iProj);
    free_real(iDProj);
    free_real(iJJtiJ);
    free_real(iJJtiJforce);
    free_real(iTMP);
    iProj       = nullptr;
    iDProj      = nullptr;
    iJJtiJ      = nullptr;
    iJJtiJforce = nullptr;
    iTMP        = nullptr;
}


#pragma mark -

/*
 Computes the projection matrix

     P = I - J' ( J J' )^-1 J

 that is associated with the length constraints:

     | point(p+1) - point(p) |^2 = lambda^2

 */
void Mecafil::makeProjection()
{
    const size_t nbc = nbSegments();             //number of constraints
    const size_t nbv = DIM * nbPoints();         //number of variables
    assert_true( nbc > 0 );
    
    //----- allocate needed temporaries:
    real* J    = new_real(nbv*nbc);
    real* JJt0 = new_real(nbc);
    real* JJt1 = new_real(nbc);
    
    //------------compute the projection matrix
    Vector v, w, dv, dw;
    zero_real(nbv*nbc, J);
    
    //set up the Jacobian matrix J and the diagonals of J * Jt
    w  = posP(0);
    dw.set(0,0,0);
    for ( size_t jj = 0; jj < nbc ; ++jj )
    {
        //set J:
        v = w;
        w = posP(jj+1);
        dv = -dw;
        dw = w - v;
        for ( size_t d = 0; d < DIM ; ++d )
        {
            J[jj+nbc*(DIM*jj+d)    ] = -dw[d];
            J[jj+nbc*(DIM*jj+DIM+d)] =  dw[d];
        }
        
        //set the diagonal and off-diagonal term of JJt:
        JJt0[jj] = 2 * dw.normSqr();  // diagonal
        JJt1[jj] = dot(dv, dw);       // off-diagonal (first term not used)
    }
    
    // JJtiJ <- J
    //blas::xcopy( nbc * nbv, J, 1, iJJtiJ, 1 );
    copy_real(nbc*nbv, J, iJJtiJ);
    
    // JJtiJ <- inv( JJt ) * J
    int info = 0;
    lapack::xptsv(nbc, nbv, JJt0, JJt1, iJJtiJ, nbc, &info);
    if ( info ) ABORT_NOW("lapack::ptsv() failed");
    
    // iProj <-  -Jt * JJtiJ
    blas::xgemm('T', 'N', nbv, nbv, nbc, -1.0, J, nbc, iJJtiJ, nbc, 0., iProj, nbv );
    
    // iProj <- iProj + I
    for ( size_t j = 0; j < nbv*nbv; j += nbv+1 )
        iProj[j] += 1.0;
    
    free_real(J);
    free_real(JJt0);
    free_real(JJt1);
}


/**
 Attention, the vector 'X' and 'Y' may point to the same address!
 */
void Mecafil::projectForces(const real* X, real* Y) const
{
    const size_t nbv = DIM * nbPoints();
    if ( X == Y )
    {
        // the BLAS::dsymv will fail if the input and output vectors are identical
        copy_real(nbv, X, iTMP);
        blas::xsymv('U', nbv, 1.0, iProj, nbv, iTMP, 1, 0.0, Y, 1);
        //blas::xgemv('N', nbv, nbv, 1.0, iProj, nbv, iTMP, 1, 0.0, Y, 1);
    }
    else
    {
        blas::xsymv('U', nbv, 1.0, iProj, nbv, X, 1, 0.0, Y, 1);
        //blas::xgemv('N', nbv, nbv, 1.0, iProj, nbv, X, 1, 0.0, Y, 1);
    }
}


void Mecafil::printProjection(std::ostream& os) const
{
    const size_t nbv = DIM * nbPoints();
    os << reference() << '\n';
    VecPrint::print(os, nbv, nbv, iProj, nbv);
}


void Mecafil::computeTensions(const real* force)
{
    const size_t nbs = nbSegments();
    const size_t nbv = DIM * nbPoints();
    
    // calculate the lagrangian multipliers:
    blas::xgemv('N', nbs, nbv, 1., iJJtiJ, nbs, force, 1, 0., iLag, 1);
}


void Mecafil::storeTensions(const real* force)
{
    computeTensions(force);
}

//------------------------------------------------------------------------------
#pragma mark -

void Mecafil::makeProjectionDiff(const real* force)
{
    const size_t nbs = nbSegments();             //number of constraints
    const size_t nbv = DIM * nbPoints();         //number of variables
    
    // calculate the lagrangian coefficients:
    blas::xgemv('N', nbs, nbv, 1., iJJtiJ, nbs, force, 1, 0., iLag, 1);
    
    //printf("Lagrange: "); VecPrint::print(std::clog, nbc, iLag);
    
    // select expensive forces ( lagrangian > 0 )
    useProjectionDiff = false;
    for ( size_t ii = 0; ii < nbs; ++ii )
    {
        if ( iLag[ii] > 0 )
        {
            iJJtiJforce[ii] = iLag[ii];
            useProjectionDiff = true;
        }
        else
            iJJtiJforce[ii] = 0.0;
    }
    
    //printf("diffP ");VecPrint::print(std::clog, nbs, iJJtiJforce);
    
    //set up the first term in the derivative of J with respect to variable x[ii]
    //set up term  P * (DJ)t (JJti) J force:
    for ( size_t jj = 0; jj < nbv; ++jj )
    {
        real* col = iDProj + nbv * jj;
        zero_real(nbv, col);
        unsigned lin = jj / DIM;
        if ( lin > 0 )
        {
            col[jj-DIM] = +iJJtiJforce[lin-1];
            col[jj    ] = -iJJtiJforce[lin-1];
        }
        if ( lin < nbs )
        {
            col[jj    ] += -iJJtiJforce[lin];
            col[jj+DIM]  = +iJJtiJforce[lin];
        }
    }

    /*
     The final matrix is symmetric, for any force,
     as can be seen from the above relations to set its columns
     */
    //printf("projectionDiff\n");
    //VecPrint::print(std::clog, nbv, nbv, iDProj, nbv);
}


void Mecafil::addProjectionDiff( const real* X, real* Y ) const
{
    assert_true(useProjectionDiff);
    size_t nbv = DIM * nbPoints();
    blas::xsymv('U', nbv, 1.0, iDProj, nbv, X, 1, 1.0, Y, 1);
}

