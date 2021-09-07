// Cytosim was created by Francois Nedelec. Copyright Cambridge University 2020
/**
 projections performed with explicit matrices
 This is a slow method, but it can be useful to compare with other methods
*/
#include "vecprint.h"


void Mecafil::initProjection()
{
    //reset all variables for the projections:
    iProj  = nullptr;
    iDProj = nullptr;
    iJJtiJ = nullptr;
    iTMP   = nullptr;
#if ADD_PROJECTION_DIFF
    iJJtJF = nullptr;
#endif
}


void Mecafil::allocateProjection(const size_t ms)
{
    //std::clog << reference() << "allocateProjection(" << nbp << ")\n";
    free_real(iProj);
    free_real(iDProj);
    free_real(iJJtiJ);
    free_real(iTMP);
    const size_t N = DIM * ms;
    iProj  = new_real(N*N);
    iDProj = new_real(N*N);
    iJJtiJ = new_real(N*ms);
    iTMP   = new_real(N);
#if ADD_PROJECTION_DIFF
    free_real(iJJtJF);
    iJJtJF = new_real(ms);
#endif
}


void Mecafil::destroyProjection()
{
    //std::clog << reference() << "destroyProjection\n";
    free_real(iProj);
    free_real(iDProj);
    free_real(iJJtiJ);
    free_real(iTMP);
    iProj  = nullptr;
    iDProj = nullptr;
    iJJtiJ = nullptr;
    iTMP   = nullptr;
#if ADD_PROJECTION_DIFF
    free_real(iJJtJF);
    iJJtJF = nullptr;
#endif
}


#pragma mark -

/**
 Compute the projection matrix

     P = I - J' ( J J' )^-1 J

 associated with the length constraints:

     | point(p+1) - point(p) |^2 = segmentation^2

 */
void Mecafil::makeProjection()
{
    const size_t nbc = nbSegments();      // number of constraints
    const size_t nbv = DIM * nbPoints();  // number of variables
    assert_true( nbc > 0 );
    
    //----- allocate temp space:
    real* J = new_real(nbc*(nbv+2));
    real* D = J + nbc*nbv;
    real* E = D + nbc;
    
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
        D[jj] = 2 * dw.normSqr();  // diagonal
        E[jj] = dot(dv, dw);       // off-diagonal (first term not used)
    }
    
    // JJtiJ <- J
    //blas::xcopy( nbc * nbv, J, 1, iJJtiJ, 1 );
    copy_real(nbc*nbv, J, iJJtiJ);
    
    // JJtiJ <- inv( JJt ) * J
    int info = 0;
    lapack::xptsv(nbc, nbv, D, E+1, iJJtiJ, nbc, &info);
    if ( info ) ABORT_NOW("lapack::ptsv() failed");
    
    // iProj <-  -Jt * JJtiJ
    blas::xgemm('T', 'N', nbv, nbv, nbc, -1.0, J, nbc, iJJtiJ, nbc, 0., iProj, nbv );
    
    // iProj <- iProj + I
    for ( size_t j = 0; j < nbv*nbv; j += nbv+1 )
        iProj[j] += 1.0;
    
    free_real(J);
}


/**
 Apply the projection associated with the length constraints
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
    VecPrint::full(os, nbv, nbv, iProj, nbv);
}


void Mecafil::computeTensions(const real* force)
{
    const size_t nbs = nbSegments();
    const size_t nbv = DIM * nbPoints();
    
    // calculate the lagrangian multipliers:
    blas::xgemv('N', nbs, nbv, 1., iJJtiJ, nbs, force, 1, 0., iLag, 1);
}

//------------------------------------------------------------------------------
#pragma mark - Correction terms to the Projection


#if ADD_PROJECTION_DIFF

void Mecafil::makeProjectionDiff(const real* force)
{
    const size_t nbs = nbSegments();             //number of constraints
    const size_t nbv = DIM * nbPoints();         //number of variables
    
    // calculate the lagrangian coefficients:
    blas::xgemv('N', nbs, nbv, 1., iJJtiJ, nbs, force, 1, 0., iLag, 1);
    
    //VecPrint::print("Lagrange: ", nbc, iLag);
    
    // select expensive forces ( lagrangian > 0 )
    useProjectionDiff = false;
    for ( size_t ii = 0; ii < nbs; ++ii )
    {
        if ( iLag[ii] > 0 )
        {
            iJJtJF[ii] = iLag[ii];
            useProjectionDiff = true;
        }
        else
            iJJtJF[ii] = 0.0;
    }
    
    //VecPrint::print("diffP ", nbs, iJJtJF);
    
    //set up the first term in the derivative of J with respect to variable x[ii]
    //set up term  P * (dJ)t * inverse(J*Jt) * J * force:
    for ( size_t jj = 0; jj < nbv; ++jj )
    {
        real* col = iDProj + nbv * jj;
        zero_real(nbv, col);
        unsigned lin = jj / DIM;
        if ( lin > 0 )
        {
            col[jj-DIM] = +iJJtJF[lin-1];
            col[jj    ] = -iJJtJF[lin-1];
        }
        if ( lin < nbs )
        {
            col[jj    ] += -iJJtJF[lin];
            col[jj+DIM]  = +iJJtJF[lin];
        }
    }

    /*
     The final matrix is symmetric, for any force,
     as can be seen from the above relations to set its columns
     */
    //printf("projectionDiff\n");
    //VecPrint::full(nbv, nbv, iDProj, nbv);
}


void Mecafil::addProjectionDiff( const real* X, real* Y ) const
{
    size_t nbv = DIM * nbPoints();
    blas::xsymv('U', nbv, 1.0, iDProj, nbv, X, 1, 1.0, Y, 1);
}

#endif
