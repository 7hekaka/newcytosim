/**
 This implements our previous idea of replacing divisions by Fused multiply-add,
 in the tridiagonal symmetric matrix factorization codes,
 which is implemented in alsatian_xpttrs(), combining with code from the LINPACK
 project where reduction is bidirectional (routine dptsl).
 
 LINPACK documentation has this note:
    Because of the two-way nature of the algorithm, systems can be solved up to
    25% faster than conventional algorithms. Although the techniques used here
    were independently discovered, they have been known for some time.
    (private communications, J. H. Wilkinson).
 
 Francois J. Nedelec, La Foret Fouesnant, 2--10 August 2022.
 */

/**
     sptsl given a positive definite tridiagonal matrix and a right
     hand side will find the solution.

     on entry

        N        integer
                 is the order of the tridiagonal matrix.

        D        real(N)
                 is the diagonal of the tridiagonal matrix.
                 on output D is destroyed.

        E        real(N)
                 is the offdiagonal of the tridiagonal matrix.
                 E(1) through E(N-1) should contain the offdiagonal.

        X        real(N)
                 is the right hand side vector.

     on return

        B        contains the solution.

     linpack. this version dated 08/14/78.
     Jack Dongarra, Argonne National Laboratory.
 
 Translated to C & transformed: Francois J Nedelec, 30 July 2022.
**/

void linpack_xptsl(int N, real D[], const real E[], real X[])
{
    if ( N == 1 )
    {
        X[0] = X[0] / D[0];
        return;
    }
    int mid = ( N - 1 ) / 2;
    int inx = N - 2;
    // zero top half of subdiagonal and bottom half of superdiagonal
    for ( int k = 0; k < mid; ++k, --inx )
    {
        real t1 = E[k] / D[k];
        D[k+1] -= t1 * E[k];
        X[k+1] -= t1 * X[k];
        real t2 = E[inx] / D[inx+1];
        D[inx] -= t2 * E[inx];
        X[inx] -= t2 * X[inx+1];
    }
    inx = mid - 1;
    // special code for even N
    if ( ( N & 1 ) == 0 )
    {
        // clean up for possible 2 x 2 block at center
        real t1 = E[mid] / D[mid];
        D[mid+1] -= t1 * E[mid];
        X[mid+1] -= t1 * X[mid];
        // back solve starting at the center, going towards the top and bottom
        X[mid+1] /= D[mid+1];
        X[mid] = ( X[mid] - E[mid] * X[mid+1] ) / D[mid];
        ++mid;
    }
    else
    {
        X[mid] = X[mid] / D[mid];
    }
    for( int kp1 = mid+1; kp1 < N; ++kp1, --inx )
    {
        X[inx] = ( X[inx] - E[inx] * X[inx+1] ) / D[inx];
        X[kp1] = ( X[kp1] - E[kp1-1] * X[kp1-1] ) / D[kp1];
    }
}


/** 
 alsatian_xptsl() was derived from linpack_xptsl(), replacing divisions by multiplications,
 and reordering some of the operations to reduce the dependency chain.
 Francois J. Nedelec, La Foret Fouesnant, 2--10 August 2022. */
void alsatian_xptsl(int N, real D[], real E[], real X[])
{
    int mid = ( N - 1 ) / 2;
    int inx = N - 2;
    // zero top half of subdiagonal and bottom half of superdiagonal
    real xk = X[0];
    real xi = X[N-1];
    real dk = D[0];
    real di = D[N-1];
    for ( int k = 0; k < mid; ++k, --inx )
    {
        dk = inverse(dk);
        D[k] = dk;
        real t1 = dk * E[k];
        X[k] = dk * xk;
        dk = D[k+1] - t1 * E[k];
        xk = X[k+1] - t1 * xk;
        E[k] = t1;

        di = inverse(di);
        D[inx+1] = di;
        real t2 = di * E[inx];
        X[inx+1] = di * xi;
        di = D[inx] - t2 * E[inx];
        xi = X[inx] - t2 * xi;
        E[inx] = t2;
    }
    if ( mid == inx ) // N is even
    {
        // clean up for possible 2 x 2 block at center
        dk = inverse(dk);
        D[mid] = dk;
        xk = xk * dk;
        real t2 = E[mid] * dk;
        di = inverse(di - t2 * E[mid]);
        D[mid+1] = di;
        xi = ( xi - E[mid] * xk ) * di;
        xk = xk - t2 * xi;
        X[mid] = xk;
        ++mid;
        --inx;
    }
    else // N is odd
    {
        di = inverse(dk+di-D[mid]);
        D[mid] = di;
        xi = ( xk + xi - X[mid] ) * di;
        xk = xi;
    }
    X[mid] = xi;
    // back solve starting at the center, going towards the top and bottom
    for ( int kp1 = mid+1; kp1 < N; ++kp1, --inx )
    {
        //X[inx] -= E[inx] * X[inx+1];  // inx from mid-1 to 0
        xk = X[inx] - E[inx] * xk;
        X[inx] = xk;
        //X[kp1] -= E[kp1-1] * X[kp1-1]; // kp1 from mid+1 to N-1
        xi = X[kp1] - E[kp1-1] * xi;
        X[kp1] = xi;
    }
}


/** 
 alsadual_factor() will modify D[] and E[] according to alsatian_xptsl() above,
 skipping the operations associated with X[] that are left for solve()
 Francois J. Nedelec, La Foret Fouesnant, 2--10 August 2022.
 */
void alsadual_factor(int N, real D[], real E[])
{
    int mid = ( N - 1 ) / 2;
    int T = N - 2;
    // zero top half of subdiagonal and bottom half of superdiagonal
    real dk = D[0];
    real di = D[N-1];
    for ( int k = 0; k < mid; ++k, --T )
    {
        dk = inverse(dk);
        D[k] = dk;
        real t1 = dk * E[k];
        dk = D[k+1] - t1 * E[k];
        E[k] = t1;

        // Attention: for odd N, k+1==T at the last iteration:
        di = inverse(di);
        D[T+1] = di;
        real t2 = di * E[T];
        di = D[T] - t2 * E[T];
        E[T] = t2;
    }
    if ( mid == T ) // N is even
    {
        // clean up for possible 2 x 2 block at center
        dk = inverse(dk);
        D[mid] = dk;
        D[mid+1] = inverse(di - dk * E[mid] * E[mid]);
    }
    else // N is odd
    {
        D[mid] = inverse(dk+di-D[mid]);
    }
}


void alsadual_xpttrf(int size, real* D, real* E, int* INFO)
{
    *INFO = 0;
    alsadual_factor(size, D, E);
}

/** 
 alsadual_xptts2() will solve for the solution X[], given
 the factorization calculated by alsadual_solve() above
Francois J. Nedelec, La Foret Fouesnant, 2--10 August 2022.
 */
void alsadual_xptts2(int N, const real D[], const real E[], real X[])
{
    int mid = ( N - 1 ) / 2;
    int T = N - 2;
    // zero top half of subdiagonal and bottom half of superdiagonal
    real xk = X[0];
    real xi = X[N-1];
    int k = 0;
    while ( k < mid )
    {
        // k going up from 0 to mid-1
        X[k] = D[k] * xk;
        xk = X[k+1] - E[k] * xk;
        // Attention: for odd N, k+1==inx at the last iteration:
        // T going down from N-1 to mid
        X[T+1] = D[T+1] * xi;
        xi = X[T] - E[T] * xi;
        ++k;
        --T;
    }
    ++k; // assert_true( k == mid+1 );
    if ( mid == T ) // N is even
    {
        // clean up for possible 2 x 2 block at center
        xk = xk * D[mid];
        real t2 = E[mid] * D[mid];
        xi = ( xi - E[mid] * xk ) * D[mid+1];
        xk = xk - t2 * xi;
        X[mid] = xk;
        X[k] = xi;
        ++k;
        --T;
    }
    else // N is odd
    {
        xi = ( xk + xi - X[mid] ) * D[mid];
        xk = xi;
        X[mid] = xi;
    }
    // back solve starting at the center, going towards the top and bottom
    while ( T >= 0 )
    {
        //X[T] -= E[T] * X[T+1];  // T from mid-1 to 0
        xk = X[T] - E[T] * xk;
        X[T] = xk;
        //X[k] -= E[k-1] * X[k-1]; // k from mid+1 to N-1
        xi = X[k] - E[k-1] * xi;
        X[k] = xi;
        ++k;
        --T;
    }
    //assert_true( k == N );
}

