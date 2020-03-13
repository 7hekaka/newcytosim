

typedef double real;
#define DIM 2

/*
void copy(real * restrict dst, real * restrict src, int n)
{
    for (size_t i = 0; i < n; ++i)
        dst[i] = src[i];
}
*/
/*
void projectForcesU(unsigned nbs, const real* dif, const real* X, real* mul)
{
    for ( size_t jj = 0; jj < nbs; ++jj )
    {
        mul[jj] = dif[DIM*jj  ] * ( X[DIM*jj+DIM  ] - X[DIM*jj  ] )
                + dif[DIM*jj+1] * ( X[DIM*jj+DIM+1] - X[DIM*jj+1] )
#if ( DIM > 2 )
                + dif[DIM*jj+2] * ( X[DIM*jj+DIM+2] - X[DIM*jj+2] )
#endif
        ;
    }
}
*/


//void projectForcesV(unsigned nbs, const real restrict* dif, const real restrict* vec, real restrict* mul)
void projectForcesV(unsigned nbs, const real* dif, const real* vec, real* mul)
{
#pragma vector unaligned
    for ( unsigned jj = 0; jj < nbs; ++jj )
    {
        const real * X = vec + DIM * jj;
        const real * d = dif + DIM * jj;
        mul[jj] = d[0] * ( X[DIM  ] - X[0] )
                + d[1] * ( X[DIM+1] - X[1] )
#if ( DIM > 2 )
                + d[2] * ( X[DIM+2] - X[2] )
#endif
        ;
    }
}


void projectForcesSOA(unsigned nbs, const real*restrict dif, const real*restrict vec, real* mul)
{
    unsigned S = nbs, T = nbs + nbs;
    const real * X = vec, * Y = vec + S, * Z = vec + T;
    const real *dX = dif, *dY = dif + S, *dZ = dif + T;
    
    for ( unsigned jj = 0; jj < nbs; ++jj )
    {
        mul[jj] = dX[jj] * ( X[jj+1] - X[jj] )
                + dY[jj] * ( Y[jj+1] - Y[jj] )
#if ( DIM > 2 )
                + dZ[jj] * ( Z[jj+1] - Z[jj] )
#endif
        ;
    }
}

