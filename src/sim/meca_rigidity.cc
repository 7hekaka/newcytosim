// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University.


/**
 Set bending elasticity terms with modulus 'R1' in diagonal and lower parts of `mat`,
 for a filament with 'cnt' points.
 */
template < typename MATRIX >
void addBendingRigidityMatrix0(MATRIX& mat, const size_t inx, const size_t cnt, const real R1)
{
    assert_true( cnt > 2 );
    const real R2 = R1 * 2;
    const real R4 = R1 * 4;

    const size_t end = inx + cnt - 1;
    
    for ( size_t i = inx+1; i < end; ++i )
    {
        mat(i-1, i-1) -= R1;
        mat(i  , i-1) += R2;
        mat(i+1, i-1) -= R1;
        mat(i  , i  ) -= R4;
        mat(i+1, i  ) += R2;
        mat(i+1, i+1) -= R1;
    }
}


/**
 Set bending elasticity terms on the diagonal and lower triangle of `mat`,
 for a filament with 'cnt' points and bending modulus 'R'.
 */
template < typename MATRIX >
void addBendingRigidityMatrix(MATRIX& mat, const size_t inx, const size_t cnt, const real R)
{
    assert_true( cnt > 2 );

    const real R1 = -1 * R;
    const real R2 = +2 * R;
    const real R4 = +4 * R;
    const real R6 = -6 * R;

    const size_t s = inx;
    const size_t e = s + ( cnt - 2 );

    mat(s  , s  ) += R1;
    mat(s+1, s  ) += R2;
    mat(s+2, s  ) += R1;
    
    mat(e+1, e+1) += R1;
    mat(e+1, e  ) += R2;

    if ( 3 < cnt )
    {
        const real R5 = -5 * R;
        mat(s+1, s+1) += R5;
        mat(s+2, s+1) += R4;
        mat(s+3, s+1) += R1;
        mat(e  , e  ) += R5;
    }
    else
    {
        mat(s+1, s+1) -= R4;
    }
    
    for ( size_t n = s+2; n < e ; n += 1 )
    {
        mat(n,   n) += R6;
        mat(n+1, n) += R4;
        mat(n+2, n) += R1;
    }
}


template < int ORD, typename MATRIX>
void addBendingRigidityBlockMatrix(MATRIX& mat, const size_t inx, const size_t cnt, const real R)
{
    assert_true( cnt > 2 );

    const real R1 = -1 * R;
    const real R2 = +2 * R;
    const real R4 = +4 * R;
    const real R6 = -6 * R;
    
    constexpr size_t U = ORD, D = ORD*2, T = ORD*3;
    const size_t s = ORD * inx;
    const size_t e = s + ORD * ( cnt - 2 );
    
    mat.block(s  , s  ).add_diag(R1);
    mat.block(s+U, s  ).add_diag(R2);
    mat.block(s+D, s  ).add_diag(R1);
    
    mat.block(e+U, e+U).add_diag(R1);
    mat.block(e+U, e  ).add_diag(R2);
    
    if ( 3 < cnt )
    {
        const real R5 = -5 * R;
        mat.block(s+U, s+U).add_diag(R5);
        mat.block(s+D, s+U).add_diag(R4);
        mat.block(s+T, s+U).add_diag(R1);
        mat.block(e  , e  ).add_diag(R5);
    }
    else
    {
        mat.block(s+U, s+U).add_diag(-R4);
    }
    
    for ( size_t n = s+U; n < e ; n += U )
    {
        mat.block(n,   n).add_diag(R6);
        mat.block(n+U, n).add_diag(R4);
        mat.block(n+D, n).add_diag(R1);
    }
}


/**
 Set elements of matrix `mat` corresponding to the elastic terms of the Fiber.
 The array `mat` must be square of dimension `ORD * cnt`
 Only terms corresponding to the first subspace are set
 */
template < size_t ORD >
void addBendingRigidity(real* mat, size_t ldd, size_t cnt, const real R)
{
    assert_true( cnt > 2 );

    const real R1 = -1 * R;
    const real R2 = +2 * R;
    const real R4 = +4 * R;
    const real R6 = -6 * R;

    constexpr size_t U = ORD, D = ORD*2, T = ORD*3;
    const size_t e = ORD * ( cnt - 2 );
    const size_t f = ORD * ( cnt - 1 );
    
    mat[0    ] += R1;
    mat[ldd*U] += R2;
    mat[ldd*D] += R1;
    mat[U    ] += R2;
    mat[D    ] += R1;

    mat[e+ldd*f] += R2;
    mat[f+ldd*f] += R1;
    mat[f+ldd*e] += R2;

    if ( 3 < cnt )
    {
        const real R5 = -5 * R;
        mat[U+ldd*U] += R5;
        mat[U+ldd*D] += R4;
        mat[U+ldd*T] += R1;
        mat[D+ldd*U] += R4;
        mat[T+ldd*U] += R1;
        mat[e+ldd*e] += R5;
    }
    else
    {
        mat[U+ldd*U] -= R4;
    }
    
    for ( size_t n = D; n < e; n += U )
    {
        mat[n+ldd* n   ] += R6;
        mat[n+ldd*(n+U)] += R4;
        mat[n+ldd*(n+D)] += R1;
        mat[n+U+ldd* n ] += R4;
        mat[n+D+ldd* n ] += R1;
    }
}


/**
 Set elements of matrix `mat` corresponding to the elastic terms of the Fiber.
 The array `mat` must be square of dimension `dim * this->nPoints`
 Only terms above the diagonal and corresponding to the first subspace are set
 */
template < size_t ORD >
void addBendingRigidityLower(real* mat, size_t ldd, size_t cnt, const real R)
{
    assert_true( cnt > 2 );

    const real R2 = +2 * R;
    const real R4 = +4 * R;
    const real R6 = -6 * R;

    constexpr size_t U = ORD, D = ORD*2, T = ORD*3;
    const size_t e = ORD * ( cnt - 2 );
    const size_t f = ORD * ( cnt - 1 );
    
    mat[0] -= R;
    mat[U] += R2;
    mat[D] -= R;
    
    mat[f+ldd*e] += R2;
    mat[f+ldd*f] -= R;
    
    if ( 3 < cnt )
    {
        const real R5 = -5 * R;
        mat[U+ldd*U] += R5;
        mat[D+ldd*U] += R4;
        mat[T+ldd*U] -= R;
        mat[e+ldd*e] += R5;
    }
    else
    {
        mat[U+ldd*U] -= R4;
    }
    
    for ( size_t n = D; n < e; n += U )
    {
        mat[  n+ldd*n] += R6;
        mat[U+n+ldd*n] += R4;
        mat[D+n+ldd*n] -= R;
    }
}


/*
 set lower triangle of `mat' corresponding to bending elasticity with parameter R
 */
template < size_t ORD >
void setBendingRigidity(real* mat, size_t ldd, size_t cnt, const real R)
{
    assert_true( cnt > 2 );
    constexpr size_t U = ORD, D = ORD*2, T = ORD*3;

    const real R1 = -1 * R;
    const real R2 = +2 * R;
    const real R4 = +4 * R;
    const real R6 = -6 * R;

    const size_t e = ORD * ( cnt - 2 );
    const size_t f = ORD * ( cnt - 1 );

    mat[0] = R1;
    mat[U] = R2;
    mat[D] = R1;
    
    mat[  f+f*ldd] = R1;
    mat[U+f+f*ldd] = 0;
    mat[D+f+f*ldd] = 0;
    mat[U+e+e*ldd] = R2;
    mat[D+e+e*ldd] = 0;

    if ( 3 < cnt )
    {
        const real R5 = -5 * R;
        mat[U+U*ldd] = R5;
        mat[D+U*ldd] = R4;
        mat[T+U*ldd] = R1;
        mat[e+e*ldd] = R5;
    }
    else
    {
        mat[U+ORD*ldd] = -R4;
    }
    
    for ( size_t n = D; n < e; n += ORD )
    {
        mat[  n+n*ldd] = R6;
        mat[U+n+n*ldd] = R4;
        mat[D+n+n*ldd] = R1;
    }
}

