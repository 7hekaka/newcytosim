// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University.


/**
 Set rigidity terms with modulus 'R1' in diagonal and lower parts of `mat`,
 for a filament with 'cnt' points.
 */
template<typename MATRIX>
void addRigidityMatrix0(MATRIX& mat, const size_t inx, const size_t cnt, const real R1)
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
 Set rigidity terms with modulus 'R1' in diagonal and lower parts of `mat`,
 for a filament with 'cnt' points.
 */
template<typename MATRIX>
void addRigidityMatrix(MATRIX& mat, const size_t inx, const size_t cnt, const real R)
{
    assert_true( cnt > 2 );

    const real R1 = -1 * R;
    const real R2 =  2 * R;
    const real R4 =  4 * R;
    const real R5 = -5 * R;
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


template<typename MATRIX>
void addRigidityBlockMatrix(MATRIX& mat, const size_t inx, const size_t cnt, const real R)
{
    assert_true( cnt > 2 );

    const real R1 = -1 * R;
    const real R2 =  2 * R;
    const real R4 =  4 * R;
    const real R5 = -5 * R;
    const real R6 = -6 * R;
    
    constexpr size_t U = DIM, D = DIM*2, T = DIM*3;
    const size_t s = DIM * inx;
    const size_t e = s + DIM * ( cnt - 2 );
    
    mat.block(s  , s  ).add_diag(R1);
    mat.block(s+U, s  ).add_diag(R2);
    mat.block(s+D, s  ).add_diag(R1);
    
    mat.block(e+U, e+U).add_diag(R1);
    mat.block(e+U, e  ).add_diag(R2);
    
    if ( 3 < cnt )
    {
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
 The array `mat` must be square of dimension `dim * this->nPoints`
 Only terms above the diagonal and corresponding to the first subspace are set
 */
void addRigidityUpper(real* mat, size_t ldd, size_t cnt, const real R)
{
    assert_true( cnt > 2 );

    const real R1 = -1 * R;
    const real R2 =  2 * R;
    const real R4 =  4 * R;
    const real R5 = -5 * R;
    const real R6 = -6 * R;

    constexpr size_t U = DIM, D = DIM*2, T = DIM*3;
    const size_t e = DIM * ( cnt - 2 );
    const size_t f = DIM * ( cnt - 1 );
    
    mat[0      ] += R1;
    mat[  ldd*U] += R2;
    mat[  ldd*D] += R1;
    
    mat[e+ldd*f] += R2;
    mat[f+ldd*f] += R1;
    
    if ( 3 < cnt )
    {
        mat[U+ldd*U] += R5;
        mat[U+ldd*D] += R4;
        mat[U+ldd*T] += R1;
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
    }
}


/**
 Set elements of matrix `mat` corresponding to the elastic terms of the Fiber.
 The array `mat` must be square of dimension `dim * this->nPoints`
 Only terms above the diagonal and corresponding to the first subspace are set
 */
void addRigidityLower(real* mat, size_t ldd, size_t cnt, const real R)
{
    assert_true( cnt > 2 );

    const real R1 = -1 * R;
    const real R2 =  2 * R;
    const real R4 =  4 * R;
    const real R5 = -5 * R;
    const real R6 = -6 * R;

    constexpr size_t U = DIM, D = DIM*2, T = DIM*3;
    const size_t e = DIM * ( cnt - 2 );
    const size_t f = DIM * ( cnt - 1 );
    
    mat[0] += R1;
    mat[U] += R2;
    mat[D] += R1;
    
    mat[f+ldd*e] += R2;
    mat[f+ldd*f] += R1;
    
    if ( 3 < cnt )
    {
        mat[U+ldd*U] += R5;
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
        mat[ n   +ldd*n] += R6;
        mat[(n+U)+ldd*n] += R4;
        mat[(n+D)+ldd*n] += R1;
    }
}


/**
 Set rigidity terms with modulus 'R1' in diagonal and lower parts of `mat`,
 for a filament with 'cnt' points.
 */
void addRigidityBanded(real* mat, size_t ldd, size_t cnt, const real R)
{
    assert_true( cnt > 2 );

    const real R1 = -1 * R;
    const real R2 =  2 * R;
    const real R4 =  4 * R;
    const real R5 = -5 * R;
    const real R6 = -6 * R;

    const size_t s = 0;
    const size_t e = cnt - 2;

    mat[  s*ldd] += R1;
    mat[1+s*ldd] += R2;
    mat[2+s*ldd] += R1;
    
    mat[  (e+1)*ldd] += R1;
    mat[1+ e   *ldd] += R2;

    if ( 3 < cnt )
    {
        mat[  (s+1)*ldd] += R5;
        mat[1+(s+1)*ldd] += R4;
        mat[2+(s+1)*ldd] += R1;
        mat[   e   *ldd] += R5;
    }
    else
    {
        mat[(s+1)*ldd] -= R4;
    }
    
    for ( size_t n = s+2; n < e; ++n )
    {
        mat[  n*ldd] += R6;
        mat[1+n*ldd] += R4;
        mat[2+n*ldd] += R1;
    }
}


void setRigidityBanded(real* mat, size_t ldd, size_t cnt, const real R)
{
    assert_true( cnt > 2 );

    const real R1 = -1 * R;
    const real R2 =  2 * R;
    const real R4 =  4 * R;
    const real R5 = -5 * R;
    const real R6 = -6 * R;

    const size_t e = cnt - 2;

    mat[0] = R1;
    mat[1] = R2;
    mat[2] = R1;
    
    mat[  (e+1)*ldd] = R1;
    mat[1+(e+1)*ldd] = 0;
    mat[2+(e+1)*ldd] = 0;
    mat[1+  e*ldd]   = R2;
    mat[2+  e*ldd]   = 0;

    if ( 3 < cnt )
    {
        mat[  ldd] = R5;
        mat[1+ldd] = R4;
        mat[2+ldd] = R1;
        mat[e*ldd] = R5;
    }
    else
    {
        mat[ldd] = -R4;
    }
    
    for ( size_t n = 2; n < e; ++n )
    {
        mat[  n*ldd] = R6;
        mat[1+n*ldd] = R4;
        mat[2+n*ldd] = R1;
    }
}

