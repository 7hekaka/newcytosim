// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University.


/**
 Fill-in matrix 'dst' as the duplicate of 'src', for each 'ORD' dimension.
 For 'ORD==1', this makes a simple copy
 For 'ORD==2', 'src' is copied twice, into odd indices, and into even indices.
 For 'ORD==3', three copies of 'src' are made into 'dst'.
 
 Both 'src' and 'dst' must be symmetrix square matrices.
 The size of 'dst' is ORD times the size of 'src'.
 Only the upper diagonal of 'src' is specified.
 Matrix Y is specified in full.
 */

template < size_t ORD >
void duplicate_matrix(size_t siz, real const* src, size_t ldd, real * dst)
{
    size_t ddd = ORD * siz;
    size_t lll = ORD * ldd;
    
    zero_real(ddd*ddd, dst);
    
    for ( size_t i = 0; i < siz; ++i )
    {
        real dia = src[i+ldd*i];
        size_t ii = ORD * i;
        
        for ( size_t d = 0; d < ORD; ++d )
            dst[(1+lll)*(ii+d)] = dia;
        
        for ( size_t j = i+1; j < siz; ++j )
        {
            real val = src[i+ldd*j];
            size_t jj = ORD * j;
            for ( size_t d = 0; d < ORD; ++d )
            {
                dst[ii+d+lll*(jj+d)] = val;
                dst[jj+d+lll*(ii+d)] = val;
            }
        }
    }
    
#if ( 0 )
    std::clog << "\nduplicate_matrix:\n";
    VecPrint::print(std::clog, siz, siz, src, ldd);
    std::clog << "Duplicated:\n";
    VecPrint::print(std::clog, ddd, ddd, dst, lll);
#endif
}


/**
 This will symmetrize matrix `mat`, by copying the upper triangle to the lower one
 It will also copy the terms that are within the first subspace `X` into the other
 dimensions.
 
 input: upper triangular matrix
 ouput: full symmetric matrix
 */
template < size_t ORD >
void expand_upper_matrix(size_t siz, real * mat, size_t ldd)
{
#if ( 0 )
    std::clog << "\nexpand_upper_matrix:\n";
    VecPrint::print(std::clog, siz, siz, mat, ldd);
#endif
    
    for ( size_t jj = 0; jj < siz; jj += ORD  )
    {
        for ( size_t ii = 0; ii < jj; ii += ORD  )
        {
            real val = mat[ii+ldd*jj];
            // expand term in other dimensions:
            for ( size_t d = 1; d < ORD; ++d )
                mat[ii+d+ldd*(jj+d)] = val;
            
            // symmetrize matrix:
            for ( size_t d = 0; d < ORD; ++d )
                mat[jj+d+ldd*(ii+d)] = val;
        }
        // expand diagonal term in other dimensions:
        real val = mat[jj+ldd*jj];
        for ( size_t d = 1; d < ORD; ++d )
            mat[jj+d+ldd*(jj+d)] = val;
    }

#if ( 0 )
    std::clog << "Expanded:\n";
    VecPrint::print(std::clog, siz, siz, mat, ldd);
#endif
}

/**
 This will symmetrize matrix `mat`, by copying the lower triangle to the upper one
 It will also copy the terms that are within the first subspace `X` into the other
 dimensions.
 
 input: lower triangular matrix
 ouput: full symmetric matrix
 */
template < size_t ORD >
void expand_lower_matrix(size_t siz, real * mat, size_t ldd)
{
#if ( 0 )
    size_t S = std::min(12UL, siz);
    std::clog << "\nexpand_lower_matrix:\n";
    VecPrint::print(std::clog, S, S, mat, ldd);
#endif
    
    for ( size_t jj = 0; jj < siz; jj += ORD  )
    {
        for ( size_t ii = jj; ii < siz; ii += ORD  )
        {
            real val = mat[ii+ldd*jj];
            // expand term in other dimensions:
            for ( size_t d = 1; d < ORD; ++d )
                mat[ii+d+ldd*(jj+d)] = val;
            
            // symmetrize matrix:
            for ( size_t d = 0; d < ORD; ++d )
                mat[jj+d+ldd*(ii+d)] = val;
        }
        // expand diagonal term in other dimensions:
        real val = mat[jj+ldd*jj];
        for ( size_t d = 1; d < ORD; ++d )
            mat[jj+d+ldd*(jj+d)] = val;
    }
    
#if ( 0 )
    std::clog << "Expanded:\n";
    VecPrint::print(std::clog, S, S, mat, ldd);
#endif
}


/*
 This averages the terms in the different subspaces
 */
template < size_t ORD >
void average_matrix(size_t siz, real* src, size_t ldd)
{
#if ( 0 )
    size_t S = std::min(12UL, siz);
    std::clog << "\naverage_matrix:\n";
    VecPrint::print(std::clog, S, S, src, ldd);
#endif
    for ( size_t jj = 0; jj < siz; jj += ORD  )
    for ( size_t ii = 0; ii < siz; ii += ORD  )
    {
        real* ptr = src + jj * ldd + ii;
        real val = ptr[0];
        for ( size_t d = 1; d < ORD; ++d )
            val += ptr[d*(ldd+1)];
        val /= (real)ORD;
        for ( size_t d = 0; d < ORD; ++d )
        for ( size_t u = 0; u < ORD; ++u )
            ptr[d*ldd+u] = 0;
        for ( size_t d = 0; d < ORD; ++d )
            ptr[d*(ldd+1)] = val;
    }
#if (0 )
    std::clog << "Averaged:\n";
    VecPrint::print(std::clog, S, S, src, ldd);
#endif
}


/*
 This averages the terms in the different subspaces
 */
template < size_t ORD >
void project_matrix(size_t siz, real const* src, size_t lll, real* dst, size_t ldd)
{
#if ( 0 )
    size_t S = std::min(12UL, siz);
    std::clog << "\nproject_matrix:\n";
    VecPrint::print(std::clog, ORD*S, ORD*S, src, lll);
#endif
    for ( size_t jj = 0; jj < siz; ++jj )
    for ( size_t ii = 0; ii < siz; ++ii )
    {
        real const* ptr = src + ORD * ( jj * lll + ii );
        real val = ptr[0];
        for ( size_t d = 1; d < ORD; ++d )
            val += ptr[d*(lll+1)];
        dst[ii+ldd*jj] = val / (real)ORD;
    }
#if ( 0 )
    std::clog << "Projected:\n";
    VecPrint::print(std::clog, S, S, dst, ldd);
#endif
}



/**
 reset terms that are below diagonal `kl`, or above diagonal `ku`
 if 'ku==0' and 'kl==0', only the diagonal is kept.
 if 'ku==1' and 'kl==1', the matrix is made tri-diagonal.
 */
void truncate_matrix(size_t siz, real* mat, size_t ldd, size_t kl, size_t ku)
{
#if ( 0 )
    std::clog << "\ntruncate_matrix:\n";
    VecPrint::print(std::clog, siz, siz, mat, ldd);
#endif

    for ( size_t j = 0; j < siz; ++j )
    {
        real * col = mat + ldd * j;
        //zero out terms above the diagonal:
        for ( size_t i = 0; i+ku < j; ++i )
            col[i] = 0;
        
        //zero out terms below the diagonal:
        for ( size_t i = j+kl+1; i < siz; ++i )
            col[i] = 0;
    }
    
#if ( 0 )
    std::clog << "Truncated:\n";
    VecPrint::print(std::clog, siz, siz, mat, ldd);
#endif
}


/// sum(element^2) / sum(diagonal^2)
real off_diagonal_norm(size_t siz, real * mat)
{
    real all = 0;
    for ( size_t k = 0; k < siz*siz; ++k )
        all += mat[k] * mat[k];

    real dia = 0;
    for ( size_t k = 0; k < siz*siz; k+=siz+1 )
        dia += mat[k] * mat[k];

    return std::sqrt( ( all - dia ) / dia );
}


/// set all values between '-val' and 'val' to zero
void threshold_matrix(size_t siz, real * mat, real val)
{
    for ( size_t k = 0; k < siz*siz; ++k )
    {
        if ( abs_real(mat[k]) < val )
            mat[k] = 0.0;
    }
}


/// set 'mat' of order `siz` with `diag` on the diagonal and 'off' elsewhere
void init_matrix(size_t siz, real * mat, real dia, real off)
{
    for ( size_t k = 0; k < siz*siz; ++k )
        mat[k] = off;
    for ( size_t k = 0; k < siz*siz; k+=siz+1 )
        mat[k] = dia;
}


/// erase all off-diagonal terms in `mat` of order `siz`
void make_diagonal(size_t siz, real * mat, size_t ldd)
{
    for ( size_t j = 0; j < siz; ++j )
    {
        real * col = mat + j * ldd;
        for ( size_t i = 0; i < j; ++i )
            col[i] = 0.0;
        for ( size_t i = j+1; i < siz; ++i )
            col[i] = 0.0;
    }
}


/// a test matrix with integer components
void test_matrix(size_t siz, real * mat, size_t ldd)
{
    for ( size_t i = 0; i < siz; ++i )
    for ( size_t j = 0; j < siz; ++j )
        mat[i+ldd*j] = j - i;
}

/**
Convert a full matrix into a LAPACK banded matrix data suitable for
Cholesky factorization by DPBTRF() or DPBTF2().

 `src` is a DOUBLE PRECISION square matrix of dimension `N * N`
 `dst` is a DOUBLE PRECISION array, dimension `ldd * N`
 
 The lower triangle of the symmetric band matrix `src`, is transferred into
 the first KD+1 rows of `dst`.
 
 The j-th column of `src` is stored in the j-th column of the array `dst`
 as follows:
 
      dst(i-j,j) = src(i,j)  for  j <= i <= min(N-1, j+KD).
 
 This should work even if 'src==dst' provided `ldd <= N`
 */
void lower_band_storage(size_t N, real const* src, size_t kd, real* dst, size_t ldd)
{
    assert_true( ldd == kd+1 );
    for ( size_t j = 0; j < N; ++j )
    {
        size_t sup = std::min(N-1, j+kd);
        real const* S = src + N * j;
        real * D = dst + ldd * j;
        for ( size_t i = j; i <= sup; ++i )
            D[i-j] = S[i];
    }
    // zero-out unused values:
    for ( size_t j = N-kd; j < N; ++j )
    {
        for ( size_t i = N-j; i < ldd; ++i )
            dst[i+ldd*j] = 0;
    }
}


/**
 Convert a full matrix into a LAPACK banded matrix data suitable for
 LU factorization by DGTRF() or DGTF2().
 
 `src` is a square matrix of dimension `N * N`
 `dst` is of dimension `ldd * N`, with `ldd > ku+2*kl+1`
 
 create matrix `dst` in band storage, in rows KL+1 to  2*KL+KU+1;
 rows 1 to KL of the array are not set.
 
 The j-th column of `src` is stored in the j-th column of `dst` as follows:
           dst(KL+KU+i-j, j) = src(i,j)
 for max(0,j-KU) <= i <= min(N-1, j+KL)
*/
void band_storage(size_t N, real const* src, size_t kl, size_t ku, real* dst, size_t ldd)
{
    assert_true( ldd == 2*kl+ku+1 );
    
    for ( size_t u = 0; u < N * ldd ; ++u )
        dst[u] = 0;
    
    for ( size_t j = 0; j < N; ++j )
    {
        size_t inf = j - std::min(j, ku);
        size_t sup = std::min(N-1, j+kl);
        real * D = dst + ldd * j + kl + ku - j;
        real const* S = src + N * j;
        for ( size_t i = inf; i <= sup; ++i )
            D[i] = S[i];
    }
}


/*
 uses power iterations to estimate the largest eigenvalue of `mat * tam + alpha * I`
 Vector `vec` is used to initialize the algorithm
 @returns an estimate of the largest eigenvalue
 The precision of the estimate is low: 10%
 */
real largest_eigenvalue(size_t siz, real const* blk, int const* piv, real const* mat, real alpha, real * vec, real * tmp)
{
    assert_true(siz > 0);
    const real TOLERANCE = 0.05;
    real oge, eig = blas::nrm2(siz, vec);
    //fprintf(stderr, "      power size %i eig %10.6f\n", siz, eig);

    int info = 0;
    for ( size_t n = 0; n < siz; n += 2 )
    {
        blas::xcopy(siz, vec, 1, tmp, 1);
        lapack::xgetrs('N', siz, 1, blk, siz, piv, tmp, siz, &info);
        assert_true(info==0);
        blas::xgemv('N', siz, siz, 1.0/eig, mat, siz, tmp, 1, alpha/eig, vec, 1);
        oge = blas::nrm2(siz, vec);
        //VecPrint::print(std::clog, std::min(16UL, siz), vec, 3);
        
        blas::xcopy(siz, vec, 1, tmp, 1);
        lapack::xgetrs('N', siz, 1, blk, siz, piv, tmp, siz, &info);
        assert_true(info==0);
        blas::xgemv('N', siz, siz, 1.0/oge, mat, siz, tmp, 1, alpha/oge, vec, 1);
        eig = blas::nrm2(siz, vec);
        //VecPrint::print(std::clog, std::min(16UL, siz), vec, 3);
        
        //fprintf(stderr, "      power iter %3i: eigen %10.6f %10.6f\n", n, eig, oge);
        
        if ( abs_real(oge-eig) < TOLERANCE * ( abs_real(eig) + abs_real(oge) ) )
            break;
    }
    //fprintf(stderr, "      power iter %3i: eigen %10.6f %10.6f\n", n, eig, oge);
    
    return std::max(eig, oge);
}


/*
 uses power iterations to estimate the largest eigenvalue of `mat * tam + alpha * I`
 Vector `vec` is used to initialize the algorithm
 @returns an estimate of the largest eigenvalue
 The precision of the estimate is low: 10%
 */
real largest_eigenvalue(size_t siz, real const* mat, real const* tam, real alpha, real * vec, real * tmp)
{
    const real TOLERANCE = 0.05;
    real oge, eig = blas::nrm2(siz, vec);
    
    for ( size_t n = 0; n < siz; n += 2 )
    {
        blas::xgemv('N', siz, siz, 1.0/eig, mat, siz, vec, 1,       0.0, tmp, 1);
        blas::xgemv('N', siz, siz, 1.0,     tam, siz, tmp, 1, alpha/eig, vec, 1);
        oge = blas::nrm2(siz, vec);
        //VecPrint::print(std::clog, std::min(16UL, siz), vec, 3);
        
        blas::xgemv('N', siz, siz, 1.0/oge, mat, siz, vec, 1,       0.0, tmp, 1);
        blas::xgemv('N', siz, siz, 1.0,     tam, siz, tmp, 1, alpha/oge, vec, 1);
        eig = blas::nrm2(siz, vec);
        //VecPrint::print(std::clog, std::min(16UL, siz), vec, 3);
        //fprintf(stderr, "      power iter %3i: eigen %10.6f %10.6f\n", n, eig, oge);
        if ( abs_real(oge-eig) < TOLERANCE * ( abs_real(eig) + abs_real(oge) ) )
            break;
    }
    //fprintf(stderr, "      power size %4i iter %3i: eigen %10.6f %10.6f\n", siz, n, eig, oge);
    
    return std::max(eig, oge);
}
