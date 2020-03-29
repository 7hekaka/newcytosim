// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University.


/**
 Fill-in matrix 'dst' as the duplicate of 'src', for each 'DIM' dimension.
 For 'DIM==1', this makes a simple copy
 For 'DIM==2', 'src' is copied twice, into odd indices, and into even indices.
 For 'DIM==3', three copies of 'src' are made into 'dst'.
 
 Both 'src' and 'dst' must be symmetrix square matrices.
 The size of 'dst' is DIM times the size of 'src'.
 Only the upper diagonal of 'src' is specified.
 Matrix Y is specified in full.
 */

void duplicate_matrix(size_t siz, real const* src, real * dst)
{
    size_t ddd = DIM * siz;
    
    zero_real(ddd*ddd, dst);
    
    for ( size_t ii = 0; ii < siz; ++ii )
    {
        real xx = src[ii + siz * ii];
        
        size_t kk = ( ddd+1 ) * DIM * ii;
        for ( size_t d = 0; d < DIM; ++d, kk += ddd+1 )
            dst[kk] = xx;
        
        for ( size_t jj = ii+1; jj < siz; ++jj )
        {
            xx = src[ii + siz * jj];
            kk = DIM * ( ii + ddd * jj );
            size_t ll = DIM * ( jj + ddd * ii );
            for ( size_t d = 0; d < DIM; ++d )
            {
                dst[kk] = xx;
                dst[ll] = xx;
                kk += ddd+1;
                ll += ddd+1;
            }
        }
    }
    
#if ( 0 )
    std::clog << "\nOriginal:\n";
    VecPrint::print(std::clog, siz, siz, src, siz);
    std::clog << "Duplicated:\n";
    VecPrint::print(std::clog, ddd, ddd, dst, ddd);
#endif
}


/**
 This will symmetrize matrix `mat`, by copying the upper triangle to the lower one
 It will also copy the terms that are within the first subspace `X` into the other
 dimensions.
 
 input: upper triangular matrix
 ouput: full symmetric matrix
 */
void expand_upper_matrix(size_t siz, real * mat)
{
#if ( 0 )
    std::clog << "\nOriginal:\n";
    VecPrint::print(std::clog, siz, siz, mat, siz);
#endif
    
    for ( size_t jj = 0; jj < siz; jj += DIM  )
    {
        for ( size_t ii = 0; ii < jj; ii += DIM  )
        {
            real val = mat[ii+siz*jj];
            // expand term in other dimensions:
            for ( size_t d = 1; d < DIM; ++d )
                mat[ii+d+siz*(jj+d)] = val;
            
            // symmetrize matrix:
            for ( size_t d = 0; d < DIM; ++d )
                mat[jj+d+siz*(ii+d)] = val;
        }
        // expand diagonal term in other dimensions:
        real val = mat[jj+siz*jj];
        for ( size_t d = 1; d < DIM; ++d )
            mat[jj+d+siz*(jj+d)] = val;
    }

#if ( 0 )
    std::clog << "Expanded:\n";
    VecPrint::print(std::clog, siz, siz, mat, siz);
#endif
}

/**
 This will symmetrize matrix `mat`, by copying the lower triangle to the upper one
 It will also copy the terms that are within the first subspace `X` into the other
 dimensions.
 
 input: lower triangular matrix
 ouput: full symmetric matrix
 */
void expand_lower_matrix(size_t siz, real * mat)
{
#if ( 0 )
    std::clog << "\nOriginal:\n";
    VecPrint::print(std::clog, siz, siz, mat, siz);
#endif
    
    for ( size_t jj = 0; jj < siz; jj += DIM  )
    {
        for ( size_t ii = jj; ii < siz; ii += DIM  )
        {
            real val = mat[ii+siz*jj];
            // expand term in other dimensions:
            for ( size_t d = 1; d < DIM; ++d )
                mat[ii+d+siz*(jj+d)] = val;
            
            // symmetrize matrix:
            for ( size_t d = 0; d < DIM; ++d )
                mat[jj+d+siz*(ii+d)] = val;
        }
        // expand diagonal term in other dimensions:
        real val = mat[jj+siz*jj];
        for ( size_t d = 1; d < DIM; ++d )
            mat[jj+d+siz*(jj+d)] = val;
    }

#if ( 0 )
    std::clog << "Expanded:\n";
    VecPrint::print(std::clog, siz, siz, mat, siz);
#endif
}


/**
 Set to zero all the terms that are not within 'diag' from the diagonal.
 With 'diag==0' the entire matrix is set to zero.
 With 'diag==1', only the diagonal is kept.
 With 'diag==2', the matrix is made tri-diagonal
 etc.
 */
void truncate_matrix(size_t siz, real* mat, size_t diag)
{
#if ( 0 )
    std::clog << "\nOriginal:\n";
    VecPrint::print(std::clog, siz, siz, mat, siz);
#endif
    
    for ( size_t ii = 0; ii < siz; ++ii )
    {
        real * col = mat + siz * ii;
        for ( size_t jj = 0; jj+diag < ii; ++jj )
            col[jj] = 0.0;
        
        for ( size_t kk = ii+diag+1; kk < siz; ++kk )
            col[kk] = 0.0;
    }
    
#if ( 0 )
    std::clog << "Truncated:\n";
    VecPrint::print(std::clog, siz, siz, mat, siz);
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

    return sqrt( ( all - dia ) / dia );
}


/// set all values between '-val' and 'val' to zero
void threshold_matrix(size_t siz, real * mat, real val)
{
    for ( size_t k = 0; k < siz*siz; ++k )
    {
        if ( fabs(mat[k]) < val )
            mat[k] = 0.0;
    }
}


/// set 'mat' of order `siz` to `diag * I`
void diagonal_matrix(size_t siz, real * mat, real val)
{
    for ( size_t k = 0; k < siz*siz; ++k )
        mat[k] = 0.0;
    for ( size_t k = 0; k < siz*siz; k+=siz+1 )
        mat[k] = val;
}


/// erase all off-diagonal terms in `mat` of order `siz`
void make_diagonal(size_t siz, real * mat)
{
    for ( size_t j = 0; j < siz; ++j )
    {
        real * col = mat + j * siz;
        size_t i;
        for ( i = 0; i < j; ++i )
            col[i] = 0.0;
        for ( i = j+1; i < siz; ++i )
            col[i] = 0.0;
    }
}


/// a test matrix with integer components
void test_matrix(int siz, real * mat)
{
    for ( int i = 0; i < siz; ++i )
    for ( int j = 0; j < siz; ++j )
        mat[i+siz*j] = j - i;
}


/**
 Convert a full matrix into a LAPACK banded matrix data suitable for LU factorization.
 `src` is a square matrix of side 'siz'
 `dst` is of size ldd * siz, with ldd > ku+2*kl+1
 
 with indices starting at 1 (LAPACK documentation):
 src(i,j) is stored in dst(ku+1+i-j, j) for max(1, j-ku) <= i <= min(siz, j+kl)

 considering that ku <- ku + kl:
 src(i,j) is stored in dst(ku+kl+1+i-j, j) for max(1, j-ku) <= i <= min(siz, j+kl)

 with indices starting at 0:
 src(i,j) is stored in dst(ku+kl+i-j, j) for max(0, j-ku) <= i <= min(siz-1, j+kl)

 */
void banded_matrix(int siz, real const* src, int kl, int ku, real * dst, int ldd)
{
    assert_true( ldd == kl+kl+ku+1 );
#if ( 0 )
    if ( siz < 64 )
    {
        std::clog << "\noriginal:\n";
        VecPrint::print(std::clog, siz, siz, src, siz);
    }
#endif
    for ( int jj = 0; jj < siz; ++jj )
    {
        // dst[ii+ldd-kl-1-jj+ldd*jj] = src[ii+siz*jj];
        int add = ldd - kl - 1 - jj + ldd * jj;
        int inf = add + std::max(0, jj-ku);
        int sup = add + std::min(siz-1, jj+kl);
        int off = siz * jj - add;

        //std::clog << " inf " << inf << " sup " << sup << " off " << off << " jj " << jj << "  shift " << shift << "\n";
        //std::clog << " dst[ " << inf << " " << sup << " ] <- src[ " << inf+off << " " << sup+off << " ] \n";
        
        if ( off > 0 )
        {
            for (int ii = inf; ii <= sup; ++ii )
                dst[ii] = src[ii+off];
        }
        else
        {
            for (int ii = sup; ii >= inf; --ii )
                dst[ii] = src[ii+off];
        }
    }
    // zero out values:
    for ( int jj = 0; jj < siz; ++jj )
    {
        int add = ldd - kl - 1 - jj;
        int inf = add + std::max(0, jj-ku);
        int sup = add + std::min(siz-1, jj+kl);

        real * col = dst + ldd * jj;
        for ( int ii = 0; ii < inf; ++ii )
            col[ii] = 0.0;
        for ( int ii = sup+1; ii < ldd; ++ii )
            col[ii] = 0.0;
    }
#if ( 0 )
    if ( siz < 64 )
    {
        std::clog << " banded_matrix (size " << siz << ") :\n";
        VecPrint::print(std::clog, ldd, siz, dst, ldd);
    }
#endif
}


/*
 uses power iterations to estimate the largest eigenvalue of `mat * tam + alpha * I`
 Vector `vec` is used to initialize the algorithm
 @returns an estimate of the largest eigenvalue
 The precision of the estimate is low: 10%
 */
real largest_eigenvalue(int siz, real const* blk, int const* piv, real const* mat, real alpha, real * vec, real * tmp)
{
    assert_true(siz > 0);
    const real TOLERANCE = 0.05;
    real oge, eig = blas::nrm2(siz, vec);
    //fprintf(stderr, "      power size %i eig %10.6f\n", siz, eig);

    int n, info = 0;
    for ( n = 0; n < siz; n += 2 )
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
        
        if ( fabs(oge-eig) < TOLERANCE * ( fabs(eig) + fabs(oge) ) )
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
real largest_eigenvalue(int siz, real const* mat, real const* tam, real alpha, real * vec, real * tmp)
{
    const real TOLERANCE = 0.05;
    real oge, eig = blas::nrm2(siz, vec);
    
    int n;
    for ( n = 0; n < siz; n += 2 )
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
        if ( fabs(oge-eig) < TOLERANCE * ( fabs(eig) + fabs(oge) ) )
            break;
    }
    //fprintf(stderr, "      power size %4i iter %3i: eigen %10.6f %10.6f\n", siz, n, eig, oge);
    
    return std::max(eig, oge);
}
