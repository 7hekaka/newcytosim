// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University.


/**
 Fill-in matrix 'dst' as the duplicate of 'src', for each 'ORD' dimension.
 For 'ORD==1', this does nothing
 For 'ORD==2', the X-subspace of 'mat' is copied to the 'Y' subspace
 For 'ORD==3', the X-subspace of 'mat' is copied to the 'Y' and 'Z' subspaces
 
 The size of 'mat' is `ORD*lin x DOR*col`
 Only the X-suspace of 'src' is used.
 */

template < size_t ORD >
static void copy_lower_subspace(size_t siz, real* mat, size_t ldd, size_t rank)
{
    //VecPrint::full("copy_subspace", siz, siz, mat, ldd);
    
    for ( size_t j = 0; j < siz; j += ORD )
    for ( size_t i = j; i <= std::min(siz-1, j+rank); i += ORD )
    {
        real * dst = mat + i + ldd * j;
        real val = dst[0];
        for ( size_t d = 1; d < ORD; ++d )
            dst[d+ldd*d] = val;
    }
    
    //VecPrint::full("copied", siz, siz, mat, ldd);
}


/**
 This will copy the terms that are within the first subspace `X` into the other
 dimensions. if 'SYMMETRIZE == true', this will also copy the upper triangle to
 the lower one to make the matrix symmetric
 
 at entry: `mat` is an upper triangular matrix
 at exit: `mat` is a full symmetric matrix if `SYMMETRIZE==true`
 */
template < size_t ORD, bool SYMMETRIZE >
static void copy_upper_subspace(size_t siz, real* mat, size_t ldd)
{
    //VecPrint::full("copy_upper_subspace", siz, siz, mat, ldd);
    
    for ( size_t jj = 0; jj < siz; jj += ORD  )
    for ( size_t ii = 0; ii <= jj; ii += ORD  )
    {
        real * dst = mat + ii + ldd * jj;
        real val = dst[0];
        // expand term in other dimensions:
        for ( size_t d = 1; d < ORD; ++d )
            dst[d+ldd*d] = val;
        
        if ( SYMMETRIZE )
        {
            real * tsd = mat + jj + ldd * ii;
            for ( size_t d = 0; d < ORD; ++d )
                tsd[d+ldd*d] = val;
        }
    }

    //VecPrint::full("Expanded", siz, siz, mat, ldd);
}

/**
 This will copy the terms that are within the first subspace `X` into the other
 dimensions. if 'SYMMETRIZE == true', this will also copy the lower triangle to
 the upper one to make the matrix symmetric

 at entry, `mat` is a lower triangular matrix
 at exit, `mat` is a full symmetric matrix if `SYMMETRIZE==true`
 */
template < size_t ORD, bool SYMMETRIZE >
static void copy_lower_subspace(size_t siz, real* mat, size_t ldd)
{
    //VecPrint::full("\ncopy_lower_subspace", siz, siz, mat, ldd);
    
    for ( size_t jj =  0; jj < siz; jj += ORD )
    for ( size_t ii = jj; ii < siz; ii += ORD )
    {
        real * dst = mat + ii + ldd * jj;
        real val = dst[0];
        // expand term in other dimensions:
        for ( size_t d = 1; d < ORD; ++d )
            dst[d+ldd*d] = val;
        
        if ( SYMMETRIZE )
        {
            real * tsd = mat + jj + ldd * ii;
            for ( size_t d = 0; d < ORD; ++d )
                tsd[d+ldd*d] = val;
        }
    }
    
    //VecPrint::full("Expanded", siz, siz, mat, ldd);
}


/*
 This averages the terms in the different subspaces
 */
template < size_t ORD >
static void average_matrix(size_t siz, real* src, size_t ldd)
{
    //VecPrint::full("\naverage_matrix", siz, siz, src, ldd);

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
    //VecPrint::full("Averaged", siz, siz, src, ldd);
}


/*
 This averages the terms in the different subspaces
 */
template < size_t ORD >
static void project_matrix(size_t siz, real const* src, size_t lll, real* dst, size_t ldd)
{
    //VecPrint::full("\nproject_matrix", ORD*siz, ORD*siz, src, lll);

    for ( size_t jj = 0; jj < siz; ++jj )
    for ( size_t ii = 0; ii < siz; ++ii )
    {
        real const* ptr = src + ORD * ( jj * lll + ii );
        real val = ptr[0];
        for ( size_t d = 1; d < ORD; ++d )
            val += ptr[d*(lll+1)];
        dst[ii+ldd*jj] = val / (real)ORD;
    }
    //VecPrint::full("Projected", siz, siz, dst, ldd);
}


/*
 This averages the terms along X
 */
template < size_t ORD > [[maybe_unused]]
static void modify_matrix(size_t lin, size_t col, real* mat, size_t ldd)
{
    //VecPrint::full("\nproject_matrix", ORD*siz, ORD*siz, src, lll);

    for ( size_t jj = 0; jj < col; jj += ORD )
    for ( size_t ii = 0; ii < lin; ++ii )
    {
        real * ptr = mat + ( ii + jj * ldd );
        real val = ptr[0];
        for ( size_t d = 1; d < ORD; ++d )
            val += ptr[d*ldd];
        val /= (real)ORD;
        for ( size_t d = 0; d < ORD; ++d )
            ptr[d*ldd] = val;
    }
    //VecPrint::full("Modified", siz, siz, dst, ldd);
}


/**
 reset terms that are below diagonal `kl`, or above diagonal `ku`
 if 'ku==0' and 'kl==0', only the diagonal is kept.
 if 'ku==1' and 'kl==1', the matrix is made tri-diagonal.
 */
[[maybe_unused]]
static void truncate_matrix(size_t siz, real* mat, size_t ldd, size_t kl, size_t ku)
{
    //VecPrint::full("\ntruncate_matrix", siz, siz, mat, ldd);

    for ( size_t j = 0; j < siz; ++j )
    {
        real* col = mat + ldd * j;
        //zero out terms above the diagonal:
        for ( size_t i = 0; i+ku < j; ++i )
            col[i] = 0;
        
        //zero out terms below the diagonal:
        for ( size_t i = j+kl+1; i < siz; ++i )
            col[i] = 0;
    }
    
    //VecPrint::full("Truncated", siz, siz, mat, ldd);
}


/// sum(element^2) / sum(diagonal^2)
[[maybe_unused]]
static real off_diagonal_norm(size_t siz, real* mat)
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
[[maybe_unused]]
static void threshold_matrix(size_t siz, real* mat, real val)
{
    for ( size_t k = 0; k < siz*siz; ++k )
    {
        if ( abs_real(mat[k]) < val )
            mat[k] = 0.0;
    }
}


/// remove 2 bytes of fraction data, rounding up the value
inline static float truncate_float(const float arg)
{
    constexpr uint32_t MASK{0xffff0000};
    union { float f; uint32_t i; } udi { arg };
    udi.i &= MASK;
    return udi.f;
}

/// reduce precision of floating point by zeroing out part of the mantissa
[[maybe_unused]]
static void truncate_floats(size_t num, float* mat)
{
    for ( size_t i = 0; i < num; ++i )
    {
        float y = truncate_float(mat[i]);
        //printf(" %.12f --> %.12f\n", mat[i], y);
        mat[i] = y;
    }
}

/// keep only two bytes from floats
[[maybe_unused]]
static void tile_doubles(size_t num, double const* src, float * dst)
{
    static_assert(sizeof(double) == 2*sizeof(uint32_t), "wrong type in tile_doubles");
    uint32_t * fast = (uint32_t*)src;
    uint32_t * slow = (uint32_t*)dst;
    uint32_t * end = (uint32_t*)(src+num);
    //printf("\n%10.7f %10.7f -->", src[0], src[1]);
    while ( fast < end )
    {
#if __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
        slow[0] = fast[0];
#else
        slow[0] = fast[1];
#endif
        slow += 1;
        fast += 2;
    }
    //vec2 v = load1d(dst), vv = load2d(dst);
    //printf("  %10.7f %10.7f | %10.7f %10.7f\n", v[0], v[1], vv[0], vv[1]);
}


/// convert doubles to floats
static void convert_to_floats(size_t cnt, double const* src, float* dst)
{
    #pragma omp simd
    for ( size_t i = 0; i < cnt; ++i )
        dst[i] = (float)src[i];
}

#if 0
typedef uint16_t bfloat16;

static inline bfloat16 convert_to_bfloat16(float val)
{
#if __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
    return *reinterpret_cast<bfloat16*>(
        reinterpret_cast<uint16_t*>(&val));
#else
    return *reinterpret_cast<bfloat16*>(
        &(reinterpret_cast<uint16_t*>(&val)[1]));
#endif
}


/// convert doubles to floats
[[maybe_unused]]
static void convert_to_bfloat16s(size_t cnt, double const* src, bfloat16* dst)
{
    #pragma omp simd
    for ( size_t i = 0; i < cnt; ++i )
        dst[i] = convert_to_bfloat16(src[i]);
}
#endif

/// set 'mat' of order `siz` with `diag` on the diagonal and 'off' elsewhere
[[maybe_unused]]
static void init_matrix(size_t siz, real* mat, real dia, real off)
{
    for ( size_t k = 0; k < siz*siz; ++k )
        mat[k] = off;
    for ( size_t k = 0; k < siz*siz; k+=siz+1 )
        mat[k] = dia;
}


/// erase all off-diagonal terms in `mat` of order `siz`
[[maybe_unused]]
static void make_diagonal(size_t siz, real* mat, size_t ldd)
{
    for ( size_t j = 0; j < siz; ++j )
    {
        real* col = mat + j * ldd;
        for ( size_t i = 0; i < j; ++i )
            col[i] = 0.0;
        for ( size_t i = j+1; i < siz; ++i )
            col[i] = 0.0;
    }
}


/// a test matrix with integer components
[[maybe_unused]]
static void build_test_matrix(size_t siz, real* mat, size_t ldd)
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
template < int KD >
static void lower_band_storage(int N, real const* src, real* dst, int ldd)
{
    assert_true( ldd == KD+1 );
    assert_true( dst != src || ldd <= N );
    const int nkd = std::max(0, N-KD);
    for ( int j = 0; j < N; ++j )
    {
        int sup = std::min(N-1, j+KD);
        real const* S = src + N * j;
        real* D = dst + ldd * j;
        for ( int i = j; i <= sup; ++i )
            D[i-j] = S[i];
    }
    // zero-out unused values:
    for ( int j = nkd; j < N; ++j )
    {
        for ( int i = std::max(0, N-j); i < ldd; ++i )
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
[[maybe_unused]]
static void band_storage(size_t N, real const* src, size_t kl, size_t ku, real* dst, size_t ldd)
{
    assert_true( ldd == 2*kl+ku+1 );
    
    for ( size_t u = 0; u < N * ldd ; ++u )
        dst[u] = 0;
    
    for ( size_t j = 0; j < N; ++j )
    {
        size_t inf = j - std::min(j, ku);
        size_t sup = std::min(N-1, j+kl);
        real* D = dst + ldd * j + kl + ku - j;
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
[[maybe_unused]]
static real largest_eigenvalue(int siz, real const* blk, int const* piv, real const* mat, real alpha, real* vec, real* tmp)
{
    assert_true(siz > 0);
    const real TOLERANCE = 0.05;
    real oge, eig = blas::nrm2(siz, vec);
    //fprintf(stderr, "      power size %i eig %10.6f\n", siz, eig);

    int info = 0;
    for ( int n = 0; n < siz; n += 2 )
    {
        blas::xcopy(siz, vec, 1, tmp, 1);
        lapack::xgetrs('N', siz, 1, blk, siz, piv, tmp, siz, &info);
        assert_true(info==0);
        blas::xgemv('N', siz, siz, 1.0/eig, mat, siz, tmp, 1, alpha/eig, vec, 1);
        oge = blas::nrm2(siz, vec);
        //VecPrint::head(siz, vec, 3);
        
        blas::xcopy(siz, vec, 1, tmp, 1);
        lapack::xgetrs('N', siz, 1, blk, siz, piv, tmp, siz, &info);
        assert_true(info==0);
        blas::xgemv('N', siz, siz, 1.0/oge, mat, siz, tmp, 1, alpha/oge, vec, 1);
        eig = blas::nrm2(siz, vec);
        //VecPrint::head(siz, vec, 3);
        
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
[[maybe_unused]]
static real largest_eigenvalue(int siz, real const* mat, real const* tam, real alpha, real* vec, real* tmp)
{
    const real TOLERANCE = 0.05;
    real oge, eig = blas::nrm2(siz, vec);
    
    for ( int n = 0; n < siz; n += 2 )
    {
        blas::xgemv('N', siz, siz, 1.0/eig, mat, siz, vec, 1,       0.0, tmp, 1);
        blas::xgemv('N', siz, siz, 1.0,     tam, siz, tmp, 1, alpha/eig, vec, 1);
        oge = blas::nrm2(siz, vec);
        //VecPrint::head(siz, vec, 3);
        
        blas::xgemv('N', siz, siz, 1.0/oge, mat, siz, vec, 1,       0.0, tmp, 1);
        blas::xgemv('N', siz, siz, 1.0,     tam, siz, tmp, 1, alpha/oge, vec, 1);
        eig = blas::nrm2(siz, vec);
        //VecPrint::head(siz, vec, 3);
        //fprintf(stderr, "      power iter %3i: eigen %10.6f %10.6f\n", n, eig, oge);
        if ( abs_real(oge-eig) < TOLERANCE * ( abs_real(eig) + abs_real(oge) ) )
            break;
    }
    //fprintf(stderr, "      power size %4i iter %3i: eigen %10.6f %10.6f\n", siz, n, eig, oge);
    
    return std::max(eig, oge);
}
