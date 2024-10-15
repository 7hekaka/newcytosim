// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

#include "../base/bitmap.cc"


//------------------------------------------------------------------------------
#pragma mark - Connectivity Analysis


/// equalize flags for any existing matrix element between Mecables
template < typename MATRIX >
static void flagConnectedMecables(Array<Mecable*> const mecables, const index_t sup,
                                  Mecable** table, MATRIX const& MAT)
{
    // process all matrix columns:
    for ( index_t j = 0; j < sup; ++j )
    {
        Mecable const* A = table[j];
        for ( index_t n = 0; n < MAT.column_size(j); ++n )
        {
            /* we do not check the value of the elements here, but just
             the fact of having a block at these indices */
            size_t i = MAT.column_index(j, n);
            assert_true( i < sup );
            Mecable const* B = table[i];
            if ( A->flag() != B->flag() )
            {
                ObjectFlag f = std::min(A->flag(), B->flag());
                ObjectFlag g = std::max(A->flag(), B->flag());
                // replace g -> f everywhere:
                for ( Mecable * mec : mecables )
                {
                    if ( mec->flag() == g )
                        mec->flag(f);
                }
            }
        }
    }
}


/** Assuming that Mecable::flag() have been set already */
void Meca::flagClusters() const
{
    const index_t MAX = nbVertices();
    Mecable ** table = new Mecable*[MAX]{nullptr};
    
    for ( Mecable * mec : mecables )
    {
        const index_t inx = mec->matIndex();
        const index_t end = mec->nbPoints() + inx;
        assert_true( end <= MAX );
        for ( size_t i = inx; i < end; ++i )
            table[i] = mec;
    }
    
    flagConnectedMecables(mecables, MAX, table, mFUL);
#if USE_ISO_MATRIX
    flagConnectedMecables(mecables, MAX, table, mISO);
#endif
    delete[] table;
}


//------------------------------------------------------------------------------
#pragma mark - Matrix Extraction

/**
 Count number of non-zero entries in the full system matrix
 */
size_t Meca::countTerms(const real threshold) const
{
    const size_t dim = dimension();
    real * src = new_real(dim);
    real * dst = new_real(dim);
    zero_real(dim, src);
    
    size_t cnt = 0;
    for ( size_t j = 0; j < dim; ++j )
    {
        src[j] = 1;
        multiply(src, dst);
        for ( size_t i = 0; i < dim; ++i )
            cnt += ( abs_real(dst[i]) >= threshold );
        src[j] = 0;
    }
    
    free_real(dst);
    free_real(src);
    return cnt;
}


/** This is the same as Meca::multiply() without the projection */
void Meca::multiplyElasticity(const real* X, real* Y) const
{
#if USE_ISO_MATRIX
    // Y <- mFUL * X
    if ( useFullMatrix )
        mFUL.vecMul(X, Y);
    else
        zero_real(dimension(), Y);
    // Y <- Y + mISO * X
    mISO.VECMULADDISO(X, Y);
#else
    mFUL.vecMul(X, Y);
#endif
    
    for ( Mecable * mec : mecables )
    {
        const size_t inx = DIM * mec->matIndex();
#if SEPARATE_RIGIDITY_TERMS
        mec->addRigidity(X+inx, Y+inx);
#endif
#if ADD_PROJECTION_DIFF
        if ( mec->hasProjectionDiff() )
            mec->addProjectionDiff(X+inx, Y+inx);
#endif
        // Y <- X + beta * Y
        const real beta = -tau_ * mec->leftoverMobility();
        blas::xpay(DIM*mec->nbPoints(), X+inx, beta, Y+inx);
    }
}


/**
 Extract the matrix associated defined by MULTIPLY().
 The array `mat[]` should be preallocated to hold `dim*lda` real scalars,
 with `dim >= Meca::dimension()`, and `lda >= dim` should be the leading
 dimension of the array.
 */
template < Meca::MultiplyFuncPtr MULTIPLY >
void Meca::getMatrix(real * mat, index_t lda) const
{
    size_t dim = dimension();
    if ( lda < dim )
        throw InvalidIO("insufficient matrix dimension");
    real * src = new_real(dim);
    zero_real(dim, src);
    
    for ( size_t j = 0; j < dim; ++j )
    {
        src[j] = 1;
        (this->*MULTIPLY)(src, mat+j*lda);
        src[j] = 0;
    }
    
    free_real(src);
}

//------------------------------------------------------------------------------
#pragma mark - Text Export

static void saveVector(FILE * fp, size_t dim, real const* VEC)
{
    fprintf(fp, "%% This is a vector produced by Cytosim\n");
    fprintf(fp, "%% author: Francois J. Nedelec\n");
    fprintf(fp, "%% kind: biological cell simulation (cytoskeleton)\n");
    
    fprintf(fp, "%lu\n", dim);
    for ( size_t i = 0; i < dim; ++i )
        fprintf(fp, "%f\n", VEC[i]);
}


void Meca::saveObjectID(FILE * fp) const
{
    int i = 1;
    for ( Mecable const* mec : mecables )
    {
        const index_t nbp = DIM * mec->nbPoints();
        for ( index_t p = 0; p < nbp; ++p )
            fprintf(fp, "%if\n", i);
        ++i;
    }
}

void Meca::saveMobility(FILE * fp) const
{
    for ( Mecable const* mec : mecables )
    {
        const index_t nbp = mec->nbPoints();
        const real val = mec->pointMobility();
        for ( index_t p = 0; p < DIM * nbp; ++p )
            fprintf(fp, "%f\n", val);
    }
}

/**
 Save a sparse matrix in Matrix Market format
 https://math.nist.gov/MatrixMarket/formats.html
 This is a Sparse text format
 */
template < Meca::MultiplyFuncPtr MULTIPLY >
void Meca::saveMatrix(FILE * fp, const size_t dim, real threshold) const
{
    fprintf(fp, "%%%%MatrixMarket matrix coordinate real general\n");
    fprintf(fp, "%% This is a matrix produced by Cytosim\n");
    fprintf(fp, "%% author: Francois J. Nedelec\n");
    fprintf(fp, "%% kind: biological cell simulation (cytoskeleton)\n");

    real * src = new_real(dim);
    real * dst = new_real(dim);
    zero_real(dim, src);
    
    fpos_t pos;
    size_t cnt = 0, top = 0;
    fprintf(fp, "%lu %lu ", dim, dim);
    fgetpos(fp, &pos);
    fprintf(fp, "%12lu\n", cnt);

    for ( size_t j = 0; j < dim; ++j )
    {
        src[j] = 1;
        (this->*MULTIPLY)(src, dst);
        for ( size_t i = 0; i < dim; ++i )
            if ( abs_real(dst[i]) > threshold )
            {
                fprintf(fp, "%3lu %3lu %f\n", i, j, dst[i]);
                ++cnt;
                top += ( i > j );
            }
        src[j] = 0;
    }
    
    fsetpos(fp, &pos);
    fprintf(fp, "%12lu\n", cnt);
    printf("saveMatrix size %lu (%lu, %lu)\n", dim, cnt-dim-top, top);

    free_real(dst);
    free_real(src);
}


/**
 Save full matrix, elasticity matrix and right-hand-side vector
 */
void Meca::saveSystem() const
{
    size_t dim = dimension();
    
    FILE * f = FilePath::open_file("matrix.mtx", "w");
    if ( f )
    {
        saveMatrix<&Meca::multiply>(f, dim, 0);
        fclose(f);
    }
    
    f = FilePath::open_file("elasticity.mtx", "w");
    if ( f )
    {
        saveMatrix<&Meca::multiplyElasticity>(f, dim, 0);
        fclose(f);
    }

    f = FilePath::open_file("vector.mtx", "w");
    if ( f )
    {
        saveVector(f, dim, vRHS);
        fclose(f);
    }
}


/**
 save vectors and matrices in a text-based sparse formats
 */
void Meca::exportSystem() const
{
#if SEPARATE_RIGIDITY_TERMS
    std::clog << "incorrect dump since SEPARATE_RIGIDITY_TERMS is defined\n";
#endif
    FILE * f = FilePath::open_file("ord.txt", "w");
    fprintf(f, "%u %u %lu\n", dimension(), DIM, sizeof(real));
    fclose(f);
    
    f = FilePath::open_file("stp.txt", "w");
    fprintf(f, "%f %f\n", tau_, tolerance_);
    fclose(f);
    
    f = FilePath::open_file("mob.txt", "w");
    saveMobility(f);
    fclose(f);
    
    f = FilePath::open_file("obj.txt", "w");
    saveObjectID(f);
    fclose(f);
    
    f = FilePath::open_file("sol.txt", "w");
    VecPrint::dump(f, dimension(), vPTS);
    fclose(f);

    f = FilePath::open_file("rhs.txt", "w");
    VecPrint::dump(f, dimension(), vRHS);
    fclose(f);
    
    std::ofstream os("full.txt");
    mFUL.printSparse(os, 0);
    os.close();

#if USE_ISO_MATRIX
    os.open("iso.txt");
    mISO.printSparse(os, 0);
    os.close();
#endif
        
    index_t alc = 0;
    for ( Mecable const* mec : mecables )
        alc = std::max(alc, mec->nbPoints());

    real * tmp1 = new_real(DIM*alc);
    real * tmp2 = new_real(DIM*DIM*alc*alc);
    
    f = FilePath::open_file("diag.txt", "w");
    
    for ( Mecable * mec : mecables )
    {
        const index_t bks = DIM * mec->nbPoints();
        extractBlock(mec, tmp2);
        VecPrint::sparse_off(f, bks, bks, tmp2, bks, DIM*mec->matIndex());
    }
    os.close();
    
    free_real(tmp1);
    free_real(tmp2);
}


//------------------------------------------------------------------------------
#pragma mark - Binary Export

static void dumpVector(FILE * fp, size_t dim, real* vec, bool nat)
{
    static float * low = nullptr;
    static size_t alc = 0;
    
    if ( !fp )
    {
        delete[] low;
        low = nullptr;
        return;
    }
    if ( !nat && std::is_same<real, double>::value )
    {
        if ( dim > alc )
        {
            delete[] low;
            low = new float[dim];
            alc = dim;
        }
        copy_real(dim, vec, low);
        fwrite(low, sizeof(float), dim, fp);
    }
    else
        fwrite(vec, sizeof(real), dim, fp);
}


void Meca::dumpObjectID(FILE * fp) const
{
    uint32_t * vec = new uint32_t[largestMecable()];
    
    uint32_t i = 1;
    for ( Mecable const* mec : mecables )
    {
        const index_t nbp = mec->nbPoints();
        for ( index_t p = 0; p < nbp; ++p )
            vec[p] = i;
        for ( int d = 0; d < DIM; ++d )
            fwrite(vec, sizeof(uint32_t), nbp, fp);
        ++i;
    }
    
    delete[](vec);
}


void Meca::dumpMobility(FILE * fp, bool nat) const
{
    real * vec = new_real(largestMecable());
    
    for ( Mecable const* mec : mecables )
    {
        const index_t nbp = mec->nbPoints();
        const real val = mec->pointMobility();
        for ( index_t p=0; p < nbp; ++p )
            vec[p] = val;
        for ( int d = 0; d < DIM; ++ d )
            dumpVector(fp, nbp, vec, nat);
    }
    
    free_real(vec);
}


/**
 Save the full matrix associated with multiply(), in binary format
 */
template < Meca::MultiplyFuncPtr MULTIPLY >
void Meca::dumpMatrix(FILE * fp, bool nat) const
{
    const size_t dim = dimension();
    real * src = new_real(dim);
    real * res = new_real(dim);
    
    zero_real(dim, src);
    
    for ( size_t ii = 0; ii < dim; ++ii )
    {
        src[ii] = 1;
        (this->*MULTIPLY)(src, res);
        dumpVector(fp, dim, res, nat);
        src[ii] = 0;
    }
    
    free_real(res);
    free_real(src);
}


/**
 Save the projection matrix multiplied by the mobility, in binary format
 */
void Meca::dumpProjection(FILE * fp, bool nat) const
{
    const size_t dim = dimension();
    real * vec = new_real(dim);
        
    for ( size_t i = 0; i < dim; ++i )
    {
        zero_real(dim, vec);
        vec[i] = 1;
        
        for ( Mecable const* mec : mecables )
        {
            const size_t inx = DIM * mec->matIndex();
            // this includes the mobility, but not the time_step:
            mec->projectForces(vec+inx, vec+inx);
            blas::xscal(DIM*mec->nbPoints(), mec->leftoverMobility(), vec+inx, 1);
        }
        // write column to fp directly:
        dumpVector(fp, dim, vec, nat);
    }
    
    free_real(vec);
}


/**
 Save matrix associated with the preconditionner, in binary format
 This relies on `Meca::precondition()`, which may apply a dummy preconditionner
 */
void Meca::dumpPreconditionner(FILE * fp, bool nat) const
{
    const size_t dim = dimension();
    real * vec = new_real(dim);
    
    for ( size_t i = 0; i < dim; ++i )
    {
        zero_real(dim, vec);
        vec[i] = 1;
        for ( Mecable const* mec : mecables )
        {
            const size_t inx = DIM * mec->matIndex();
            applyPreconditionner(mec, vec+inx);
        }
        dumpVector(fp, dim, vec, nat);
    }
    
    free_real(vec);
}


/**
 This dump the total matrix and some vectors in binary files.
 
 This MATLAB code should read the output:
 
     ord = load('ord.txt');
     time_step = load('stp.txt');
     precision = 'double' % or float?
     obj = fread(fopen('obj.bin'), ord, 'uint32');
     drg = fread(fopen('drg.bin'), ord, precision);
     sys = fread(fopen('sys.bin'), [ord, ord], precision);
     ela = fread(fopen('ela.bin'), [ord, ord], precision);
     mob = fread(fopen('mob.bin'), [ord, ord], precision);
     con = fread(fopen('con.bin'), [ord, ord], precision);
     pts = fread(fopen('pts.bin'), ord, precision);
     rhs = fread(fopen('rhs.bin'), ord, precision);
     sol = fread(fopen('sol.bin'), ord, precision);
 
 To display the matrices:

     imshow(abs(sys))
     imshow(abs(ela))
 
 You can then compare the results with matlab's own iterative method,
 and compare the result using a scatter plot:
 
     x = bicgstab(sys, rhs, 0.001, ord);
     plot(x, sol, '.');
 
 */
void Meca::dumpSystem(bool nat) const
{
    FILE * f = FilePath::open_file("ord.txt", "w");
    fprintf(f, "%u %u %lu\n", dimension(), DIM, sizeof(real));
    fclose(f);
    
    f = FilePath::open_file("stp.txt", "w");
    fprintf(f, "%.12f %.12f\n", tau_, tolerance_);
    fclose(f);
    
    f = FilePath::open_file("mob.bin", "wb");
    dumpMobility(f, nat);
    fclose(f);
    
    f = FilePath::open_file("obj.bin", "wb");
    dumpObjectID(f);
    fclose(f);
    
    f = FilePath::open_file("rhs.bin", "wb");
    dumpVector(f, dimension(), vRHS, nat);
    fclose(f);
    
    f = FilePath::open_file("sol.bin", "wb");
    dumpVector(f, dimension(), vSOL, nat);
    fclose(f);
    
    f = FilePath::open_file("pts.bin", "wb");
    dumpVector(f, dimension(), vPTS, nat);
    fclose(f);
    
    f = FilePath::open_file("sys.bin", "wb");
    dumpMatrix<&Meca::multiply>(f, nat);
    fclose(f);
    
    f = FilePath::open_file("ela.bin", "wb");
    dumpMatrix<&Meca::multiplyElasticity>(f, nat);
    fclose(f);
    
    f = FilePath::open_file("prj.bin", "wb");
    dumpProjection(f, nat);
    fclose(f);
    
    f = FilePath::open_file("con.bin", "wb");
    dumpPreconditionner(f, nat);
    fclose(f);
    
    dumpVector(nullptr, 0, nullptr, nat);
}


//------------------------------------------------------------------------------
#pragma mark - Bitmap Export


// Just considering Couple between Fibers here:
void markConnectivity(BitMap<1>& bmap, Array<Mecable*> const& mecs)
{
    bmap.clear();
    ObjectFlag i = 0;
    for ( Mecable * mec : mecs )
        mec->flag(i++);
    if ( (size_t)i != mecs.size() )
        throw InvalidParameter("ObjectFlag overflow in markConnectivity()");

    for ( Mecable const* mec : mecs )
    {
        i = mec->flag();

        Fiber const* fib = Fiber::toFiber(mec);
        for ( Hand const* h = fib->firstHand(); h; h = h->next() )
        {
            Hand const* g = h->otherHand();
            if ( g > h  &&  g->attached() )
            {
                ObjectFlag j = g->fiber()->flag();
                bmap.set(i, j, 1);
            }
            else if ( g )
            {
                
            }
        }
    }
}


void Meca::saveConnectivityBitmap() const
{
    static unsigned cnt = 0;
    const index_t nbv = mecables.size();
    BitMap<1> bmap(nbv, nbv);
    char str[32] = { 0 };
    
    snprintf(str, sizeof(str), "net%08u.bmp", cnt++);
    FILE * f = fopen(str, "w");
    if ( f ) {
        if ( !ferror(f) ) {
            markConnectivity(bmap, mecables);
            bmap.save(f);
        }
        fclose(f);
    }
}


template < typename MATRIX >
static void markMatrix(BitMap<1>& bmap, index_t sup, MATRIX const& mat)
{
    for ( index_t j = 0; j < sup; ++j )
    {
        for ( index_t n = 0; n < mat.column_size(j); ++n )
        {
            index_t i = mat.column_index(j, n);
            if ( i < sup )
                bmap.set(i, j, 1);
        }
    }
}

/// add diagonal elements
[[maybe_unused]]
static void markDiagonal(BitMap<1>& bmap, index_t mag, Array<Mecable*> const& mecs)
{
    for ( Mecable * mec : mecs )
    {
        index_t i = mag * mec->matIndex();
        index_t k = mag * mec->nbPoints();
        for ( index_t j = 0; j < k; ++j )
            bmap.set_if(i+j, i+j, 1);
    }
}

/// add vertical and horizontal lines to indicate mecables indices
static void markMecables(BitMap<1>& bmap, index_t mag, Array<Mecable*> const& mecs)
{
    for ( Mecable * mec : mecs )
    {
        index_t i = mag * mec->matIndex();
        for ( index_t j = 1; j < 6; ++j )
            bmap.set_if(i-j-1, i+j, 1);
    }
}

template < typename MATRIX >
static void saveMatrixBitmap(MATRIX& mat, Array<Mecable*> mecs, index_t nbv, const char str[])
{
    BitMap<1> bmap(nbv, nbv);
    FILE * f = fopen(str, "w");
    if ( f ) {
        if ( !ferror(f) ) {
            bmap.clear();
            markMatrix(bmap, nbv, mat);
            index_t mag = DIM*nbv/mat.size();
            //markDiagonal(bmap, mag, mecs);
            markMecables(bmap, mag, mecs);
            bmap.save(f);
        }
        fclose(f);
    }
}

void Meca::saveMatrixBitmaps(const char prefix[], unsigned inc) const
{
    static unsigned cnt = 0;
    char str[64] = { 0 };
    
#if USE_ISO_MATRIX
    snprintf(str, sizeof(str), "%siso%08u.bmp", prefix, cnt);
    saveMatrixBitmap(mISO, mecables, mISO.num_columns(), str);
#endif
    snprintf(str, sizeof(str), "%sful%08u.bmp", prefix, cnt);
    saveMatrixBitmap(mFUL, mecables, mFUL.num_columns(), str);
    cnt += inc;
}

