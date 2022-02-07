// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University


//------------------------------------------------------------------------------
#pragma mark - Connectivity Analysis


/// equalize flags for any existing matrix element between Mecables
template < typename MatrixClass >
void computeClusters(Array<Mecable*> const mecables, const size_t MAX, Mecable** table,
                     MatrixClass const& MAT, size_t ORD)
{
    for ( size_t j = 0; j < MAX; ++j )
    {
        size_t jj = ORD * j;
        Mecable const* A = table[j];
        for ( size_t n = 0; n < MAT.column_size(jj); ++n )
        {
            // we do not check the value here, but just having a block
            size_t i = MAT.column_index(jj, n) / ORD;
            assert_true( i < MAX );
            Mecable const* B = table[i];
            if ( A->flag() != B->flag() )
            {
                ObjectFlag a = std::min(A->flag(), B->flag());
                ObjectFlag b = std::max(A->flag(), B->flag());
                // replace b -> a everywhere:
                for ( Mecable * mec : mecables )
                {
                    if ( mec->flag() == b )
                        mec->flag(a);
                }
            }
        }
    }
}


/** Assuming that Mecable::flag() have been set already */
void Meca::flagClusters() const
{
    const size_t MAX = nbVertices();
    Mecable ** table = new Mecable*[MAX]{nullptr};
    
    for ( Mecable * mec : mecables )
    {
        const size_t inx = mec->matIndex();
        const size_t end = mec->nbPoints() + inx;
        assert_true( end <= MAX );
        for ( size_t i = inx; i < end; ++i )
            table[i] = mec;
    }
    
#if USE_MATRIX_BLOCK
    computeClusters(mecables, MAX, table, mFUL, DIM);
#endif
#if USE_ISO_MATRIX
    computeClusters(mecables, MAX, table, mISO, 1);
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

/**
 Extract the full matrix associated with Meca::multiply().
 The array `mat[]` should be preallocated to hold `dim*lda` real scalars,
 with `dim >= Meca::dimension()`, and `lda >= dim` the leading dimension of
 the array.
 */
void Meca::getMatrix(real * mat, size_t lda) const
{
    size_t dim = dimension();
    if ( lda < dim )
        throw InvalidIO("invalid matrix dimensions");
    real * src = new_real(dim);
    zero_real(dim, src);
    
    for ( size_t j = 0; j < dim; ++j )
    {
        src[j] = 1;
        multiply(src, mat+j*lda);
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
        const size_t nbp = DIM * mec->nbPoints();
        for ( size_t p = 0; p < nbp; ++p )
            fprintf(fp, "%if\n", i);
        ++i;
    }
}

void Meca::saveMobility(FILE * fp) const
{
    for ( Mecable const* mec : mecables )
    {
        const size_t nbp = mec->nbPoints();
        const real val = mec->pointMobility();
        for ( size_t p = 0; p < DIM * nbp; ++p )
            fprintf(fp, "%f\n", val);
    }
}

/**
 Save a sparse matrix in Matrix Market format
 https://math.nist.gov/MatrixMarket/formats.html
 This is a Sparse text format
 */
void Meca::saveMatrix(FILE * fp, real threshold) const
{
    fprintf(fp, "%%%%MatrixMarket matrix coordinate real general\n");
    fprintf(fp, "%% This is a matrix produced by Cytosim\n");
    fprintf(fp, "%% author: Francois J. Nedelec\n");
    fprintf(fp, "%% kind: biological cell simulation (cytoskeleton)\n");

    const size_t dim = dimension();
    real * src = new_real(dim);
    real * dst = new_real(dim);
    zero_real(dim, src);
    
    fprintf(fp, "%lu %lu ", dim, dim);

    fpos_t pos;
    fgetpos(fp, &pos);
    size_t cnt = 0;
    fprintf(fp, "%10lu\n", cnt);

    for ( size_t j = 0; j < dim; ++j )
    {
        src[j] = 1;
        multiply(src, dst);
        for ( size_t i = 0; i < dim; ++i )
            if ( abs_real(dst[i]) > threshold )
            {
                fprintf(fp, "%3lu %3lu %f\n", i, j, dst[i]);
                ++cnt;
            }
        src[j] = 0;
    }
    
    fsetpos(fp, &pos);
    fprintf(fp, "%10lu\n", cnt);

    free_real(dst);
    free_real(src);
}


/**
 Save Matrix and Right-hand-side Vector
 */
void Meca::saveSystem() const
{
    FILE * f = FilePath::open_file("matrix.mtx", "w");
    saveMatrix(f, 0);
    fclose(f);
    
    f = FilePath::open_file("vector.mtx", "w");
    saveVector(f, dimension(), vRHS);
    fclose(f);
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
    fprintf(f, "%lu %i %lu\n", dimension(), DIM, sizeof(real));
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
    
    std::ofstream os("sol.txt");
    VecPrint::dump(os, dimension(), vPTS);
    os.close();
    
    os.open("rhs.txt");
    VecPrint::dump(os, dimension(), vRHS);
    os.close();
    
#if USE_ISO_MATRIX
    os.open("iso.txt");
    mISO.printSparse(os, 0);
    os.close();
#endif
    
    os.open("full.txt");
    mFUL.printSparse(os, 0);
    os.close();
        
    size_t alc = 0;
    for ( Mecable const* mec : mecables )
        alc = std::max(alc, mec->nbPoints());

    real * tmp1 = new_real(DIM*alc);
    real * tmp2 = new_real(DIM*DIM*alc*alc);
    
    os.open("diag.txt");
    
    for ( Mecable * mec : mecables )
    {
        const size_t bks = DIM * mec->nbPoints();
        extractBlock(mec, tmp2);
        VecPrint::sparse_off(os, bks, bks, tmp2, bks, DIM*mec->matIndex());
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
        const size_t nbp = mec->nbPoints();
        for ( size_t p = 0; p < nbp; ++p )
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
        const size_t nbp = mec->nbPoints();
        const real val = mec->pointMobility();
        for ( size_t p=0; p < nbp; ++p )
            vec[p] = val;
        for ( int d = 0; d < DIM; ++ d )
            dumpVector(fp, nbp, vec, nat);
    }
    
    free_real(vec);
}


/**
 Save the full matrix associated with multiply(), in binary format
 */
void Meca::dumpMatrix(FILE * fp, bool nat) const
{
    const size_t dim = dimension();
    real * src = new_real(dim);
    real * res = new_real(dim);
    
    zero_real(dim, src);
    
    for ( size_t ii = 0; ii < dim; ++ii )
    {
        src[ii] = 1;
        multiply(src, res);
        dumpVector(fp, dim, res, nat);
        src[ii] = 0;
    }
    
    free_real(res);
    free_real(src);
}


/**
 Save the elasticity matrix, in binary format
 */
void Meca::dumpElasticity(FILE * fp, bool nat) const
{
    const size_t dim = dimension();
    real * src = new_real(dim);
    real * res = new_real(dim);
    
    zero_real(dim, src);
    
    for ( size_t ii = 0; ii < dim; ++ii )
    {
        src[ii] = 1;
        
        mFUL.vecMul(src, res);
#if USE_ISO_MATRIX
        mISO.VECMULADDISO(src, res);
#endif
#if SEPARATE_RIGIDITY_TERMS
        addAllRigidity(src, res);
#endif

#if ADD_PROJECTION_DIFF
        for ( Mecable const* mec : mecables )
        {
            if ( mec->hasProjectionDiff() )
            {
                const size_t inx = DIM * mec->matIndex();
                mec->addProjectionDiff(src+inx, res+inx);
            }
        }
#endif
        
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
    fprintf(f, "%lu %i %lu\n", dimension(), DIM, sizeof(real));
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
    dumpMatrix(f, nat);
    fclose(f);
    
    f = FilePath::open_file("ela.bin", "wb");
    dumpElasticity(f, nat);
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
#pragma mark - Matrix Bitmap Export

#include "../base/save_bitmap.cc"


// Just considering Couple between Fibers here:
void setConnectivity(BitMap<1>& bmap, Array<Mecable*> const& mecables)
{
    bmap.clear();
    unsigned i = 0;
    for ( Mecable * mec : mecables )
        mec->flag(i++);
    
    for ( Mecable const* mec : mecables )
    {
        i = mec->flag();

        Fiber const* fib = Fiber::toFiber(mec);
        for ( Hand * h = fib->firstHand(); h; h = h->next() )
        {
            HandMonitor const* m = h->monitor();
            Hand const* oh = m->otherHand(h);
            if ( oh > h  &&  oh->attached() )
            {
                unsigned j = oh->fiber()->flag();
                bmap.set(i, j, 15);
            }
            else if ( oh )
            {
                
            }
        }
    }
}


void Meca::saveConnectivityBitmap() const
{
    static size_t cnt = 0;
    const size_t nbv = mecables.size();
    BitMap<1> bmap(nbv, nbv);
    char str[32] = { 0 };
    
    snprintf(str, sizeof(str), "net%08lu.bmp", cnt++);
    FILE * f = fopen(str, "w");
    if ( f ) {
        if ( !ferror(f) ) {
            setConnectivity(bmap, mecables);
            bmap.save(f);
        }
        fclose(f);
    }
}


template < typename MatrixClass >
static void setMatrixBitmap(BitMap<1>& bmap, size_t nbv, MatrixClass const& MAT, size_t ORD)
{
    bmap.clear();
    for ( size_t j = 0; j < nbv; ++j )
    {
        size_t jj = j * ORD;
        for ( size_t n = 0; n < MAT.column_size(jj); ++n )
        {
            size_t i = MAT.column_index(jj, n) / ORD;
            // swap i,j and flip i to display matrix properly on screen
            bmap.set(i, j, 1);
        }
    }
}


/// add vertical and horizontal lines to indicate mecables indices
static void markMecables(BitMap<1>& bmap, Array<Mecable*> const& mecables)
{
    for ( Mecable * mec : mecables )
    {
        size_t i = mec->matIndex();
        size_t s = i + mec->nbPoints() - 1;
        for ( size_t j = i+1; j < s; ++j )
        {
            bmap.set(i, j, 1);
            bmap.set(j, s, 1);
        }
        bmap.set(i, s, 1);
    }
}


void Meca::saveMatrixBitmaps() const
{
    static size_t cnt = 0;
    const size_t nbv = nbVertices();
    BitMap<1> bmap(nbv, nbv);
    char str[32] = { 0 };
    
#if USE_ISO_MATRIX
    snprintf(str, sizeof(str), "iso%08lu.bmp", cnt);
    FILE * f = fopen(str, "w");
    if ( f ) {
        if ( !ferror(f) ) {
            setMatrixBitmap(bmap, nbv, mISO, 1);
            markMecables(bmap, mecables);
            bmap.save(f);
        }
        fclose(f);
    }
#endif
    snprintf(str, sizeof(str), "ful%08lu.bmp", cnt++);
    FILE * g = fopen(str, "w");
    if ( g ) {
        if ( !ferror(g) ) {
            setMatrixBitmap(bmap, nbv, mFUL, DIM);
            markMecables(bmap, mecables);
            bmap.save(g);
        }
        fclose(g);
    }
}

