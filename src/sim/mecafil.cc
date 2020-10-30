// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "dim.h"
#include "sim.h"
#include "mecafil.h"
#include "blas.h"
#include "lapack.h"
#include "random.h"
#include "vecprint.h"
//#include "cytoblas.h"

//------------------------------------------------------------------------------
Mecafil::Mecafil()
{
    buildProjection();
    iPointMobility = 0;
    iRigidity = 0;
    iDir = nullptr;
    iLag = nullptr;
    iLLG = nullptr;
#if NEW_ANISOTROPIC_FIBER_DRAG
    iAni = nullptr;
#endif
    useProjectionDiff = false;
}

void Mecafil::release()
{
    destroyProjection();
    free_real(iDir);
    iDir = nullptr;
    iLag = nullptr;
    iLLG = nullptr;
#if NEW_ANISOTROPIC_FIBER_DRAG
    iAni = nullptr;
#endif
}

Mecafil::~Mecafil()
{
    release();
}


//------------------------------------------------------------------------------
Mecafil::Mecafil(Mecafil const&)
{
    ABORT_NOW("unfinished: cannot copy a Mecafil");
}


Mecafil& Mecafil::operator=(Mecafil const&)
{
    ABORT_NOW("unfinished: cannot copy a Mecafil");
}


//------------------------------------------------------------------------------
size_t Mecafil::allocateMecable(const size_t nbp)
{
    size_t ms = Mecable::allocateMecable(nbp);
    /*
     if Mecable::allocateMecable() allocated memory, it will return the 
     size of the new array, and we allocate the same size for other arrays.
     */
    if ( ms )
    {
        //std::clog << "Mecafil::allocatePoints " << ms << '\n';
        allocateProjection(ms);
        
        // allocate memory:
        free_real(iDir);
        
#if NEW_ANISOTROPIC_FIBER_DRAG
        // allocations: iDir=DIM*N  iLag=N  iLLG=DIM*N  iAni=DIM*N
        iDir = new_real(ms*(3*DIM+1));
        iLag = iDir + ms*DIM;
        iLLG = iLag + ms;
        iAni = iLLG + ms*DIM;
#else
        // allocations: iDir=DIM*N  iLag=N  iLLG=N
        iDir = new_real(ms*(DIM+2));
        iLag = iDir + ms*DIM;
        iLLG = iLag + ms;
#endif
        
        // reset Lagrange multipliers
        zero_real(ms, iLag);
        zero_real(ms, iLLG); // generally not needed
    }
    return ms;
}


//------------------------------------------------------------------------------
#pragma mark -

/**
 The argument should be: sc = kT / timestep;
 */
real Mecafil::addBrownianForces(real const* rnd, real alpha, real* rhs) const
{
    real b = std::sqrt( 2 * alpha / iPointMobility );

    for ( size_t jj = 0; jj < DIM*nPoints; ++jj )
        rhs[jj] += b * rnd[jj];
    
    return b * iPointMobility;
}


//------------------------------------------------------------------------------

/**
 Calculate the normalized difference between successive vertices of the fiber:

     const real alpha = 1.0 / segmentation();
     for ( int n = 0; n < DIM*lastPoint(); ++n )
         iDir[n] = alpha * ( pPos[n+DIM] - pPos[n] );

 */

void Mecafil::storeDirections()
{
    //iDirValid = true;
#if ( 1 )
    /*
     we assume here that successive points are correctly separated by 'segmentation',
     such that we can normalize the vector simply by dividing by 'segmentation'
     */
    //checkSegmentation(0.01);
    const real val = 1.0 / segmentation();
    const size_t end = DIM * lastPoint();
    #pragma ivdep
    for ( size_t i = 0; i < end; ++i )
        iDir[i] = val * ( pPos[i+DIM] - pPos[i] );
#else
    for ( size_t p = 0; p < lastPoint(); ++p )
        normalize(diffPoints(p)).store(iDir+DIM*p);
#endif
    
#if NEW_ANISOTROPIC_FIBER_DRAG
    /*
     Calculate the average filament direction at each vertex
     */
    const real* dir = iDir;
    // for the extremities, copy direction of the nearby segment:
    for ( size_t d = 0; d < DIM; ++d )
    {
        iAni[d]     = dir[d];
        iAni[d+end] = dir[d+end-DIM];
    }
    
    // for inner vertices, average directions of the flanking segments:
    for ( size_t p = DIM ; p < end; ++p )
        iAni[p] = 0.5 * ( dir[p-DIM] + dir[p] );

    //VecPrint::print(std::clog, last+DIM, iAni);
#endif
}

//------------------------------------------------------------------------------
#pragma mark - Project

#if ( DIM > 1 )

#  if PROJECT_WITH_MATRIX
#     include "mecafil_projectmat.cc"
#     warning "projection matrices are built explicitly"
#  else
#     include "mecafil_project.cc"
#  endif

#else

void Mecafil::buildProjection()   {}  //DIM == 1
void Mecafil::makeProjection()    {}  //DIM == 1
void Mecafil::destroyProjection() {}  //DIM == 1
void Mecafil::allocateProjection(size_t) {}  //DIM == 1

void Mecafil::projectForces(const real* X, real* Y) const
{
    real sum = X[0];
    for ( size_t ii = 1; ii < nPoints; ++ii )
        sum += X[ii];
    
    sum = sum / (real)nPoints;
#if NEW_ANISOTROPIC_FIBER_DRAG
    sum *= 2;
#endif
    for ( size_t ii = 0; ii < nPoints; ++ii )
        Y[ii] = sum;
}

void Mecafil::computeTensions(const real*) {} //DIM == 1
void Mecafil::makeProjectionDiff(const real*) {} //DIM == 1
void Mecafil::addProjectionDiff(const real*, real*) const {} //DIM == 1

#endif


void Mecafil::printTensions(FILE * out, char c) const
{
    fprintf(out, "\n%c%s ", c, reference().c_str());
    VecPrint::print(out, nbSegments(), iLag, 2);
    fprintf(out, "  fM"); netForceEndM().print(out);
    fprintf(out, "  fP"); netForceEndP().print(out);
}


void Mecafil::getForces(const real* ptr)
{
    Mecable::getForces(ptr);
    //fprintf(stderr, "\nF "); VecPrint::print(stderr, DIM*nbPoints(), ptr, 2, DIM);
    computeTensions(ptr);
    //printTensions(stderr);
}

//-----------------------------------------------------------------------
#pragma mark -

/*
 This is the reference implementation
 */
void add_rigidity0(const size_t nbt, const real* X, const real rigid, real* Y)
{
    assert_true( X != Y );
    for ( size_t jj = 0; jj < nbt; ++jj )
    {
        real f = rigid * (( X[jj+DIM*2] - X[jj+DIM] ) - ( X[jj+DIM] - X[jj] ));
        Y[jj      ] -= f;
        Y[jj+DIM  ] += f*2.0;
        Y[jj+DIM*2] -= f;
    }
}

/*
 This is an optimized implementation
 */
void add_rigidityF(const size_t nbt, const real* X, const real R1, real* Y)
{
    assert_true(nbt > DIM);
    const real R2 = R1 * 2;
    const real R4 = R1 * 4;
    const real R6 = R1 * 6;
    
    const size_t end = nbt;
    // in the general case all values can be computed independently:
    #pragma ivdep
    for ( size_t i = DIM*2; i < end; ++i )
        Y[i] += R4 * ((X-DIM)[i]+(X+DIM)[i]) - R1 * ((X-DIM*2)[i]+(X+DIM*2)[i]) - R6 * X[i];

    // special cases at the edges:
    real      * Z = Y + nbt;
    real const* E = X + nbt + DIM;
    #pragma ivdep
    for ( size_t d = 0; d < DIM; ++d )
    {
        Y[d+DIM] -= R1 * ((X+DIM)[d]+(X+DIM*3)[d]) + R4 * ((X+DIM)[d]-(X+DIM*2)[d]) - R2 * X[d];
        Z[d    ] -= R1 * ((E-DIM)[d]+(E-DIM*3)[d]) + R4 * ((E-DIM)[d]-(E-DIM*2)[d]) - R2 * E[d];
        Y[d    ] -= R1 * ((X+DIM*2)[d]+X[d]) - R2 * (X+DIM)[d];
        Z[d+DIM] -= R1 * ((E-DIM*2)[d]+E[d]) - R2 * (E-DIM)[d];
    }
}

/**
 Add rigidity terms between three points {A, B, C}
 Done with Serge DMITRIEFF, 2015
 */
void add_rigidity(size_t A, size_t B, size_t C, const real* X, const real R1, real* Y)
{
#if ( DIM > 1 )
    const real R2 = 2.0 * R1;
    for ( size_t d = 0; d < DIM; ++ d )
    {
        real f = 2.0 * X[B*DIM+d] - ( X[A*DIM+d] + X[C*DIM+d] );
        Y[A*DIM+d] += f * R1;
        Y[B*DIM+d] -= f * R2;
        Y[C*DIM+d] += f * R1;
    }
#endif
}


/**
 This is the bending elasticity terms, as obtained by derivation of the
 Hamiltonian representing bending elasticity.

     F1 = k * ( t1 * dot(t1, t2) - t2 )
     F3 = k * ( t1 - dot(t1, t2) * t2 )
     F2 = -F1 -F3
 
 These forces are normal to the segments: dot(F1, t1) = dot(F3, t2) = 0
 The cosinus are obtained here from the normalized difference vector 'dir'.

 Ivan Hornak & Heiko Rieger in:
     Stochastic Model of T Cell Repolarization during Target Elimination
     https://doi.org/10.1016/j.bpj.2020.01.045
 claimed that this would lead to a better estimation of bending elasticity.
 However, this is not true, and using these formula makes strictly no difference,
 because compared to our standard implementation:
 
     F1 = k * ( t1 - t2 )
     F3 = k * ( t1 - t2 )
     F2 = -F1 -F3

 the forces only differ by a vector that is tangent to the segments, and any such
 tangent force is fully absorbed by the constraints imposed on the lengths of the
 segments. Thus there is no advantage in using these (more exact) formula.
 It makes no difference.
 */
void add_rigidityN(const size_t nbt, const real* X, const real rigid, real* Y, real const* dir)
{
    assert_true( X != Y );
    for ( size_t jj = 0; jj < nbt; jj+=DIM )
    {
        // cosine of the angle between two consecutive segments:
        const real C = dot(Vector(dir+jj), Vector(dir+jj+DIM));
        for ( int d = 0; d < DIM; ++d )
        {
            int i = jj + d;
            real f1 = rigid * (C * ( X[i+DIM] - X[i] ) - ( X[i+DIM*2] - X[i+DIM] ));
            real f3 = rigid * (( X[i+DIM] - X[i] ) - C * ( X[i+DIM*2] - X[i+DIM] ));
            Y[i      ] += f1;
            Y[i+DIM  ] -= f1+f3;
            Y[i+DIM*2] += f3;
        }
    }
}

//------------------------------------------------------------------------------

#define CHECK_RIGIDITY 0

/**
 calculates the second-derivative of point's coordinates,
 scale by the rigidity term, and add to vector Y
*/
void Mecafil::addRigidity(const real* X, real* Y) const
{
#if CHECK_RIGIDITY
    // compare to default implementation:
    real * tmp = new_real(DIM*nPoints);
    copy_real(DIM*nPoints, Y, tmp);
    add_rigidity0(DIM*(nPoints-2), X, iRigidity, tmp);
#endif
    if ( nPoints > 3 )
    {
        const size_t nbt = DIM * ( nPoints - 2 );  // number of triplet values

#if ( DIM == 2 ) && REAL_IS_DOUBLE && defined(__AVX__)
        add_rigidityF(nbt, X, iRigidity, Y);
#elif ( DIM == 2 ) && REAL_IS_DOUBLE && defined(__SSE3__)
        add_rigiditySSE(nbt, X, iRigidity, Y);
#elif ( DIM > 1 )
        add_rigidityF(nbt, X, iRigidity, Y);
#endif
    
#if NEW_FIBER_LOOP
        if ( iRigidityLoop )
        {
            /*
             With Serge DMITRIEFF:
             Link first and last point in the same way as all other points,
             making the fiber mechanically homogeneous and all points equivalent
             */
            const size_t L = lastPoint();
            add_rigidity(L,   0, 1, X, iRigidity, Y);
            add_rigidity(L-1, L, 0, X, iRigidity, Y);
        }
#endif
    }
    else if ( nPoints > 2 )
    {
        //add_rigidityN(1, X, iRigidity, Y, iDir);
        add_rigidity(0, 1, 2, X, iRigidity, Y);
    }
    
#if CHECK_RIGIDITY
    static size_t cnt = 0;
    real err = blas::max_diff(DIM*nPoints, tmp, Y);
    if ( err > 1.0e-6 || ++cnt > 100 )
    {
        cnt = 0;
        printf("addRigidity(%lu) error %e\n", nPoints, err);
    }
    free_real(tmp);
#endif
}

