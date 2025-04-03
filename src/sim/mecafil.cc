// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "dim.h"
#include "mecafil.h"
#include "blas.h"
#include "lapack.h"
#include "random.h"
#include "vecprint.h"
//#include "cytoblas.h"

//------------------------------------------------------------------------------
Mecafil::Mecafil()
{
    initProjection();
    iPointMobility = 0;
    iRigidity = 0;
    iDir = nullptr;
    iLag = nullptr;
    iLLG = nullptr;
}

void Mecafil::release()
{
    destroyProjection();
    iDir = nullptr;
    iLag = nullptr;
    iLLG = nullptr;
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


Mecafil& Mecafil::operator = (Mecafil const&)
{
    ABORT_NOW("unfinished: cannot copy a Mecafil");
}


//------------------------------------------------------------------------------
void Mecafil::allocateMecable(const index_t nbp)
{
    index_t add = ADD_PROJECTION_DIFF ? 3 : 2;
    index_t top = DIM+2;
    real * ptr = Mecable::allocateMemory(nbp, add+top);
    /*
     if Mecable::allocateMecable() allocated memory, it will return a
     non-zero pointer, with extra memory usable for pointers here.
     */
    if ( ptr )
    {
        size_t all = allocated();
        //std::clog << "Mecafil::allocateMecable " << ms << '\n';
        allocateProjection(all, ptr);
        
        // allocations: iDir=DIM*N  iLag=N  iLLG=N
        iDir = ptr + all*add;
        iLag = iDir + all*DIM;
        iLLG = iLag + all;
        
        // reset Lagrange multipliers
        zero_real(all, iLag);
        zero_real(all, iLLG); // generally not needed
    }
}


//------------------------------------------------------------------------------
#pragma mark -

/**
 On entry, rhs[] should contain unit Brownian random numbers, and fce[] the forces
 On exit, rhs[] is updated to include Brownian components: rhs <- fce + b * random()
 The argument `alpha = 2 * kT / timestep` is used to calculate `b`
 */
real Mecafil::addBrownianForces(real const* fce, real alpha, real* rhs) const
{
    real b = std::sqrt( alpha / iPointMobility );

    for ( size_t j = 0; j < DIM*nPoints; ++j )
        rhs[j] = fce[j] + b * rhs[j];
    
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
#if ( 1 )
    /*
     we assume here that successive points are correctly separated by 'segmentation',
     such that we can normalize the vector simply by dividing by 'segmentation'
     */
    const real alpha = 1.0 / segmentation();
    const size_t end = DIM * lastPoint();
    #pragma omp simd
    for ( size_t i = 0; i < end; ++i )
        iDir[i] = alpha * ( pPos[i+DIM] - pPos[i] );
#else
    for ( index_t p = 0; p < lastPoint(); ++p )
        normalize(diffPoints(p)).store(iDir+DIM*p);
#endif
    //VecPrint::print("iDir", end, iDir);
}

//------------------------------------------------------------------------------
#pragma mark - Project

#if ( DIM > 1 )

#include "mecafil_project.cc"

#else

void Mecafil::initProjection() {}  //DIM == 1
void Mecafil::makeProjection() {}  //DIM == 1
void Mecafil::destroyProjection() {}  //DIM == 1
void Mecafil::allocateProjection(size_t, real*) {}  //DIM == 1

void Mecafil::projectForces(const real* X, real* Y) const
{
    real sum = X[0];
    for ( size_t ii = 1; ii < nPoints; ++ii )
        sum += X[ii];
    
    sum = sum / (real)nPoints;
    for ( size_t ii = 0; ii < nPoints; ++ii )
        Y[ii] = sum;
}

void Mecafil::computeTensions(const real*) {} //DIM == 1

#if ADD_PROJECTION_DIFF
void Mecafil::makeProjectionDiff(const real*) {} //DIM == 1
void Mecafil::addProjectionDiff(const real*, real*) const {} //DIM == 1
void Mecafil::addProjectionDiff(real*) const {} //DIM == 1
#endif
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
    if ( ptr ) computeTensions(ptr);
    //printTensions(stderr);
}

//-----------------------------------------------------------------------
#pragma mark -

/*
 This is the reference implementation
 */
void add_rigidity0(const int nbp, const real* X, const real R1, real* Y)
{
    assert_true(R1 >= 0);
    const real R2 = 2.0 * R1;
    const real two = 2.0;
    const unsigned end = DIM * ( nbp - 2 );
    #pragma omp simd
    for ( unsigned jj = 0; jj < end; ++jj )
    {
        real f = ( X[jj] + X[jj+DIM*2] ) - two * X[jj+DIM];
        Y[jj      ] -= f * R1;
        Y[jj+DIM  ] += f * R2;
        Y[jj+DIM*2] -= f * R1;
    }
}

/*
 This is a different implementation
 */
inline void add_rigidityF(const int nbp, const real* X, const real R1, real* Y)
{
    assert_true(R1 >= 0);
    const real six = 6.0;
    const real R4 = R1 * 4.0;
    const real R2 = R1 * 2.0;

    const size_t end = DIM * ( nbp - 2 );
    #pragma omp simd
    for ( size_t i = DIM*2; i < end; ++i )
        Y[i] = Y[i] + R4 * (X[i-DIM]+X[i+DIM]) - R1 * (six*X[i]+(X[i-DIM*2]+X[i+DIM*2]));
    
    // special cases near the edges:
    real      * Z = Y + DIM * ( nbp - 1 );
    real const* E = X + DIM * ( nbp - 1 );
    if ( nbp > 3 )
    {
        #pragma omp simd
        for ( int d = 0; d < DIM; ++d )
        {
            Y[d+DIM] -= R1 * (X[d+DIM]+X[d+DIM*3]) + R4 * (X[d+DIM]-X[d+DIM*2]) - R2 * X[d];
            Z[d-DIM] -= R1 * (E[d-DIM]+E[d-DIM*3]) + R4 * (E[d-DIM]-E[d-DIM*2]) - R2 * E[d];
        }
    }
    else
    {
        for ( int d = 0; d < DIM; ++d )
            Y[d+DIM] += R2 * (X[d+DIM*2]+X[d]) - R4 * X[d+DIM];
    }
    for ( int d = 0; d < DIM; ++d )
    {
        Y[d] -= R1 * (X[d+DIM*2]+X[d]) - R2 * X[d+DIM];
        Z[d] -= R1 * (E[d-DIM*2]+E[d]) - R2 * E[d-DIM];
    }
}

/// In works in 2D, but the loop has dependencies preventing unrolling
inline void add_rigidity2D(const int nbp, const real* X, const real R1, real* Y)
{
    assert_true(R1 >= 0);
    real fx = 0;
    real fy = 0;
    real y0 = Y[0];
    real y1 = Y[1];
    const real two = -2.0;
    real const*const end = X + DIM * ( nbp - 2 );
    while ( X < end )
    {
        real gx = X[2] * two + ( X[4] + X[0] );
        real gy = X[3] * two + ( X[5] + X[1] );
        real rx = fx - gx;
        real ry = fy - gy;
        fx = gx;
        fy = gy;
        Y[0] = y0 + R1 * rx;
        Y[1] = y1 + R1 * ry;
        y0 = Y[2] - R1 * rx;
        y1 = Y[3] - R1 * ry;
        X += DIM;
        Y += DIM;
    }
    Y[0] = R1 * fx + y0;
    Y[1] = R1 * fy + y1;
    Y[2] -= R1 * fx;
    Y[3] -= R1 * fy;
}

/// In works in 3D, but the loop has dependencies preventing unrolling
inline void add_rigidity3D(const int nbp, const real* X, const real R1, real* Y)
{
    assert_true(R1 >= 0);
    real fx = 0;
    real fy = 0;
    real fz = 0;
    real y0 = Y[0];
    real y1 = Y[1];
    real y2 = Y[2];
    const real two = -2.0;
    real const*const end = X + DIM * ( nbp - 2 );
    while ( X < end )
    {
        real gx = X[3] * two + ( X[0] + X[6] );
        real gy = X[4] * two + ( X[1] + X[7] );
        real gz = X[5] * two + ( X[2] + X[8] );
        real rx = fx - gx;
        real ry = fy - gy;
        real rz = fz - gz;
        fx = gx;
        fy = gy;
        fz = gz;
        Y[0] = y0 + R1 * rx;
        Y[1] = y1 + R1 * ry;
        Y[2] = y2 + R1 * rz;
        y0 = Y[3] - R1 * rx;
        y1 = Y[4] - R1 * ry;
        y2 = Y[5] - R1 * rz;
        X += DIM;
        Y += DIM;
    }
    Y[0] = R1 * fx + y0;
    Y[1] = R1 * fy + y1;
    Y[2] = R1 * fz + y2;
    Y[3] -= R1 * fx;
    Y[4] -= R1 * fy;
    Y[5] -= R1 * fz;
}

/**
 Add bending elasticity terms between three points {A, B, C}
 Done with Serge DMITRIEFF, 2015
 */
void add_rigidityP(const int A, const int B, const int C, const real* X, const real R1, real* Y)
{
    assert_true(R1 >= 0);
#if ( DIM > 1 )
    const real R2 = 2 * R1;
    for ( int d = 0; d < DIM; ++ d )
    {
        real f = 2 * X[B*DIM+d] - ( X[A*DIM+d] + X[C*DIM+d] );
        Y[A*DIM+d] += f * R1;
        Y[B*DIM+d] -= f * R2;
        Y[C*DIM+d] += f * R1;
    }
#endif
}


/**
 This adds bending elasticity terms, as obtained by derivation of the
 Hamiltonian representing bending elasticity:

     F1 = k * ( t1 * dot(t1, t2) - t2 )
     F3 = k * ( t1 - dot(t1, t2) * t2 )
     F2 = -F1 -F3
 
 These forces are normal to the segments: dot(F1, t1) = dot(F3, t2) = 0
 The cosine are obtained here from the normalized difference vector 'dir'.

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
void add_rigidityN(const int nbp, const real* X, const real R1, real* Y, real const* dir)
{
    assert_true(R1 >= 0);
    assert_true( X != Y );
    const int end = DIM * ( nbp - 2 );
    for ( int jj = 0; jj < end; jj += DIM )
    {
        // cosine of the angle between two consecutive segments:
        const real C = dot(Vector(dir+jj), Vector(dir+jj+DIM));
        for ( int d = 0; d < DIM; ++d )
        {
            int i = jj + d;
            real f1 = R1 * (C * ( X[i+DIM] - X[i] ) - ( X[i+DIM*2] - X[i+DIM] ));
            real f3 = R1 * (( X[i+DIM] - X[i] ) - C * ( X[i+DIM*2] - X[i+DIM] ));
            Y[i      ] += f1;
            Y[i+DIM  ] -= f1+f3;
            Y[i+DIM*2] += f3;
        }
    }
}

//------------------------------------------------------------------------------

/**
 calculates the second-derivative of point's coordinates,
 scale by the bending elasticity scalar, and add to vector Y
*/
void Mecafil::addRigidity(const real* X, real* Y) const
{
#if ( DIM >= 3 )
    //add_rigidityF(nPoints, X, iRigidity, Y);
    add_rigidity3D(nPoints, X, iRigidity, Y);
#elif ( DIM == 2 )
    //add_rigidity0(nPoints, X, iRigidity, Y);
    add_rigidity2D(nPoints, X, iRigidity, Y);
#endif
    
#if NEW_FIBER_LOOP
    if ( iRigidityLoop && ( nPoints > 3 ))
    {
        /*
         With Serge DMITRIEFF:
         Link fiber end points in the same way as consecutive points triplets,
         making the fiber mechanically homogeneous for bending elasticity
         */
        const index_t L = nbPoints() - 2;
        add_rigidityP(L+1, 0, 1, X, iRigidity, Y);
        add_rigidityP(L, L+1, 0, X, iRigidity, Y);
    }
#endif
}

