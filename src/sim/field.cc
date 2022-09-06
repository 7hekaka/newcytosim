// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University
// Created by Francois Nedelec on 18/12/07.

#include "field.h"
#include "simul_part.h"
#include "fiber_site.h"
#include "fiber_set.h"
#include "blas.h"
#include "cymdef.h"


/**
 Initialize the diffusion matrix using periodic boundary conditions
 if the underlying space is peridic
 */
void Field::prepareDiffusion(real theta)
{
    const size_t nbc = mGrid.nbCells();
    fiDiffusionMatrix.resize(nbc);
    fiDiffusionMatrix.reset();
    
    for ( size_t c = 0; c < nbc; ++c )
    {
        for ( int d = 0; d < DIM; ++d )
        {
            size_t n = mGrid.next(c, d);
            
            if ( n != c )
            {
                fiDiffusionMatrix(c, n) += theta;
                fiDiffusionMatrix(c, c) -= theta;
                fiDiffusionMatrix(n, n) -= theta;
            }
        }
    }
    fiDiffusionMatrix.prepareForMultiply(1);
    
    /*
    std::clog << "tight Field has diffusion matrix with ";
    std::clog << fiDiffusionMatrix.nbElements() << " elements" << '\n';
     */
}


/**
 Initialize the diffusion matrix.
 Diffusion is allowed between neighboring cells that are in the same domain:

     ( domain[c] > 0 ) && ( domain[c] == domain[n] )

 */
void Field::prepareDiffusion(real theta, unsigned char * domain)
{
    const size_t nbc = mGrid.nbCells();
    fiDiffusionMatrix.resize(nbc);
    fiDiffusionMatrix.reset();
    
    for ( size_t c = 0; c < nbc; ++c )
    {
        if ( domain[c] )
        {
            for ( size_t d = 0; d < DIM; ++d )
            {
                size_t n = c + mGrid.stride(d);
                
                if ( n < nbc  &&  domain[c] == domain[n] )
                {
                    fiDiffusionMatrix(c, n) += theta;
                    fiDiffusionMatrix(c, c) -= theta;
                    fiDiffusionMatrix(n, n) -= theta;
                }
            }
        }
    }
    fiDiffusionMatrix.prepareForMultiply(1);
    
    /*
     std::clog << "Field has diffusion matrix with ";
     std::clog << fiDiffusionMatrix.nbElements() << " elements" << '\n';
     */
}


/**
 Initialize Field to be ready for step()
 */
void Field::prepare()
{
    Space const* spc = prop->field_space_ptr;

    if ( !spc )
        throw InvalidParameter("A Space must be created before the field");

    const size_t nbc = mGrid.nbCells();
    assert_true( nbc > 0 );
    
    free_real(fiTMP);
    fiTMP = new_real(nbc);
    fiTMPSize = nbc;

    if ( prop->slow_diffusion > 0 )
    {
        const real tau = time_step(simul());
        real theta = prop->slow_diffusion * tau / ( prop->step * prop->step );

        if ( DIM == 1 || prop->field_periodic )
            prepareDiffusion(theta);
        else
        {
            unsigned char * domain = new unsigned char[nbc];
            
            // determine which cell is inside the space:
#if ( 1 )
            for ( size_t c = 0; c < nbc; ++c )
            {
                Vector pos;
                mGrid.setPositionFromIndex(pos, c, 0.5);
                domain[c] = spc->inside(pos);
            }
#else
            // extended covered area:
            const real range = 2 * cellWidth();
            for ( size_t c = 0; c < nbc; ++c )
            {
                Vector pos;
                mGrid.setPositionFromIndex(pos, c, 0.5);
                domain[c] = ! spc->allOutside(pos, range);
            }
#endif
            
            prepareDiffusion(theta, domain);
            
            delete[] domain;
        }
    }
}


//------------------------------------------------------------------------------
#pragma mark - Diffusion


void Field::diffuseX(real * field, real c)
{
    const auto nbc = mGrid.nbCells();

    size_t nx = mGrid.breadth(0);
    size_t nyz = nbc / nx;
    size_t ide = mGrid.stride(1);
    
    real * a = new_real(nyz);
    real * b = new_real(nyz);
    
    // diffusion in X-direction:
    zero_real(nyz, a);
    for ( size_t x = 1; x < nx; ++x )
    {
        real * h = field + x - 1;
        real * n = field + x;
        // b = n - h
        blas::xcopy(nyz,  n, ide, b, 1);
        blas::xaxpy(nyz, -1, h, ide, b, 1);
        // a = a - b
        blas::xaxpy(nyz, -1, b, 1, a, 1);
        // h = h - c * a
        blas::xaxpy(nyz, -c, a, 1, h, ide);
        // swap a and b
        real * t = a;
        a = b;
        b = t;
    }
    real * h = field + nx - 1;
    blas::xaxpy(nyz, -c, a, 1, h, ide);
    
    if ( prop->field_periodic )
    {
        real * n = field;
        blas::xcopy(nyz,  n, ide, b, 1);
        blas::xaxpy(nyz, -1, h, ide, b, 1);
        blas::xaxpy(nyz,  c, b, 1, h, ide);
        blas::xaxpy(nyz, -c, b, 1, n, ide);
    }
    
    free_real(a);
    free_real(b);
}


void Field::laplacian(const real* field, real * mat) const
{
    const size_t nbc = mGrid.nbCells();
    const real sc = 2 * DIM;

    for ( size_t c = 0; c < nbc; ++c )
        mat[c] = sc * field[c];
    
    const size_t nx = mGrid.breadth(0);
#if ( 1 )
    // derivative in the X-direction:
    const size_t nyz = nbc / nx;
    for ( size_t xx = 1; xx < nx; ++xx )
    {
        blas::xaxpy(nyz, -1, field+xx-1, nx, mat+xx  , nx);
        blas::xaxpy(nyz, -1, field+xx  , nx, mat+xx-1, nx);
    }
    // index of last valid X index:
    int xx = mGrid.breadth(0) - 1;
    
    if ( prop->field_periodic )
    {
        blas::xaxpy(nyz, -1, field+xx, nx, mat   , nx);
        blas::xaxpy(nyz, -1, field   , nx, mat+xx, nx);
    }
    else
    {
        blas::xaxpy(nyz, -1, field   , nx, mat   , nx);
        blas::xaxpy(nyz, -1, field+xx, nx, mat+xx, nx);
    }
#endif
    
#if ( DIM == 2 )
    // derivative in the Y-direction:
    blas::xaxpy(nbc-nx, -1, field,    1, mat+nx, 1);
    blas::xaxpy(nbc-nx, -1, field+nx, 1, mat,    1);
    
    size_t yy = mGrid.breadth(1) - 1;
    if ( prop->field_periodic )
    {
        blas::xaxpy(nx, -1, field+nx*yy, 1, mat      , 1);
        blas::xaxpy(nx, -1, field      , 1, mat+nx*yy, 1);
    }
    else
    {
        blas::xaxpy(nx, -1, field      , 1, mat      , 1);
        blas::xaxpy(nx, -1, field+nx*yy, 1, mat+nx*yy, 1);
    }
#endif

#if ( DIM >= 3 )
    // derivative in the Y-direction:
    const size_t ss = mGrid.stride(2);
    for ( size_t yy = 1; yy < mGrid.breadth(1); ++yy )
    for ( size_t zz = 0; zz < mGrid.breadth(2); ++zz )
    {
        blas::xaxpy(nx, -1, field+nx*(yy-1)+ss*zz, 1, mat+nx*(yy  )+ss*zz, 1);
        blas::xaxpy(nx, -1, field+nx*(yy  )+ss*zz, 1, mat+nx*(yy-1)+ss*zz, 1);
    }
    size_t yy = mGrid.breadth(1) - 1;
    
    if ( prop->field_periodic )
    {
        for ( size_t zz = 0; zz < mGrid.breadth(2); ++zz )
        {
            blas::xaxpy(nx, -1, field+nx*yy+ss*zz, 1, mat      +ss*zz, 1);
            blas::xaxpy(nx, -1, field      +ss*zz, 1, mat+nx*yy+ss*zz, 1);
        }
    }
    else
    {
        for ( size_t zz = 0; zz < mGrid.breadth(2); ++zz )
        {
            blas::xaxpy(nx, -1, field      +ss*zz, 1, mat      +ss*zz, 1);
            blas::xaxpy(nx, -1, field+nx*yy+ss*zz, 1, mat+nx*yy+ss*zz, 1);
        }
    }
#endif

#if ( DIM >= 3 )
    // derivative in the Z-direction:
    const size_t nxy = nbc / mGrid.breadth(2);
    assert_true( nxy == ss );
    blas::xaxpy(nbc-nxy, -1, field,     1, mat+nxy, 1);
    blas::xaxpy(nbc-nxy, -1, field+nxy, 1, mat,     1);
    size_t zz = mGrid.breadth(2) - 1;
    
    if ( prop->field_periodic )
    {
        blas::xaxpy(nxy, -1, field+ss*zz, 1, mat      , 1);
        blas::xaxpy(nxy, -1, field      , 1, mat+ss*zz, 1);
    }
    else
    {
        blas::xaxpy(nxy, -1, field      , 1, mat      , 1);
        blas::xaxpy(nxy, -1, field+ss*zz, 1, mat+ss*zz, 1);
    }
#endif
}


void Field::setEdgesX(real * field, real val)
{
    const size_t nbc = mGrid.nbCells();
    
    // set X-edges:
    const size_t nx = mGrid.breadth(0);
    
    real * lastf = field + nx - 1;
    for ( size_t xx = 0; xx < nbc; xx += nx )
    {
        field[xx] = val;
        lastf[xx] = val;
    }
}


void Field::setEdgesY(real * field, real val)
{
#if ( DIM > 1 )
    const size_t nbc = mGrid.nbCells();
    const size_t nx = mGrid.breadth(0);
#endif
    
#if ( DIM == 2 )
    // set Y-edges:
    real * lastf = field + nbc - nx;
    for ( size_t xx = 0; xx < nx; ++xx )
    {
        field[xx] = val;
        lastf[xx] = val;
    }
#endif
    
#if ( DIM >= 3 )
    // set Y-edges:
    const size_t nz = mGrid.breadth(2);
    const size_t nxy = nbc / nz;
    
    real * lastf = field + nxy - nx;
    for ( size_t zz = 0; zz < nz; ++zz )
    {
        for ( size_t xx = 0; xx < nx; ++xx )
        {
            field[xx+zz*nxy] = val;
            lastf[xx+zz*nxy] = val;
        }
    }
#endif
}


void Field::setEdgesZ(real * field, real val)
{
#if ( DIM >= 3 )
    const size_t nbc = mGrid.nbCells();
    
    const size_t nz = mGrid.breadth(2);
    const size_t nxy = nbc / nz;
    
    real * lastf = field + nxy * ( nz - 1 );
    for ( size_t xy = 0; xy < nxy; ++xy )
    {
        field[xy] = val;
        lastf[xy] = val;
    }
#endif
}


/**
 //\todo implement Crank-Nicholson for diffusion
 */
void Field::step(FiberSet& fibers)
{
    assert_true( prop );
    const real tau = time_step(simul());

    // we cast FieldScalar to floating-point type :
    static_assert(sizeof(FieldScalar) == sizeof(real), "unexpected FieldScalar type");
    real * field = reinterpret_cast<real*>(mGrid.data());
    const auto nbc = mGrid.nbCells();
    
    real * dup = fiTMP;

    // decay:
    if ( prop->decay_rate > 0 )
    {
        // calculate decay coefficient during interval `time-step`
        real frac = std::exp( -tau * prop->decay_rate );
        // field = field * frac:
        blas::xscal(nbc, frac, field, 1);
    }

    // full grid diffusion:
    if ( prop->full_diffusion > 0 )
    {
        real c = prop->full_diffusion * tau / ( prop->step * prop->step );

#if ( DIM > 1 )
        laplacian(field, dup);
        blas::xaxpy(nbc, -c, dup, 1, field, 1);
#else
        diffuseX(field, c);
#endif
    }

    // diffusion:
    if ( prop->slow_diffusion > 0 )
    {
        assert_true( fiTMP );
        assert_true( fiTMPSize == nbc );
        assert_true( fiDiffusionMatrix.size() == nbc );

        // dup = field:
        blas::xcopy(nbc, field, 1, dup, 1);
        
        // field = field + fiDiffusionMatrix * dup:
        fiDiffusionMatrix.vecMulAdd(dup, field);
    }

    if ( prop->boundary_condition & 1 )
        setEdgesX(field, prop->boundary_value * cellVolume());
    
#if ( DIM > 1 )
    if ( prop->boundary_condition & 2 )
        setEdgesY(field, prop->boundary_value * cellVolume());
#endif
    
#if ( DIM >= 3 )
    if ( prop->boundary_condition & 4 )
        setEdgesZ(field, prop->boundary_value * cellVolume());
#endif
    
    
#if ( 0 ) // disabled features below

    Array<FiberSite> loc(1024);
    
    // instantaneous transport along Fibers
    if ( prop->transport_strength > 0 )
    {
        const real spread = 0.5 * cellWidth();
        const real rate = prop->transport_strength * spread / cellVolume();
        const real frac = -std::expm1( -rate * tau );
        
        if ( frac >= 0.5 )
            throw InvalidParameter("field:transport_strength is too high");
        
        fibers.uniFiberSites(loc, spread);
        for ( FiberSite & i : loc )
        {
            // abscissa for exit point of transport:
            real abs = i.abscissa() + RNG.exponential(prop->transport_length);

            // find index of cell:
            FieldCell cell = mGrid.cell(i.pos());
            
            // amount to be transferred:
            real mass = cell * frac;
            
            // transport:
            cell -= mass;
            field[mGrid.index(i.fiber()->pos(abs))] += mass;
        }
    }
    
    // direct cutting of fiber
    // this is deprecated in favor of fiber:mesh_cut_fiber
    if ( prop->cut_fibers )
    {
        LOG_ONCE("!!!! Field severs fibers\n");
        const real spread = 0.5 / tau;
        const real fac = spread * tau / cellVolume();
        
        fibers.uniFiberSites(loc, spread);
        for ( FiberSite & i : loc )
        {
            real val = field[mGrid.index(i.pos())];
            if ( prop->cut_fibers == 2 )
                val = val * val / cellVolume();
            if ( RNG.test(-std::expm1(-fac*val)) )
                i.fiber()->sever(i.abscissa(), STATE_RED, STATE_GREEN);
        }
    }
    
    if ( prop->chew_fibers )
    {
        LOG_ONCE("!!!! Field chews PLUS_END\n");
        const real fac = -tau / cellVolume();
        for ( Fiber * fib = fibers.first(); fib ; fib = fib->next() )
            fib->growP(fac*cell(fib->posEndP()));
    }
#endif
}

//------------------------------------------------------------------------------
#pragma mark - Display

#ifdef DISPLAY

#include "gle.h"
#include "grid_display.h"
#include "gym_view.h"
#include "gym_cap.h"

class FieldDisplayParameters
{
public:
    FieldDisplayParameters()
    {
        amp = 0;
        spc = nullptr;
    }
    
    /// amplification for color
    real amp;
    
    /// Space for cropping
    Space const* spc;
};


// set color for scalar field
static gym_color field_color(FieldDisplayParameters fdp, FieldScalar const& cell, Vector const& pos)
{
    if ( fdp.spc && ! fdp.spc->inside(pos) )
        return gym_color(0, 0, 0);
    gym_color::COLOF x(fdp.amp * cell.val);
    if ( x > 0 )
        return gym_color::jet_color_dark(x, 1.0);
    return gym_color(-x, 0, -x);
}


/// set color for Vector field
template < size_t N >
static gym_color field_color(FieldDisplayParameters fdp, FieldVector<N> const& cell, Vector const& pos)
{
    if ( fdp.spc && ! fdp.spc->inside(pos) )
        return gym_color(0, 0, 0);
    //this maps val[0] to the red channel, val[1] to green and val[2] to blue
    gym_color::COLOF rgb[3] = { 0, 0, 0 };
    const int sup = std::min(3UL, N);
    for ( int c = 0; c < sup; ++c )
        rgb[c] = fdp.amp * cell[c];
    return gym_color(rgb[0], rgb[1], rgb[2]);
}


/// display all cells
void Field::draw() const
{
    FieldDisplayParameters fdp;
    fdp.amp = 1.0 / ( prop->display_scale * mGrid.cellVolume() );
    fdp.spc = nullptr;
    
    gym::disableDepthTest();
    gym::disableCullFace();
    gym::disableLighting();
    gym::ref_view();
#if ( DIM <= 3 )
    drawValues(mGrid, field_color, fdp);
#endif
#if ( 0 )
    float col[] = {1,0,0,1};
    drawBoundaries(mGrid, col, 0.5f);
#endif
    gym::restoreDepthTest();
    gym::restoreCullFace();
    gym::restoreLighting();
}


/// display only cells that are inside `spc`
void Field::draw(Space const* spc, Vector3 const& dir, const real pos) const
{
    FieldDisplayParameters fdp;
    fdp.amp = 1.0 / ( prop->display_scale * mGrid.cellVolume() );
    fdp.spc = spc;
    
    gym::disableDepthTest();
    gym::disableCullFace();
    gym::disableLighting();
    gym::ref_view();
#if ( DIM == 3 )
    drawValues(mGrid, field_color, fdp, dir, pos);
#elif ( DIM <= 2 )
    drawValues(mGrid, field_color, fdp);
#endif
    gym::restoreDepthTest();
    gym::restoreCullFace();
    gym::restoreLighting();
}

#else

void Field::draw() const
{
    LOG_ONCE("no field:draw()\n");
}

void Field::draw(Space const*, Vector3 const& dir, const real pos) const
{
    LOG_ONCE("no field:draw()\n");
}

#endif

