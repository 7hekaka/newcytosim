// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "spherical_code.h"
#include "random.h"

SphericalCode::SphericalCode()
: num_points_(0), coord_(nullptr)
{
}


SphericalCode::SphericalCode(size_t nbp, real precision, size_t mx_nb_iterations)
: num_points_(0), coord_(nullptr)
{
    distributePoints(nbp, precision, mx_nb_iterations);
}


SphericalCode::SphericalCode(size_t nbp)
: num_points_(0), coord_(nullptr)
{
    distributePoints(nbp, 1e-4, 1<<14);
}


SphericalCode::~SphericalCode()
{
    free_real(coord_);
    coord_ = nullptr;
}

//------------------------------------------------------------------------------

void SphericalCode::putPoint(real ptr[3], const size_t i)
{
    if ( i < num_points_ )
    {
        ptr[0] = coord_[3*i+0];
        ptr[1] = coord_[3*i+1];
        ptr[2] = coord_[3*i+2];
    }
}


void SphericalCode::putPoint(double* x, double* y, double* z, const size_t i)
{
    if ( i < num_points_ )
    {
        *x = coord_[3*i+0];
        *y = coord_[3*i+1];
        *z = coord_[3*i+2];
    }
}

void SphericalCode::putPoint(float* x, float* y, float* z, const size_t i)
{
    if ( i < num_points_ )
    {
        *x = coord_[3*i+0];
        *y = coord_[3*i+1];
        *z = coord_[3*i+2];
    }
}


void SphericalCode::putPoints(real ptr[], const size_t sup)
{
    size_t end = std::min(3*num_points_, sup);
    for ( size_t i = 0; i < end; ++i )
        ptr[i] = coord_[i];
}


void SphericalCode::scale(const real factor)
{
    for ( size_t ii = 0; ii < 3*num_points_; ++ii )
        coord_[ii] *= factor;
}


void SphericalCode::printAllPositions( FILE* file )
{
    for ( size_t ii = 0; ii < num_points_; ++ii )
        fprintf( file, "%f %f %f\n", coord_[3*ii], coord_[3*ii+1], coord_[3*ii+2]);
}


real SphericalCode::distance3( const real P[], const real Q[] )
{
    return std::sqrt( (P[0]-Q[0])*(P[0]-Q[0]) + (P[1]-Q[1])*(P[1]-Q[1]) + (P[2]-Q[2])*(P[2]-Q[2]) );
}


real SphericalCode::distance3Sqr( const real P[], const real Q[] )
{
    return (P[0]-Q[0])*(P[0]-Q[0]) + (P[1]-Q[1])*(P[1]-Q[1]) + (P[2]-Q[2])*(P[2]-Q[2]);
}

//------------------------------------------------------------------------------


bool SphericalCode::project(real P[3], const real S[3])
{
    real n = S[0]*S[0] + S[1]*S[1] + S[2]*S[2];
    if ( n > 0 )
    {
        n = 1.0 / std::sqrt(n);
        P[0] = S[0] * n;
        P[1] = S[1] * n;
        P[2] = S[2] * n;
        return false;
    }
    return true;
}


/**
 hypercube rejection method
 */
void SphericalCode::randomize(real P[3])
{
    real n;
    do {
        P[0] = RNG.sreal();
        P[1] = RNG.sreal();
        P[2] = RNG.sreal();
        n = P[0]*P[0] + P[1]*P[1] + P[2]*P[2];
        if ( n == 0 )
        {
            fprintf(stderr, "The Random Number Generator may not be properly initialized");
        }
    } while ( n > 1.0 );
    
    n = std::sqrt(n);
    P[0] /= n;
    P[1] /= n;
    P[2] /= n;
}


//------------------------------------------------------------------------------
/**
 With N points on the sphere according to a triagular lattice, 
 each of ~2N triangles should occupy an area of S = 4*PI/2*N, 
 and the distance between points should be ~2 * std::sqrt(S/std::sqrt(3)).
 */
real SphericalCode::expectedDistance(size_t n)
{
    real surface = 2 * M_PI / (real)n;
    return 2 * std::sqrt( surface / std::sqrt(3) );
}


real SphericalCode::minimumDistance()
{
    real res = INFINITY;
    for ( size_t ii = 1; ii < num_points_; ++ii )
    {
        for ( size_t jj = 0; jj < ii; ++jj )
        {
            real dis = distance3Sqr(&coord_[3*ii], &coord_[3*jj]);
            if ( dis < res )
                res = dis;
        }
    }
    return std::sqrt(res);
}


real SphericalCode::coulombEnergy( const real P[] )
{
    real dist, result = 0;
    for ( size_t ii = 1; ii < num_points_; ++ii )
    {
        for ( size_t jj = 0; jj < ii; ++jj )
        {
            dist = distance3( P + 3 * ii, P + 3 * jj );
            if ( dist > 0 ) result += 1.0 / dist;
        }
    }
    return result;
}

//------------------------------------------------------------------------------

void SphericalCode::setForces( real forces[], real threshold )
{
    real dx[3];
    real dist;
    
    //--------- reset forces:
    for ( size_t ii = 0; ii < 3 * num_points_; ++ii )
        forces[ii] = 0.0;
    
    //--------- calculate Coulomb pair interactions:
    // first particle is ii, second one is jj:
    for ( size_t ii = 1; ii < num_points_; ++ii )
    {
        for ( size_t jj = 0; jj < ii; ++jj )
        {
            //calculate vector and distance^2 between from jj to ii
            dist = 0;
            for ( size_t d = 0; d < 3 ; ++d )
            {
                dx[d] = coord_[3*ii+d] - coord_[3*jj+d];
                dist += dx[d] * dx[d];
            }
            
            if ( dist == 0 )
            {   //if ii and jj overlap, we use a random force
                for ( size_t d = 0 ; d < 3; ++d )
                {
                    dx[d] = 0.1 * RNG.sreal();
                    forces[3*ii+d] += dx[d];
                    forces[3*jj+d] -= dx[d];
                }
            }
            else if ( dist < threshold )
            {
                // points do not overlap:
                //force = vector / r^3, but here dist = r^2
                dist = std::sqrt(dist) / ( dist * dist );
                //update forces for jj and ii:
                for ( size_t d = 0 ; d < 3; ++d )
                {
                    dx[d] *= dist;
                    forces[3*ii+d] += dx[d];
                    forces[3*jj+d] -= dx[d];
                }
            }
        }
    }
    

#if ( 1 )
    /*
     Remove centripetal contribution of forces:
     assuming here that points are already on the sphere (norm=1)
     ( the algorithm converge even without this, but slower )
     */
    for ( size_t ii = 0; ii < num_points_; ++ii )
    {
        dist = 0;
        for ( size_t d = 0; d < 3; ++d )
            dist += coord_[3*ii+d] * forces[3*ii+d];
        
        for ( size_t d = 0; d < 3; ++d )
            forces[3*ii+d] -= dist * coord_[3*ii+d];
    }
#endif
}


/**
 Move the points in the direction of the forces, with scaling factor S
 
 @todo: use Newton's method to optimize the energy, rather than a fixed step
 
     x = x - F'(x) / F''(x)
 
 where F(x) = energy of configuration
 */
void SphericalCode::movePoints(real Pnew[], const real Pold[], real forces[], real S)
{
    for ( size_t ii = 0; ii < num_points_; ++ii )
    {
        real W[3];
        
        W[0] = Pold[3*ii+0] + S * forces[3*ii+0];
        W[1] = Pold[3*ii+1] + S * forces[3*ii+1];
        W[2] = Pold[3*ii+2] + S * forces[3*ii+2];
        
        if ( project(Pnew+3*ii, W) )
            randomize(Pnew+3*ii);
    }
}


void SphericalCode::allocate(size_t nbp)
{
    //reallocate the array if needed:
    if ( num_points_ != nbp || ! coord_ )
    {
        free_real(coord_);
        coord_ = new_real(3*nbp);
        num_points_ = nbp;
    }
}


void SphericalCode::initializePoints()
{
#if 0
    // distribute the points randomly on the sphere:
    for ( size_t i = 0; i < num_points_; ++i )
        randomize(coord_+i*3);
#else
    // use Fibonacci's spiral on the sphere:
    real gold = 4 * M_PI / ( 1.0 + sqrt(5) );
    real beta = 2.0 / num_points_;
    for ( size_t i = 0; i < num_points_; ++i )
    {
        real a = i * gold;
        real p = acos(1-(i+0.5)*beta);
        coord_[  i*3] = cos(a)*sin(p);
        coord_[1+i*3] = sin(a)*sin(p);
        coord_[2+i*3] = cos(p);
    }
#endif
}


/**
 create a relatively even distribution of nbp points on the sphere
 the coordinates are stored in real array coord_[]
 */
size_t SphericalCode::distributePoints(size_t nbp, real precision, size_t mx_iterations)
{
    allocate(nbp);
    initializePoints();

    //--------- for one point only, we return:
    if ( nbp > 1 )
        return refinePoints(precision, mx_iterations);
    return 0;
}


/**
 create a relatively even distribution of nbp points on the sphere
 the coordinates are stored in real array coord_[]
 */
size_t SphericalCode::refinePoints(real precision, size_t mx_iterations)
{
    // the precision is rescaled with the expected distance:
    real len = expectedDistance(num_points_);

    /* 
     Threshold cut-off for repulsive force:
     The best results are obtained for threshold > 2
     */
    real threshold = 10 * len;
    real mag = 0.1 * len * len * len * len / (real)num_points_;
    precision *= mag;

    //------------ calculate the initial energy:
    
    // allocate forces and new coordinates:
    real * coord = new_real(3*num_points_);
    real * force = new_real(3*num_points_);
    
    //make an initial guess for the step size:
    unsigned history = 0;
    energy_ = coulombEnergy(coord_);

    size_t step = 0;
    for ( step = 0; step < mx_iterations; ++step )
    {
        setForces(force, threshold);
        
        while ( 1 ) 
        {
            movePoints(coord, coord_, force, mag);
            
            // energy of new configuration:
            real energy = coulombEnergy(coord);
            
            //printf("%3lu   %5lu : energy %18.10f   %18.10f\n", num_points_, step, energy, energy-energy_);

            if ( energy < energy_ )
            {
                // swap pointers to accept configuration:
                real* m = coord_;
                coord_  = coord;
                coord   = m;
                energy_ = energy;
                
                /*
                 After 'SEVEN' successful moves at a given step size, we increase
                 the step size. Values for 'magic_seven' were tested in term of
                 convergence, and 7 seems to work well.
                 */
                if ( ++history >= SEVEN )
                {
                    mag *= 1.4147;   //this value is somewhat arbitrary
                    history = 0;
                }
                break;
            }
            else
            {
                /*
                 If the new configuration has higher energy,
                 we try a smaller step size with the same forces:
                 */
                history = 0;
                mag /= 2;
                
                //exit when the desired precision is reached
                if ( mag < precision )
                    goto done;
            }
        }
    }
done:
    free_real(coord);
    free_real(force);
    return step;
}

