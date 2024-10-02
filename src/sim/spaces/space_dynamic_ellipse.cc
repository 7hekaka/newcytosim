// Cytosim 3.0 - F. Nedelec and Laboratory, Copyright 2020 Cambridge University.

#include "dim.h"
#include "space_dynamic_ellipse.h"
#include "space_dynamic_prop.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "glossary.h"


/// prefactor for volume computation: Pi in 2D and 4/3 Pi in 3D
constexpr real VPREF = (DIM+1)*M_PI/3.0;

/// power for ellipsoid surface calculation
constexpr real POW = 1.6075;

/// building block for area of an ellipsoid in 3D
inline real surf_block(const real a, const real b)
{
    return std::pow(a*b,POW);
}

/// building block for area of an ellipsoid in 3D
inline real surf_block(const real a, const real b, const real c)
{
    return std::pow(a*b,POW) + std::pow(b*c,POW) + std::pow(a*c,POW);
}


SpaceDynamicEllipse::SpaceDynamicEllipse(SpaceDynamicProp const* p)
: SpaceEllipse(p)
{
    if ( DIM == 1 )
        throw InvalidParameter("dynamic_ellipse is not usable in 1D");
    
    pressure = 0;
    mat = MatrixD::one();
    inv = MatrixD::one();
    
    reset_forces();
    rad_forces.set(0,0,0);
}


void SpaceDynamicEllipse::dump(std::ostream& out) const
{
    char str[1024];
    
    real S = surfaceEllipse(radius_);
    snprintf(str, sizeof(str), " dEllipse %7.3f %7.3f %7.3f  S %7.3f  V %7.3f  P %7.3f",
             radius_[0], radius_[1], radius_[2], S, volumeEllipse(radius_), pressure);
    out << str;
    
#if ( DIM > 1 )
    snprintf(str, sizeof(str), " rF %7.2f %7.2f", rad_forces.XX, rad_forces.YY);
    out << str;
    
    Vector tf = tension_forces();
    Vector pf = pressure_forces(pressure);
    snprintf(str, sizeof(str), " tF %7.2f %7.2f pF %7.2f %7.2f\n", tf.XX, tf.YY, pf.XX, pf.YY);
    out << str;
#endif
}

//-------------------------------------------------------------------------------------
//        Set interactions and update forces felt by ellipse.
//-------------------------------------------------------------------------------------

void SpaceDynamicEllipse::setInteractions(Meca&, Simul const&) const
{
    reset_forces();
}


void SpaceDynamicEllipse::setConfinement(Vector const&pos, Mecapoint const& mp,
                                         Meca& meca, real stiff) const
{
    Vector prj;
    prj = project(pos);
    Vector dir = pos - prj;
    real n = dir.normSqr();
    if ( n > 0 )
    {
        // Register the force to the ellipse
        decompose_force(stiff * dir, prj, dir);
        // And to the meca
        meca.addPlaneClamp(mp, prj, dir, stiff/n);
    }
}



//------------------------------------------------------------------------------
//        Computing forces
//------------------------------------------------------------------------------

void SpaceDynamicEllipse::reset_forces() const
{
    Torques = nullTorque;
    Rforces.set(0,0,0);
}

/**
 register forces applied to the space
 */
void SpaceDynamicEllipse::decompose_force(Vector const& forces, Vector const& pos,Vector const& dir) const
{
#if ( 0 )
    // keep only force component in the normal direction:
    Vector nfo = dir * ((forces*dir) / dir.normSqr());
    add_radial_force(nfo, proj);
#else
    add_radial_force(forces, pos);
#endif
    Torques += cross(pos, forces);
}

/**
 Add a point-like force acting on the ellipse
 */
void SpaceDynamicEllipse::add_radial_force(Vector const& forces, Vector const& pos) const
{
    Vector U = director(0);
    Rforces.XX += dot(U, forces) * dot(U, pos) / radius_[0];
#if ( DIM >= 2 )
    Vector V = director(1);
    Rforces.YY += dot(V, forces) * dot(V, pos) / radius_[1];
#endif
#if ( DIM > 2 )
    Vector W = director(2);
    Rforces.ZZ += dot(W, forces) * dot(W, pos) / radius_[2];
#endif
}


// ----------------------------------------------
//  Internal forces
// ----------------------------------------------

/**
 The derivative of pressure energy with respect to each ellipse parameter:
 EP = Volume * Pressure
 dEP/da = dV/da * Pressure
*/
Vector SpaceDynamicEllipse::pressure_forces(const real P) const
{
    const real S = VPREF * P;
#if ( DIM == 1 )
    return Vector(0, 0, 0);
#elif ( DIM == 2 )
    return Vector(S*radius_[1], S*radius_[0]);
#elif ( DIM > 2 )
    return Vector(S*radius_[1]*radius_[2], S*radius_[2]*radius_[0], S*radius_[0]*radius_[1]);
#endif
}

/*
 Pressure is a Lagrange multiplier associated with volume conservation
 We follow Newtons's method to minimize:
 
     F = Volume(next_time_step) - prop()->volume
 
 Hence we iterate:
 
     P = P - F / dF/dP
 
 until the machine precision is exhausted
*/
real SpaceDynamicEllipse::compute_pressure(Vector const& sizes, Vector const& radif) const
{
    real P = pressure;
    real err = INFINITY, last_err;
    
    if ( prop()->mobility_dt <= 0 )
        return 0;
    
    size_t cnt = 0;
    do {
        last_err = err;
        
        // the objective is to reach desired volume at the next time-step:
        Vector dim = sizes + ( radif + pressure_forces(P) ) * prop()->mobility_dt;
        
        real der = VPREF * VPREF * prop()->mobility_dt;
        
        real r0 = dim.XX;
#if ( DIM == 2 )
        // r0 = A0 + P * B0;   B0 = VPREF * r1 * mobility_dt
        // r1 = A1 + P * B1;   B1 = VPREF * r0 * mobility_dt
        // vol = VPREF * ( A0 + P * B0 ) * ( A1 + P * B1 );
        // der = VPREF * ( B0 * r1 + r0 * B1 )
        real r1 = dim.YY;
        err = VPREF * r0 * r1 - prop()->volume;
        der *= square(r0) + square(r1);
#elif ( DIM > 2 )
        // r0 = A0 + P * B0;   B0 = VPREF * r1 * r2 * mobility_dt
        // r1 = A1 + P * B1;   B1 = VPREF * r2 * r0 * mobility_dt
        // r2 = A2 + P * B2;   B2 = VPREF * r0 * r1 * mobility_dt
        // vol = VPREF * ( A0 + P * B0 ) * ( A1 + P * B1 ) * ( A2 + P * B2 );
        // der = VPREF * ( B0 * r1 * r2 + r0 * B1 * r2 + r0 * r1 * B2 )
        real r1 = dim.YY;
        real r2 = dim.ZZ;
        err = VPREF * r0 * r1 * r2 - prop()->volume;
        der *= square(r0*r1) + square(r0*r2) + square(r1*r2);
#endif

        P -= err / der;
        
        if ( ++cnt > 256 )
        {
            std::clog << "pressure calculation failed at " << err << '\n';
            return 0;
        }

    } while ( abs_real(err) < abs_real(last_err) );
    //std::clog << "volume error " <<  err << '\n';

    return P;
}


/**
 The derivative of surface energy with respect to each ellipse parameter:
 ES = Surface * tension
 dES/da = dS/da * tension
 */
Vector SpaceDynamicEllipse::tension_forces() const
{
    Vector res(0,0,0);

#if ( DIM == 2 )

    real S = -prop()->tension * M_PI;
    real N = std::sqrt( (3.0*radius_[0]+radius_[1])*(radius_[0]+3.0*radius_[1]) );

    res.XX = S * (3.0 - ( 3.0*radius_[0] + 5.0*radius_[1] ) / N );
    res.YY = S * (3.0 - ( 3.0*radius_[1] + 5.0*radius_[0] ) / N );
    
#elif ( DIM > 2 )
    
    real S = -prop()->tension * surfaceEllipse(radius_);
    
    real pXY = surf_block(radius_[0], radius_[1]);
    real pXZ = surf_block(radius_[0], radius_[2]);
    real pYZ = surf_block(radius_[1], radius_[2]);
    real XYZ = surf_block(radius_[0], radius_[1], radius_[2]);

    res.XX = S * ( pXY + pXZ ) / ( radius_[0] * XYZ );
    res.YY = S * ( pXY + pYZ ) / ( radius_[1] * XYZ );
    res.ZZ = S * ( pXZ + pYZ ) / ( radius_[2] * XYZ );

#endif
    return res;
}

//-------------------------------------------------------------------------------------
///    Update ellipse shape
//-------------------------------------------------------------------------------------
    

void SpaceDynamicEllipse::step()
{
    if ( prop()->volume > 0 )
    {
        rad_forces = Rforces;
        
        // calculate forces:
        Rforces += tension_forces();
        pressure = compute_pressure(radius_, Rforces);
        Rforces += pressure_forces(pressure);

        // implement changes in shape:
        if ( prop()->mobility_dt > 0 )
        {
            Vector delta = prop()->mobility_dt * Rforces;
            for ( int i=0; i<DIM; ++i )
            {
                assert_true(delta[i] == delta[i]);
                radius_[i] += delta[i];
            }
            //dump(std::clog);
            //std::clog << "%  balance " << Rforces << "\n";
        }
        
        // implement rotation:
        if ( prop()->mobility_rot_dt > 0 )
        {
            //std::clog << "% DynamicEllipse torque " << Torques << "\n";
#if ( DIM == 2 )
            real theta = prop()->mobility_rot_dt * Torques;
            if ( theta > REAL_EPSILON )
            {
                real c = std::cos(theta);
                real s = std::sin(theta);
                mat = Matrix22(c, s, -s, c) * mat;
            }
#elif ( DIM > 2 )
            real n = Torques.norm();
            real theta = prop()->mobility_rot_dt * n;
            if ( theta > REAL_EPSILON )
            {
                MatrixD rot = MatrixD::rotationAroundAxis(Torques/n, std::cos(theta), std::sin(theta));
                mat = rot * mat;
            }
#endif
            // update rotation:
            inv = mat.transposed();
        }
    }
    
    reset_forces();
}




/// Checking consistency of ellipse sizes
void SpaceDynamicEllipse::resize(Glossary& opt)
{
    SpaceEllipse::resize(opt);
    if ( prop()->volume <= 0 )
    {
        prop()->volume = volume();
        //std::clog << "DynamicEllipse:volume <- " << prop()->volume << '\n';
    }
}


//-------------------------------------------------------------
// Utilities
//-------------------------------------------------------------


Vector SpaceDynamicEllipse::director(index_t ix) const
{
    assert_true(ix < DIM);
    return mat.column(ix);
}


real SpaceDynamicEllipse::surfaceEllipse(Vector const& sizes)
{
#if ( DIM > 2 )
    real a = sizes.XX;
    real b = sizes.YY;
    real c = sizes.ZZ;
    return (4.0*M_PI)*std::pow(surf_block(a,b,c)/3.0, 1.0/POW);
#elif ( DIM == 2 )
    // In 2D, the 'surface' is a line
    real a = sizes.XX;
    real b = sizes.YY;
    real r = square( (a-b)/(a+b) );
    return M_PI*(a+b)*(1.0+3.0*r/(10.0+std::sqrt(4.0-3.0*r)));
#else
    return 0;
#endif
}


real SpaceDynamicEllipse::volumeEllipse(Vector const& sizes)
{
#if ( DIM > 2 )
    return VPREF*sizes[0]*sizes[1]*sizes[2];
#elif ( DIM == 2 )
    return VPREF*sizes[0]*sizes[1];
#endif
    return 0;
}


void SpaceDynamicEllipse::read(Inputter& in, Simul& sim, ObjectTag tag)
{
    SpaceEllipse::read(in, sim, tag);
    unsigned n = in.readUInt16();
    if ( n != 10 )
        throw InvalidIO("unexpected data in SpaceDynamicEllipse");
    // desired volume
    prop()->volume = in.readFloat();
#if ( DIM > 2 )
    // read 3x3 orientation matrix:
    for ( unsigned j = 0; j < 3; ++j )
    for ( unsigned i = 0; i < 3; ++i )
        mat(i,j) = in.readFloat();
#else
    // read 3x3 orientation matrix:
    real m[9] = { 0 };
    for ( unsigned i = 0; i < 9; ++i )
        m[i] = in.readFloat();
#  if ( DIM == 2 )
    mat(0,0) = m[0];
    mat(1,0) = m[1];
    mat(0,1) = m[3];
    mat(1,1) = m[4];
#  endif
#endif
}


void SpaceDynamicEllipse::write(Outputter& out) const
{
    SpaceEllipse::write(out);
    out.writeUInt16(10);
    out.writeFloat(prop()->volume);
#if ( DIM > 2 )
    // write matrix in column-major format
    for ( unsigned j = 0; j < 3; ++j )
    for ( unsigned i = 0; i < 3; ++i )
        out.writeFloat(mat(i,j));
#else
    real m[9] = { 0 };
#if ( DIM == 2 )
    m[0] = mat(0,0);
    m[1] = mat(1,0);
    m[3] = mat(0,1);
    m[4] = mat(1,1);
#endif
    for ( unsigned i = 0; i < 9; ++i )
        out.writeFloat(m[i]);
#endif
}


/**
 Reports the radii of the ellipsoid
 */
void SpaceDynamicEllipse::report(std::ostream& os) const
{
    os << "  radius_0" << radius_[0];
    os << "  radius_1" << radius_[1];
#if ( DIM == 3 )
    os << "  radius_2" << radius_[2];
#endif
}

//------------------------------------------------------------------------------
#pragma mark - OpenGL display

#ifdef DISPLAY

#include "gle.h"
#include "gym_view.h"

void SpaceDynamicEllipse::draw3D() const
{
#if ( 0 )
    const float rad = 0.01;
    for ( unsigned n = 0; n < DIM; ++n )
    {
        gym::transScale(radius_[n]*director(n), rad);
        gle::sphere2();
    }
#endif

    float MM[16] = { 0 };
    MM[ 5]=1.0f;
    MM[10]=1.0f;
    MM[15]=1.0f;
    for ( unsigned i=0; i<DIM; ++i )
    for ( unsigned j=0; j<DIM; ++j )
        MM[i+4*j] = mat(i,j);
    gym::apply(MM);
    SpaceEllipse::draw3D();
}

#else

void SpaceDynamicEllipse::draw3D() const {}

#endif

