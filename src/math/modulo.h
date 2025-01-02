// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#ifndef MODULO_H
#define MODULO_H

#include "real.h"
#include "vector.h"

///\def compile switch to enable/disable support for periodic boundary conditions
#define ENABLE_PERIODIC_BOUNDARIES 1


/// Modulo is a helper class used to implement periodic boundary conditions
/**
 This class implements periodic boundaries conditions for Vector,
 in all dimensions and selectively to only some of the coordinates: X, Y or Z.
 
 We follow the method B (Do not restrict the particle coordinates) described in
 [wikipedia](https://en.wikipedia.org/wiki/Periodic_boundary_conditions)

 In this way, we can tell the movements of the filaments, and track all paths
 without any loss of information. Attention: our origin is in the box center,
 unlike on Wikipedia’s description.
 */
class Modulo
{
public:
    
    /// the period in each dimension
    real period_[4];
    
    /// inverse of the period in each dimension
    real inv_period_[4];

    /// bitfield indicating the dimensions that are periodic
    int mMode;

    
    /// adjust 'x' to canonical image in dimension i
    inline real fold_(const real x, int i) const
    {
        const real P = period_[i];
        //return std::remainder(x, P);
        long q = (long)std::nearbyint(x * inv_period_[i]);
        /*
        real y = std::remainder(x, P);
        if ( std::fabs(y-x+w*P) > 0.001 )
            std::clog << "   " << x-q*P << "   " << y << "\n";
         */
        return x - q * P;
    }
    
    /// adjust 'x' to canonical image in dimension i
    inline float foldf(const float x, int i) const
    {
        const real P = period_[i];
        long q = (long)std::nearbyint(x * inv_period_[i]);
        return float(x - q * P);
    }

public:
    
    /// set as non periodic
    void reset() { mMode = 0; for (int d=0; d<4; ++d) period_[d] = 0; }
    
    /// constructor
    Modulo() { reset(); }

    /// destructor
    ~Modulo() {}
    
    /// disable periodicity in all dimensions
    void disable() { mMode = 0; }
    
    /// enable periodicity in dimension 'd'
    void enablePeriodic(int d, real size);
    
    /// true if at least one direction has periodic boundaries
    bool isPeriodic() const { return mMode; }

    /// true if direction `d` has periodic boundaries
    bool isPeriodic(int d) const { return mMode & (1<<d); }
    
    /// return the d-th direction of periodicity
    Vector3 period(int d) const;
    
    /// shift `pos` to its canonical image, which is the one closest to the origin
    void fold(Vector& pos) const;
    
    /// shift `pos` to its image which is closest to `ref`
    void fold(Vector& pos, Vector const& ref) const;
    
    /// return translation necessary to bring `pos` to its canonical image
    Vector offset(Vector const& pos) const;
    
    /// set `pos` to its canonical image, and return offset = pos - fold(pos)
    Vector foldOffset(Vector& pos) const;
    
    /// shift `pos` to its canonical image
    void fold_float(float* pos) const;

    /// shift `pos` to its image which is closest to `ref`
    void fold_float(float* pos, float const* ref) const;

};


#if ENABLE_PERIODIC_BOUNDARIES

/**
 This is a global variable that is initialized in Simul
 It is used to implement periodic boundary conditions
 */
extern Modulo const* modulo;

#else

/**
 Any code under 'if ( modulo )' will not be executed, and should even
 be discarded during compilation if optimizations are enabled.
 */
constexpr Modulo * modulo = nullptr;

#endif

#endif


