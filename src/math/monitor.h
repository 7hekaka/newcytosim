// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#ifndef MONITOR_H
#define MONITOR_H

#include "real.h"
#include "blas.h"
#include "cytoblas.h"

/// Iterative methods to solve a system of linear equations
namespace LinearSolvers
{
    /// records the number of iterations, and the convergence
    class Monitor
    {
    private:
        
        /// desired residual
        real target_;

        /// residual achieved in last step
        real residual_;
        
        /// maximum allowed number of iterations
        unsigned limit_;

        /// counter for number of matrix-vector operations
        unsigned count_;
        
        /// exit flag
        unsigned flag_;
        
        /// true if achieved desired residual
        bool converged_;

    public:
        
        /// set the maximum number of iterations, and the residual threshold
        Monitor(unsigned i, real r) { reset(); limit_=i; target_=r; }
        
        /// reset state variables (counters, flags and residual)
        void reset() { flag_=0; count_=0; residual_=INFINITY; converged_=0; }
        
        /// increment counter
        void operator ++() { ++count_; }
        
        /// increment counter by `i`
        void operator +=(size_t i) { count_ += i; }
       
        /// value of return flag
        unsigned flag() const { return flag_; }
        
        /// set flag to `f`
        void flag(unsigned f) { flag_ = f; }

        /// iteration count
        unsigned count() const { return count_; }
        
        /// last achieved residual
        real residual()  const { return residual_; }
        
        /// true if achieve residual < residual threshold
        bool converged() const { return converged_; }
        
        /// true if converged or number of multiplications exceeded
        bool finished() const { return converged_ | ( count_ > limit_ ); }

        /// register given residual and return true if target is achieved
        bool finished(real res)
        {
            residual_ = res;
            converged_ = ( res < target_ );
            return converged_ | ( count_ > limit_ );
        }

        /// calculate residual from `vec` and return true if threshold is achieved
        bool finished(size_t size, const real* vec)
        {
            //fprintf(stderr, "   Monitor %4u residual %12.6f %9.6f\n", count_, blas::nrm2(size, vec), blas::nrm8(size, vec));
#if ( 1 )
            // use the 'infinite' norm (i.e. the largest element)
            real res = blas::nrm8(size, vec);
#else
            // use the standard Euclidian norm:
            real res = blas::nrm2(size, vec);
#endif
            //fprintf(stderr, "   Monitor %4u  residual %12.6f\n", count_, res);
            if ( res != res )
            {
                fprintf(stderr, "Invalid solution at step %3u : %lu NaNs)\n", count_, has_nan(size, vec));
                return true;
            }
            return finished(res);
        }

        /// calculate residual from `vec` and return `residual < margin * target`
        bool finished(size_t size, const real* vec, real margin)
        {
            finished(size, vec);
            return ( residual_ < margin * target_ ) | ( count_ > limit_ );
        }
        
        /// calculate residual from `vec` and set flag to `f`
        void finish(unsigned f, size_t size, const real* vec) { flag_ = f; finished(size, vec); }
    };
}

#endif

