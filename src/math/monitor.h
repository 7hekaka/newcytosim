// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef MONITOR_H
#define MONITOR_H

#include "real.h"
#include "blas.h"
#include "cytoblas.h"

/// Iterative Solvers
namespace LinearSolvers
{
    /// records the number of iterations, and the convergence
    class Monitor
    {
    private:
        
        /// flag
        int      flag_;
        
        /// counter for iterations or number of matrix-vector operations
        size_t   cnt_, cntMax_, cntOld_;
        
        /// desired residual
        real     target_;

        /// achieved residual
        real     res_;
        
        ///
        
    public:
        
        /// set the maximum number of iterations, and the residual threshold
        Monitor(size_t i, real r) { reset(); cntMax_ = i; target_ = r; }
        
        /// reset state variables (counters, flags and residual)
        void reset() { flag_ = 0; cnt_ = 0; res_ = INFINITY; cntOld_ = 32; }
        
        /// increment counter
        void operator ++() { ++cnt_; }
        
        /// increment counter by `i`
        void operator +=(size_t i) { cnt_ += i; }
       
        /// value of return flag
        int flag()       const { return flag_; }
        
        /// set flag to `f`
        void flag(const int f) { flag_ = f; }

        /// iteration count
        size_t count() const { return cnt_; }
        
        /// last achieved residual
        real residual()  const { return res_; }
        
        /// true if achieve residual < residual threshold
        bool converged() const { return res_ < target_; }
        
        /// check given residual and return true if threshold is achieved
        bool finished(real res)
        {
            res_ = res;
            
            if ( cnt_ > cntMax_ )
                return true;
            
            return ( res < target_ );
        }

        /// calculate residual from `x` and return true if threshold is achieved
        bool finished(size_t size, const real* x)
        {
            //fprintf(stderr, "Solver %3u residual %9.6f %9.6f\n", cnt_, blas::nrm2(size, x), blas::nrm8(size, x));
            
#if ( 1 )
            // use the 'infinite' norm (i.e. the largest element)
            real res = blas::nrm8(size, x);
#else
            // use the standard Euclidian norm:
            real res = blas::nrm2(size, x);
#endif
#if ( 1 )
            if ( cnt_ > cntOld_+128 )
            {
                if ( res > 2*res_ )
                {
                    printf("Warning: slow convergence (reduce time_step?)");
                    printf(" residual %.3e at iteration %lu, %.3e at %lu\n", res_, cntOld_, res, cnt_);
                }
                cntOld_ = cnt_;
            }
#endif
            
            if ( res != res )
            {
                fprintf(stderr, "Solver diverged at step %3lu (residual is not a number)\n", cnt_);
                return true;
            }
            
            return finished(res);
        }

        /// calculate residual from `x` and set flag to `f`
        void finish(int f, size_t size, const real* x) { flag_ = f; finished(size, x); }
    };
}

#endif

