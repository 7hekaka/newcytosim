// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
// A test for linear iterative solver (BCGS, GMRES, etc).
// FJN, 19.08.2019

#include "real.h"
#include "allocator.h"
#include "monitor.h"
#include "bicgstab.h"
#include "gmres.h"
#include "cblas.h"

/// interface for a linear system
class System
{
    /// dimension
    int   dim;
    
    /// data matrix
    real* mat;
    
public:
    
    /// initialize matrix
    System(int d)
    {
        dim = d;
        mat = new_real(d*d);
        // set all values to 0:
        zero_real(d*d, mat);
        
        // set diagonal terms to 1:
        for ( int i = 0; i < d; ++i )
            mat[i+dim*i] = 1.0;
    }
    
    /// destructor
    ~System()
    {
        free_real(mat);
    }
    
    /// size of the matrix M
    int dimension() const { return dim; };
    
    /// multiply a vector ( Y <- M * X )
    void multiply(const real* X, real* Y) const
    {
        blas::xgemv('N', dim, dim, 1.0, mat, dim, X, 1, 0.0, Y, 1);
    }
    
    /// apply preconditionner ( Y <- P * X )
    void precondition(const real* X, real* Y) const
    {
        copy_real(dim, X, Y);
    }
    
    /// multiply vector and apply preconditionner ( Y <- P * M * X )
    void multiplyP(real* X, real* Y) const
    {
        multiply(X, Y);       // Y <- M*X
        precondition(Y, X);   // X <- P*Y
    }
};


int main(int argc, char* argv[])
{
    System sys(8);
    const int dim = sys.dimension();

    LinearSolvers::Allocator alc, tmp;      // memory allocation class
    LinearSolvers::Matrix mH, mV;           // Matrices used for GMRES
    LinearSolvers::Monitor mon(dim, 0.001); // max_iteration, absolute_tolerance

    // create vectors
    real * rhs = new_real(dim);
    real * sol = new_real(dim);
    
    for ( int i; i < dim; ++i )
        rhs[i] = ( i % 4 );
    
    print_real(stdout, dim, rhs, " rhs\n");

    for ( int RS : {4, 8, 16, 32} )
    {
        mon.reset();
        zero_real(dim, sol);
        LinearSolvers::GMRES(sys, rhs, sol, RS, mon, alc, mH, mV, tmp);
        print_real(stdout, dim, sol, " sol |");
        fprintf(stdout, " GMRES-%02i  count %4i  residual %10.6f\n", RS, mon.count(), mon.residual());
    }

    free_real(sol);
    free_real(rhs);
    return EXIT_SUCCESS;
}

