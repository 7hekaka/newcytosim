// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
// A test for linear iterative solver (BCGS, GMRES, etc).
// FJN, 19.08.2019

#include <cstdio>
#include "real.h"
#include "cblas.h"
#include "monitor.h"
#include "allocator.h"
#include "bicgstab.h"
#include "gmres.h"

/// interface for a linear system
class System
{
    /// dimension
    int   dim;
    
    /// data matrix
    real* mat;
    
public:
    
    /// initialize matrix
    System()
    {
        mat = nullptr;
    }
    
    /// destructor
    ~System()
    {
        free_real(mat);
    }
    
    /// allocate for given size
    void allocate(size_t d)
    {
        dim = d;
        free_real(mat);
        mat = new_real(d*d);
        zero_real(d*d, mat);
    }
    
    /// size of the matrix M
    size_t dimension() const { return dim; };
    
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

    
    /// read MatrixMarket format
    int read(FILE * file)
    {
        constexpr size_t MAX = 1024;
        char str[MAX], * ptr;
        do {
            if ( 0 == fgets(str, MAX, file) ) return 1;
            // skip comments:
        } while ( str[0] == '%' );
        // parse dimension line:
        printf(" reading matrix: %s", str);
        size_t lin = strtoul(str, &ptr, 10);
        size_t col = strtoul(ptr, &ptr, 10);
        size_t cnt = strtoul(ptr, &ptr, 10);
        if ( lin != col )
            return 2;
        allocate(lin);
        for ( size_t i = 0; i < cnt; ++i )
        {
            if ( 0 == fgets(str, MAX, file) ) return 3;
            lin = strtoul(str, &ptr, 10);
            col = strtoul(ptr, &ptr, 10);
            real val = strtof(ptr, &ptr);
            mat[lin+dim*col] = val;
        }
        return 0;
    }
    
    /// read matrix from file
    int read(const char filename[])
    {
        int err = 0;
        FILE * f = fopen(filename, "r");
        if ( f && ~ferror(f) )
        {
            err = read(f);
            if ( err )
                fprintf(stderr, "failed to read matrix (error %i)\n", err);
            fclose(f);
        }
        return err;
    }
};


int readVector(FILE * file, size_t dim, real * vec)
{
    constexpr size_t MAX = 1024;
    char str[MAX];
    do {
        if ( 0 == fgets(str, MAX, file) ) return 1;
        // skip comments:
    } while ( str[0] == '%' );
    // parse dimension line:
    printf(" reading vector: %s", str);
    size_t cnt = strtoul(str, 0, 10);
    for ( size_t i = 0; i < cnt; ++i )
    {
        if ( 0 == fgets(str, MAX, file) ) return 3;
        real val = strtof(str, 0);
        if ( i < dim ) vec[i] = val;
    }
    return 0;
}


int readVector(const char filename[], size_t dim, real * vec)
{
    int err = 0;
    FILE * f = fopen(filename, "r");
    if ( f && ~ferror(f) )
    {
        err = readVector(f, dim, vec);
        if ( err )
            fprintf(stderr, "failed to read vector (error %i)\n", err);
    }
    fclose(f);
    return err;
}


int main(int argc, char* argv[])
{
    System sys;
    sys.read("matrix.mtx");
    const int dim = sys.dimension();

    LinearSolvers::Matrix mH, mV;           // Matrices used for GMRES
    LinearSolvers::Allocator alc, tmp;      // memory allocation class
    LinearSolvers::Monitor mon(2*dim, 0.001); // max_iteration, absolute_tolerance

    // create vectors
    real * rhs = new_real(dim);
    real * sol = new_real(dim);
    real * vec = new_real(dim);

    // get system's right-hand-side
    zero_real(dim, rhs);
    readVector("rhs.mtx", dim, rhs);
    print_real(stdout, std::min(16, dim), rhs, " rhs\n");

    if ( 1 )
    {
        mon.reset();
        zero_real(dim, sol);
        LinearSolvers::BCGS(sys, rhs, sol, mon, alc);
        print_real(stdout, std::min(16, dim), sol, " sol |");
        
        // calculate true residual:
        sys.multiply(sol, vec);
        blas::xaxpy(dim, -1.0, rhs, 1, vec, 1);
        real res = blas::nrm2(dim, vec);
        fprintf(stdout, " BiCGStab count %4i  residual %10.6f\n", mon.count(), res);
    }

    for ( int RS : {2, 4, 8, 16, 32, 64, 128} )
    {
        mon.reset();
        zero_real(dim, sol);
        LinearSolvers::GMRES(sys, rhs, sol, RS, mon, alc, mH, mV, tmp);
        print_real(stdout, std::min(16, dim), sol, " sol |");
        
        // calculate true residual:
        sys.multiply(sol, vec);
        blas::xaxpy(dim, -1.0, rhs, 1, vec, 1);
        real res = blas::nrm2(dim, vec);
        fprintf(stdout, " GMRES%03i count %4i  residual %10.6f\n", RS, mon.count(), res);
    }

    free_real(sol);
    free_real(rhs);
    free_real(vec);
    return EXIT_SUCCESS;
}

