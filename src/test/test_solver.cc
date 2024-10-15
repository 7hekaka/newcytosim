// Cytosim was created by Francois Nedelec. Copyright 2023 Cambridge University.
/*
 A test for linear iterative solver (BCGS, etc):
 Read 'A' from "matrix.mtx" and 'b' from "vector.mtx" and solves the linear system
     A.x = b
 using various iterative solvers.
 FJN, 19.08.2019
*/

#include <cstdio>

#include "real.h"
#include "real_print.h"
#include "blas.h"
#include "cytoblas.h"
#include "monitor.h"
#include "allocator.h"
#include "bicgstab.h"

/// interface for a linear system
class System
{
    /// dimension
    int dim;
    
    /// data matrix
    real* mat;
    
public:
    
    /// initialize matrix
    System()
    {
        dim = 0;
        mat = nullptr;
    }
    
    /// destructor
    ~System()
    {
        free_real(mat);
    }
    
    /// allocate for given size
    void allocate(int d)
    {
        dim = d;
        free_real(mat);
        mat = new_real(d*d);
        zero_real(d*d, mat);
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

    
    /// read MatrixMarket format
    int readMatrix(FILE * file)
    {
        constexpr size_t MAX = 1024;
        char str[MAX], * ptr;
        do {
            if ( !fgets(str, MAX, file) ) return 1;
            // skip comments:
        } while ( str[0] == '%' );
        // parse dimension line:
        unsigned long lin = strtoul(str, &ptr, 10);
        unsigned long col = strtoul(ptr, &ptr, 10);
        unsigned long cnt = strtoul(ptr, &ptr, 10);
        if ( lin != col )
            return 2;
        allocate(lin);
        for ( size_t i = 0; i < cnt; ++i )
        {
            if ( !fgets(str, MAX, file) ) return 3;
            lin = strtoul(str, &ptr, 10);
            col = strtoul(ptr, &ptr, 10);
            real val = strtof(ptr, &ptr);
            mat[lin+dim*col] = val;
        }
        return 0;
    }
    
    /// read matrix from file
    int readMatrix(const char filename[])
    {
        int err = 1;
        FILE * f = fopen(filename, "r");
        if ( f )
        {
            if ( ~ferror(f) )
            {
                printf(" reading `%s`", filename);
                err = readMatrix(f);
                if ( err )
                    printf(" ... failed (error %i)", err);
                printf("\n");
            }
            fclose(f);
        }
        else
            printf("file not found `%s`\n", filename);
        return err;
    }
};


int readVector(FILE * file, size_t dim, real * vec)
{
    constexpr size_t MAX = 1024;
    char str[MAX];
    do {
        if ( !fgets(str, MAX, file) ) return 1;
        // skip comments:
    } while ( str[0] == '%' );
    // parse dimension line:
    unsigned long cnt = strtoul(str, nullptr, 10);
    for ( size_t i = 0; i < cnt; ++i )
    {
        if ( !fgets(str, MAX, file) ) return 3;
        real val = strtof(str, nullptr);
        if ( i < dim ) vec[i] = val;
    }
    return 0;
}


int readVector(const char filename[], size_t dim, real * vec)
{
    int err = 1;
    FILE * f = fopen(filename, "r");
    if ( f )
    {
        if ( ~ferror(f) )
        {
            printf(" reading `%s`", filename);
            err = readVector(f, dim, vec);
            if ( err )
                printf("... failed (error %i)", err);
            printf("\n");
        }
        fclose(f);
    }
    else
        printf("file not found `%s`\n", filename);
    return err;
}


int main(int argc, char* argv[])
{
    System sys;
    if ( sys.readMatrix("matrix.mtx") )
        return 1;
    const int dim = sys.dimension();

    LinearSolvers::Allocator alc, tmp;      // memory allocation class
    LinearSolvers::Monitor mon(2*dim, 0.001); // max_iteration, absolute_tolerance

    // create vectors
    real * rhs = new_real(dim);

    // get system's right-hand-side
    zero_real(dim, rhs);
    if ( readVector("vector.mtx", dim, rhs) )
    {
        free_real(rhs);
        return 2;
    }
    
    real * sol = new_real(dim);
    real * vec = new_real(dim);

    print_real(stdout, std::min(16, dim), rhs, " rhs |");
    fprintf(stdout, "\n");

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
        fprintf(stdout, " BiCGStab count %4u  residual %10.6f\n", mon.count(), res);
    }

    free_real(sol);
    free_real(rhs);
    free_real(vec);
}

