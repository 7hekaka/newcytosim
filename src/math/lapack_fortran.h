// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

/*
 This contains C-wrappers around some of the Fortran routines of LAPACK.
 see http://www.netlib.org/lapack
 
 LAPACK contains more than 1000 linear algebra functions in LAPACK,
 but here we just propagate the functions needed in Cytosim,
 allowing direct link with the fortran LAPACK library
 */


#ifndef LAPACK_H
#define LAPACK_H

#include "real.h"

/// macro will expand to the FORTRAN function name
#if REAL_IS_DOUBLE
#   define LAPACK(x) d##x##_
#else
#   define LAPACK(x) s##x##_
#endif

extern "C"
{
int ilaenv_(int *ispec, char const*name, char const*opts, int *n1, int *n2, int *n3, int *n4);
void LAPACK(ptsv)(int*, int*, real*, real*, real*, int*, int*);
void LAPACK(ptsvx)(char*, int*, int*, const real*, const real*, real*, real*, const real*, int*, real*, int*, real*, real*, real*, real*, int*);
void LAPACK(pttrf)(int*, real*, real*, int*);
void LAPACK(pttrs)(int*, int*, const real*, const real*, real*, int*, int*);
void LAPACK(ptts2)(int*, int*, const real*, const real*, real*, int*);
void LAPACK(posv)(char*, int*, int*, real*, int*, real*, int*, int*);
void LAPACK(potrf)(char*, int*, real*, int*, int*);
void LAPACK(potrs)(char*, int*, int*, const real*, int*, real*, int*, int*);
void LAPACK(potf2)(char*, int*, real*, int*, int*);
void LAPACK(potri)(char*, int*, real*, int*, int*);
void LAPACK(pptrf)(char*, int*, real*, int*);
void LAPACK(pptrs)(char*, int*, int*, const real*, real*, int*, int*);
void LAPACK(pptri)(char*, int*, real*, int*);
void LAPACK(trtrs)(char*, char*, char*, int*, int*, const real*, int*, real*, int*, int*);
void LAPACK(gbtrf)(int*, int*, int*, int*, real*, int*, int*, int*);
void LAPACK(gbtf2)(int*, int*, int*, int*, real*, int*, int*, int*);
void LAPACK(gbtrs)(char*, int*, int*, int*, int*, const real*, int*, const int*, real*, int*, int*);
void LAPACK(gesv)(int*, int*, real*, int*, int*, real*, int*, int*);
void LAPACK(getrf)(int*, int*, real*, int*, int*, int*);
void LAPACK(getf2)(int*, int*, real*, int*, int*, int*);
void LAPACK(getrf2)(int*, int*, real*, int*, int*, int*);
void LAPACK(getri)(int*, real*, int*, const int*, real*, const int*, int*);
void LAPACK(getrs)(char*, int*, int*, const real*, int*, const int*, real*, int*, int*);
void LAPACK(gerfs)(char*, int*, int*, const real*, int*, const real*, int*, const int*, real*, int*, real*, int*, real*, real*, real*, int*, int*);
void LAPACK(laswp)(int*, const real*, int*, int*, int*, const int*, int*);
void LAPACK(sysv)(char*, int*, int*, real*, int*, int*, real*, int*, real*, int*, int*);
void LAPACK(sytrf)(char*, int*, real*, int*, int*, real*, int*, int*);
void LAPACK(sytrs)(char*, int*, int*, real*, int*, int*, real*, int*, int*);
void LAPACK(syev)(char*, char*, int*, real*, int*, real*, real*, int*, int*);
void LAPACK(syevd)(char*, char*, int*, real*, int*, real*, real*, int*, int*, int*, int*);
void LAPACK(syevx)(char*, char*, char*, int*, real*, int*, real*, real*, int*, int*, real*, int*, real*, real*, int*, real*, int*, const int*, int*, int*);
void LAPACK(gtsv)(int*, int*, real*, real*, real*, real*, int*, int*);
void LAPACK(gttrf)(int*, real*, real*, real*, real*, int*, int*);
void LAPACK(gttrs)(char*, int*, int*, real*, real*, real*, real*, int*, real*, int*, int*);
void LAPACK(geqrf)(int*, int*, real*, int*, real*, real*, const int*, int*);
void LAPACK(ormqr)(char*, char*, int*, int*, int*, const real*, int*, const real*, real*, int*, real*, const int*, int*);
void LAPACK(gels)(char*, int*, int*, int*, const real*, int*, const real*, int*, real*, const int*, int*);
}


namespace lapack
{

inline int ilaenv(int ispec, char const*name, char const*opts, int n1, int n2, int n3, int n4)
{
    return ilaenv_(&ispec, name, opts, &n1, &n2, &n3, &n4);
}
    
inline void xptsv(int N, int NRHS, real* D, real* E, real* B, int LDB, int* INFO)
{
    LAPACK(ptsv)(&N, &NRHS, D, E, B, &LDB, INFO);
}

inline void xptsvx(char fact, int N, int NRHS, const real* D, const real* E, real* DF, real* EF,
                   const real* B, int LDB, real* X, int LDX, real* RCOND, real* FERR, real* BERR, real* work, int* INFO)
{
    LAPACK(ptsvx)(&fact, &N, &NRHS, D, E, DF, EF, B, &LDB, X, &LDX, RCOND, FERR, BERR, work, INFO);
}

inline void xpttrf(int N, real* D, real* E, int* INFO)
{
    LAPACK(pttrf)(&N, D, E, INFO);
}

inline void xpttrs(int N, int NRHS, const real* D, const real* E, real* B, int LDB, int* INFO)
{
    LAPACK(pttrs)(&N, &NRHS, D, E, B, &LDB, INFO);
}

inline void xptts2(int N, int NRHS, const real* D, const real* E, real* B, int LDB)
{
    LAPACK(ptts2)(&N, &NRHS, D, E, B, &LDB);
}

inline void xptts2(int N, const real* D, const real* E, real* B)
{
    int NRHS = 1;
    int LDB = 1;
    LAPACK(ptts2)(&N, &NRHS, D, E, B, &LDB);
}

inline void xposv(char UPLO, int N, int NRHS, real* A, int LDA, real* B, int LDB, int* INFO)
{
    LAPACK(posv)(&UPLO, &N, &NRHS, A, &LDA, B, &LDB, INFO);
}

inline void xpotrf(char UPLO, int N, real* A, int LDA, int* INFO)
{
    LAPACK(potrf)(&UPLO, &N, A, &LDA, INFO);
}

inline void xpotrs(char UPLO, int N, int NRHS, const real* A, int LDA, real* B, int LDB, int* INFO)
{
    LAPACK(potrs)(&UPLO, &N, &NRHS, A, &LDA, B, &LDB, INFO);
}

inline void xpotf2(char UPLO, int N, real* A, int LDA, int* INFO)
{
    LAPACK(potf2)(&UPLO, &N, A, &LDA, INFO);
}

inline void xpotri(char UPLO, int N, real* A, int LDA, int* INFO)
{
    LAPACK(potri)(&UPLO, &N, A, &LDA, INFO);
}

inline void xpptrf(char UPLO, int N, real* A, int* INFO)
{
    LAPACK(pptrf)(&UPLO, &N, A, INFO);
}

inline void xpptrs(char UPLO, int N, int NRHS, const real* A, real* B, int LDB, int* INFO)
{
    LAPACK(pptrs)(&UPLO, &N, &NRHS, A, B, &LDB, INFO);
}

inline void xpptri(char UPLO, int N, real* A, int* INFO)
{
    LAPACK(pptri)(&UPLO, &N, A, INFO);
}

inline void xtrtrs(char UPLO, char trans, char diag, int N, int NRHS, const real* A, int LDA, real* B, int LDB, int* INFO)
{
    LAPACK(trtrs)(&UPLO, &trans, &diag, &N, &NRHS, A, &LDA, B, &LDB, INFO);
}

inline void xgbtrf(int M, int N, int KL, int KU, real* AB, int LDAB, int* IPIV, int* INFO)
{
    LAPACK(gbtrf)(&M, &N, &KL, &KU, AB, &LDAB, IPIV, INFO);
}

inline void xgbtf2(int M, int N, int KL, int KU, real* AB, int LDAB, int* IPIV, int* INFO)
{
    LAPACK(gbtf2)(&M, &N, &KL, &KU, AB, &LDAB, IPIV, INFO);
}

inline void xgbtrs(char trans, int N, int KL, int KU, int NRHS, const real* AB, int LDAB, const int* IPIV, real* B, int LDB, int* INFO)
{
    LAPACK(gbtrs)(&trans, &N, &KL, &KU, &NRHS, AB, &LDAB, IPIV, B, &LDB, INFO);
}

inline void xgesv(int N, int NRHS, real* A, int LDA, int* IPIV, real* B, int LDB, int* INFO)
{
    LAPACK(gesv)(&N, &NRHS, A, &LDA, IPIV, B, &LDB, INFO);
}

inline void xgetrf(int M, int N, real* A, int LDA, int* IPIV, int* INFO)
{
    LAPACK(getrf)(&M, &N, A, &LDA, IPIV, INFO);
}

inline void xgetf2(int M, int N, real* A, int LDA, int* IPIV, int* INFO)
{
    LAPACK(getf2)(&M, &N, A, &LDA, IPIV, INFO);
}

inline void xgetrf2(int M, int N, real* A, int LDA, int* IPIV, int* INFO)
{
    LAPACK(getrf2)(&M, &N, A, &LDA, IPIV, INFO);
}

inline void xgetri(int N, real* A, int LDA, const int* IPIV, real* WORK, int LWORK, int* INFO)
{
    LAPACK(getri)(&N, A, &LDA, IPIV, WORK, &LWORK, INFO);
}

inline void xgetrs(char trans, int N, int NRHS, const real* A, int LDA, const int* IPIV, real* B, int LDB, int* INFO)
{
    LAPACK(getrs)(&trans, &N, &NRHS, A, &LDA, IPIV, B, &LDB, INFO);
}

inline void xgerfs(char trans, int N, int NRHS, const real* A, int LDA, const real* AF, int LDAF, const int* IPIV, real* B, int LDB, real* X, int LDX, real* FERR, real* BERR, real* WORK, int* IWORK, int* INFO)
{
    LAPACK(gerfs)(&trans, &N, &NRHS, A, &LDA, AF, &LDAF, IPIV, B, &LDB, X, &LDX, FERR, BERR, WORK, IWORK, INFO);
}

inline void xlaswp(int N, const real* A, int LDA, int K1, int K2, const int* IPIV, int INCX)
{
    LAPACK(laswp)(&N, A, &LDA, &K1, &K2, IPIV, &INCX);
}

inline void xsysv(char UPLO, int N, int NRHS, real* A, int LDA, int* IPIV, real*B, int LDB, real* WORK, int LWORK, int* INFO)
{
    LAPACK(sysv)(&UPLO, &N, &NRHS, A, &LDA, IPIV, B, &LDB, WORK, &LWORK, INFO);
}

inline void xsytrf(char UPLO, int N, real* A, int LDA, int* IPIV, real* WORK, int LWORK, int* INFO)
{
    LAPACK(sytrf)(&UPLO, &N, A, &LDA, IPIV, WORK, &LWORK, INFO);
}

inline void xsytrs(char UPLO, int N, int NRHS, real* A, int LDA, int* IPIV, real* B, int LDB, int* INFO)
{
    LAPACK(sytrs)(&UPLO, &N, &NRHS, A, &LDA, IPIV, B, &LDB, INFO);
}

inline void xsyev(char JOBZ, char UPLO, int N, real* A, int LDA, real* W, real* WORK, int LWORK, int* INFO)
{
    LAPACK(syev)(&JOBZ, &UPLO, &N, A, &LDA, W, WORK, &LWORK, INFO);
}

inline void xsyevd(char JOBZ, char UPLO, int N, real* A, int LDA, real* W, real* WORK, int LWORK, int* IWORK, int LIWORK, int* INFO)
{
    LAPACK(syevd)(&JOBZ, &UPLO, &N, A, &LDA, W, WORK, &LWORK, IWORK, &LIWORK, INFO);
}

inline void xsyevx(char JOBZ, char RANGE, char UPLO, int N, real* A, int LDA, real VL, real VU, int IL, int IU,
                   real ABSTOL, int* M, real* W, real* Z, int LDZ, real* WORK, int LWORK, const int* IWORK, int* IFAIL, int* INFO)
{
    LAPACK(syevx)(&JOBZ, &RANGE, &UPLO, &N, A, &LDA, &VL, &VU, &IL, &IU, &ABSTOL, M, W, Z, &LDZ, WORK, &LWORK, IWORK, IFAIL, INFO);
}

inline void xgtsv(int N, int NRHS, real* DL, real* D, real* DU, real* B, int LDB, int* INFO)
{
    LAPACK(gtsv)(&N, &NRHS, DL, D, DU, B, &LDB, INFO);
}

inline void xgttrf(int N, real* DL, real* D, real* DU, real* DU2, int* IPIV, int* INFO)
{
    LAPACK(gttrf)(&N, DL, D, DU, DU2, IPIV, INFO);
}

inline void xgttrs(char TRANS, int N, int NRHS, real* DL, real* D, real* DU, real* DU2, int* IPIV, real* B, int LDB, int* INFO)
{
    LAPACK(gttrs)(&TRANS, &N, &NRHS, DL, D, DU, DU2, IPIV, B, &LDB, INFO);
}

inline void xgeqrf(int M, int N, real* A, int LDA,  real* TAU, real* WORK, int LWORK, int* INFO)
{
    LAPACK(geqrf)(&M, &N, A, &LDA, TAU, WORK, &LWORK, INFO);
}

inline void xormqr(char side, char trans, int M, int N, int K, const real* A, int LDA, const real* TAU,
                   real* C, int LDC, real* WORK, int LWORK, int* INFO)
{
    LAPACK(ormqr)(&side, &trans, &M, &N, &K, A, &LDA, TAU, C, &LDC, WORK, &LWORK, INFO);
}

inline void xgels(char trans, int M, int N, int NRHS, const real* A, int LDA, const real* B, int LDB, real* WORK, int LWORK, int* INFO)
{
    LAPACK(gels)(&trans, &M, &N, &NRHS, A, &LDA, B, &LDB, WORK, &LWORK, INFO);
}

}

#endif
