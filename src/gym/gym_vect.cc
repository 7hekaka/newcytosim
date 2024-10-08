// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

#include "gym_vect.h"

/**
 This calculates vectors X and Y of norm N, assuming that Z is given with norm N
 Derived from `Building an Orthonormal Basis, Revisited`,
 Tom Duff et al. Journal of Computer Graphics Techniques Vol. 6 N.1, 2017
 optimized by Marc B. Reynolds
 */
void gym::orthonormal(const float z[3], float N, float x[3], float y[3])
{
    const float s = std::copysign(1.f, z[2]);
    const float a = z[1] / ( z[2] + s * N );
    const float b = z[1] * a;
    const float c = z[0] * a;
    x[0] = -z[2] - b;
    x[1] = c;
    x[2] = z[0];
    y[0] = s * c;
    y[1] = s * b - N;
    y[2] = s * z[1];
}

/**
 `R` is the transverse scaling done in the XY plane after rotation
 */
void gym::stretchAlignZ(Vector1 const& A, Vector1 const& B, float R)
{
    float X = B.XX-A.XX;
    float Y = std::copysign(R, X);
    //warning! this matrix appears here transposed
    float mat[16] = {
        0, Y, 0, 0,
        0, 0, R, 0,
        X, 0, 0, 0,
        float(A.XX), 0, 0, 1 };
    apply(mat);
}

/**
 Translate and rotate to place A in (0,0,0) and B at (0,0,1).
 Scale XY plane by `R' and Z axis by 1/|AB|
 `R` is the transverse scaling done in the XY plane after rotation
 */
void gym::stretchAlignZ(Vector2 const& A, Vector2 const& B, float R)
{
    float X(B.XX-A.XX);
    float Y(B.YY-A.YY);
    float r = R / std::sqrt(X*X+Y*Y);
    //warning! this matrix appears here transposed
    float mat[16] = {
        -r*Y, r*X, 0, 0,
        0, 0, R, 0,
        X, Y, 0, 0,
        float(A.XX), float(A.YY), 0, 1 };
    apply(mat);
}

/**
 Translate and rotate to place A in (0,0,0) and B at (0,0,1).
 Scale XY plane by `R' and Z axis by 1/|AB|
 */
void gym::stretchAlignZ(Vector3 const& A, Vector3 const& B, float R)
{
    float X(B.XX-A.XX);
    float Y(B.YY-A.YY);
    float Z(B.ZZ-A.ZZ);
    float N = R / std::sqrt(X*X+Y*Y+Z*Z);
    float vec[3] = { N*X, N*Y, N*Z };
    float mat[16] = {
        0, 0, 0, 0,
        0, 0, 0, 0,
        X, Y, Z, 0,
        float(A.XX), float(A.YY), float(A.ZZ), 1};
    gym::orthonormal(vec, std::fabs(R), mat, mat+4);
    apply(mat);
}


void gym::rotate(Vector3 const& A, Vector3 const& B, Vector3 const& C)
{
    float mat[16];
    for ( int i = 0; i < 3; ++i )
    {
        mat[i  ] = A[i];
        mat[i+4] = B[i];
        mat[i+8] = C[i];
        mat[i+12]  = 0;
        mat[i*4+3] = 0;
    }
    mat[15] = 1;
    apply(mat);
}

void gym::rotateInverse(Vector3 const& A, Vector3 const& B, Vector3 const& C)
{
    float mat[16];
    for ( int i = 0; i < 3; ++i )
    {
        mat[4*i  ] = A[i];
        mat[4*i+1] = B[i];
        mat[4*i+2] = C[i];
        mat[4*i+3] = 0;
        mat[i+12]  = 0;
    }
    mat[15] = 1;
    apply(mat);
}


// translate to center T and rotate to align with axes A, B, C
void gym::transRotate(Vector2 const& T, Vector2 const& A,
                      Vector2 const& B)
{
    //warning! this matrix is displayed here transposed
    float mat[16] = {
        (float)A.XX, (float)A.YY, 0, 0,
        (float)B.XX, (float)B.YY, 0, 0,
        0, 0, 1, 0,
        (float)T.XX, (float)T.YY, 0, 1 };
    apply(mat);
}

// translate to center T and rotate to align with axes A, B, C
void gym::transRotate(Vector3 const& T, Vector3 const& A,
                      Vector3 const& B, Vector3 const& C)
{
    //warning! this matrix is displayed here transposed
    float mat[16] = {
        (float)A.XX, (float)A.YY, (float)A.ZZ, 0,
        (float)B.XX, (float)B.YY, (float)B.ZZ, 0,
        (float)C.XX, (float)C.YY, (float)C.ZZ, 0,
        (float)T.XX, (float)T.YY, (float)T.ZZ, 1 };
    apply(mat);
}

#pragma mark - Vector1

// rotate to align Z with 'D' and translate to center 'P'
void gym::transAlignZ(Vector1 const& P, float R, Vector1 const& D)
{
    float X = std::copysign(R, float(D.XX));
    float mat[16] = {
        0, -X,  0,  0,
        0,  0, -R,  0,
        X,  0,  0,  0,
        float(P.XX), 0, 0, 1};
    apply(mat);
}

// rotate to align Z with 'D', assuming norm(D)==1, and translate to center 'P'
void gym::stretchAlignZ1(Vector1 const& P, float R, Vector1 const& D, float S)
{
    float X = R * D.XX;
    float Z = S * D.XX;
    float mat[16] = {
        0, -X,  0,  0,
        0,  0, -R,  0,
        Z,  0,  0,  0,
        float(P.XX), 0, 0, 1};
    apply(mat);
}

#pragma mark - Vector2

/**
 Rotate to align Z with 'D' and translate to center 'P',
 scale uniformly by 'R'
 */
void gym::transAlignZ(Vector2 const& P, float R, Vector2 const& D)
{
    float X(D.XX);
    float Y(D.YY);
    float n = R / std::sqrt(X*X+Y*Y);
    X *= n;
    Y *= n;
    float mat[16] = {
        Y, -X,  0,  0,
        0,  0, -R,  0,
        X,  Y,  0,  0,
        float(P.XX), float(P.YY), 0, 1};
    apply(mat);
}

/**
 Assuming norm(D)==1,
 Rotate to align Z with 'D' and translate to center 'P',
 scale XY by 'R' and Z by 'S'
 */
void gym::stretchAlignZ1(Vector2 const& P, float R, Vector2 const& D, float S)
{
    float X(D.XX);
    float Y(D.YY);
    float mat[16] = {
        R*Y, -R*X,  0,  0,
        0,      0, -R,  0,
        S*X,  S*Y,  0,  0,
        float(P.XX), float(P.YY), 0, 1};
    apply(mat);
}

#pragma mark - Vector3

/**
 Rotate to align Z with 'D' and translate to center 'P',
 scale XY by 'R'
 */
void gym::transAlignZ(Vector3 const& P, float R, Vector3 const& D)
{
    float X(D.XX);
    float Y(D.YY);
    float Z(D.ZZ);
    float N = R / std::sqrt(X*X+Y*Y+Z*Z);
    float vec[3] = { N*X, N*Y, N*Z };
    float mat[16] = {
        0, 0, 0, 0,
        0, 0, 0, 0,
        vec[0], vec[1], vec[2], 0,
        float(P.XX), float(P.YY), float(P.ZZ), 1};
    gym::orthonormal(vec, std::fabs(R), mat, mat+4);
    apply(mat);
}

/**
 Assuming norm(D)==1,
 Rotate to align Z with 'D' and translate to center 'P',
 scale XY by 'R' and Z by 'S'
 */
void gym::stretchAlignZ1(Vector3 const& P, float R, Vector3 const& D, float S)
{
    float X(D.XX);
    float Y(D.YY);
    float Z(D.ZZ);
    float vec[3] = { R*X, R*Y, R*Z };
    float mat[16] = {
        0, 0, 0, 0,
        0, 0, 0, 0,
        S*X, S*Y, S*Z, 0,
        float(P.XX), float(P.YY), float(P.ZZ), 1};
    gym::orthonormal(vec, std::fabs(R), mat, mat+4);
    apply(mat);
}
