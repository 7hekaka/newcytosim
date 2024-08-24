// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

#include "gym_matrix.h"
#include <stdio.h>
#include <math.h>

/*
mat4x4_inverse is derived from GLFW linmath.h by Camilla Löwy <elmindreda@elmindreda.org>
*/

void mat_print(FILE * f, const float mat[16])
{
    char se[] = "/||\\";
    for ( int i = 0; i < 4; ++i )
    {
        fprintf(f, "%c", se[i]);
        for ( int j = 0; j < 4; ++j )
            fprintf(f, "%6.2f ", mat[i+4*j]);
        fprintf(f, "%c\n", se[3-i]);
    }
}

void gym::mat_rotation(float M[16], float X, float Y, float Z, float C, float S)
{
    float dX = X - C * X, dY = Y - C * Y, dZ = Z - C * Z;
    float sX =     S * X, sY =     S * Y, sZ =     S * Z;

    M[ 0] = dX * X + C;
    M[ 1] = dY * X + sZ;
    M[ 2] = dZ * X - sY;
    M[ 3] = 0.f;

    M[ 4] = dX * Y - sZ;
    M[ 5] = dY * Y + C;
    M[ 6] = dZ * Y + sX;
    M[ 7] = 0.f;

    M[ 8] = dX * Z + sY;
    M[ 9] = dY * Z - sX;
    M[10] = dZ * Z + C;
    M[11] = 0.f;

    M[12] = 0.f;
    M[13] = 0.f;
    M[14] = 0.f;
    M[15] = 1.f;
}

void gym::mat_ortho(float M[16], float L, float R, float B, float T, float N, float F)
{
    mat_zero(M);
    
    M[ 0] = -2.f / ( L - R );
    M[ 5] = -2.f / ( B - T );
    M[10] =  2.f / ( N - F );
    
    M[12] = ( R + L ) / ( L - R );
    M[13] = ( T + B ) / ( B - T );
    M[14] = ( F + N ) / ( N - F );
    M[15] = 1.f;
}

void gym::mat_frustum(float M[16], float L, float R, float B, float T, float N, float F)
{
    mat_zero(M);

    M[ 0] = 2.f * N / ( R - L );
    M[ 5] = 2.f * N / ( T - B );
    M[ 8] = ( R + L ) / ( R - L );
    M[ 9] = ( T + B ) / ( T - B );
    M[10] = ( F + N ) / ( N - F );
    M[11] = -1.f;
    M[14] = ( F * N ) * (2.f / ( N - F ));
}

void gym::mat_perspective(float M[16], float y_fov, float aspect, float N, float F)
{
    float const a = 1.f / tanf(y_fov / 2.f);

    mat_zero(M);

    M[ 0] = a / aspect;
    M[ 5] = a;
    M[10] = ( F + N ) / ( N - F );
    M[11] = -1.f;
    M[14] = ( F * N ) * (2.f / ( N - F ));
}

void gym::mat_pick(float M[16], float cx, float cy, float dx, float dy, const int vp[4])
{
    mat_diagonal(M, 1);
    float tx = ( vp[2] - 2.f * ( cx - vp[0] )) / dx;
    float ty = ( vp[3] - 2.f * ( cy - vp[1] )) / dy;
    mat_translate(M, tx, ty, 0.f);
    mat_scale(M, vp[2]/dx, vp[3]/dy, 1.f);
}


/** out <- M x in */
void gym::mat_mulvec(float out[4], const float M[16], const float in[4])
{
    for ( int i = 0; i < 4; ++i )
        out[i] = in[0] * M[i] + in[1] * M[4+i] + in[2] * M[8+i] + in[3] * M[12+i];
}

/** Attention: M <- A x B */
void gym::mat_multiply(float M[16], const float A[16], const float B[16])
{
    for ( int i = 0; i < 4; ++i )
        mat_mulvec(M+4*i, A, B+4*i);
}

/** Attention: M <- M x B */
void gym::mat_multiply(float M[16], const float B[16])
{
    float tmp[16];
    mat_copy(tmp, M);
    for ( int i = 0; i < 4; ++i )
        mat_mulvec(M+4*i, tmp, B+4*i);
}

void gym::mat_rotateX(float M[16], float C, float S)
{
    for ( int i = 0; i < 4; ++i )
    {
        float X = M[4+i];
        float Y = M[8+i];
        M[4+i] = C * X + S * Y;
        M[8+i] = C * Y - S * X;
    }
}

void gym::mat_rotateY(float M[16], float C, float S)
{
    for ( int i = 0; i < 4; ++i )
    {
        float X = M[8+i];
        float Y = M[0+i];
        M[8+i] = C * X + S * Y;
        M[0+i] = C * Y - S * X;
    }
}

void gym::mat_rotateZ(float M[16], float C, float S)
{
    for ( int i = 0; i < 4; ++i )
    {
        float X = M[0+i];
        float Y = M[4+i];
        M[0+i] = C * X + S * Y;
        M[4+i] = C * Y - S * X;
    }
}

void gym::mat_scale(float M[16], float S)
{
    for ( int i = 0; i < 12; ++i )
        M[i] *= S;
}

void gym::mat_scale(float M[16], float X, float Y, float Z)
{
    for ( int i = 0; i < 4; ++i )
        M[i] *= X;
    for ( int i = 0; i < 4; ++i )
        M[4+i] *= Y;
    for ( int i = 0; i < 4; ++i )
        M[8+i] *= Z;
}

void gym::mat_translate(float M[16], float X, float Y, float Z)
{
    for ( int i = 0; i < 4; ++i )
        M[12+i] += X * M[i] + Y * M[4+i] + Z * M[8+i];
}

void gym::mat_copytrans(float M[16], const float R[16], float X, float Y, float Z)
{
    for ( int i = 0; i < 12; ++i )
        M[i] = R[i];
    for ( int i = 0; i < 4; ++i )
        M[12+i] = R[12+i] + X * R[i] + Y * R[4+i] + Z * R[8+i];
}

void gym::mat_translate(float M[16], const float R[16], float X, float Y, float Z)
{
    for ( int i = 0; i < 4; ++i )
        M[12+i] = R[12+i] + X * R[i] + Y * R[4+i] + Z * R[8+i];
}

void gym::mat_transscale(float M[16], float X, float Y, float Z, float S)
{
    for ( int i = 0; i < 4; ++i )
    {
        M[12+i] += X * M[i] + Y * M[4+i] + Z * M[8+i];
        M[i  ] *= S;
        M[4+i] *= S;
        M[8+i] *= S;
    }
}

void gym::mat_transscale(float M[16], const float R[16], float X, float Y, float Z, float S)
{
    for ( int i = 0; i < 12; i += 4 )
    {
        M[i  ] = R[i  ] * S;
        M[i+1] = R[i+1] * S;
        M[i+2] = R[i+2] * S;
        M[i+3] = R[i+3] * S;
    }
    for ( int i = 0; i < 4; ++i )
        M[12+i] = R[12+i] + X * R[i] + Y * R[i+4] + Z * R[i+8];
}


void gym::mat_unrotate(float T[16], const float M[16])
{
    float SX = sqrtf(M[0] * M[0] + M[1] * M[1] + M[2] * M[2]);
    float SY = sqrtf(M[4] * M[4] + M[5] * M[5] + M[6] * M[6]);
    float SZ = sqrtf(M[8] * M[8] + M[9] * M[9] + M[10] * M[10]);
    float X = M[12];
    float Y = M[13];
    float Z = M[14];
    mat_zero(T);
    T[ 0] = SX;
    T[ 5] = SY;
    T[10] = SZ;
    T[14] = 0.f;
    T[15] = 1.f;
    T[12] = X;
    T[13] = Y;
    T[14] = Z;
}


int gym::mat4x4_inverse(float T[16], const float M[16])
{
    float s[6];
    float c[6];
    s[0] = M[0+4*0]*M[1+4*1] - M[1+4*0]*M[0+4*1];
    s[1] = M[0+4*0]*M[1+4*2] - M[1+4*0]*M[0+4*2];
    s[2] = M[0+4*0]*M[1+4*3] - M[1+4*0]*M[0+4*3];
    s[3] = M[0+4*1]*M[1+4*2] - M[1+4*1]*M[0+4*2];
    s[4] = M[0+4*1]*M[1+4*3] - M[1+4*1]*M[0+4*3];
    s[5] = M[0+4*2]*M[1+4*3] - M[1+4*2]*M[0+4*3];
    
    c[0] = M[2+4*0]*M[3+4*1] - M[3+4*0]*M[2+4*1];
    c[1] = M[2+4*0]*M[3+4*2] - M[3+4*0]*M[2+4*2];
    c[2] = M[2+4*0]*M[3+4*3] - M[3+4*0]*M[2+4*3];
    c[3] = M[2+4*1]*M[3+4*2] - M[3+4*1]*M[2+4*2];
    c[4] = M[2+4*1]*M[3+4*3] - M[3+4*1]*M[2+4*3];
    c[5] = M[2+4*2]*M[3+4*3] - M[3+4*2]*M[2+4*3];
    
    float det = s[0]*c[5] - s[1]*c[4] + s[2]*c[3] + s[3]*c[2] - s[4]*c[1] + s[5]*c[0];
    
    if ( det == 0 )
        return 1;
    det = 1.f / det;

    T[0+4*0] = M[1+4*1] * c[5] - M[1+4*2] * c[4] + M[1+4*3] * c[3];
    T[1+4*0] =-M[1+4*0] * c[5] + M[1+4*2] * c[2] - M[1+4*3] * c[1];
    T[2+4*0] = M[1+4*0] * c[4] - M[1+4*1] * c[2] + M[1+4*3] * c[0];
    T[3+4*0] =-M[1+4*0] * c[3] + M[1+4*1] * c[1] - M[1+4*2] * c[0];

    T[0+4*1] =-M[0+4*1] * c[5] + M[0+4*2] * c[4] - M[0+4*3] * c[3];
    T[1+4*1] = M[0+4*0] * c[5] - M[0+4*2] * c[2] + M[0+4*3] * c[1];
    T[2+4*1] =-M[0+4*0] * c[4] + M[0+4*1] * c[2] - M[0+4*3] * c[0];
    T[3+4*1] = M[0+4*0] * c[3] - M[0+4*1] * c[1] + M[0+4*2] * c[0];

    T[0+4*2] = M[3+4*1] * s[5] - M[3+4*2] * s[4] + M[3+4*3] * s[3];
    T[1+4*2] =-M[3+4*0] * s[5] + M[3+4*2] * s[2] - M[3+4*3] * s[1];
    T[2+4*2] = M[3+4*0] * s[4] - M[3+4*1] * s[2] + M[3+4*3] * s[0];
    T[3+4*2] =-M[3+4*0] * s[3] + M[3+4*1] * s[1] - M[3+4*2] * s[0];

    T[0+4*3] =-M[2+4*1] * s[5] + M[2+4*2] * s[4] - M[2+4*3] * s[3];
    T[1+4*3] = M[2+4*0] * s[5] - M[2+4*2] * s[2] + M[2+4*3] * s[1];
    T[2+4*3] =-M[2+4*0] * s[4] + M[2+4*1] * s[2] - M[2+4*3] * s[0];
    T[3+4*3] = M[2+4*0] * s[3] - M[2+4*1] * s[1] + M[2+4*2] * s[0];
    
    for (int i = 0; i < 16; i++)
        T[i] = T[i] * det;
    
    //float R[16]; mat_multiply(R, T, M); mat_print(stdout, R);
    return 0;
}


/// Invert 3x3 matrix.
int gym::mat3x3_inverse(float inv[9], const float m[9])
{
    float det = m[0] * ( m[4]*m[8] - m[5]*m[7] ) +
                m[1] * ( m[5]*m[6] - m[3]*m[8] ) +
                m[2] * ( m[3]*m[7] - m[4]*m[6] ) ;
    
    if ( det != 0 )
    {
        det = 1 / det;
        inv[0] = ( m[4]*m[8] - m[5]*m[7] ) * det;
        inv[3] = ( m[5]*m[6] - m[3]*m[8] ) * det;
        inv[6] = ( m[3]*m[7] - m[4]*m[6] ) * det;
        inv[1] = ( m[2]*m[7] - m[1]*m[8] ) * det;
        inv[4] = ( m[0]*m[8] - m[2]*m[6] ) * det;
        inv[7] = ( m[1]*m[6] - m[0]*m[7] ) * det;
        inv[2] = ( m[1]*m[5] - m[2]*m[4] ) * det;
        inv[5] = ( m[2]*m[3] - m[0]*m[5] ) * det;
        inv[8] = ( m[0]*m[4] - m[1]*m[3] ) * det;
        return 0;
    }
    return 1;
}


/*
 The vector `vec` is both input and output. Note that vec[3] should be 1
 Assuming GLint == int for 'viewport'
 */
int gym::unproject(float vec[4],
                   const float modelMatrix[16],
                   const float projMatrix[16],
                   const int viewport[4])
{
    float mat[16];
    float inv[16];
    
    gym::mat_multiply(mat, projMatrix, modelMatrix);
    
    if ( gym::mat4x4_inverse(inv, mat) )
        return 1;
    
    /* Map x and y from window coordinates */
    vec[0] = (vec[0] - viewport[0]) / viewport[2];
    vec[1] = (vec[1] - viewport[1]) / viewport[3];
    
    /* Map to range -1 to 1 */
    vec[0] = vec[0] * 2.f - 1.f;
    vec[1] = vec[1] * 2.f - 1.f;
    vec[2] = vec[2] * 2.f - 1.f;
    
    float out[4];
    gym::mat_mulvec(out, inv, vec);

    if ( out[3] == 0 )
        return 2;

    vec[0] = out[0] / out[3];
    vec[1] = out[1] / out[3];
    vec[2] = out[2] / out[3];
    vec[3] = 1;
    return 0;
}

