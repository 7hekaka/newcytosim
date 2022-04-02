// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University

#include "gym_matrix.h"
#include <stdio.h>
#include <math.h>

/*
 Some of this is derived from GLFW linmath.h by Camilla Löwy <elmindreda@elmindreda.org>
 */


namespace gym
{
    static void mat_zero(float M[16])
    {
        for ( int d = 0; d < 16; ++d )
            M[d] = 0.f;
    }
    
    void mat_print(FILE * f, float M[16])
    {
        for ( int i = 0; i < 4; ++i )
        {
            for ( int j = 0; j < 4; ++j )
                fprintf(f, "%5.2f ", M[i+4*j]);
            fprintf(f, "\n");
        }
    }
    
    /**
     This set a matrix like glOrtho()
     */
    void mat_ortho(float M[16], float L, float R, float B, float T, float N, float F)
    {
        mat_zero(M);
        
        M[ 0] =  2.f / ( R - L );
        M[ 5] =  2.f / ( T - B );
        M[10] = -2.f / ( F - N );
        
        M[12] = -( R + L ) / ( R - L );
        M[13] = -( T + B ) / ( T - B );
        M[14] = -( F + N ) / ( F - N );
        M[15] = 1.f;
    }


    void mat_frustum(float M[16], float l, float r, float b, float t, float n, float f)
    {
        mat_zero(M);

        M[ 0] = 2.f*n/(r-l);
        M[ 5] = 2.f*n/(t-b);
        M[ 8] = (r+l)/(r-l);
        M[ 9] = (t+b)/(t-b);
        M[10] = -(f+n)/(f-n);
        M[11] = -1.f;
        M[14] = -2.f*(f*n)/(f-n);
    }

    void mat_perspective(float M[16], float y_fov, float aspect, float n, float f)
    {
        float const a = 1.f / tanf(y_fov / 2.f);

        mat_zero(M);

        M[ 0] = a / aspect;
        M[ 5] = a;
        M[10] = -((f + n) / (f - n));
        M[11] = -1.f;
        M[14] = -((2.f * f * n) / (f - n));
    }
    


    void mat_mulvec(float out[4], const float M[16], const float in[4])
    {
        for ( int i = 0; i < 4; ++i )
            out[i] = in[0] * M[i] + in[1] * M[4+i] + in[2] * M[8+i] + in[3] * M[12+i];
    }


    void mat_mul(float out[16], const float A[16], const float B[16])
    {
        for ( int i = 0; i < 4; ++i )
            mat_mulvec(out+4*i, A, B+4*i);
    }
    

    int mat4x4_inverse(float T[16], const float M[16])
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
        
        float det = s[0]*c[5]-s[1]*c[4]+s[2]*c[3]+s[3]*c[2]-s[4]*c[1]+s[5]*c[0];
        
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
        
        //float R[16]; mat_mul(R, T, M); mat_print(stdout, R);
        return 0;
    }


    // Invert 3x3 matrix.
    int mat3x3_inverse(float inv[9], const float m[9])
    {
        float det = m[0]*m[4]*m[8] + m[2]*m[3]*m[7] + m[1]*m[5]*m[6]
                  - m[2]*m[4]*m[6] - m[1]*m[3]*m[8] - m[0]*m[5]*m[7];
        
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

    
}
