// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University

#ifndef GYM_MATRIX_H
#define GYM_MATRIX_H


namespace gym
{

    // glOrtho()
    void mat_ortho(float M[16], float L, float R, float B, float T, float N, float F);


    // glFrustum()
    void mat_frustum(float M[16], float l, float r, float b, float t, float n, float f);

    
    // glPerspective()
    void mat_perspective(float M[16], float y_fov, float aspect, float n, float f);
    
    
    // multiply Matrix Vector
    void mat_mulvec(float out[4], const float M[16], const float in[4]);

    
    // multiply matrices
    void mat_mul(float out[16], const float A[16], const float B[16]);
    
    
    // inverse matrix
    int mat4x4_inverse(float T[16], const float M[16]);
    
    // inverse matrix
    int mat3x3_inverse(float inv[9], const float m[9]);

}


#endif
