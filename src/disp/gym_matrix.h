// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

#ifndef GYM_MATRIX_H
#define GYM_MATRIX_H


namespace gym
{

    // glOrtho()
    void mat_ortho(float[16], float L, float R, float B, float T, float N, float F);

    // glFrustum()
    void mat_frustum(float[16], float l, float r, float b, float t, float n, float f);

    // glPerspective()
    void mat_perspective(float [16], float y_fov, float aspect, float n, float f);
    
    
    // multiply Matrix Vector
    void mat_mulvec(float[4], const float[16], const float[4]);

    // multiply matrices
    void mat_mul(float[16], const float[16], const float[16]);

    // scale matrix, like glScale()
    void mat_scale(float[16], float X, float Y, float Z);
    
    // translate matrix, like glTranslate()
    void mat_translate(float[16], float X, float Y, float Z);
    
    // translate matrix and then scale, like glTranslate() followed by glScale()
    void mat_transscale(float[16], float X, float Y, float Z, float S);

    
    // inverse matrix
    int mat4x4_inverse(float[16], const float[16]);
    
    // inverse matrix
    int mat3x3_inverse(float[9], const float[9]);

}


#endif
