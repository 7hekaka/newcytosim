// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

#ifndef GYM_MATRIX_H
#define GYM_MATRIX_H


namespace gym
{
    inline void mat_zero(float M[16])
    {
        for ( int d = 0; d < 16; ++d )
            M[d] = 0.f;
    }
    
    inline void mat_diagonal(float M[16], float diag)
    {
        mat_zero(M);
        M[ 0] = diag;
        M[ 5] = diag;
        M[10] = diag;
        M[15] = diag;
    }

    inline void mat_copy(float M[16], const float Q[16])
    {
        for ( int d = 0; d < 16; ++d )
            M[d] = Q[d];
    }
    
    /// rotation with axis (X, Y, Z) by angle defined by Cos, Sin
    void mat_rotation(float M[16], float X, float Y, float Z, float C, float S);
    
    /// rotate matrix around axis X, by angle defined by Cos, Sin
    void mat_rotateX(float M[16], float C, float S);
    
    /// rotate matrix around axis Y, by angle defined by Cos, Sin
    void mat_rotateY(float M[16], float C, float S);
    
    /// rotate matrix around axis Z, by angle defined by Cos, Sin
    void mat_rotateZ(float M[16], float C, float S);

    /// glOrtho()
    void mat_ortho(float[16], float L, float R, float B, float T, float N, float F);

    /// glFrustum()
    void mat_frustum(float[16], float l, float r, float b, float t, float n, float f);

    /// glPerspective()
    void mat_perspective(float [16], float y_fov, float aspect, float n, float f);
    
    
    /// multiply Matrix Vector
    void mat_mulvec(float[4], const float[16], const float[4]);

    /// multiply matrices
    void mat_mul(float[16], const float[16], const float[16]);
    
    /// multiply matrices
    void mat_mul(float[16], const float[16]);

    /// scale matrix, like glScale()
    void mat_scale(float[16], float X, float Y, float Z);
    
    /// translate matrix, like glTranslate()
    void mat_translate(float[16], float X, float Y, float Z);
    
    /// translate matrix and then scale, like glTranslate() followed by glScale()
    void mat_transscale(float[16], float X, float Y, float Z, float S);

    
    /// inverse matrix
    int mat4x4_inverse(float[16], const float[16]);
    
    /// inverse matrix
    int mat3x3_inverse(float[9], const float[9]);

}


#endif
