// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

#ifndef GYM_MATRIX_H
#define GYM_MATRIX_H


namespace gym
{
    // FILE <- M
    //void mat_print(FILE*, const float M[16]);

    /// M <- zero
    inline void mat_zero(float M[16])
    {
        for ( int d = 0; d < 16; ++d )
            M[d] = 0.f;
    }
    
    // M <- Identity * diag
    inline void mat_diagonal(float M[16], float diag, float depth = 0)
    {
        mat_zero(M);
        M[ 0] = diag;
        M[ 5] = diag;
        M[10] = diag;
        M[14] = depth;
        M[15] = 1.f;
    }

    /// M <- Q
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
    
    /// gluPickMatrix()
    void mat_pick(float [16], float X, float Y, float W, float H, const int[4]);

    
    /// multiply Matrix and Vector
    void mat_mulvec(float[4], const float[16], const float[4]);

    /// multiply matrices
    void mat_multiply(float[16], const float[16], const float[16]);
    
    /// multiply matrices: out = out x B */
    void mat_multiply(float[16], const float[16]);

    /// scale matrix, like glScale()
    void mat_scale(float[16], float S);

    /// scale matrix, like glScale()
    void mat_scale(float[16], float X, float Y, float Z);
    
    /// translate matrix, like glTranslate()
    void mat_translate(float[16], float X, float Y, float Z);
    
    /// copy and translate matrix
    void mat_copytrans(float[16], const float[16], float X, float Y, float Z);
    
    /// translate matrix according to second view
    void mat_translate(float[16], const float[16], float X, float Y, float Z);

    /// translate matrix and then scale, like glTranslate() followed by glScale()
    void mat_transscale(float[16], float X, float Y, float Z, float S);

    /// translate matrix and then scale, like glTranslate() followed by glScale()
    void mat_transscale(float[16], const float[16], float X, float Y, float Z, float S);

    /// keep translation and scale but not rotation component
    void mat_unrotate(float[16], const float[16]);

    /// inverse matrix
    int mat4x4_inverse(float[16], const float[16]);
    
    /// inverse matrix
    int mat3x3_inverse(float[9], const float[9]);

    /// unproject a 3D point from window coordinates to world coordinates
    int unproject(float vec[4], const float modelMatrix[16],
                  const float projMatrix[16], const int viewport[4]);

}


#endif
