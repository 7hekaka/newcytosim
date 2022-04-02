// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University



/*
 We use a different name to avoid a possible name collision with GLUT
 */
int myUnproject(float winx, float winy, float winz,
                const float modelMatrix[16],
                const float projMatrix[16],
                const GLint viewport[4],
                float XYZ[4])
{
    float mat[16];
    float inv[16];
    float in[4];
    float out[4];
    
    gym::mat_mul(mat, projMatrix, modelMatrix);
    
    if ( gym::mat4x4_inverse(inv, mat) )
        return 1;
    
    in[0] = winx;
    in[1] = winy;
    in[2] = winz;
    in[3] = 1;
    
    /* Map x and y from window coordinates */
    in[0] = (in[0] - viewport[0]) / viewport[2];
    in[1] = (in[1] - viewport[1]) / viewport[3];
    
    /* Map to range -1 to 1 */
    in[0] = in[0] * 2 - 1;
    in[1] = in[1] * 2 - 1;
    in[2] = in[2] * 2 - 1;
    
    gym::mat_mulvec(out, inv, in);

    if ( out[3] == 0 )
        return 2;

    XYZ[0] = out[0] / out[3];
    XYZ[1] = out[1] / out[3];
    XYZ[2] = out[2] / out[3];
    return 0;
}

