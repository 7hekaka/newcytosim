// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University

#include <cmath>
#include <cctype>
#include <cstdlib>
#include "assert_macro.h"
#include "gle.h"
#include "glut.h"
#include "tesselator.h"
#include "simd.h"
#include "simd_float.h"
#include "simd_math.h"
#include "gle_flute.h"

namespace gle
{
    /// values of cosine, sine over a full circle
    float cir_[2*pi_twice+8] = { 0 };
    
    /// non modifiable values of cosine, sine
    const float* circle_ = const_cast<float*>(cir_);
        
    /// vertex buffer objects for static draw
    GLuint buf_[4] = { 0 };

    /// offset for objects data stored in buffers
    GLsizei tubes_[24] = { 0 };

    /// offset for objects data stored in buffers
    GLsizei cubes_[12] = { 0 };
    
    /// offset for objects data stored in buffers
    GLsizei blobs_[4] = { 0 };
    
    /// offset for objects data stored in buffers
    GLsizei discs_[2] = { 0 };

    /// vertex buffer objects for icosahedrons
    GLsizei ico_pts_[8] = { 0 };
    GLsizei ico_idx_[8] = { 0 };
    /// number of faces in icosahedrons
    GLsizei ico_cnt_[8] = { 0 };

    void initBuffers()
    {
        glGenBuffers(2, buf_);
        glBindBuffer(GL_ARRAY_BUFFER, buf_[0]);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buf_[1]);
        setBuffers();
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }
    
    void initialize()
    {
        CHECK_GL_ERROR("before gle:initialize()");
#ifndef __APPLE__
        //need to initialize GLEW on Linux
        const GLenum err = glewInit();
        if ( GLEW_OK != err )
        {
            /* Problem: glewInit failed, something is seriously wrong. */
            fprintf(stderr, "Error: %s\n", glewGetErrorString(err));
            exit(1);
        }
#endif
        // circle starts at index 4
        compute_circle(pi_twice, cir_+4, 1, 0);
        
        if ( !glIsBuffer(buf_[0]) )
        {
            initBuffers();
            initStreams();
            CHECK_GL_ERROR("gle:initBuffers()");
            std::atexit(quit);
        }
    }
    
    void quit()
    {
        // The system will release all GPU resources, so this is not necessary
        /*
        glDeleteBuffers(2, buf_);
        for (int i=0; i<4; ++i) buf_[i] = 0;
        releaseStreams();
        */
    }
    
    //-----------------------------------------------------------------------
    #pragma mark - Compute Arc and Circle
    
    /// Calculates coordinates over an arc of circle
    /**
    Set ptr[] to coordinates around a circle:
     delta = 2 * PI / cnt;
     for i = 0 : cnt
        ptr[  2*i] = rad * cos(start+i*delta) + cX
        ptr[1+2*i] = rad * sin(start+i*delta) + cY
     ptr[] should be allocated to hold `2*cnt+2' values
    */
    void set_arc(size_t cnt, float ptr[], double rad,
                 double start, double delta, float cX, float cY)
    {
#ifdef __SSE3__
        return set_arc_SEE(cnt, ptr, rad, start, delta, cX, cY);
#else
        const double c = std::cos(delta);
        const double s = std::sin(delta);

        double t;
        double x = rad * std::cos(start);
        double y = rad * std::sin(start);
        
        for( size_t n = 0; n < cnt; ++n )
        {
            ptr[  2*n] = float(x) + cX;
            ptr[1+2*n] = float(y) + cY;
            //apply the rotation matrix
            t = x;
            x = c * x - s * y;
            y = s * t + c * y;
            //std::clog << n << " " << x << " " << y << "\n";
        }
        ptr[  2*cnt] = float(x);
        ptr[1+2*cnt] = float(y);
#endif
    }
    
    void compute_circle(size_t cnt, float ptr[], double rad, double start)
    {
        set_arc(cnt, ptr, rad, start, 2*M_PI/cnt, 0, 0);
        ptr[  2*cnt] = ptr[0];
        ptr[1+2*cnt] = ptr[1];
    }
    
    void compute_arc(size_t cnt, float ptr[], double rad, double start,
                     double angle, float cX, float cY)
    {
        set_arc(cnt, ptr, rad, start, angle/(cnt-1), cX, cY);
    }
    
    //-----------------------------------------------------------------------
    #pragma mark - Rotations
    
    /**
     This works if norm(z[]) == n!
     Derived from `Building an Orthonormal Basis, Revisited`,
     Tom Duff et al. Journal of Computer Graphics Techniques Vol. 6 N.1, 2017
     optimized by Marc B. Reynolds
    */
    void orthonormal(const float z[3], float n, float x[3], float y[3])
    {
        const float s = std::copysign(1.f, z[2]);
        const float a = z[1] / ( z[2] + s * n );
        const float b = z[1] * a;
        const float c = z[0] * a;
        x[0] = -z[2] - b;
        x[1] = c;
        x[2] = z[0];
        y[0] = s * c;
        y[1] = s * b - n;
        y[2] = s * z[1];
    }
    
    /**
     `R` is the transverse scaling done in the XY plane after rotation
     */
    void stretchAlignZ(Vector1 const& A, Vector1 const& B, float R)
    {
        float X = B.XX-A.XX;
        float Y = std::copysign(R, X);
        //warning! this matrix appears here transposed
        float mat[16] = {
            0, Y, 0, 0,
            0, 0, R, 0,
            X, 0, 0, 0,
            float(A.XX), 0, 0, 1 };
        glMultMatrixf(mat);
    }

    /**
      Translate and rotate to place A in (0,0,0) and B at (0,0,1).
      Scale XY plane by `rad' and Z axis by 1/|AB|
     `R` is the transverse scaling done in the XY plane after rotation
     */
    void stretchAlignZ(Vector2 const& A, Vector2 const& B, float R)
    {
        float X(B.XX-A.XX);
        float Y(B.YY-A.YY);
        float r = R * invsqrt(X*X+Y*Y);
        //warning! this matrix appears here transposed
        float mat[16] = {
            -r*Y, r*X, 0, 0,
            0, 0, R, 0,
            X, Y, 0, 0,
            float(A.XX), float(A.YY), 0, 1 };
        glMultMatrixf(mat);
    }
    
    /**
     Translate and rotate to place A in (0,0,0) and B at (0,0,1).
     Scale XY plane by `rad' and Z axis by 1/|AB|
     */
    void stretchAlignZ(Vector3 const& A, Vector3 const& B, float R)
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
        orthonormal(vec, R, mat, mat+4);
        glMultMatrixf(mat);
    }

    
    void rotate(Vector3 const& A, Vector3 const& B, Vector3 const& C)
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
        glMultMatrixf(mat);
    }
    
    void rotateInverse(Vector3 const& A, Vector3 const& B, Vector3 const& C)
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
        glMultMatrixf(mat);
    }

    
    void transRotate(Vector3 const& O, Vector3 const& A,
                     Vector3 const& B, Vector3 const& C)
    {
        //warning! this matrix is displayed here transposed
        float mat[16] = {
            (float)A.XX, (float)A.YY, (float)A.ZZ, 0,
            (float)B.XX, (float)B.YY, (float)B.ZZ, 0,
            (float)C.XX, (float)C.YY, (float)C.ZZ, 0,
            (float)O.XX, (float)O.YY, (float)O.ZZ, 1 };
        glMultMatrixf(mat);
    }
    
    // rotate to align Z with 'D' and translate to center 'P'
    void transAlignZ(Vector1 const& P, float R, Vector1 const& D)
    {
        float X = std::copysign(R, float(D.XX));
        float mat[16] = {
            0, -X,  0,  0,
            0,  0, -R,  0,
            X,  0,  0,  0,
            float(P.XX), 0, 0, 1};
        glMultMatrixf(mat);
    }
    
    // rotate to align Z with 'D' and translate to center 'P'
    void transAlignZ(Vector2 const& P, float R, Vector2 const& D)
    {
        float X(D.XX);
        float Y(D.YY);
        float n = R * invsqrt(X*X+Y*Y);
        X *= n;
        Y *= n;
        float mat[16] = {
            Y, -X,  0,  0,
            0,  0, -R,  0,
            X,  Y,  0,  0,
            float(P.XX), float(P.YY), 0, 1};
        glMultMatrixf(mat);
    }
    
    // rotate to align Z with 'D' and translate to center 'P', scale uniformly by `R`
    void transAlignZ(Vector3 const& P, float R, Vector3 const& D)
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
        orthonormal(vec, R, mat, mat+4);
        glMultMatrixf(mat);
    }

    // rotate to align Z with 'D', assuming norm(D)==1, and translate to center 'P'
    void stretchAlignZ1(Vector1 const& P, float R, Vector1 const& D, float S)
    {
        float X = std::copysign(R, float(D.XX));
        float Z = std::copysign(S, float(D.XX));
        float mat[16] = {
            0, -X,  0,  0,
            0,  0, -R,  0,
            Z,  0,  0,  0,
            float(P.XX), 0, 0, 1};
        glMultMatrixf(mat);
    }
    
    // rotate to align Z with 'D', assuming norm(D)==1, and translate to center 'P'
    void stretchAlignZ1(Vector2 const& P, float R, Vector2 const& D, float S)
    {
        float X(D.XX);
        float Y(D.YY);
        float mat[16] = {
            R*Y, -R*X,  0,  0,
            0,      0, -R,  0,
            S*X,  S*Y,  0,  0,
            float(P.XX), float(P.YY), 0, 1};
        glMultMatrixf(mat);
    }
    
    // rotate to align Z with 'D', assuming norm(D)==1, and translate to center 'P', scale Z axis by S
    void stretchAlignZ1(Vector3 const& P, float R, Vector3 const& D, float S)
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
        orthonormal(vec, R, mat, mat+4);
        glMultMatrixf(mat);
    }

    // rotate to align Z with X and translate to center 'P'
    void transAlignZX(float P, float R, float X)
    {
        float Y = std::copysign(R, X);
        float mat[16] = {
            0, Y, 0, 0,
            0, 0, R, 0,
            Y, 0, 0, 0,
            P, 0, 0, 1};
        glMultMatrixf(mat);
    }
    
    // rotate to align Z with X and translate to center 'A'
    void stretchAlignZX(float A, float B, float R)
    {
        float X = B - A;
        float Y = std::copysign(R, X);
        //warning! this matrix appears here transposed
        float mat[16] = {
            0, Y, 0, 0,
            0, 0, R, 0,
            X, 0, 0, 0,
            A, 0, 0, 1 };
        glMultMatrixf(mat);
    }
    
    
    void setClipPlane(GLenum glp, Vector1 const& dir, Vector1 const& pos)
    {
        GLdouble eq[4] = { dir.XX, 0, 0, -dot(dir, pos) };
        glClipPlane(glp, eq);
    }
    
    void setClipPlane(GLenum glp, Vector2 const& dir, Vector2 const& pos)
    {
        GLdouble eq[4] = { dir.XX, dir.YY, 0, -dot(dir, pos) };
        glClipPlane(glp, eq);
    }
    
    void setClipPlane(GLenum glp, Vector3 const& dir, Vector3 const& pos)
    {
        GLdouble eq[4] = { dir.XX, dir.YY, dir.ZZ, -dot(dir, pos) };
        glClipPlane(glp, eq);
    }
    
    
    /**
     draw back first, and then front of object,
     GL_CULL_FACE is temporarily enabled for this
     */
    void dualPass(void primitive())
    {
        GLboolean cull = glIsEnabled(GL_CULL_FACE);
        glEnable(GL_CULL_FACE);
        glCullFace(GL_FRONT);
        primitive();
        glCullFace(GL_BACK);
        primitive();
        if ( !cull ) glDisable(GL_CULL_FACE);
    }

    //-----------------------------------------------------------------------
    #pragma mark - Some 3D objects

    /// function callback
    using drawCall = flute6* (*)(flute6*, size_t, float const*, float const*);

    /// code to calculate and print normals
    void printNormals(size_t cnt, float const* pts)
    {
        for ( size_t i = 0; i < cnt; ++i )
        {
            float const* x = pts + 9 * i;
            Vector3 a(x[0], x[1], x[2]);
            Vector3 b(x[3], x[4], x[5]);
            Vector3 c(x[6], x[7], x[8]);
            Vector3 N = cross(b-a, c-b).normalized();
            printf("%lu ", i);
            for ( int n = 0; n < 3; ++n )
                printf("%+9.5f,%+9.5f,%+9.5f,", N.XX, N.YY, N.ZZ);
            printf("\n");
        }
    }
    
    /// set triangle strip for tube with hexagonal crosssection
    size_t setHexTube(flute6* flu, float B, float T, float rad)
    {
        // the hexagon has the same surface as a disc of radius rad.
        constexpr float C = 0.8660254037844386f; //std::sqrt(3)/2;
        constexpr float S = 0.5f;
        const float R = rad * 1.0996361107912678f; //std::sqrt( 2 * M_PI / ( 3 * std::sqrt(3) ));
        const float Y = R * C, X = R * S;
        
        flu[0] = { R, 0, T, 1, 0, 0};
        flu[1] = { R, 0, B, 1, 0, 0};
        flu[2] = { X, Y, T, S, C, 0};
        flu[3] = { X, Y, B, S, C, 0};
        flu[4] = {-X, Y, T,-S, C, 0};
        flu[5] = {-X, Y, B,-S, C, 0};
        flu[6] = {-R, 0, T,-1, 0, 0};
        flu[7] = {-R, 0, B,-1, 0, 0};
        flu[8] = {-X,-Y, T,-S, -C, 0};
        flu[9] = {-X,-Y, B,-S, -C, 0};
        flu[10] = { X,-Y, T, S,-C, 0};
        flu[11] = { X,-Y, B, S,-C, 0};
        flu[12] = { R, 0, T, 1, 0, 0};
        flu[13] = { R, 0, B, 1, 0, 0};
        return 14;
    }

    /// Tetrahedron is make of 4 triangles = 12 vertices
    size_t setTetrahedron(flute6* flt, float R=1.2f)
    {
        const float S = R / M_SQRT3;
        const float Y = 2 * S;
        const float B = -M_SQRT1_2 * S;
        const float Z = -3 * B;

        // -R,-S, B
        // +R,-S, B
        //  0, Y, B
        //  0, 0, Z
        flt[0] = { R,-S, B, 0.00000, 0.00000,-1.00000};
        flt[1] = {-R,-S, B, 0.00000, 0.00000,-1.00000};
        flt[2] = { 0, Y, B, 0.00000, 0.00000,-1.00000};
        flt[3] = { R,-S, B, 0.81650, 0.47140, 0.33333};
        flt[4] = { 0, Y, B, 0.81650, 0.47140, 0.33333};
        flt[5] = { 0, 0, Z, 0.81650, 0.47140, 0.33333};
        flt[6] = { 0, Y, B,-0.81650, 0.47140, 0.33333};
        flt[7] = {-R,-S, B,-0.81650, 0.47140, 0.33333};
        flt[8] = { 0, 0, Z,-0.81650, 0.47140, 0.33333};
        flt[9] = {-R,-S, B, 0.00000,-0.94281, 0.33333};
        flt[10] = { R,-S, B, 0.00000,-0.94281, 0.33333};
        flt[11] = { 0, 0, Z, 0.00000,-0.94281, 0.33333};
        return 12;
    }
    
    /// inversed Tetrahedrons by central symmetry
    size_t invTetrahedron(flute6* flt, float R=1.2f)
    {
        const float S = R / M_SQRT3;
        const float Y = 2 * S;
        const float B = -M_SQRT1_2 * S;
        const float Z = -3.0 * B;

        // reversed tetrahedron by central symmetry
        flt[0] = { R, S,-B, 0.00000, 0.00000, 1.00000};
        flt[1] = {-R, S,-B, 0.00000, 0.00000, 1.00000};
        flt[2] = { 0,-Y,-B, 0.00000, 0.00000, 1.00000};
        flt[3] = { 0,-Y,-B,-0.81650,-0.47140,-0.33333};
        flt[4] = {-R, S,-B,-0.81650,-0.47140,-0.33333};
        flt[5] = { 0, 0,-Z,-0.81650,-0.47140,-0.33333};
        flt[6] = { R, S,-B, 0.81650,-0.47140,-0.33333};
        flt[7] = { 0,-Y,-B, 0.81650,-0.47140,-0.33333};
        flt[8] = { 0, 0,-Z, 0.81650,-0.47140,-0.33333};
        flt[9] = {-R, S,-B, 0.00000, 0.94281,-0.33333};
        flt[10] = { R, S,-B, 0.00000, 0.94281,-0.33333};
        flt[11] = { 0, 0,-Z, 0.00000, 0.94281,-0.33333};
        return 12;
    }

    
    /// Cube is made of 12 triangles = 36 vertices
    size_t setCube(flute6* flt, float R)
    {
        flt[0] = { R, R, R, 1, 0, 0};
        flt[1] = { R,-R,-R, 1, 0, 0};
        flt[2] = { R, R,-R, 1, 0, 0};
        flt[3] = { R,-R,-R, 1, 0, 0};
        flt[4] = { R, R, R, 1, 0, 0};
        flt[5] = { R,-R, R, 1, 0, 0};
        flt[6] = { R, R, R, 0, 1, 0};
        flt[7] = { R, R,-R, 0, 1, 0};
        flt[8] = {-R, R,-R, 0, 1, 0};
        flt[9] = { R, R, R, 0, 1, 0};
        flt[10] = {-R, R,-R, 0, 1, 0};
        flt[11] = {-R, R, R, 0, 1, 0};
        flt[12] = {-R, R, R, 0, 0, 1};
        flt[13] = {-R,-R, R, 0, 0, 1};
        flt[14] = { R,-R, R, 0, 0, 1};
        flt[15] = { R, R, R, 0, 0, 1};
        flt[16] = {-R, R, R, 0, 0, 1};
        flt[17] = { R,-R, R, 0, 0, 1};
        flt[18] = {-R,-R,-R,-1, 0, 0};
        flt[19] = {-R,-R, R,-1, 0, 0};
        flt[20] = {-R, R, R,-1, 0, 0};
        flt[21] = {-R,-R,-R,-1, 0, 0};
        flt[22] = {-R, R, R,-1, 0, 0};
        flt[23] = {-R, R,-R,-1, 0, 0};
        flt[24] = { R,-R, R, 0,-1, 0};
        flt[25] = {-R,-R,-R, 0,-1, 0};
        flt[26] = { R,-R,-R, 0,-1, 0};
        flt[27] = { R,-R, R, 0,-1, 0};
        flt[28] = {-R,-R, R, 0,-1, 0};
        flt[29] = {-R,-R,-R, 0,-1, 0};
        flt[30] = { R, R,-R, 0, 0,-1};
        flt[31] = {-R,-R,-R, 0, 0,-1};
        flt[32] = {-R, R,-R, 0, 0,-1};
        flt[33] = { R, R,-R, 0, 0,-1};
        flt[34] = { R,-R,-R, 0, 0,-1};
        flt[35] = {-R,-R,-R, 0, 0,-1};
        return 36;
    }
    
    /// Octahedron is make of 8 triangles = 24 vertices
    size_t setOctahedron(flute6* flt, float R=1.46459188756f)
    {
        // the default size is set to match the volume of the unit sphere
        // Eight triangles, ordered counterclockwise
        const float N = 1 / M_SQRT3;
        flt[0] = { R, 0, 0,  N,-N, N};
        flt[1] = { 0, 0, R,  N,-N, N};
        flt[2] = { 0,-R, 0,  N,-N, N};
        flt[3] = { 0, 0,-R, -N, N,-N};
        flt[4] = {-R, 0, 0, -N, N,-N};
        flt[5] = { 0, R, 0, -N, N,-N};
        flt[6] = { 0, 0, R, -N,-N, N};
        flt[7] = {-R, 0, 0, -N,-N, N};
        flt[8] = { 0,-R, 0, -N,-N, N};
        flt[9] = { 0, 0,-R,  N, N,-N};
        flt[10] = { 0, R, 0,  N, N,-N};
        flt[11] = { R, 0, 0,  N, N,-N};
        flt[12] = {-R, 0, 0, -N, N, N};
        flt[13] = { 0, 0, R, -N, N, N};
        flt[14] = { 0, R, 0, -N, N, N};
        flt[15] = { 0, 0,-R,  N,-N,-N};
        flt[16] = { R, 0, 0,  N,-N,-N};
        flt[17] = { 0,-R, 0,  N,-N,-N};
        flt[18] = { 0, 0, R,  N, N, N};
        flt[19] = { R, 0, 0,  N, N, N};
        flt[20] = { 0, R, 0,  N, N, N};
        flt[21] = { 0, 0,-R, -N,-N,-N};
        flt[22] = { 0,-R, 0, -N,-N,-N};
        flt[23] = {-R, 0, 0, -N,-N,-N};
        return 24;
    }

    
#if ( 0 )
    // this is used to calculate the vertices of the icosahedron
    void icoFace(float* a, float* b, float* c)
    {
        float pts[9] = {a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2]};
        float nx = (a[0]+b[0]+c[0]) / 3.0f;
        float ny = (a[1]+b[1]+c[1]) / 3.0f;
        float nz = (a[2]+b[2]+c[2]) / 3.0f;
        float n = std::sqrt( nx * nx + ny * ny + nz * nz );
        nx /= n;
        ny /= n;
        nz /= n;
        float nor[9] = {nx, ny, nz, nx, ny, nz, nx, ny, nz};
        glVertexPointer(3, GL_FLOAT, 0, pts);
        glNormalPointer(GL_FLOAT, 0, nor);
        glDrawArrays(GL_TRIANGLES, 0, 3);
        if ( 0 ) {
            printf("%2.0f, %2.0f, %2.0f,  ", a[0], a[1], a[2]);
            printf("%2.0f, %2.0f, %2.0f,  ", b[0], b[1], b[2]);
            printf("%2.0f, %2.0f, %2.0f,\n", c[0], c[1], c[2]);
        }
        if ( 0 ) {
            printf("%+9.7f,%+9.7f,%+9.7f, ", nx, ny, nz);
            printf("%+9.7f,%+9.7f,%+9.7f, ", nx, ny, nz);
            printf("%+9.7f,%+9.7f,%+9.7f,\n", nx, ny, nz);
        }
    }
    
    void icoFace(float* pts, size_t a, size_t b, size_t c)
    {
        icoFace(pts+3*a, pts+3*b, pts+3*c);
    }
    
    void icosahedron()
    {
        const float G = 0.5+0.5*std::sqrt(5.0);
        const float O = 1/std::sqrt(G*G+1); //0.5257311121f;
        const float T = G * O;   //0.8506508084f;
        
        // Twelve vertices of icosahedron on unit sphere
        float pts[3*12] = {
            +T,  O,  0, // 0
            -T, -O,  0, // 1
            -T,  O,  0, // 2
            +T, -O,  0, // 3
            +O,  0,  T, // 4
            -O,  0, -T, // 5
            +O,  0, -T, // 6
            -O,  0,  T, // 7
            +0,  T,  O, // 8
            +0, -T, -O, // 9
            +0, -T,  O, // 10
            +0,  T, -O  // 11
        };
        
        glEnableClientState(GL_NORMAL_ARRAY);
        /* The faces are ordered with increasing Z */
        icoFace(pts, 5,  6, 9);
        icoFace(pts, 5, 11, 6);
        
        icoFace(pts, 6, 3,  9);
        icoFace(pts, 2, 11, 5);
        icoFace(pts, 1, 5,  9);
        icoFace(pts, 0, 6, 11);
        
        icoFace(pts, 0, 3,  6);
        icoFace(pts, 1, 2,  5);
        
        icoFace(pts, 1, 9, 10);
        icoFace(pts, 0, 11, 8);
        icoFace(pts, 8, 11, 2);
        icoFace(pts, 9, 3, 10);
        
        icoFace(pts, 0, 4,  3);
        icoFace(pts, 1, 7,  2);
        
        icoFace(pts, 0, 8,  4);
        icoFace(pts, 1, 10, 7);
        icoFace(pts, 3, 4, 10);
        icoFace(pts, 7, 8,  2);
        
        icoFace(pts, 4, 8,  7);
        icoFace(pts, 4, 7, 10);
        glDisableClientState(GL_NORMAL_ARRAY);
    }
#endif
    
    /// Icosahedrong with 20 triangles = 60 vertices
    size_t setIcosahedron(flute6* flt, float R=1.0f)
    {
        const float T = R * 0.8506508084f;      // (1 + sqrt(5))/2
        const float O = R * 0.5257311121f;      // 1 / sqrt(1+T^2)
        const float N = 1 / M_SQRT3; // 0.5773503
        const float X = 0.3568221;
        const float Y = 0.9341724;

        size_t i = 0;
        flt[i++] = {-O,  0, -T, +0,-X,-Y};
        flt[i++] = { O,  0, -T, +0,-X,-Y};
        flt[i++] = { 0, -T, -O, +0,-X,-Y};
        flt[i++] = {-O,  0, -T, +0,+X,-Y};
        flt[i++] = { 0,  T, -O, +0,+X,-Y};
        flt[i++] = { O,  0, -T, +0,+X,-Y};
        flt[i++] = { O,  0, -T, +N,-N,-N};
        flt[i++] = { T, -O,  0, +N,-N,-N};
        flt[i++] = { 0, -T, -O, +N,-N,-N};
        flt[i++] = {-T,  O,  0, -N,+N,-N};
        flt[i++] = { 0,  T, -O, -N,+N,-N};
        flt[i++] = {-O,  0, -T, -N,+N,-N};
        flt[i++] = {-T, -O,  0, -N,-N,-N};
        flt[i++] = { O,  0, -T, -N,-N,-N};
        flt[i++] = { 0, -T, -O, -N,-N,-N};
        flt[i++] = { T,  O,  0, +N,+N,-N};
        flt[i++] = { O,  0, -T, +N,+N,-N};
        flt[i++] = { 0,  T, -O, +N,+N,-N};
        flt[i++] = { T,  O,  0, +Y,+0,-X};
        flt[i++] = { T, -O,  0, +Y,+0,-X};
        flt[i++] = { O,  0, -T, +Y,+0,-X};
        flt[i++] = {-T, -O,  0, -Y,+0,-X};
        flt[i++] = { T,  O,  0, -Y,+0,-X};
        flt[i++] = {-O,  0, -T, -Y,+0,-X};
        flt[i++] = {-T, -O,  0, -X,-Y,+0};
        flt[i++] = { 0, -T, -O, -X,-Y,+0};
        flt[i++] = { 0, -T,  O, -X,-Y,+0};
        flt[i++] = { T,  O,  0, +X,+Y,+0};
        flt[i++] = { 0,  T, -O, +X,+Y,+0};
        flt[i++] = { 0,  T,  O, +X,+Y,+0};
        flt[i++] = { 0,  T,  O, -X,+Y,+0};
        flt[i++] = { 0,  T, -O, -X,+Y,+0};
        flt[i++] = {-T,  O,  0, -X,+Y,+0};
        flt[i++] = { 0, -T, -O, +X,-Y,+0};
        flt[i++] = { T, -O,  0, +X,-Y,+0};
        flt[i++] = { 0, -T,  O, +X,-Y,+0};
        flt[i++] = { T,  O,  0, +Y,+0,+X};
        flt[i++] = { O,  0,  T, +Y,+0,+X};
        flt[i++] = { T, -O,  0, +Y,+0,+X};
        flt[i++] = {-T, -O,  0, -Y,+0,+X};
        flt[i++] = { O,  0,  T, -Y,+0,+X};
        flt[i++] = {-T,  O,  0, -Y,+0,+X};
        flt[i++] = { T,  O,  0, +N,+N,+N};
        flt[i++] = { 0,  T,  O, +N,+N,+N};
        flt[i++] = { O,  0,  T, +N,+N,+N};
        flt[i++] = {-T, -O,  0, -N,-N,+N};
        flt[i++] = { 0, -T,  O, -N,-N,+N};
        flt[i++] = {-O,  0,  T, -N,-N,+N};
        flt[i++] = { T, -O,  0, +N,-N,+N};
        flt[i++] = { O,  0,  T, +N,-N,+N};
        flt[i++] = { 0, -T,  O, +N,-N,+N};
        flt[i++] = {-O,  0,  T, -N,+N,+N};
        flt[i++] = { 0,  T,  O, -N,+N,+N};
        flt[i++] = {-T,  O,  0, -N,+N,+N};
        flt[i++] = { O,  0,  T, +0,+X,+Y};
        flt[i++] = { 0,  T,  O, +0,+X,+Y};
        flt[i++] = {-O,  0,  T, +0,+X,+Y};
        flt[i++] = { O,  0,  T, +0,-X,+Y};
        flt[i++] = { O,  0,  T, +0,-X,+Y};
        flt[i++] = { 0, -T,  O, +0,-X,+Y};
        assert_true( i == 60 );
        return i;
    }
    
    /// Three fins similar to the tail of a V2 rocket
    size_t setArrowTail(flute6* flt, float R=0.1f, float B=-0.5f,
                         float H=-1.5f, float L=2.0f)
    {
        const float T = B + L;
        const float U = H + L;
        const float C = 0.5f;
        const float S = 0.5f * M_SQRT3;
        const float cR = R * C;
        const float sR = R * S;
        size_t i = 0;
        flt[i++] = { cR,-sR, B,  0, -1, 0};
        flt[i++] = {  1,  0, H,  0, -1, 0};
        flt[i++] = {  1,  0, U,  0, -1, 0};
        flt[i++] = { cR,-sR, B,  0, -1, 0};
        flt[i++] = {  1,  0, U,  0, -1, 0};
        flt[i++] = {  0,  0, T,  0, -1, 0};
        flt[i++] = { cR, sR, B,  0, +1, 0};
        flt[i++] = {  0,  0, T,  0, +1, 0};
        flt[i++] = {  1,  0, U,  0, +1, 0};
        flt[i++] = { cR, sR, B,  0, +1, 0};
        flt[i++] = {  1,  0, U,  0, +1, 0};
        flt[i++] = {  1,  0, H,  0, +1, 0};
        flt[i++] = { cR,-sR, B,  S, -C, 0};
        flt[i++] = {  0,  0, T,  S, -C, 0};
        flt[i++] = { -C, -S, U,  S, -C, 0};
        flt[i++] = { cR,-sR, B,  S, -C, 0};
        flt[i++] = { -C, -S, U,  S, -C, 0};
        flt[i++] = { -C, -S, H,  S, -C, 0};
        flt[i++] = { -R,  0, B, -S,  C, 0};
        flt[i++] = { -C, -S, H, -S,  C, 0};
        flt[i++] = { -C, -S, U, -S,  C, 0};
        flt[i++] = { -R,  0, B, -S,  C, 0};
        flt[i++] = { -C, -S, U, -S,  C, 0};
        flt[i++] = {  0,  0, T, -S,  C, 0};
        flt[i++] = { cR, sR, B,  S,  C, 0};
        flt[i++] = { -C,  S, H,  S,  C, 0};
        flt[i++] = { -C,  S, U,  S,  C, 0};
        flt[i++] = { cR, sR, B,  S,  C, 0};
        flt[i++] = { -C,  S, U,  S,  C, 0};
        flt[i++] = {  0,  0, T,  S,  C, 0};
        flt[i++] = { -R,  0, B, -S, -C, 0};
        flt[i++] = {  0,  0, T, -S, -C, 0};
        flt[i++] = { -C,  S, U, -S, -C, 0};
        flt[i++] = { -R,  0, B, -S, -C, 0};
        flt[i++] = { -C,  S, U, -S, -C, 0};
        flt[i++] = { -C,  S, H, -S, -C, 0};
        flt[i++] = { cR, sR, B,  C, -S,-1};
        flt[i++] = { -R,  0, B,  C, -S,-1};
        flt[i++] = { -C,  S, H,  C, -S,-1};
        flt[i++] = { -R,  0, B,  C,  S,-1};
        flt[i++] = { cR,-sR, B,  C,  S,-1};
        flt[i++] = { -C, -S, H,  C,  S,-1};
        flt[i++] = { cR,-sR, B, -1,  0,-1};
        flt[i++] = { cR, sR, B, -1,  0,-1};
        flt[i++] = {  1,  0, H, -1,  0,-1};
        assert_true( i == 45 );
        return i;
    }
    
    //-----------------------------------------------------------------------
    
    /// this only sets vertices, skipping normals
    size_t setCuboid(flute3* flu, float R)
    {
        const float U = -R;
        flu[0] = {R, U, U};
        flu[1] = {U, U, U};
        flu[2] = {R, R, U};
        flu[3] = {U, R, U};
        flu[4] = {U, R, R};
        flu[5] = {U, U, U};
        flu[6] = {U, U, R};
        flu[7] = {R, U, U};
        flu[8] = {R, U, R};
        flu[9] = {R, R, U};
        flu[10] = {R, R, R};
        flu[11] = {U, R, R};
        flu[12] = {R, U, R};
        flu[13] = {U, U, R};
        return 14;
    }

    size_t setBlob(flute3* flu)
    {
        constexpr GLfloat R = 1.f, U = -1.f, H(M_SQRT2);
        /* start from a centerred cube, rotated appropriately
         with vertices ordered to draw all surfaces of the cube
         with a single triangle strip. */
        const float pts[] = {
             H, R, 0, 0, R, H,
             H, U, 0, 0, U, H,
            -H, U, 0, 0, R, H,
            -H, R, 0, H, R, 0,
             0, R,-H, H, U, 0,
             0, U,-H,-H, U, 0,
             0, R,-H,-H, R, 0 };
        // refine the triangle strip:
        size_t i = 0, n = 0;
        flute3 a, b, c;
        for ( n = 0; n < 11; n += 2 )
        {
            a = {pts[3*n  ], pts[3*n+1], pts[3*n+2]};
            b = {pts[3*n+3], pts[3*n+4], pts[3*n+5]};
            c = {pts[3*n+6], pts[3*n+7], pts[3*n+8]};
            flu[i++] = normalize(a+a);
            flu[i++] = normalize(a+b);
            flu[i++] = normalize(a+c);
            flu[i++] = normalize(b+c);
        }
        a = {pts[3*n+3], pts[3*n+4], pts[3*n+5]};
        flu[i++] = normalize(c+c);
        flu[i++] = normalize(a+c);
        for ( n = 13; n > 2; n -= 2 )
        {
            a = {pts[3*n  ], pts[3*n+1], pts[3*n+2]};
            b = {pts[3*n-3], pts[3*n-2], pts[3*n-1]};
            c = {pts[3*n-6], pts[3*n-5], pts[3*n-4]};
            flu[i++] = normalize(a+a);
            flu[i++] = normalize(a+b);
            flu[i++] = normalize(a+c);
            flu[i++] = normalize(b+c);
        }
        a = {pts[0], pts[1], pts[2]};
        flu[i++] = normalize(c+c);
        flu[i++] = normalize(a+c);
        assert_true(i==52);
        return i;
    }
    
    /* This moves some vertices to smoothen the blob */
    void refineBlob(flute3 * flu, size_t cnt)
    {
    }

    /* This moves some vertices to add an hexagonal needle to the blob */
    void modifyBlob(flute3 * flu)
    {
        const float R = 0.6f;
        const float Y = R * 0.8660254037844386f; // sqrtf(3)/2;
        const float X = R * 0.5f;
        const float Z = 0.6f;
        //{ 1, 3, 5, 7, 9, 11, 40, 41, 43, 45, 47, 49, 51 }
        // pull some vertices far away in Z
        for ( int u : { 42, 44, 46, 48, 50 } ) flu[u] = { 0, 0, 6 };
        // set the 6 vertex of an hexagon:
        for ( int u : { 3, 49 } ) flu[u] = { R, 0, Z};
        for ( int u : { 5, 47 } ) flu[u] = { X,-Y, Z};
        for ( int u : { 7, 45 } ) flu[u] = {-X,-Y, Z};
        for ( int u : { 9, 43 } ) flu[u] = {-R, 0, Z};
        for ( int u : {11, 41 } ) flu[u] = {-X, Y, Z};
        for ( int u : {1,40,51} ) flu[u] = { X, Y, Z};
    }
    
    size_t setPin(flute3* flu)
    {
        size_t i = setBlob(flu);
        modifyBlob(flu);
        return i;
    }

    void thing()
    {
        flute3 flu[52];
        setBlob(flu);
        //modifyBlob(flu);
        glVertexPointer(3, GL_FLOAT, 0, flu);
        glEnableClientState(GL_NORMAL_ARRAY);
        glNormalPointer(GL_FLOAT, 0, flu);

        gle_color(1,1,1).load_both();
        glPointSize(8);
        for ( unsigned u : { 1, 3, 5, 7, 9, 11, 40, 41, 43, 45, 47, 49, 51 } )
            glDrawArrays(GL_POINTS, u, 1);
        glLineWidth(1);
        glDrawArrays(GL_LINE_STRIP, 0, 52);

        gle_color(0,1,1,0.5).load_both();
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 26);
        gle_color(1,1,0,0.5).load_both();
        glDrawArrays(GL_TRIANGLE_STRIP, 26, 26);
        glDisableClientState(GL_NORMAL_ARRAY);
    }
    
    void nothing()
    {
    }
    
    //-----------------------------------------------------------------------
    
    /// draw triangles with interleaved Vertex + Normal data
    static void drawVNTriangles(flute6* ptr, GLsizei cnt)
    {
        glEnableClientState(GL_NORMAL_ARRAY);
        glVertexPointer(3, GL_FLOAT, 6*sizeof(float), ptr);
        glNormalPointer(GL_FLOAT, 6*sizeof(float), ptr);
        glDrawArrays(GL_TRIANGLES, 0, cnt);
        glDisableClientState(GL_NORMAL_ARRAY);
    }
    
    static size_t sizeCubeBuffers()
    {
        return ( 12 + 12 + 60 + 45 + 36 + 24 + 3 * 14 );
    }
    
    flute6* setCubeBuffers(flute6* ptr, flute6* const ori, GLsizei idx[])
    {
        size_t i = 0, s = ptr - ori;
        idx[0] = i+s; i += setTetrahedron(ptr);
        idx[1] = i+s; i += invTetrahedron(ptr+i);
        idx[2] = i+s; i += setOctahedron(ptr+i);
        idx[3] = i+s; i += setIcosahedron(ptr+i);
        idx[4] = i+s; i += setArrowTail(ptr+i);
        idx[5] = i+s; i += setCube(ptr+i, 0.5773502692f);
        idx[6] = i+s; i += setHexTube(ptr+i, 0, 1, 1.0f);
        idx[7] = i+s; i += setHexTube(ptr+i, 0, 1, 0.5f);
        idx[8] = i+s; i += setHexTube(ptr+i, 0, 256.f, 0.5f);
        return ptr + i;
    }
    
    static size_t sizeBlobBuffers()
    {
        return ( 52 + 52 + 14 );
    }

    flute3* setBlobBuffers(flute3* ptr, flute3* const ori, GLsizei idx[])
    {
        size_t i = 0, s = ptr - ori;
        idx[0] = i+s; i += setBlob(ptr+i);
        idx[1] = i+s; i += setPin(ptr+i);
        idx[2] = i+s; i += setCuboid(ptr+i, 1.0);
        return ptr + i;
    }

    void drawVNBuffer(GLenum mode, GLint start, size_t cnt)
    {
        glEnableClientState(GL_NORMAL_ARRAY);
        glVertexPointer(3, GL_FLOAT, sizeof(flute6), nullptr);
        glNormalPointer(GL_FLOAT, sizeof(flute6), (void*)(3*sizeof(float)));
        glDrawArrays(mode, start, cnt);
        glDisableClientState(GL_NORMAL_ARRAY);
    }

    void drawVNTriangles(GLsizei start, GLsizei cnt)
    {
        glBindBuffer(GL_ARRAY_BUFFER, buf_[0]);
        drawVNBuffer(GL_TRIANGLES, start, cnt);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }
    
    void drawVNStrip(GLsizei start, GLsizei cnt)
    {
        glBindBuffer(GL_ARRAY_BUFFER, buf_[0]);
        drawVNBuffer(GL_TRIANGLE_STRIP, start, cnt);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

    void drawTriangleStrip(GLsizei start, GLsizei cnt)
    {
        //glEnableClientState(GL_VERTEX_ARRAY);
        glEnableClientState(GL_NORMAL_ARRAY);
        glBindBuffer(GL_ARRAY_BUFFER, buf_[0]);
        glVertexPointer(3, GL_FLOAT, 0, nullptr);
        glNormalPointer(GL_FLOAT, 0, nullptr);
        glDrawArrays(GL_TRIANGLE_STRIP, start, cnt);
        glDisableClientState(GL_NORMAL_ARRAY);
        //glPointSize(7); glDrawArrays(GL_POINTS, start, cnt);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

    void tetrahedron() { drawVNTriangles(cubes_[0], 12); }
    void octahedron()  { drawVNTriangles(cubes_[2], 24); }
    void icosahedron() { drawVNTriangles(cubes_[3], 60); }
    void ICOSAHEDRON() { flute6 tmp[60]; setIcosahedron(tmp); drawVNTriangles(tmp, 60); }
    
    void arrowTail() { drawVNTriangles(cubes_[4], 45); }
    void cube()      { drawVNTriangles(cubes_[5], 36); }
    void star()      { drawVNTriangles(cubes_[0], 24); }
    
    void hexTube()      { drawVNStrip(cubes_[6], 14); }
    void thinTube()     { drawVNStrip(cubes_[7], 14); }
    void thinLongTube() { drawVNStrip(cubes_[8], 14); }

    void blob()   { drawTriangleStrip(blobs_[0], 52); }
    void needle() { drawTriangleStrip(blobs_[1], 52); }
    void cuboid() { drawTriangleStrip(blobs_[2], 14); }

    //-----------------------------------------------------------------------
    #pragma mark - 2D Circle
    
    size_t setCircle(flute2* flu, size_t inc, float Z, float R, float N)
    {
        size_t i = 0;
        for ( size_t n = 0; n <= pi_twice; n += inc )
        {
            float C = cos_(n), S = sin_(n);
            flu[i++] = {R*C, R*S};
        }
        return i;
    }
    
    size_t sizeCircBuffers()
    {
        return 4 + 2 * pi_twice;
    }
    
    flute2* setCircBuffers(flute2* ptr, flute2* const ori, GLsizei idx[])
    {
        size_t i = 0, s = ptr - ori;
        idx[0] = i+s; i += setCircle(ptr, 1, 0, 1, 1);
        idx[1] = i+s; i += setCircle(ptr, 1, 0, 1, 2);
        return ptr + i;
    }

    //-----------------------------------------------------------------------
    #pragma mark - Tubes

    /// set triangle strip for a tube of constant radius 1 with Z in [B, T]
    size_t setTube(flute6* flu, size_t inc, float B, float T)
    {
        assert_true( B <= T );
        size_t i = 0;
        for( size_t p = 0; p <= pi_twice; p += inc )
        {
            float C = cos_(p), S = sin_(p);
            flu[i++] = { C, S, T, C, S, 0 };
            flu[i++] = { C, S, B, C, S, 0 };
        }
        return i;
    }
    
    /// set triangle strip for a tube of constant radius R with Z in [B, T]
    size_t setTube(flute6* flu, size_t inc, float B, float T, float R)
    {
        assert_true( B <= T );
        size_t i = 0;
        for( size_t p = 0; p <= pi_twice; p += inc )
        {
            float C = cos_(p), S = sin_(p);
            flu[i++] = { R*C, R*S, T, C, S, 0 };
            flu[i++] = { R*C, R*S, B, C, S, 0 };
        }
        return i;
    }
    
    /// set triangle strip for a double-walled tube of inner radius R with Z in [B, T]
    size_t setDoubleTube(flute6* flu, size_t inc, float B, float T, float R)
    {
        size_t i = setTube(flu, inc, B, T);
        i += setTube(flu+i, inc, B, T, R);
        return i;
    }

    /// set triangle strip for a cone of radius rB at Z=B, and rT at Z=T
    size_t setCone(flute6* flu, size_t inc, float B, float T, float rB, float rT)
    {
        assert_true( B <= T );
        const float H(T-B), R(rT-rB);
        const float k(std::copysign(1.f,H)/sqrtf(H*H+R*R));
        const float tC(k*H);
        const float tS(k*R);
        size_t i = 0;
        for( size_t n = 0; n <= pi_twice; n += inc )
        {
            float S = sin_(n), C = cos_(n);
            flu[i++] = { rT*C, rT*S, T, tC*C, tC*S, tS };
            flu[i++] = { rB*C, rB*S, B, tC*C, tC*S, tS };
        }
        return i;
    }

    /// set triangle strip for a disc at Z, with given normal in Z
    size_t setDisc(flute6* flu, size_t inc, float Z, float N)
    {
        size_t i = 0;
        flu[i++] = { 1, 0, Z, 0, 0, N };
        for( size_t n = inc; n < pi_once; n += inc )
        {
            float S = std::copysign(sin_(n), N), C = cos_(n);
            flu[i++] = { C,  S, Z, 0, 0, N };
            flu[i++] = { C, -S, Z, 0, 0, N };
        }
        flu[i++] = {-1, 0, Z, 0, 0, N };
        return i;
    }
    
    size_t sizeTubeBuffers()
    {
        return 22 * pi_twice;  // this is empirical!
    }
    
    flute6* setTubeBuffers(flute6* ptr, flute6* const ori, GLsizei idx[])
    {
        /* The value of T limits the aspect ratio of tubes that can be drawn */
        const float B = -32.f, T = 256.f, E = 0.03125;
        size_t i = 0, s = ptr - ori;
        idx[0] = i+s; i += setTube(ptr+i, 1, 0, 1);
        idx[1] = i+s; i += setTube(ptr+i, 2, 0, 1);
        idx[2] = i+s; i += setTube(ptr+i, 4, 0, 1);
        idx[3] = i+s; i += setTube(ptr+i, 4, 0, 1+E);
        idx[4] = i+s; i += setTube(ptr+i, 4,-E, 1+E);
        idx[5] = i+s; i += setTube(ptr+i, 4,-E, 1);
        idx[6] = i+s; i += setTube(ptr+i, 1, B, T);
        idx[8] = i+s; i += setTube(ptr+i, 2, B, T);
        idx[9] = i+s; i += setTube(ptr+i, 4, B, T);
        idx[10] = i+s; i += setTube(ptr+i, 1, 0, T);
        idx[11] = i+s; i += setTube(ptr+i, 2, 0, T);
        idx[12] = i+s; i += setTube(ptr+i, 4, 0, T);
        idx[13] = i+s; i += setCone(ptr+i, 1, 0, 1, 1, 0);
        idx[14] = i+s; i += setCone(ptr+i, 2, 0, 1, 1, 0);
        idx[15] = i+s; i += setCone(ptr+i, 2, 0, 1, 1, 0.25);
        idx[16] = i+s; i += setCone(ptr+i, 2, 0, 9, 1, 0);
        idx[17] = i+s; i += setDisc(ptr+i, 1, 0, 1);
        idx[18] = i+s; i += setDisc(ptr+i, 2, 0, 1);
        idx[19] = i+s; i += setDisc(ptr+i, 1, 1, 1);
        idx[20] = i+s; i += setDisc(ptr+i, 2, 1, 1);
        idx[21] = i+s; i += setDisc(ptr+i, 1, 0, -1);
        idx[22] = i+s; i += setDisc(ptr+i, 2, 0, -1);
        return ptr + i;
    }

    inline size_t nbTrianglesTube(size_t inc)
    {
        return 2 * ( 1 + pi_twice / inc );
    }

    inline void drawTubeStrip(GLint start, size_t cnt)
    {
        glBindBuffer(GL_ARRAY_BUFFER, buf_[0]);
        drawVNBuffer(GL_TRIANGLE_STRIP, start, cnt);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }
    
    inline void drawLineStrip(GLint start, size_t cnt, size_t skip, GLenum mode)
    {
        glBindBuffer(GL_ARRAY_BUFFER, buf_[0]);
        glVertexPointer(2, GL_FLOAT, skip*sizeof(flute2), nullptr);
        glDrawArrays(mode, start, cnt/skip);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }
    
    // using Vertex Buffer Objects
    void tube1()         { drawTubeStrip(tubes_[0], nbTrianglesTube(1)); }
    void tube2()         { drawTubeStrip(tubes_[1], nbTrianglesTube(2)); }
    void tube4()         { drawTubeStrip(tubes_[2], nbTrianglesTube(4)); }
    void tubeS()         { drawTubeStrip(tubes_[3], nbTrianglesTube(4)); }
    void tubeM()         { drawTubeStrip(tubes_[4], nbTrianglesTube(4)); }
    void tubeE()         { drawTubeStrip(tubes_[5], nbTrianglesTube(4)); }
    void longTube1()     { drawTubeStrip(tubes_[7], nbTrianglesTube(1)); }
    void longTube2()     { drawTubeStrip(tubes_[8], nbTrianglesTube(2)); }
    void longTube4()     { drawTubeStrip(tubes_[9], nbTrianglesTube(4)); }
    void halfTube1()     { drawTubeStrip(tubes_[10], nbTrianglesTube(1)); }
    void halfTube2()     { drawTubeStrip(tubes_[11], nbTrianglesTube(2)); }
    void halfTube4()     { drawTubeStrip(tubes_[12], nbTrianglesTube(4)); }
    void cone1()         { drawTubeStrip(tubes_[13], nbTrianglesTube(1)); }
    void cone2()         { drawTubeStrip(tubes_[14], nbTrianglesTube(2)); }
    void truncatedCone() { drawTubeStrip(tubes_[15], nbTrianglesTube(2)); }
    void cone3()         { drawTubeStrip(tubes_[16], nbTrianglesTube(2)); }
    
    void disc1()         { drawTubeStrip(tubes_[17], pi_twice); }
    void disc2()         { drawTubeStrip(tubes_[18], pi_twice/2); }
    void discTop1()      { drawTubeStrip(tubes_[19], pi_twice); }
    void discTop2()      { drawTubeStrip(tubes_[20], pi_twice/2); }
    void discBottom1()   { drawTubeStrip(tubes_[21], pi_twice); }
    void discBottom2()   { drawTubeStrip(tubes_[22], pi_twice/2); }
    
    void circle()        { drawLineStrip(discs_[0], 1+pi_twice, 1, GL_LINE_STRIP); }
    void circle2()       { drawLineStrip(discs_[1], 1+pi_twice, 2, GL_LINE_STRIP); }
    void circle_dotted() { drawLineStrip(discs_[1], 1+pi_twice, 1, GL_LINES); }

    void disc() { disc1(); }
    void cone() { cone2(); discBottom2(); }
    void cylinder1() { tube2(); discBottom2(); discTop2(); }
    // these primitices do not preserve the modelview transformation
    void cylinder2() { glTranslatef(0,0,-0.5f); discBottom2(); tube2(); discTop2(); }
    void longCone() { glTranslatef(0,0,-1); glScalef(1,1,3); cone2(); discBottom2(); }
    void shortCone() { glTranslatef(0,0,-0.333f); glScalef(1,1,0.5); cone2(); discBottom2(); }

    //-----------------------------------------------------------------------
#pragma mark - Spheres made from refined Icosahedrons
    
    /// using icosahedrons to render the sphere:
    void setIcoBuffer(Tesselator& ico, int i, float*& ptr, float* const ptr0, unsigned*& idx, unsigned* const idx0)
    {
        /*
         check pointer alignment, which is required to get indices
         by substracting pointer values below */
        assert_true( 0 == ( ptr - ptr0 ) % 3 );
        //fprintf(stderr, "setIcoBuffer %i: %u %u\n", i, ico.max_vertices(), ico.num_vertices());
        ico_pts_[i] = ptr - ptr0;
        ico_idx_[i] = idx - idx0;
        ico_cnt_[i] = 3 * ico.num_faces();
        
        ico.store_vertices(ptr);
        ptr += 3 * ico.num_vertices();
        
        size_t cnt = 3 * ico.num_faces();
        memcpy(idx, ico.face_data(), cnt*sizeof(unsigned));
        idx += cnt;
    }
    
    void drawIcoBuffer(GLsizei pts, GLsizei inx, GLsizei cnt)
    {
        glEnableClientState(GL_NORMAL_ARRAY);
        glBindBuffer(GL_ARRAY_BUFFER, buf_[0]);
        glVertexPointer(3, GL_FLOAT, 0, (void*)(pts*sizeof(GLfloat)));
        // the normal in each vertex is equal to the vertex!
        glNormalPointer(GL_FLOAT, 0, (void*)(pts*sizeof(GLfloat)));
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buf_[1]);
        glDrawElements(GL_TRIANGLES, cnt, GL_UNSIGNED_INT, (void*)(inx*sizeof(GLuint)));
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glDisableClientState(GL_NORMAL_ARRAY);
    }
    
    void dualPassIcoBuffer(GLsizei pts, GLsizei inx, GLsizei cnt)
    {
        assert_true(glIsEnabled(GL_CULL_FACE));
        glEnableClientState(GL_NORMAL_ARRAY);
        glBindBuffer(GL_ARRAY_BUFFER, buf_[0]);
        glVertexPointer(3, GL_FLOAT, 0, (void*)(pts*sizeof(GLfloat)));
        // the normal in each vertex is equal to the vertex!
        glNormalPointer(GL_FLOAT, 0, (void*)(pts*sizeof(GLfloat)));
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buf_[1]);
        glCullFace(GL_FRONT);
        glDrawElements(GL_TRIANGLES, cnt, GL_UNSIGNED_INT, (void*)(inx*sizeof(GLuint)));
        glCullFace(GL_BACK);
        glDrawElements(GL_TRIANGLES, cnt, GL_UNSIGNED_INT, (void*)(inx*sizeof(GLuint)));
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glDisableClientState(GL_NORMAL_ARRAY);
    }

    void sphere1()     { drawIcoBuffer(ico_pts_[0], ico_idx_[0], ico_cnt_[0]); }
    void sphere2()     { drawIcoBuffer(ico_pts_[1], ico_idx_[1], ico_cnt_[1]); }
    void sphere4()     { drawIcoBuffer(ico_pts_[2], ico_idx_[2], ico_cnt_[2]); }
    void sphere8()     { drawIcoBuffer(ico_pts_[3], ico_idx_[3], ico_cnt_[3]); }
    void hemisphere1() { drawIcoBuffer(ico_pts_[4], ico_idx_[4], ico_cnt_[4]); }
    void hemisphere2() { drawIcoBuffer(ico_pts_[5], ico_idx_[5], ico_cnt_[5]); }
    void hemisphere4() { drawIcoBuffer(ico_pts_[6], ico_idx_[6], ico_cnt_[6]); }
    
    void dualPassSphere1() { dualPassIcoBuffer(ico_pts_[0], ico_idx_[0], ico_cnt_[0]); }
    void dualPassSphere2() { dualPassIcoBuffer(ico_pts_[1], ico_idx_[1], ico_cnt_[1]); }
    void dualPassSphere4() { dualPassIcoBuffer(ico_pts_[2], ico_idx_[2], ico_cnt_[2]); }
    void dualPassSphere8() { dualPassIcoBuffer(ico_pts_[3], ico_idx_[3], ico_cnt_[3]); }
    
    void setBuffers()
    {
        Tesselator ico[7];
        ico[0].buildIcosahedron(std::max(1UL, gle::finesse/2));
        ico[1].buildIcosahedron(gle::finesse);
        ico[2].buildIcosahedron(gle::finesse*2);
        ico[3].buildIcosahedron(gle::finesse*4);
        ico[4].buildHemisphere(std::max(1UL, gle::finesse/2));
        ico[5].buildHemisphere(gle::finesse);
        ico[6].buildHemisphere(gle::finesse*2);

        size_t f = 0;
        size_t s = 0;
        for ( int i = 0; i < 7; ++i )
        {
            f += 3 * ico[i].max_faces();
            s += 3 * ico[i].max_vertices();
        }

        // required buffer size:
        size_t t = 6 * sizeTubeBuffers();
        size_t c = 6 * sizeCubeBuffers();
        size_t b = 3 * sizeBlobBuffers();
        size_t o = 2 * sizeCircBuffers();

        glBufferData(GL_ARRAY_BUFFER, (s+t+c+b+o)*sizeof(float), nullptr, GL_STATIC_DRAW);
        float* ptr = (float*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
        
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, 3*f*sizeof(unsigned), nullptr, GL_STATIC_DRAW);
        unsigned* idx = (unsigned*)glMapBuffer(GL_ELEMENT_ARRAY_BUFFER, GL_WRITE_ONLY);

        float*const ptr0 = ptr;
        unsigned*const idx0 = idx;
        
        ptr = (float*)setTubeBuffers((flute6*)ptr, (flute6*)ptr0, tubes_);
        //fprintf(stderr, "setTubeBuffers : %li %li\n", ptr-ptr0, t); float* sub=ptr;
        assert_true( ptr < ptr0 + t );

        ptr = (float*)setCubeBuffers((flute6*)ptr, (flute6*)ptr0, cubes_);
        //fprintf(stderr, "setCubeBuffer : %li %li\n", ptr-sub, c); sub=ptr;
        assert_true( ptr < ptr0 + c + t );

        for ( int i = 0; i < 7; ++i )
            setIcoBuffer(ico[i], i, ptr, ptr0, idx, idx0);
        //fprintf(stderr, "setIcosBuffers : %li %li -- %li\n", ptr-sub, s, idx-idx0); sub=ptr;
        assert_true( ptr < ptr0 + s + c + t );

        ptr = (float*)setBlobBuffers((flute3*)ptr, (flute3*)ptr0, blobs_);
        //fprintf(stderr, "setBlobBuffers : %li %li\n", ptr-sub, b); sub=ptr;
        assert_true( ptr < ptr0 + s + t + c + b );

        // align pointer:
        ptr += (ptr-ptr0) & 1;
        
        ptr = (float*)setCircBuffers((flute2*)ptr, (flute2*)ptr0, discs_);
        //fprintf(stderr, "setCircBuffers : %li %li\n", ptr-sub, o);

        assert_true( ptr < ptr0 + (s+t+c+b+o) );
        glUnmapBuffer(GL_ELEMENT_ARRAY_BUFFER);
        glUnmapBuffer(GL_ARRAY_BUFFER);
    }

    //-----------------------------------------------------------------------
    #pragma mark - 3D primitives
    
    /// set Triangle strip for a Torus of radius R and a thickness 2*T
    void torusZ(float R, float T, size_t inc)
    {
        for ( size_t n = 0; n < pi_twice; n += inc )
        {
            flute6 * flu = mapBufferV3N3(2+pi_twice);
            float X0 = cos_(n), X1 = cos_(n+inc);
            float Y0 = sin_(n), Y1 = sin_(n+inc);
            size_t i = 0;
            for ( size_t p = 0; p <= pi_twice; p += 2*inc )
            {
                float S = sin_(p), C = cos_(p);
                flu[i++] = {X0*(R+T*C), Y0*(R+T*C), T*S, X0*C, Y0*C, S};
                flu[i++] = {X1*(R+T*C), Y1*(R+T*C), T*S, X1*C, Y1*C, S};
            }
            unmapBufferV3N3();
            glDrawArrays(GL_TRIANGLE_STRIP, 0, i);
        }
        glDisableClientState(GL_NORMAL_ARRAY);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

    /**
     Draw a cylindrical band on the equator of a sphere of radius 1.
     The band is in the XY plane. The axis of the cylinder is Z.
     The band is made of triangles indicating the clockwise direction.
     */
    void arrowStrip(float width, const size_t inc)
    {
        float A(M_PI * inc / pi_half);
        float W(width * A / M_SQRT3);
        float R(1.0f / cosf(A*0.5f));
        
        flute6 * flu = mapBufferV3N3(3*(1+pi_twice/(2*inc)));
        size_t i = 0;
        flu[i++] = {R, 0, W, 1, 0, 0};
        flu[i++] = {R, 0,-W, 1, 0, 0};
        for ( size_t n = 0; n < pi_twice; n += 2*inc )
        {
            float c = R * cos_(n);
            float s = R * sin_(n);
            flu[i++] = {c, s, 0, c, s, 0};
            flu[i++] = {c, s, W, c, s, 0};
            flu[i++] = {c, s,-W, c, s, 0};
        }
        flu[i++] = {R, 0, 0, 1, 0, 0};
        unmapBufferV3N3();
        glDrawArrays(GL_TRIANGLES, 0, i);
        glDisableClientState(GL_NORMAL_ARRAY);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }
    
    
    void threeArrowStrip(float w, const size_t inc)
    {
        arrowStrip(w, inc);
        glRotated(-90,1,0,0);
        arrowStrip(w, inc);
        glRotated(90,0,1,0);
        arrowStrip(w, inc);
    }
    
    /**
     Draw a surface of revolution around the Z-axis.
     The surface goes from Z to Z_max, and its radius is
     given by the function `radius`(z) provided as argument.
     */
    void drawRevolution(float (*radius)(float), float Z, float T, float dZ)
    {
        float R = radius(Z);
        while ( Z < T )
        {
            flute6 * flu = mapBufferV3N3(2+2*pi_twice);
            float Y = Z;
            float Q = R;
            Z += dZ;
            R = radius(Z);
            
            float dR = ( R - Q ) / dZ;
            float dN = 1.0f / sqrtf( 1 + dR * dR );
            dR = dR * dN;
            
            size_t i = 0;
            for ( size_t n = 0; n <= pi_twice; ++n )
            {
                float S = sin_(n), C = cos_(n);
                flu[i++] = {R*C, R*S, Z, dN*C, dN*S,-dR};
                flu[i++] = {Q*C, Q*S, Y, dN*C, dN*S,-dR};
            }
            unmapBufferV3N3();
            glDrawArrays(GL_TRIANGLE_STRIP, 0, i);
        }
        glDisableClientState(GL_NORMAL_ARRAY);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

    /// some volume of revolution with axis along Z
    float dumbbellRadius(float z) { return sinf(3.14159265f*z) * (1.3f+cosf(6.28318530f*z)); }
    float barrelRadius(float z) { return sinf(3.14159265f*z); }
    void dumbbell() { drawRevolution(dumbbellRadius, 0, 1, 0.0625); }
    void barrel() { drawRevolution(barrelRadius, 0, 1, 0.0625); }

    //-----------------------------------------------------------------------
#pragma mark - Tubes
    
    
    void stretchTube(Vector1 const& A, Vector1 const& B, float R, void (*obj)())
    {
        glPushMatrix();
        stretchAlignZ(A, B, R);
        obj();
        glPopMatrix();
    }
    
    void stretchTube(Vector2 const& A, Vector2 const& B, float R, void (*obj)())
    {
        glPushMatrix();
        stretchAlignZ(A, B, R);
        obj();
        glPopMatrix();
    }
    
    void stretchTube(Vector3 const& A, Vector3 const& B, float R, void (*obj)())
    {
        glPushMatrix();
        stretchAlignZ(A, B, R);
        obj();
        glPopMatrix();
    }
    
    //-----------------------------------------------------------------------

    void drawTube(Vector1 const& A, float R, Vector1 const& B, void (*obj)())
    {
        glPushMatrix();
        transAlignZ(A, R, B-A);
        obj();
        glPopMatrix();
    }
    
    void drawTube(Vector2 const& A, float R, Vector2 const& B, void (*obj)())
    {
        glPushMatrix();
        transAlignZ(A, R, B-A);
        obj();
        glPopMatrix();
    }
    
    void drawTube(Vector3 const& A, float R, Vector3 const& B, void (*obj)())
    {
        glPushMatrix();
        transAlignZ(A, R, B-A);
        obj();
        glPopMatrix();
    }

    //-----------------------------------------------------------------------
    
    void drawBand(Vector2 const& A, Vector2 const& B, real rad)
    {
        Vector2 d = ( B - A ).orthogonal();
        real n = d.norm();
        if ( n > 0 )
        {
            rad /= n;
            flute2 pts[4] = { A+rad*d, A-rad*d, B+rad*d, B-rad*d };
            glVertexPointer(2, GL_FLOAT, 0, pts);
            glDrawArrays(GL_TRIANGLE_STRIP, 0, 6);
        }
    }
    
    
    void drawBand(Vector1 const& A, GLfloat rA,
                  Vector1 const& B, GLfloat rB)
    {
        GLfloat AX(A.XX);
        GLfloat BX(B.XX);
        GLfloat pts[8] = { AX, rA, AX, -rA, BX, rB, BX, -rB };
        glVertexPointer(2, GL_FLOAT, 0, pts);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
    }
    
    void drawBand(Vector2 const& A, GLfloat rA,
                  Vector2 const& B, GLfloat rB)
    {
        Vector2 d = ( B - A ).orthogonal();
        real n = d.norm();
        if ( n > 0 )
        {
            GLfloat dX(d.XX/n), dY(d.YY/n);
            GLfloat AX(A.XX), AY(A.YY);
            GLfloat BX(B.XX), BY(B.YY);
            GLfloat pts[8] = { AX+rA*dX, AY+rA*dY, AX-rA*dX, AY-rA*dY,
                               BX+rB*dX, BY+rB*dY, BX-rB*dX, BY-rB*dY };
            glVertexPointer(2, GL_FLOAT, 0, pts);
            glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
        }
    }
    
    void drawBand(Vector1 const& A, GLfloat rA, gle_color cA,
                  Vector1 const& B, GLfloat rB, gle_color cB)
    {
        GLfloat AX(A.XX);
        GLfloat BX(B.XX);
        flute6 * flu = mapBufferC4V2(4);
        flu[0] = { cA, AX, -rA };
        flu[1] = { cA, AX,  rA };
        flu[2] = { cB, BX, -rB };
        flu[3] = { cB, BX,  rB };
        unmapBufferC4V2();
        glEnableClientState(GL_COLOR_ARRAY);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
        glDisableClientState(GL_COLOR_ARRAY);
    }
    
    void drawBand(Vector2 const& A, GLfloat rA, gle_color cA,
                  Vector2 const& B, GLfloat rB, gle_color cB)
    {
        Vector2 d = ( B - A ).orthogonal();
        real n = d.norm();
        if ( n > 0 )
        {
            GLfloat dX(d.XX/n), dY(d.YY/n);
            GLfloat AX(A.XX), AY(A.YY);
            GLfloat BX(B.XX), BY(B.YY);
            flute6 * flu = mapBufferC4V2(6);
            flu[0] = { cA, AX+rA*dX, AY+rA*dY };
            flu[1] = { cA, AX-rA*dX, AY-rA*dY };
            flu[2] = { cB, BX+rB*dX, BY+rB*dY };
            flu[3] = { cB, BX-rB*dX, BY-rB*dY };
            unmapBufferC4V2();
            glEnableClientState(GL_COLOR_ARRAY);
            glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
            glDisableClientState(GL_COLOR_ARRAY);
        }
    }
    
    
    /**
     This will displays a rectangle if the connection is parallel,
     and a hourglass if the connection is antiparallel
     */
    void drawHourglass(Vector2 const& a, Vector2 const& da,
                       Vector2 const& b, Vector2 const& db)
    {
        flute2 pts[6] = { b-db, b, a-da, a+da, b, b+db };
        glVertexPointer(2, GL_FLOAT, 0, pts);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 6);
    }
    
    /**
     This will displays a rectangle if the connection is parallel,
     and a hourglass if the connection is antiparallel
     */
    void drawHourglass(Vector2 const& a, Vector2 const& da, gle_color cA,
                       Vector2 const& b, Vector2 const& db, gle_color cB)
    {
        flute6 * flu = mapBufferC4V2(6);
        flu[0] = { cB, b-db };
        flu[1] = { cB, b };
        flu[2] = { cA, a-da };
        flu[3] = { cA, a+da };
        flu[4] = { cB, b };
        flu[5] = { cB, b+db };
        unmapBufferC4V2();
        glEnableClientState(GL_COLOR_ARRAY);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 6);
        glDisableClientState(GL_COLOR_ARRAY);
    }
    
    
    /**
     This will displays a rectangle if the connection is antiparallel,
     and a hourglass if the connection is parallel
     */
    void drawCross(Vector2 const& A, Vector2 const& dA,
                   Vector2 const& B, Vector2 const& dB, real rad)
    {
        glLineWidth(0.5);
        flute2 pts[8] = { A-rad*dA, A, B, B-rad*dB,
                          A+rad*dA, A, B, B+rad*dB };
        glVertexPointer(2, GL_FLOAT, 0, pts);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
        glDrawArrays(GL_TRIANGLE_STRIP, 4, 4);
    }
    
    void drawBar(Vector3 const& A, Vector3 const& dA,
                 Vector3 const& B, Vector3 const& dB, real rad)
    {
        Vector3 ab = normalize( A - B );
        Vector3 ea = cross(ab, dA);
        Vector3 eb = cross(ab, dB);
        flute3 pts[16] = { A-rad*(dA-ea), A-rad*(dA+ea), B-rad*(dB-eb), B-rad*(dB+eb),
                           A+rad*(dA-ea), A+rad*(dA+ea), B+rad*(dB-eb), B+rad*(dB+eb),
                           A-rad*dA, A+rad*dA, B-rad*dB, B+rad*dB,
                           A-rad*dA, A+rad*dA, B-rad*dB, B+rad*dB };
        glVertexPointer(3, GL_FLOAT, 0, pts);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
        glDrawArrays(GL_TRIANGLE_STRIP, 4, 4);
        glDrawArrays(GL_TRIANGLE_STRIP, 8, 4);
        glDrawArrays(GL_TRIANGLE_STRIP,12, 4);
    }
    
    
    /**
     Two hexagons linked by a rectangle
     hexagons have the same surface as a disc of radius 1.
     */
    void drawDumbbell(Vector2 const& A, Vector2 const& B, GLfloat diameter)
    {
        const float S(1.0996361107912678f); //sqrt( 2 * M_PI / ( 3 * sqrt(3) ));
        const float R(diameter * S);
        Vector2 x = ( B - A ).normalized(R*0.8660254037844386f);
        Vector2 y = x.orthogonal(R*0.5f);
        flute2 pts[20] = {
            A-x-y, A-x+y, A-y-y, A+y+y, A+x-y, A+x+y,
            B-x-y, B-x+y, B-y-y, B+y+y, B+x-y, B+x+y };
        glVertexPointer(2, GL_FLOAT, 0, pts);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 12);
    }
    
    //-----------------------------------------------------------------------
#pragma mark - Arrows
    
    void drawCone(Vector1 const& pos, Vector1 const& dir, const GLfloat rad)
    {
        GLfloat dx = rad*dir.XX, cx = pos.XX;
        GLfloat pts[8] = {cx-dx, dx, cx-dx/2, 0, cx+dx+dx, 0, cx-dx ,-dx };
        glVertexPointer(2, GL_FLOAT, 0, pts);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
    }
    
    void drawCone(Vector2 const& pos, Vector2 const& dir, const GLfloat rad)
    {
        GLfloat dx(rad*dir.XX),  cx(pos.XX);
        GLfloat dy(rad*dir.YY),  cy(pos.YY);
        GLfloat dxy = dx + dy, dyx = dy - dx;
        GLfloat pts[8] = {cx-dxy, cy-dyx, cx-dx/2, cy-dy/2, cx+2*dx, cy+2*dy,cx+dyx, cy-dxy};
        glVertexPointer(2, GL_FLOAT, 0, pts);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
    }
    
    void drawCone(Vector3 const& pos, Vector3 const& dir, const GLfloat rad)
    {
        glPushMatrix();
        transAlignZ(pos, rad, dir);
        longCone();
        glPopMatrix();
    }
    
    //-----------------------------------------------------------------------
    
    void drawCylinder(Vector1 const& pos, Vector1 const& dir, float rad)
    {
        GLfloat cx(pos.XX);
        GLfloat dx(rad * dir.XX * 0.5);
        GLfloat pts[8] = {cx-dx, -rad, cx-dx, rad, cx+dx, -rad, cx+dx, rad };
        glVertexPointer(2, GL_FLOAT, 0, pts);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
    }
    
    void drawCylinder(Vector2 const& pos, Vector2 const& dir, float rad)
    {
        GLfloat dx(rad * dir.XX), cx(pos.XX - dx * 0.5);
        GLfloat dy(rad * dir.YY), cy(pos.YY - dy * 0.5);
        GLfloat pts[8] = {cx+dy, cy-dx, cx-dy, cy+dx, cx+dx+dy, cy+dy-dx, cx+dx-dy, cy+dy+dx};
        glVertexPointer(2, GL_FLOAT, 0, pts);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
    }
    
    void drawCylinder(Vector3 const& pos, Vector3 const& dir, float rad)
    {
        glPushMatrix();
        transAlignZ(pos, rad, dir);
        cylinder2();
        glPopMatrix();
    }
    
    
    //-----------------------------------------------------------------------
    
    void drawArrowTail(Vector1 const& pos, Vector1 const& dir, float rad)
    {
        GLfloat dx(rad * dir.XX);
        GLfloat cx(pos.XX - dx * 0.5 );
        GLfloat pts[12] = {cx-dx, -dx, cx+dx, -dx, cx, 0,
                           cx+2*dx, 0, cx-dx, dx, cx+dx, dx};
        glVertexPointer(2, GL_FLOAT, 0, pts);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 6);
    }
    
    void drawArrowTail(Vector2 const& pos, Vector2 const& dir, float rad)
    {
        GLfloat dx(rad * dir.XX);
        GLfloat dy(rad * dir.YY);
        GLfloat cx(pos.XX - 1.5f * dx);
        GLfloat cy(pos.YY - 1.5f * dy);
        GLfloat ex(cx + 2 * dx);
        GLfloat ey(cy + 2 * dy);
        GLfloat pts[12] = {cx+dy, cy-dx, ex+dy, ey-dx, cx+dx, cy+dy,
                           ex+dx, ey+dy, cx-dy, cy+dx, ex-dy, ey+dx};
        glVertexPointer(2, GL_FLOAT, 0, pts);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 6);
    }
    
    void drawArrowTail(Vector3 const& pos, Vector3 const& dir, float rad)
    {
        glPushMatrix();
        transAlignZ(pos, rad, dir);
        arrowTail();
        glPopMatrix();
    }
    
    //-----------------------------------------------------------------------
    void drawArrow(Vector1 const& A, Vector1 const& B, float R)
    {
        glPushMatrix();
        stretchAlignZ(A, B, R);
        tube1();
        glTranslatef(0, 0, 1);
        glScalef(3.0, 3.0, 3*R);
        longCone();
        glPopMatrix();
    }
    
    void drawArrow(Vector2 const& A, Vector2 const& B, float R)
    {
        glPushMatrix();
        stretchAlignZ(A, B, R);
        tube1();
        glTranslatef(0, 0, 1);
        glScalef(3.0, 3.0, 3*R);
        longCone();
        glPopMatrix();
    }
    
    void drawArrow(Vector3 const& A, Vector3 const& B, float R)
    {
        glPushMatrix();
        stretchAlignZ(A, B, R);
        tube1();
        glTranslatef(0, 0, 1);
        glScalef(3.0, 3.0, 3*R);
        longCone();
        glPopMatrix();
    }
    
    //-----------------------------------------------------------------------
#pragma mark - Text
    
    
    int fontHeight(void* font)
    {
        if ( font == GLUT_BITMAP_8_BY_13 )        return 13;
        if ( font == GLUT_BITMAP_9_BY_15 )        return 15;
        if ( font == GLUT_BITMAP_TIMES_ROMAN_10 ) return 11;
        if ( font == GLUT_BITMAP_TIMES_ROMAN_24 ) return 26;
        if ( font == GLUT_BITMAP_HELVETICA_10 )   return 11;
        if ( font == GLUT_BITMAP_HELVETICA_12 )   return 15;
        if ( font == GLUT_BITMAP_HELVETICA_18 )   return 22;
        return 13;
    }
    
    
    /**
     Compute the max width of all the lines in the given text
     This uses GLUT, which should be initialized.
     */
    int maxTextWidth(const char text[], void* font, int& lines)
    {
        int res = 0;
        lines = 0;
        int w = 0;
        for (const char* c = text; *c != '\0' ; ++c)
        {
            if ( *c == '\n' )
            {
                res = std::max(res, w);
                ++lines;
                w = 0;
            }
            else if ( isspace(*c))
            {
                w += glutBitmapWidth(font, ' ');
            }
            else if ( isprint(*c))
            {
                w += glutBitmapWidth(font, *c);
            }
        }
        res = std::max(res, w);
        if ( res > 0 )
            lines = std::max(1, lines);
        return res;
    }
    
    //-----------------------------------------------------------------------
    /**
     draw the string character per character using:
     glutBitmapCharacter()
     */
    void bitmapString(const char text[], void* font, GLfloat vshift)
    {
        if ( !font )
        {
            font = GLUT_BITMAP_HELVETICA_12;
            vshift = sign_real(vshift) * fontHeight(font);
        }
        if ( vshift == 0 )
            vshift = -fontHeight(font);
        
        GLfloat ori[4], pos[4];
        glGetFloatv(GL_CURRENT_RASTER_POSITION, ori);
        
        for (const char* p = text; *p; ++p)
        {
            if ( *p == '\n' )
            {
                glGetFloatv(GL_CURRENT_RASTER_POSITION, pos);
                glBitmap(0, 0, 0, 0, ori[0]-pos[0], vshift, nullptr);
            }
            else if ( isspace(*p) )
            {
                glutBitmapCharacter(font, ' ');
            }
            else if ( isprint(*p) )
            {
                glutBitmapCharacter(font, *p);
            }
        }
    }
    
    
    /**
     draw text at position `vec`, if this corresponds to a valid raster position
     */
    void drawText(Vector3 const& vec, const char text[], void* font, float dx)
    {
        GLboolean valid = false;
        glRasterPos3f(vec.x(), vec.y(), vec.z());
        glGetBooleanv(GL_CURRENT_RASTER_POSITION_VALID, &valid);

        if ( valid == GL_TRUE )
        {
            GLboolean depth = glIsEnabled(GL_DEPTH_TEST);
            GLboolean alpha = glIsEnabled(GL_ALPHA_TEST);
            GLboolean light = glIsEnabled(GL_LIGHTING);
            glDisable(GL_DEPTH_TEST);
            glDisable(GL_ALPHA_TEST);
            glDisable(GL_LIGHTING);
            int L = 1;
            int H = fontHeight(font);
            int W = maxTextWidth(text, font, L);
            // center text:
            glBitmap(0,0,0,0,-W*dx,-H/3,nullptr);
            bitmapString(text, font, H);
            if ( depth ) glEnable(GL_DEPTH_TEST);
            if ( alpha ) glEnable(GL_ALPHA_TEST);
            if ( light ) glEnable(GL_LIGHTING);
            CHECK_GL_ERROR("in drawText()");
        }
    }
    
    void drawText(Vector2 const& w, const char text[], void* font, float dx)
    {
        drawText(Vector3(w.XX, w.YY, 0), text, font, dx);
    }
    
    void drawText(Vector1 const& w, const char text[], void* font, float dx)
    {
        drawText(Vector3(w.XX, 0, 0), text, font, dx);
    }
    
    //-----------------------------------------------------------------------
    
    void paintOctogon(const int rec[4], const int rad)
    {
        GLfloat L(rec[0]), B(rec[1]), R(rec[2]), T(rec[3]);
        GLfloat D(rad);
        GLfloat pts[16] = {R,B+D, R,T-D, R-D,B, R-D,T, L+D,B, L+D,T, L,B+D, L,T-D};
        glVertexPointer(2, GL_FLOAT, 0, pts);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 8);
    }
    
    void drawOctogon(const int rec[4], const int rad)
    {
        GLfloat L(rec[0]), B(rec[1]), R(rec[2]), T(rec[3]);
        GLfloat D(rad);
        GLfloat pts[18] = {L,B+D, L+D,B, R-D,B, R,B+D, R,T-D, R-D,T, L+D,T, L,T-D, L,B+D};
        glVertexPointer(2, GL_FLOAT, 0, pts);
        glDrawArrays(GL_LINE_STRIP, 0, 9);
    }

    /**
     The text is displayed in the current color.
     A background rectangle is displayed only if `bcol` is visible.
     
         glColor3f(1,1,1);
         drawText(fKeyString, GLUT_BITMAP_8_BY_13, 0x0, 1);
     
     Possible values for `position`:
     - 0: bottom-left, text going up
     - 1: bottom-right, text going up
     - 2: top-right, text going down
     - 3: top-left, text going down
     - 4: center, text going down
     .
     
     Note: width and height are the current size of the viewport (window)
     */
    void drawText(const char text[], void* font, const gle_color back,
                  const int position, int width, int height)
    {
        assert_true( width > 0 );
        assert_true( height > 0 );
        
        if ( !font )
            font = GLUT_BITMAP_9_BY_15;
        
        int n_lines = 1;
        int lineHeight = fontHeight(font);
        int textWidth = maxTextWidth(text, font, n_lines);
        
        GLint px, py;
        switch( position )
        {
            case 0: {
                //bottom-left, text going up
                px = lineHeight/2;
                py = lineHeight/2;
            } break;
            case 1: {
                //bottom-right, text going up
                px = width - textWidth - lineHeight/2;
                if ( px < 0 ) px = 0;
                py = lineHeight/2;
            } break;
            case 2: {
                //top-right, text going down
                px = width - textWidth - lineHeight/2;
                if ( px < 0 ) px = 0;
                py = height - lineHeight;
                lineHeight *= -1;
            } break;
            default:
            case 3: {
                //top-left, text going down
                px = lineHeight/2;
                py = height - lineHeight;
                lineHeight *= -1;
            } break;
            case 4: {
                //center, text going down
                px = ( width - textWidth ) / 2;
                if ( px < 0 ) px = 0;
                py = ( height + n_lines*lineHeight ) / 2;
                lineHeight *= -1;
            } break;
        }
        
        //set pixel coordinate system:
        GLboolean depth = glIsEnabled(GL_DEPTH_TEST);
        GLboolean alpha = glIsEnabled(GL_ALPHA_TEST);
        GLboolean light = glIsEnabled(GL_LIGHTING);
        glDisable(GL_DEPTH_TEST);
        glDisable(GL_ALPHA_TEST);
        glDisable(GL_LIGHTING);

        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();
        glLoadIdentity();
        
        glMatrixMode(GL_PROJECTION);
        glPushMatrix();
        glLoadIdentity();
        glOrtho(0, width, 0, height, 0, 1);
        
        glRasterPos2i(0, 0);
        glBitmap(0, 0, 0, 0, px, py, nullptr);
        
        if ( back.visible() )
        {
            float col[4];
            glGetFloatv(GL_CURRENT_COLOR, col);
            int R = abs(lineHeight);
            int B = std::min(py, py + n_lines * lineHeight);
            int T = std::max(py, py + n_lines * lineHeight);
            int rec[4] = { px-R, B, px+textWidth+R, T+R+R/2+R/4 };
            back.load();
            paintOctogon(rec, 3);
            if ( position == 4 )
            {
                glLineWidth(0.5);
                drawOctogon(rec, 3);
            }
            glColor4fv(col);
        }
        
        bitmapString(text, font, lineHeight);
        
        glMatrixMode(GL_PROJECTION);
        glPopMatrix();
        glMatrixMode(GL_MODELVIEW);
        glPopMatrix();
        if ( depth ) glEnable(GL_DEPTH_TEST);
        if ( alpha ) glEnable(GL_ALPHA_TEST);
        if ( light ) glEnable(GL_LIGHTING);
    }
    
    
    //-----------------------------------------------------------------------
    #pragma mark - Misc
    
    /**
     rectangle should be specified as [ left, bottom, right, top ]
     The rectangle will be drawn counter-clockwise
     */
    void drawRectangle(const int rec[4])
    {
        GLfloat L(rec[0]), B(rec[1]), R(rec[2]), T(rec[3]);
        GLfloat pts[10] = { L, B, R, B, R, T, L, T, L, B };
        glVertexPointer(2, GL_FLOAT, 0, pts);
        glDrawArrays(GL_LINE_LOOP, 0, 5);
    }
    
    void drawRectangle(float L, float B, float R, float T, float Z)
    {
        GLfloat pts[15] = { L, B, Z, R, B, Z, R, T, Z, L, T, Z, L, B, Z };
        glVertexPointer(3, GL_FLOAT, 0, pts);
        glDrawArrays(GL_LINE_LOOP, 0, 5);
    }
    
    void fillRectangle(float L, float B, float R, float T, float Z)
    {
        GLfloat pts[15] = { L, B, Z, R, B, Z, L, T, Z, R, T, Z };
        glVertexPointer(3, GL_FLOAT, 0, pts);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
    }

    
    void drawRectangle(const int rec[4], int width, int height)
    {
        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();
        glLoadIdentity();
        
        glMatrixMode(GL_PROJECTION);
        glPushMatrix();
        glLoadIdentity();
        glOrtho(0, width, 0, height, 0, 1);
        
        //disable advanced features
        GLboolean depth = glIsEnabled(GL_DEPTH_TEST);
        GLboolean light = glIsEnabled(GL_LIGHTING);
        glDisable(GL_DEPTH_TEST);
        glDisable(GL_LIGHTING);
        
        drawRectangle(rec);
        
        if ( depth ) glEnable(GL_DEPTH_TEST);
        if ( light ) glEnable(GL_LIGHTING);
        glPopMatrix();
        glMatrixMode(GL_MODELVIEW);
        glPopMatrix();
    }
    
    
    void drawResizeBox(int width, int height)
    {
        //set the matrices
        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();
        glLoadIdentity();
        
        glMatrixMode(GL_PROJECTION);
        glPushMatrix();
        glLoadIdentity();
        glOrtho(width, 0, 0, height, 0, 1 );
        
        //draw lines at 45 degrees
        GLfloat pts[16] = {16, 1, 1, 16, 12, 1, 1, 12, 8, 1, 1, 8, 4, 1, 1, 4};
        glVertexPointer(2, GL_FLOAT, 0, pts);
        glDrawArrays(GL_LINES, 0, 8);

        glPopMatrix();
        glMatrixMode(GL_MODELVIEW);
        glPopMatrix();
    }
    

    //-----------------------------------------------------------------------
    
    int copy_parity(const int a, const int b)
    {
        return a + (( std::abs(a) + b ) & 1 );
    }

    void drawTiledFloor(int R, float T, float Z, gle_color col, gle_color back)
    {
        float H = T * 0.5;
        int Q = std::floor( double(R) * M_SQRT1_2 );
        
        if ( back.visible() )
        {
            float U = R * T;
            back.load_load();
            fillRectangle(-U, -U, U, U, Z);
        }
        
        col.load_load();
        int x = R;
        int RX = 2 * x - 3;
        int RY = 0;
        for ( int y = 0; y <= x; ++y )
        {
            /*
             using the Midpoint circle algorithm
             https://en.wikipedia.org/wiki/Midpoint_circle_algorithm
            */
            if ( RY > RX )
            {
                RX += 4 * ( x - 1 );
                --x;
            }
            RY += 4 * y + 8;
            for ( int i = copy_parity(-x,y); i <= x; i+=2 )
            {
                float X = i * T;
                float Y = y * T;
                fillRectangle( X-H, Y-H, X+H, Y+H, Z);
                fillRectangle(-X+H,-Y+H,-X-H,-Y-H, Z);
            }
            for ( int i = copy_parity(Q,y); i <= x; i+=2 )
            {
                float X = y * T;
                float Y = i * T;
                fillRectangle( X-H, Y-H, X+H, Y+H, Z);
                fillRectangle(-X+H,-Y+H,-X-H,-Y-H, Z);
            }
            for ( int i = copy_parity(Q,y); i <= x; i+=2 )
            {
                float X = y * T;
                float Y = i * T;
                fillRectangle(-X-H, Y-H,-X+H, Y+H, Z);
                fillRectangle( X+H,-Y+H, X-H,-Y-H, Z);
            }
        }
    }

    
    //-----------------------------------------------------------------------
    void drawAxes(const float size, int dim)
    {
        const GLfloat S(size);
        const GLfloat R(S * 0.1f);
        
        glMatrixMode(GL_MODELVIEW);
        
        for (int d = 0; d < dim; ++d)
        {
            glPushMatrix();
            switch(d)
            {
                case 0:
                    gle_color(1, 0, 0, 1).load_load();
                    glRotatef( 90, 0, 1, 0);
                    break;
                case 1:
                    gle_color(0, 1, 0, 1).load_load();
                    glRotatef(-90, 1, 0, 0);
                    glRotatef(180, 0, 0, 1);
                    break;
                case 2:
                    gle_color(0, 0, 1, 1).load_load();
                    glRotatef(-90, 0, 0, 1);
                    break;
            }
            glScalef(R/2, R/2, S-R);
            tube1();
            glTranslatef(0, 0, 1);
            glScalef(3, 3, R/(S-R));
            longCone();
            glPopMatrix();
        }
        // display a white ball at the origin
        gle_color(1.0, 1.0, 1.0, 1.0).load_load();
        glPushMatrix();
        scale(R);
        sphere4();
        glPopMatrix();
    }
    
    //--------------------------------------------------------------------------
    #pragma mark - Debugging
    
    char const* errorString(GLenum code)
    {
        switch ( code )
        {
            case GL_NO_ERROR:          return "GL_NO_ERROR";
            case GL_INVALID_ENUM:      return "GL_INVALID_ENUM";
            case GL_INVALID_VALUE:     return "GL_INVALID_VALUE";
            case GL_INVALID_OPERATION: return "GL_INVALID_OPERATION";
            case GL_STACK_OVERFLOW:    return "GL_STACK_OVERFLOW";
            case GL_STACK_UNDERFLOW:   return "GL_STACK_UNDERFLOW";
            case GL_OUT_OF_MEMORY:     return "GL_OUT_OF_MEMORY";
            case GL_TABLE_TOO_LARGE:   return "GL_TABLE_TOO_LARGE";
            default:                   return "GL_UNKNOWN_ERROR";
        }
    }
    
    /**
     This is similart to glutReportError,
     but the additional argument can provide useful feedback for debugging
     */
    void reportErrors(FILE * out, const char* msg)
    {
        GLenum e = glGetError();
        while ( e != GL_NO_ERROR )
        {
            fprintf(out, "OpenGL error `%s' %s\n", errorString(e), msg);
            e = glGetError();
        }
    }
    
    void print_cap(GLenum cap, const char * str)
    {
        GLint i = glIsEnabled(cap);
        std::clog << str << " " << i << "   ";
    }
    
    void dump()
    {
        GLfloat c[4] = { 0 };
        glGetFloatv(GL_CURRENT_COLOR, c);
        std::clog << "color = " << c[0] << " " << c[1] << " " << c[2] << " " << c[3] << '\n';
        
        print_cap(GL_LIGHTING, "light");
        print_cap(GL_BLEND, "blend");
        print_cap(GL_FOG, "fog");
        print_cap(GL_DEPTH_TEST, "depth");
        print_cap(GL_ALPHA_TEST, "alpha");
        print_cap(GL_STENCIL_TEST, "stencil");
        print_cap(GL_CULL_FACE, "cull");
        print_cap(GL_COLOR_LOGIC_OP, "logic");
        print_cap(GL_COLOR_ARRAY, "array");
        print_cap(GL_COLOR_MATERIAL, "material");
        print_cap(GL_LINE_STIPPLE, "stipple");
        
        std::clog << '\n';
        
#if ( 0 )
        GLint vp[4] = { 0 };
        glGetIntegerv(GL_VIEWPORT, vp);
        std::clog << "viewport = " << vp[0] << " " << vp[1] << " " << vp[2] << " " << vp[3] << '\n';
#endif
    }
    
    /// print current color properties of OpenGL context
    void print_color_materials(std::ostream& os)
    {
        GLfloat mat[4] = { 0 };
        glGetMaterialfv(GL_FRONT, GL_AMBIENT, mat);
        os << "front  amb" << gle_color::components(mat) << '\n';
        glGetMaterialfv(GL_FRONT, GL_DIFFUSE, mat);
        os << "front  dif" << gle_color::components(mat) << '\n';
        glGetMaterialfv(GL_FRONT, GL_EMISSION, mat);
        os << "front  emi" << gle_color::components(mat) << '\n';
        os << '\n';

        glGetMaterialfv(GL_BACK, GL_AMBIENT, mat);
        os << "back  amb" << gle_color::components(mat) << '\n';
        glGetMaterialfv(GL_BACK, GL_DIFFUSE, mat);
        os << "back  dif" << gle_color::components(mat) << '\n';
        glGetMaterialfv(GL_BACK, GL_EMISSION, mat);
        os << "back  emi" << gle_color::components(mat) << '\n';
        os << '\n';
    }

}

