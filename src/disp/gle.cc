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
#include "flute.h"

namespace gle
{
    /// values of cosinus, sinus over a full circle
    float cir_[2*pi_twice+8] = { 0 };
    
    /// non modifiable values of cosinus, sinus
    const float* circle_ = const_cast<float*>(cir_);
    
    /// OpenGL buffers objects for streaming
    GLuint stream_[4] = { 0 };
    
    /// vertex buffer objects for static draw
    GLuint buf_[4] = { 0 };

    /// offset for tube objects
    GLsizei start_[32] = { 0 };
    
    /// offset for other objects
    GLsizei buf_pos_[12] = { 0 };

    /// vertex buffer objects for icosahedrons
    GLsizei ico_pts_[8] = { 0 };
    GLsizei ico_idx_[8] = { 0 };
    /// number of faces in icosahedrons
    GLsizei ico_cnt_[8] = { 0 };

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
            glGenBuffers(4, buf_);
            glBindBuffer(GL_ARRAY_BUFFER, buf_[0]);
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buf_[1]);
            setSphereBuffers();
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
            glBindBuffer(GL_ARRAY_BUFFER, buf_[2]);
            setTubeBuffers();
            glBindBuffer(GL_ARRAY_BUFFER, buf_[3]);
            setBuffers();
            CHECK_GL_ERROR("gle:setBuffers()");
            glBindBuffer(GL_ARRAY_BUFFER, 0);
        }
        if ( !glIsBuffer(stream_[0]) )
            glGenBuffers(4, stream_);
        CHECK_GL_ERROR("glGenBuffers(4, stream_)");
        std::atexit(release);
    }
    
    void release()
    {
        glDeleteBuffers(4, buf_);
        buf_[0] = 0;
        glDeleteBuffers(4, stream_);
        stream_[0] = 0;
    }

    //-----------------------------------------------------------------------
    #pragma mark - maps
    
    GLuint boundBuffer()
    {
        GLint i = 0;
        glGetIntegerv(GL_ARRAY_BUFFER_BINDING, &i);
        return (GLuint)i;
    }
    
    float* mapFloatBuffer(size_t cnt)
    {
        //assert_true(glIsBuffer(stream_[1]));
        glBindBuffer(GL_ARRAY_BUFFER, stream_[1]);
        glBufferData(GL_ARRAY_BUFFER, cnt*sizeof(float), nullptr, GL_STREAM_DRAW);
        return (float*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
    }
    
    void unmapFloatBuffer(size_t N)
    {
        assert_true(stream_[1] == boundBuffer());
        //glBindBuffer(GL_ARRAY_BUFFER, stream_[1]);
        glUnmapBuffer(GL_ARRAY_BUFFER);
        glVertexPointer(N, GL_FLOAT, N*sizeof(float), nullptr);
        //glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

    flute6* mapVertexNormalBuffer(size_t cnt)
    {
        //assert_true(glIsBuffer(stream_[1]));
        glBindBuffer(GL_ARRAY_BUFFER, stream_[1]);
        glBufferData(GL_ARRAY_BUFFER, cnt*sizeof(flute6), nullptr, GL_STREAM_DRAW);
        return (flute6*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
    }
    
    void unmapVertexNormalBuffer()
    {
        assert_true(stream_[1] == boundBuffer());
        //glBindBuffer(GL_ARRAY_BUFFER, stream_[1]);
        glUnmapBuffer(GL_ARRAY_BUFFER);
        glVertexPointer((DIM>2?3:2), GL_FLOAT, sizeof(flute6), nullptr);
        glEnableClientState(GL_NORMAL_ARRAY);
        glNormalPointer(GL_FLOAT, sizeof(flute6), (void*)0xC);
        //glBindBuffer(GL_ARRAY_BUFFER, 0);
    }
    
    void bindVertexNormalBuffer(size_t skip)
    {
        glBindBuffer(GL_ARRAY_BUFFER, stream_[1]);
        glVertexPointer((DIM>2?3:2), GL_FLOAT, skip*sizeof(flute6), nullptr);
        glNormalPointer(GL_FLOAT, skip*sizeof(flute6), (void*)0xC);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
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
    void set_arc(size_t cnt, GLfloat ptr[], double rad,
                 double start, double delta, GLfloat cX, GLfloat cY)
    {
#ifdef __AVX__
        // This assumes that 'GLfloat == float'
        return set_arcAVX(cnt, ptr, rad, start, delta, cX, cY);
#endif
        const double c = std::cos(delta);
        const double s = std::sin(delta);

        double t;
        double x = rad * std::cos(start);
        double y = rad * std::sin(start);
        
        for( size_t n = 0; n < cnt; ++n )
        {
            ptr[  2*n] = GLfloat(x) + cX;
            ptr[1+2*n] = GLfloat(y) + cY;
            //apply the rotation matrix
            t = x;
            x = c * x - s * y;
            y = s * t + c * y;
            //std::clog << n << " " << x << " " << y << "\n";
        }
        ptr[  2*cnt] = GLfloat(x);
        ptr[1+2*cnt] = GLfloat(y);
    }
    
    void compute_circle(size_t cnt, GLfloat ptr[], double rad, double start)
    {
        set_arc(cnt, ptr, rad, start, 2*M_PI/cnt, 0, 0);
        ptr[  2*cnt] = ptr[0];
        ptr[1+2*cnt] = ptr[1];
    }
    
    void compute_arc(size_t cnt, GLfloat ptr[], double rad, double start,
                     double angle, GLfloat cX, GLfloat cY)
    {
        set_arc(cnt, ptr, rad, start, angle/(cnt-1), cX, cY);
    }
    
    //-----------------------------------------------------------------------
    #pragma mark - Rotation
    
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
        float X = std::copysign(R, B.XX-A.XX);
        //warning! this matrix appears here transposed
        float mat[16] = {
            0, -X,  0,  0,
            0,  0, -R,  0,
            X,  0,  0,  0,
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
        float X = float(B.XX-A.XX);
        float Y = float(B.YY-A.YY);
        float r = R * invsqrt(X*X+Y*Y);
        //warning! this matrix appears here transposed
        float mat[16] = {
            r*Y, -r*X,  0,  0,
            0,      0, -R,  0,
            X,      Y,  0,  0,
            float(A.XX), float(A.YY), 0, 1 };
        glMultMatrixf(mat);
    }
    
    /**
     Translate and rotate to place A in (0,0,0) and B at (0,0,1).
     Scale XY plane by `rad' and Z axis by 1/|AB|
     */
    void stretchAlignZ(Vector3 const& A, Vector3 const& B, float R)
    {
        float X = float(B.XX-A.XX);
        float Y = float(B.YY-A.YY);
        float Z = float(B.ZZ-A.ZZ);
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
        float X = float(D.XX);
        float Y = float(D.YY);
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
    
    // rotate to align Z with 'D' and translate to center 'P'
    void transAlignZ(Vector3 const& P, float R, Vector3 const& D)
    {
        float X = float(D.XX);
        float Y = float(D.YY);
        float Z = float(D.ZZ);
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
    void transAlignZ1(Vector1 const& P, float R, Vector1 const& D, float S)
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
    void transAlignZ1(Vector2 const& P, float R, Vector2 const& D, float S)
    {
        float X = float(D.XX);
        float Y = float(D.YY);
        float mat[16] = {
            R*Y, -R*X,  0,  0,
            0,      0, -R,  0,
            S*X,  S*Y,  0,  0,
            float(P.XX), float(P.YY), 0, 1};
        glMultMatrixf(mat);
    }
    
    // rotate to align Z with 'D', assuming norm(D)==1, and translate to center 'P', scale Z axis by S
    void transAlignZ1(Vector3 const& P, float R, Vector3 const& D, float S)
    {
        float X = float(D.XX);
        float Y = float(D.YY);
        float Z = float(D.ZZ);
        float vec[3] = { R*X, R*Y, R*Z };
        float mat[16] = {
            0, 0, 0, 0,
            0, 0, 0, 0,
            S*X, S*Y, S*Z, 0,
            float(P.XX), float(P.YY), float(P.ZZ), 1};
        orthonormal(vec, R, mat, mat+4);
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
    
    /**
     return axis orthogonal to the display plane, and corresponding to depth
     obtained from the current modelview transformation
     */
    Vector3 depthAxis()
    {
        GLfloat mat[16];
        glGetFloatv(GL_MODELVIEW_MATRIX, mat);
        return normalize(Vector3(mat[2], mat[6], mat[10]));
    }

    //-----------------------------------------------------------------------
    #pragma mark - Some 3D objects

    /// function callback
    using drawCall = void (*)(flute6*&, size_t, float const*, float const*);

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
    void setHexTube(flute6*& flu, float B, float T, float rad)
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
        flu += 14;
    }

    /// Tetrahedron is make of 4 triangles = 12 vertices
    void setTetrahedron(drawCall func, flute6*& flt, float R=1.2f)
    {
        const float S = R / M_SQRT3;
        const float Y = 2 * S;
        const float B = -M_SQRT1_2 * S;
        const float Z = -3 * B;

        // -R,-S, B
        // +R,-S, B
        //  0, Y, B
        //  0, 0, Z
        float pts[9*4] = {
             R,-S, B,-R,-S, B, 0, Y, B,
             R,-S, B, 0, Y, B, 0, 0, Z,
             0, Y, B,-R,-S, B, 0, 0, Z,
            -R,-S, B, R,-S, B, 0, 0, Z,
        };
        //printNormals(4, pts);
        float dir[9*4] = {
            +0.00000, 0.00000,-1.00000,+0.00000, 0.00000,-1.00000,+0.00000, 0.00000,-1.00000,
            +0.81650, 0.47140, 0.33333, 0.81650, 0.47140, 0.33333, 0.81650, 0.47140, 0.33333,
            -0.81650, 0.47140, 0.33333,-0.81650, 0.47140, 0.33333,-0.81650, 0.47140, 0.33333,
            +0.00000,-0.94281, 0.33333, 0.00000,-0.94281, 0.33333, 0.00000,-0.94281, 0.33333
        };
        assert_true( sizeof(pts) == sizeof(dir) );
        func(flt, 12, pts, dir);
    }
    
    /// The star is made of two Tetrahedrons: 8 triangles = 24 vertices
    void setStar(drawCall func, flute6*& flt, float R=1.2f)
    {
        const float S = R / M_SQRT3;
        const float Y = 2 * S;
        const float B = -M_SQRT1_2 * S;
        const float Z = -3.0 * B;

        float pts[9*8] = {
             R,-S, B,-R,-S, B, 0, Y, B,
             R,-S, B, 0, Y, B, 0, 0, Z,
             0, Y, B,-R,-S, B, 0, 0, Z,
            -R,-S, B, R,-S, B, 0, 0, Z,
            // reversed tetrahedron by central symmetry
             R, S,-B,-R, S,-B, 0,-Y,-B,
             0,-Y,-B,-R, S,-B, 0, 0,-Z,
             R, S,-B, 0,-Y,-B, 0, 0,-Z,
            -R, S,-B, R, S,-B, 0, 0,-Z,
        };
        
        float dir[9*8] = {
            +0.00000, 0.00000,-1.00000, 0.00000, 0.00000,-1.00000,+0.00000, 0.00000,-1.00000,
            +0.81650, 0.47140, 0.33333, 0.81650, 0.47140, 0.33333, 0.81650, 0.47140, 0.33333,
            -0.81650, 0.47140, 0.33333,-0.81650, 0.47140, 0.33333,-0.81650, 0.47140, 0.33333,
            +0.00000,-0.94281, 0.33333, 0.00000,-0.94281, 0.33333, 0.00000,-0.94281, 0.33333,
            // reversed tetrahedron by central symmetry
            -0.00000,-0.00000, 1.00000, 0.00000, 0.00000, 1.00000, 0.00000, 0.00000, 1.00000,
            -0.81650,-0.47140,-0.33333,-0.81650,-0.47140,-0.33333,-0.81650,-0.47140,-0.33333,
            +0.81650,-0.47140,-0.33333, 0.81650,-0.47140,-0.33333, 0.81650,-0.47140,-0.33333,
            -0.00000, 0.94281,-0.33333, 0.00000, 0.94281,-0.33333, 0.00000, 0.94281,-0.33333
        };
        assert_true( sizeof(pts) == sizeof(dir) );
        func(flt, 24, pts, dir);
    }

    
    /// Cube is make of 12 triangles = 36 vertices
    void setCube(drawCall func, flute6*& flt, float R=0.5773502692f)
    {
        float pts[9*12] = {
            +R, R, R, R,-R,-R, R, R,-R,
            +R,-R,-R, R, R, R, R,-R, R,
            +R, R, R, R, R,-R,-R, R,-R,
            +R, R, R,-R, R,-R,-R, R, R,
            -R, R, R,-R,-R, R, R,-R, R,
            +R, R, R,-R, R, R, R,-R, R,
            -R,-R,-R,-R,-R, R,-R, R, R,
            -R,-R,-R,-R, R, R,-R, R,-R,
            +R,-R, R,-R,-R,-R, R,-R,-R,
            +R,-R, R,-R,-R, R,-R,-R,-R,
            +R, R,-R,-R,-R,-R,-R, R,-R,
            +R, R,-R, R,-R,-R,-R,-R,-R
        };
        
        float dir[9*12] = {
            +1, 0, 0, 1, 0, 0, 1, 0, 0,
            +1, 0, 0, 1, 0, 0, 1, 0, 0,
            +0, 1, 0, 0, 1, 0, 0, 1, 0,
            +0, 1, 0, 0, 1, 0, 0, 1, 0,
            +0, 0, 1, 0, 0, 1, 0, 0, 1,
            +0, 0, 1, 0, 0, 1, 0, 0, 1,
            -1, 0, 0,-1, 0, 0,-1, 0, 0,
            -1, 0, 0,-1, 0, 0,-1, 0, 0,
            +0,-1, 0, 0,-1, 0, 0,-1, 0,
            +0,-1, 0, 0,-1, 0, 0,-1, 0,
            +0, 0,-1, 0, 0,-1, 0, 0,-1,
            +0, 0,-1, 0, 0,-1, 0, 0,-1
        };
        assert_true( sizeof(pts) == sizeof(dir) );
        func(flt, 36, pts, dir);
    }
    
    /// Octahedron is make of 8 triangles = 24 vertices
    void setOctahedron(drawCall func, flute6*& flt, float R=1.46459188756f)
    {
        // Eight triangles, ordered counterclockwise
        // set size to match the volume of the unit sphere
        float pts[9*8] = {
            +R, 0, 0, 0, 0, R, 0,-R, 0,
            +0, 0,-R,-R, 0, 0, 0, R, 0,
            +0, 0, R,-R, 0, 0, 0,-R, 0,
            +0, 0,-R, 0, R, 0, R, 0, 0,
            -R, 0, 0, 0, 0, R, 0, R, 0,
            +0, 0,-R, R, 0, 0, 0,-R, 0,
            +0, 0, R, R, 0, 0, 0, R, 0,
            +0, 0,-R, 0,-R, 0,-R, 0, 0
        };
        
        // normals
        const float N = 1 / M_SQRT3;
        float dir[9*8] = {
            +N,-N, N, N,-N, N, N,-N, N,
            -N, N,-N,-N, N,-N,-N, N,-N,
            -N,-N, N,-N,-N, N,-N,-N, N,
            +N, N,-N, N, N,-N, N, N,-N,
            -N, N, N,-N, N, N,-N, N, N,
            +N,-N,-N, N,-N,-N, N,-N,-N,
            +N, N, N, N, N, N, N, N, N,
            -N,-N,-N,-N,-N,-N,-N,-N,-N
        };
        assert_true( sizeof(pts) == sizeof(dir) );
        func(flt, 24, pts, dir);
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
    void setIcosahedron(drawCall func, flute6*& flt, float R=1.0f)
    {
        const float T = R * 0.8506508084f;      // (1 + sqrt(5))/2
        const float O = R * 0.5257311121f;      // 1 / sqrt(1+T^2)

        float pts[9*20] = {
            -O,  0, -T,   O,  0, -T,   0, -T, -O,
            -O,  0, -T,   0,  T, -O,   O,  0, -T,
            +O,  0, -T,   T, -O,  0,   0, -T, -O,
            -T,  O,  0,   0,  T, -O,  -O,  0, -T,
            -T, -O,  0,  -O,  0, -T,   0, -T, -O,
            +T,  O,  0,   O,  0, -T,   0,  T, -O,
            +T,  O,  0,   T, -O,  0,   O,  0, -T,
            -T, -O,  0,  -T,  O,  0,  -O,  0, -T,
            -T, -O,  0,   0, -T, -O,   0, -T,  O,
            +T,  O,  0,   0,  T, -O,   0,  T,  O,
            +0,  T,  O,   0,  T, -O,  -T,  O,  0,
            +0, -T, -O,   T, -O,  0,   0, -T,  O,
            +T,  O,  0,   O,  0,  T,   T, -O,  0,
            -T, -O,  0,  -O,  0,  T,  -T,  O,  0,
            +T,  O,  0,   0,  T,  O,   O,  0,  T,
            -T, -O,  0,   0, -T,  O,  -O,  0,  T,
            +T, -O,  0,   O,  0,  T,   0, -T,  O,
            -O,  0,  T,   0,  T,  O,  -T,  O,  0,
            +O,  0,  T,   0,  T,  O,  -O,  0,  T,
            +O,  0,  T,  -O,  0,  T,   0, -T,  O,
        };
        
        float dir[9*20] = {
            +0.0000000,-0.3568221,-0.9341724, +0.0000000,-0.3568221,-0.9341724, +0.0000000,-0.3568221,-0.9341724,
            +0.0000000,+0.3568221,-0.9341724, +0.0000000,+0.3568221,-0.9341724, +0.0000000,+0.3568221,-0.9341724,
            +0.5773503,-0.5773503,-0.5773503, +0.5773503,-0.5773503,-0.5773503, +0.5773503,-0.5773503,-0.5773503,
            -0.5773503,+0.5773503,-0.5773503, -0.5773503,+0.5773503,-0.5773503, -0.5773503,+0.5773503,-0.5773503,
            -0.5773503,-0.5773503,-0.5773503, -0.5773503,-0.5773503,-0.5773503, -0.5773503,-0.5773503,-0.5773503,
            +0.5773503,+0.5773503,-0.5773503, +0.5773503,+0.5773503,-0.5773503, +0.5773503,+0.5773503,-0.5773503,
            +0.9341724,+0.0000000,-0.3568221, +0.9341724,+0.0000000,-0.3568221, +0.9341724,+0.0000000,-0.3568221,
            -0.9341724,+0.0000000,-0.3568221, -0.9341724,+0.0000000,-0.3568221, -0.9341724,+0.0000000,-0.3568221,
            -0.3568221,-0.9341724,+0.0000000, -0.3568221,-0.9341724,+0.0000000, -0.3568221,-0.9341724,+0.0000000,
            +0.3568221,+0.9341724,+0.0000000, +0.3568221,+0.9341724,+0.0000000, +0.3568221,+0.9341724,+0.0000000,
            -0.3568221,+0.9341724,+0.0000000, -0.3568221,+0.9341724,+0.0000000, -0.3568221,+0.9341724,+0.0000000,
            +0.3568221,-0.9341724,+0.0000000, +0.3568221,-0.9341724,+0.0000000, +0.3568221,-0.9341724,+0.0000000,
            +0.9341724,+0.0000000,+0.3568221, +0.9341724,+0.0000000,+0.3568221, +0.9341724,+0.0000000,+0.3568221,
            -0.9341724,+0.0000000,+0.3568221, -0.9341724,+0.0000000,+0.3568221, -0.9341724,+0.0000000,+0.3568221,
            +0.5773503,+0.5773503,+0.5773503, +0.5773503,+0.5773503,+0.5773503, +0.5773503,+0.5773503,+0.5773503,
            -0.5773503,-0.5773503,+0.5773503, -0.5773503,-0.5773503,+0.5773503, -0.5773503,-0.5773503,+0.5773503,
            +0.5773503,-0.5773503,+0.5773503, +0.5773503,-0.5773503,+0.5773503, +0.5773503,-0.5773503,+0.5773503,
            -0.5773503,+0.5773503,+0.5773503, -0.5773503,+0.5773503,+0.5773503, -0.5773503,+0.5773503,+0.5773503,
            +0.0000000,+0.3568221,+0.9341724, +0.0000000,+0.3568221,+0.9341724, +0.0000000,+0.3568221,+0.9341724,
            +0.0000000,-0.3568221,+0.9341724, +0.0000000,-0.3568221,+0.9341724, +0.0000000,-0.3568221,+0.9341724,
        };
        assert_true( sizeof(pts) == sizeof(dir) );
        func(flt, 60, pts, dir);
    }
    
    /// Three fins similar to the tail of a V2 rocket
    void setArrowTail(drawCall func, flute6*& flt, float R=0.1f,
                       float B=-0.5f, float H=-1.5f, float L=2.0f)
    {
        const float T = B + L;
        const float U = H + L;
        const float C = 0.5f;
        const float S = 0.5f * M_SQRT3;
        const float cR = R * C;
        const float sR = R * S;

        float pts[9*15] = {
            cR,-sR, B,  1,  0, H,  1,  0, U,
            cR,-sR, B,  1,  0, U,  0,  0, T,
            cR, sR, B,  0,  0, T,  1,  0, U,
            cR, sR, B,  1,  0, U,  1,  0, H,
            cR,-sR, B,  0,  0, T, -C, -S, U,
            cR,-sR, B, -C, -S, U, -C, -S, H,
            -R,  0, B, -C, -S, H, -C, -S, U,
            -R,  0, B, -C, -S, U,  0,  0, T,
            cR, sR, B, -C,  S, H, -C,  S, U,
            cR, sR, B, -C,  S, U,  0,  0, T,
            -R,  0, B,  0,  0, T, -C,  S, U,
            -R,  0, B, -C,  S, U, -C,  S, H,
            cR, sR, B, -R,  0, B, -C,  S, H,
            -R,  0, B, cR,-sR, B, -C, -S, H,
            cR,-sR, B, cR, sR, B,  1,  0, H
        };
        
        float dir[9*15] = {
            0, -1, 0,  0, -1, 0,  0, -1, 0,
            0, -1, 0,  0, -1, 0,  0, -1, 0,
            0, +1, 0,  0, +1, 0,  0, +1, 0,
            0, +1, 0,  0, +1, 0,  0, +1, 0,
            S, -C, 0,  S, -C, 0,  S, -C, 0,
            S, -C, 0,  S, -C, 0,  S, -C, 0,
           -S,  C, 0, -S,  C, 0, -S,  C, 0,
           -S,  C, 0, -S,  C, 0, -S,  C, 0,
            S,  C, 0,  S,  C, 0,  S,  C, 0,
            S,  C, 0,  S,  C, 0,  S,  C, 0,
           -S, -C, 0, -S, -C, 0, -S, -C, 0,
           -S, -C, 0, -S, -C, 0, -S, -C, 0,
            C, -S,-1,  C, -S,-1,  C, -S,-1,
            C,  S,-1,  C,  S,-1,  C,  S,-1,
            -1, 0,-1, -1,  0,-1, -1,  0,-1
        };
        assert_true( sizeof(pts) == sizeof(dir) );
        func(flt, 45, pts, dir);
    }
    
    //-----------------------------------------------------------------------

    void cube1()
    {
        constexpr GLfloat R = 1.f, U = -1.f;
        const GLfloat pts[] = {
             R, U, U, U, U, U,
             R, R, U, U, R, U,
             U, R, R, U, U, U,
             U, U, R, R, U, U,
             R, U, R, R, R, U,
             R, R, R, U, R, R,
             R, U, R, U, U, R };
        glVertexPointer(3, GL_FLOAT, 0, pts);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 14);
    }
    
    size_t setCube(flute3*& flu, float R)
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
        flu += 14;
        return 14;
    }

    /// this does not really work
    void tetrahedron1()
    {
        constexpr GLfloat R = 1.f, U = -1.f;
        GLfloat pts[3*8] = {
            0, 0, R, R, 0, 0, 0, R, 0,
            0, 0, U, U, 0, 0, 0, U, 0,
            0, 0, R, R, 0, 0 };
        glVertexPointer(3, GL_FLOAT, 0, pts);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 8);
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
        return i;
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

    void setBlob(flute3*& flu, bool modified)
    {
        size_t n = setBlob(flu);
        assert_true(n==52);
        if ( modified )
            modifyBlob(flu);
        flu += n;
    }
    
    void loadBuffer(GLuint buf)
    {
        assert_true(glIsBuffer(buf));
        glBindBuffer(GL_ARRAY_BUFFER, buf);
        glVertexPointer(3, GL_FLOAT, 0, nullptr);
        glEnableClientState(GL_NORMAL_ARRAY);
        glNormalPointer(GL_FLOAT, 0, nullptr);
    }
    void unloadBuffer()
    {
        glDisableClientState(GL_NORMAL_ARRAY);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

    void pinS()
    {
        flute3 flu[52];
        size_t i = setBlob(flu);
        assert_true(i == 52);
        modifyBlob(flu);
        glVertexPointer(3, GL_FLOAT, 0, flu);
        glEnableClientState(GL_NORMAL_ARRAY);
        glNormalPointer(GL_FLOAT, 0, flu);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 26);
        glColor4f(1,1,1,0.5);
        glDrawArrays(GL_TRIANGLE_STRIP, 26, 26);
        glDisableClientState(GL_NORMAL_ARRAY);
        glColor3f(1,1,1);
        for ( unsigned u : { 1, 3, 5, 7, 9, 11, 40, 41, 43, 45, 47, 49, 51 } )
            glDrawArrays(GL_POINTS, u, 1);
        glDrawArrays(GL_LINE_STRIP, 0, i);
    }

    //-----------------------------------------------------------------------
    
    /// draw
    void drawVNBuffer(flute6*&, size_t cnt, float const* pts, float const* dir)
    {
        glEnableClientState(GL_NORMAL_ARRAY);
        glVertexPointer(3, GL_FLOAT, 0, pts);
        glNormalPointer(GL_FLOAT, 0, dir);
        glDrawArrays(GL_TRIANGLES, 0, cnt);
        glDisableClientState(GL_NORMAL_ARRAY);
    }

    /// load data to GPU
    void setBuffer(flute6*& dst, size_t cnt, float const* pts, float const* dir)
    {
        // intertwine the data:
        for ( size_t i = 0; i < cnt; ++i )
            dst[i] = { pts[3*i], pts[1+3*i], pts[2+3*i], dir[3*i], dir[1+3*i], dir[2+3*i] };
        dst += cnt;
    }
    
    void setBuffers()
    {
        size_t cnt = 24 * pi_twice;  // this is empirical!
        glBufferData(GL_ARRAY_BUFFER, cnt*sizeof(flute6), nullptr, GL_STATIC_DRAW);
        flute6 * flu = (flute6*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
        flute6 * ptr = flu;
        buf_pos_[0] = ptr-flu; setTetrahedron(setBuffer, ptr);
        buf_pos_[1] = ptr-flu; setOctahedron(setBuffer, ptr);
        buf_pos_[2] = ptr-flu; setIcosahedron(setBuffer, ptr);
        buf_pos_[3] = ptr-flu; setArrowTail(setBuffer, ptr);
        buf_pos_[4] = ptr-flu; setCube(setBuffer, ptr);
        buf_pos_[5] = ptr-flu; setStar(setBuffer, ptr);
        buf_pos_[6] = ptr-flu; setHexTube(ptr, 0, 1, 1.0f);
        buf_pos_[7] = ptr-flu; setHexTube(ptr, 0, 1, 0.5f);
        buf_pos_[8] = ptr-flu; setHexTube(ptr, 0, 256.f, 0.5f);
        flute3 * fl3 = (flute3*)ptr;
        buf_pos_[9] = fl3-(flute3*)flu; setBlob(fl3, false);
        buf_pos_[10] = fl3-(flute3*)flu; setBlob(fl3, true);
        buf_pos_[11] = fl3-(flute3*)flu; setCube(fl3, true);
        glUnmapBuffer(GL_ARRAY_BUFFER);
    }
    
    inline void drawVNBuffer(GLenum mode, GLint start, size_t cnt)
    {
        glEnableClientState(GL_NORMAL_ARRAY);
        glVertexPointer(3, GL_FLOAT, sizeof(flute6), nullptr);
        glNormalPointer(GL_FLOAT, sizeof(flute6), (void*)(3*sizeof(float)));
        glDrawArrays(mode, start, cnt);
        glDisableClientState(GL_NORMAL_ARRAY);
    }

    void drawVNTriangles(GLsizei start, GLsizei cnt)
    {
        glBindBuffer(GL_ARRAY_BUFFER, buf_[3]);
        drawVNBuffer(GL_TRIANGLES, start, cnt);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }
    
    void drawVNStrip(GLsizei start, GLsizei cnt)
    {
        glBindBuffer(GL_ARRAY_BUFFER, buf_[3]);
        drawVNBuffer(GL_TRIANGLE_STRIP, start, cnt);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

    void drawTriangleStrip(GLsizei start, GLsizei cnt)
    {
        glBindBuffer(GL_ARRAY_BUFFER, buf_[3]);
        glVertexPointer(3, GL_FLOAT, 0, nullptr);
        glEnableClientState(GL_NORMAL_ARRAY);
        glNormalPointer(GL_FLOAT, 0, nullptr);
        glDrawArrays(GL_TRIANGLE_STRIP, start, cnt);
        glDisableClientState(GL_NORMAL_ARRAY);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

    void tetrahedron() { drawVNTriangles(buf_pos_[0], 12); }
    void octahedron()  { drawVNTriangles(buf_pos_[1], 24); }
    void icosahedron() { drawVNTriangles(buf_pos_[2], 60); }
    void icosahedron1() { flute6* ptr; setIcosahedron(drawVNBuffer, ptr); }
    
    void arrowTail() { drawVNTriangles(buf_pos_[3], 45); }
    void cube()      { drawVNTriangles(buf_pos_[4], 36); }
    void star()      { drawVNTriangles(buf_pos_[5], 24); }
    
    void hexTube()      { drawVNStrip(buf_pos_[6], 14); }
    void thinTube()     { drawVNStrip(buf_pos_[7], 14); }
    void thinLongTube() { drawVNStrip(buf_pos_[8], 14); }

    void blob()      { drawTriangleStrip(buf_pos_[9], 52); }
    void needle()    { drawTriangleStrip(buf_pos_[10], 52); }
    void smallCube() { drawTriangleStrip(buf_pos_[11], 14); }

    void loadBlobBuffer() { loadBuffer(buf_[3]); }
    void blobf()     { glDrawArrays(GL_TRIANGLE_STRIP, 0, 52); }
    void needlef()   { glDrawArrays(GL_TRIANGLE_STRIP,64, 52); }

    //-----------------------------------------------------------------------
    #pragma mark - 2D circle

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
    
    void setCircle(flute6*& flu, size_t inc, float Z, float R, float N)
    {
        size_t i = setCircle((flute2*)flu, inc, Z, R, N);
        flu += 1 + i / 3;
    }
    
    /// set triangle strip for a tube of constant radius 1 with Z in [B, T]
    void setTube(flute6*& flu, size_t inc, float B, float T)
    {
        assert_true( B <= T );
        size_t i = 0;
        for( size_t p = 0; p <= pi_twice; p += inc )
        {
            float C = cos_(p), S = sin_(p);
            flu[i++] = { C, S, T, C, S, 0 };
            flu[i++] = { C, S, B, C, S, 0 };
        }
        flu += i;
    }
    
    /// set triangle strip for a tube of constant radius R with Z in [B, T]
    void setTube(flute6*& flu, size_t inc, float B, float T, float R)
    {
        assert_true( B <= T );
        size_t i = 0;
        for( size_t p = 0; p <= pi_twice; p += inc )
        {
            float C = cos_(p), S = sin_(p);
            flu[i++] = { R*C, R*S, T, C, S, 0 };
            flu[i++] = { R*C, R*S, B, C, S, 0 };
        }
        flu += i;
    }
    
    /// set triangle strip for a double-walled tube of inner radius R with Z in [B, T]
    void setDoubleTube(flute6*& flu, size_t inc, float B, float T, float R)
    {
        setTube(flu, inc, B, T);
        setTube(flu, inc, B, T, R);
    }

    /// set triangle strip for a cone of radius rB at Z=B, and rT at Z=T
    void setCone(flute6*& flu, size_t inc, float B, float T, float rB, float rT)
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
        flu += i;
    }

    /// set triangle strip for a disc at Z, with given normal in Z
    void setDisc(flute6*& flu, size_t inc, float Z, float N)
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
        flu += i;
    }
    
    size_t setTubeBuffers(flute6* flu)
    {
        flute6* ptr = flu;
        /* The value of T limits the aspect ratio of tubes that can be drawn */
        const float B = -32.f, T = 256.f, E = 0.03125;
        start_[0] = ptr-flu; setTube(ptr, 1, 0, 1);
        start_[1] = ptr-flu; setTube(ptr, 2, 0, 1);
        start_[2] = ptr-flu; setTube(ptr, 4, 0, 1);
        start_[3] = ptr-flu; setTube(ptr, 4,-E, 1+E);
        start_[4] = ptr-flu; setTube(ptr, 4, 0, 1+E);
        //start_[5] = ptr-flu; setDoubleTube(ptr, 4,-E, 1+E, 18/25.0);
        //start_[6] = ptr-flu; setDoubleTube(ptr, 4, 0, 1+E, 18/25.0);
        start_[7] = ptr-flu; setTube(ptr, 1, B, T);
        start_[8] = ptr-flu; setTube(ptr, 2, B, T);
        start_[9] = ptr-flu; setTube(ptr, 4, B, T);
        start_[10] = ptr-flu; setTube(ptr, 1, 0, T);
        start_[11] = ptr-flu; setTube(ptr, 2, 0, T);
        start_[12] = ptr-flu; setTube(ptr, 4, 0, T);
        start_[13] = ptr-flu; setCone(ptr, 1, 0, 1, 1, 0);
        start_[14] = ptr-flu; setCone(ptr, 2, 0, 1, 1, 0);
        start_[15] = ptr-flu; setCone(ptr, 2, 0, 1, 1, 0.25);
        start_[16] = ptr-flu; setDisc(ptr, 1, 0, 1);
        start_[17] = ptr-flu; setDisc(ptr, 2, 0, 1);
        start_[18] = ptr-flu; setDisc(ptr, 1, 1, 1);
        start_[19] = ptr-flu; setDisc(ptr, 2, 1, 1);
        start_[20] = ptr-flu; setDisc(ptr, 1, 0, -1);
        start_[21] = ptr-flu; setDisc(ptr, 2, 0, -1);
        start_[22] = ptr-flu; setCircle(ptr, 1, 0, 1, 1);
        start_[23] = ptr-flu; setCircle(ptr, 1, 0, 1, 2);
        return ptr-flu;
    }
    
    inline size_t nbTrianglesTube(size_t inc)
    {
        return 2 * ( 1 + pi_twice / inc );
    }

    void setTubeBuffers()
    {
        size_t cnt = 24 * pi_twice;  // this is empirical!
        glBufferData(GL_ARRAY_BUFFER, cnt*sizeof(flute6), nullptr, GL_STATIC_DRAW);
        flute6 * flu = (flute6*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
        size_t n = setTubeBuffers(flu);
        assert_true( n <= cnt );
        glUnmapBuffer(GL_ARRAY_BUFFER);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

    inline void drawTubeStrip(GLint start, size_t cnt)
    {
        glBindBuffer(GL_ARRAY_BUFFER, buf_[2]);
        drawVNBuffer(GL_TRIANGLE_STRIP, start, cnt);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }
    
    inline void drawLineStrip(GLint start, size_t cnt, size_t skip)
    {
        glBindBuffer(GL_ARRAY_BUFFER, buf_[2]);
        glVertexPointer(2, GL_FLOAT, skip*sizeof(flute2), nullptr);
        glDrawArrays(GL_LINE_STRIP, 3*start, cnt/skip);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }
    
    // using Vertex Buffer Objects
    void tube1()         { drawTubeStrip(start_[0], nbTrianglesTube(1)); }
    void tube2()         { drawTubeStrip(start_[1], nbTrianglesTube(2)); }
    void tube4()         { drawTubeStrip(start_[2], nbTrianglesTube(4)); }
    void tubeE()         { drawTubeStrip(start_[3], nbTrianglesTube(4)); }
    void tubeF()         { drawTubeStrip(start_[4], nbTrianglesTube(4)); }
    void doubleTubeE()   { drawTubeStrip(start_[5], 2*nbTrianglesTube(4)); }
    void doubleTubeF()   { drawTubeStrip(start_[6], 2*nbTrianglesTube(4)); }
    void longTube1()     { drawTubeStrip(start_[7], nbTrianglesTube(1)); }
    void longTube2()     { drawTubeStrip(start_[8], nbTrianglesTube(2)); }
    void longTube4()     { drawTubeStrip(start_[9], nbTrianglesTube(4)); }
    void halfTube1()     { drawTubeStrip(start_[10], nbTrianglesTube(1)); }
    void halfTube2()     { drawTubeStrip(start_[11], nbTrianglesTube(2)); }
    void halfTube4()     { drawTubeStrip(start_[12], nbTrianglesTube(4)); }
    void cone1()         { drawTubeStrip(start_[13], nbTrianglesTube(1)); }
    void cone2()         { drawTubeStrip(start_[14], nbTrianglesTube(2)); }
    void truncatedCone() { drawTubeStrip(start_[15], nbTrianglesTube(2)); }
    void disc1()         { drawTubeStrip(start_[16], pi_twice); }
    void disc2()         { drawTubeStrip(start_[17], pi_twice/2); }
    void discTop1()      { drawTubeStrip(start_[18], pi_twice); }
    void discTop2()      { drawTubeStrip(start_[19], pi_twice/2); }
    void discBottom1()   { drawTubeStrip(start_[20], pi_twice); }
    void discBottom2()   { drawTubeStrip(start_[21], pi_twice/2); }
    void circle()        { drawLineStrip(start_[22], 1+pi_twice, 1); }
    void circle2()       { drawLineStrip(start_[23], 1+pi_twice, 2); }

    void disc() { disc1(); }
    void cone() { cone2(); discBottom2(); }
    void cylinder1() { tube2(); discBottom2(); discTop2(); }
    // these primitices do not preserve the modelview transformation
    void cylinder2() { glTranslatef(0,0,-0.5f); discBottom2(); tube2(); discTop2(); }
    void longCone() { glTranslatef(0,0,-1); glScalef(1,1,3); cone2(); discBottom2(); }
    void shortCone() { glTranslatef(0,0,-0.333f); glScalef(1,1,0.5); cone2(); discBottom2(); }

    //-----------------------------------------------------------------------
#pragma mark - Spheres
    
    void addSphereData(size_t i, Tesselator& ico, float* PTR, float*& ptr, unsigned* IDX, unsigned*& idx)
    {
        ico_pts_[i] = ptr - PTR;
        ico_idx_[i] = idx - IDX;
        ico_cnt_[i] = ico.nb_faces() * 3;

        ico.store_vertices(ptr);
        ptr += ico.nb_vertices() * 3;
        
        size_t cnt = ico.nb_faces() * 3;
        memcpy(idx, ico.face_data(), cnt * sizeof(unsigned));
        idx += cnt;
        //fprintf(stderr, "icosahedron %lu has %u faces\n", i, ico.nb_faces());
    }

    /// using icosahedrons to render the sphere:
    void setSphereBuffers()
    {
        Tesselator ico[8];
        ico[0].build(Tesselator::ICOSAHEDRON, std::max(1UL, gle::finesse/2));
        ico[1].build(Tesselator::ICOSAHEDRON, gle::finesse);
        ico[2].build(Tesselator::ICOSAHEDRON, gle::finesse*2);
        ico[3].build(Tesselator::ICOSAHEDRON, gle::finesse*4);
        ico[4].build(Tesselator::HEMISPHERE, std::max(1UL, gle::finesse/2));
        ico[5].build(Tesselator::HEMISPHERE, gle::finesse);
        ico[6].build(Tesselator::HEMISPHERE, gle::finesse*2);
        
        //std::clog << "initializeIco ico " << ico.nb_faces() << '\n';
        size_t n_pts = 0, n_idx = 0;
        for ( int i = 0; i < 7; ++i )
        {
            n_pts += ico[i].nb_vertices();
            n_idx += ico[i].nb_faces();
        }
        
        glBufferData(GL_ARRAY_BUFFER, n_pts*3*sizeof(float), nullptr, GL_STATIC_DRAW);
        float *const PTR = (float*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);

        glBufferData(GL_ELEMENT_ARRAY_BUFFER, n_idx*3*sizeof(unsigned), nullptr, GL_STATIC_DRAW);
        unsigned *const IDX = (unsigned*)glMapBuffer(GL_ELEMENT_ARRAY_BUFFER, GL_WRITE_ONLY);
        
        float* ptr = PTR;
        unsigned* idx = IDX;
        
        for ( int i = 0; i < 7; ++i )
            addSphereData(i, ico[i], PTR, ptr, IDX, idx);

        //fprintf(stderr, "icosahedron point buffer: %li\n", ptr-PTR);
        //fprintf(stderr, "icosahedron index buffer: %li\n", idx-IDX);

        glUnmapBuffer(GL_ELEMENT_ARRAY_BUFFER);
        glUnmapBuffer(GL_ARRAY_BUFFER);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

    
    void drawIcoBuffer(GLsizei start, GLsizei start_inx, GLsizei cnt)
    {
        glEnableClientState(GL_NORMAL_ARRAY);
        glBindBuffer(GL_ARRAY_BUFFER, buf_[0]);
        glVertexPointer(3, GL_FLOAT, 0, (void*)(start*sizeof(GLfloat)));
        // the normal in each vertex is equal to the vertex!
        glNormalPointer(GL_FLOAT, 0, (void*)(start*sizeof(GLfloat)));
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buf_[1]);
        glDrawElements(GL_TRIANGLES, cnt, GL_UNSIGNED_INT, (void*)(start_inx*sizeof(GLuint)));
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
    
    void dualPassSphere1()
    {
        assert_true(glIsEnabled(GL_CULL_FACE));
        glCullFace(GL_FRONT);
        sphere1();
        glCullFace(GL_BACK);
        sphere1();
    }
    void dualPassSphere2()
    {
        assert_true(glIsEnabled(GL_CULL_FACE));
        glCullFace(GL_FRONT);
        sphere2();
        glCullFace(GL_BACK);
        sphere2();
    }
    void dualPassSphere4()
    {
        assert_true(glIsEnabled(GL_CULL_FACE));
        glCullFace(GL_FRONT);
        sphere4();
        glCullFace(GL_BACK);
        sphere4();
    }
    void dualPassSphere8()
    {
        assert_true(glIsEnabled(GL_CULL_FACE));
        glCullFace(GL_FRONT);
        sphere8();
        glCullFace(GL_BACK);
        sphere8();
    }
    
    //-----------------------------------------------------------------------
    #pragma mark - 3D primitives
    
    /// set Triangle strip for a Torus of radius R and a thickness 2*T
    void torusZ(float R, float T, size_t inc)
    {
        for ( size_t n = 0; n < pi_twice; n += inc )
        {
            flute6 * flu = mapVertexNormalBuffer(2+pi_twice);
            float X0 = cos_(n), X1 = cos_(n+inc);
            float Y0 = sin_(n), Y1 = sin_(n+inc);
            size_t i = 0;
            for ( size_t p = 0; p <= pi_twice; p += 2*inc )
            {
                float S = sin_(p), C = cos_(p);
                flu[i++] = {X0*(R+T*C), Y0*(R+T*C), T*S, X0*C, Y0*C, S};
                flu[i++] = {X1*(R+T*C), Y1*(R+T*C), T*S, X1*C, Y1*C, S};
            }
            unmapVertexNormalBuffer();
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
        
        flute6 * flu = mapVertexNormalBuffer(3*(1+pi_twice/(2*inc)));
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
        unmapVertexNormalBuffer();
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
            flute6 * flu = mapVertexNormalBuffer(2+2*pi_twice);
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
            unmapVertexNormalBuffer();
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
        GLfloat pts[8] = { AX,-rA, AX, rA, BX, -rB, BX, rB };
        flute4 col[4] = { cA, cA, cB, cB };
        glEnableClientState(GL_COLOR_ARRAY);
        glColorPointer(4, GL_FLOAT, 0, col);
        glVertexPointer(2, GL_FLOAT, 0, pts);
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
            GLfloat pts[8] = { AX+rA*dX, AY+rA*dY, AX-rA*dX, AY-rA*dY,
                               BX+rB*dX, BY+rB*dY, BX-rB*dX, BY-rB*dY };
            flute4 col[] = { cA, cA, cB, cB };
            glEnableClientState(GL_COLOR_ARRAY);
            glColorPointer(4, GL_FLOAT, 0, col);
            glVertexPointer(2, GL_FLOAT, 0, pts);
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
    void drawHourglass(Vector2 const& a, Vector2 const& da, gle_color ca,
                       Vector2 const& b, Vector2 const& db, gle_color cb)
    {
        flute2 pts[6] = { b-db, b, a-da, a+da, b, b+db };
        flute4 col[6] = { cb, cb, ca, ca, cb, cb };
        glEnableClientState(GL_COLOR_ARRAY);
        glVertexPointer(2, GL_FLOAT, 0, pts);
        glColorPointer(4, GL_FLOAT, 0, col);
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
    
    int copyparity(int a, int b)
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
            for ( int i = copyparity(-x,y); i <= x; i+=2 )
            {
                float X = i * T;
                float Y = y * T;
                fillRectangle( X-H, Y-H, X+H, Y+H, Z);
                fillRectangle(-X+H,-Y+H,-X-H,-Y-H, Z);
            }
            for ( int i = copyparity(Q,y); i <= x; i+=2 )
            {
                float X = y * T;
                float Y = i * T;
                fillRectangle( X-H, Y-H, X+H, Y+H, Z);
                fillRectangle(-X+H,-Y+H,-X-H,-Y-H, Z);
            }
            for ( int i = copyparity(Q,y); i <= x; i+=2 )
            {
                float X = y * T;
                float Y = i * T;
                fillRectangle(-X-H, Y-H,-X+H, Y+H, Z);
                fillRectangle( X+H,-Y+H, X-H,-Y-H, Z);
            }
        }
    }

    
    //-----------------------------------------------------------------------
    void drawAxes(const real size, int dim)
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

