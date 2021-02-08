// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University

#include <cmath>
#include <cctype>
#include <cstdlib>
#include "assert_macro.h"
#include "gle.h"
#include "glut.h"
#include "platonic.h"
#include "simd.h"
#include "simd_float.h"
#include "simd_math.h"
#include "vector_float.h"

namespace gle
{
    /// values of cosinus, sinus over a full circle
    GLfloat cir_[2*ncircle+8] = { 0 };
    
    inline GLfloat cos_(size_t n) { return cir_[4+2*n]; }
    inline GLfloat sin_(size_t n) { return cir_[5+2*n]; }

    /// vertex buffer objects for tubes
    GLuint tub_[24] = { 0 };
    
    /// vertex buffer objects for hex tubes
    GLuint buf_[16] = { 0 };

    /// vertex buffer objects for icosahedrons
    GLuint ico_[16] = { 0 };
    
    /// number of faces in icosahedrons
    GLuint ico_nfaces[8] = { 0 };
    
    /// inverse square root
#if defined(__SSE3__)
    inline float invsqrt(float x) { return _mm_cvtss_f32(_mm_rsqrt_ss(_mm_set_ss(x))); }
#else
    inline float invsqrt(float x) { return 1.0f / sqrtf(x); }
#endif

    /// Calculates cosine and sinus over an arc of circle
    /**
     CS[] should be preallocated to hold `2x(cnt+1)' values.
     Fill in with a counter-clockwise circle starting at angle `start`.
     with radius `rad` and center {cenX, cenY}
    */
    void compute_arc(size_t cnt, GLfloat CS[], double rad,
                     double start, double end, GLfloat cX, GLfloat cY)
    {
        const double theta = ( end - start ) / (double)cnt;
        const double c = std::cos(theta);
        const double s = std::sin(theta);
    
        double t;
        double x = rad * std::cos(start);
        double y = rad * std::sin(start);
        
        for( size_t n = 0; n <= cnt; ++n )
        {
            CS[  2*n] = GLfloat(x) + cX;
            CS[1+2*n] = GLfloat(y) + cY;
            //apply the rotation matrix
            t = x;
            x = c * x - s * y;
            y = s * t + c * y;
            //std::clog << n << " " << x << " " << y << "\n";
        }
    }

    /// Calculates cosine and sinus over the entire circle
    /**
    Set CS[] to coordinates around a circle:
     delta = 2 * PI / cnt;
     for i = 0 : cnt
        CS[  2*i] = rad * cos(start+i*delta)
        CS[1+2*i] = rad * sin(start+i*delta)
     CS should be allocated to hold `2*cnt+2' values
    */
    void compute_circle(size_t cnt, GLfloat CS[], double rad, double start)
    {
#ifdef __AVX__
        // This assumes that 'GLfloat == float'
        return circleAVX(cnt, CS, rad, start);
#endif
        const double theta = 2 * M_PI / (double)cnt;
        const double c = std::cos(theta);
        const double s = std::sin(theta);

        double t;
        double x = rad * std::cos(start);
        double y = rad * std::sin(start);
        GLfloat x0(x), y0(y);
        
        for( size_t n = 0; n < cnt; ++n )
        {
            CS[  2*n] = GLfloat(x);
            CS[1+2*n] = GLfloat(y);
            //apply the rotation matrix
            t = x;
            x = c * x - s * y;
            y = s * t + c * y;
            //std::clog << n << " " << x << " " << y << "\n";
        }
        CS[  2*cnt] = x0;
        CS[1+2*cnt] = y0;
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
        // store { 0, 0 } at the start of circle
        for ( size_t i = 0; i < 4; ++i )
            cir_[i] = 0;
        // circle starts at index 4
        compute_circle(ncircle, cir_+4, 1);
        
        initSphereBuffers();
        CHECK_GL_ERROR("gle:initSphereBuffers()");
        initTubeBuffers();
        CHECK_GL_ERROR("gle:initTubeBuffers()");
        initBuffers();
        CHECK_GL_ERROR("gle:initBuffers()");
        std::atexit(release);
    }
    
    void release()
    {
        glDeleteBuffers(16, ico_);
        ico_[0] = 0;
        glDeleteBuffers(24, tub_);
        tub_[0] = 0;
        glDeleteBuffers(16, buf_);
        buf_[0] = 0;
    }
    
    //-----------------------------------------------------------------------
    #pragma mark - Rotation

    /** extracts axis orthogonal to the display plane, and corresponding to depth
     from the current modelview transformation: */
    Vector3 directionDepth()
    {
        GLfloat mat[16];
        glGetFloatv(GL_MODELVIEW_MATRIX, mat);
        return normalize(Vector3(mat[2], mat[6], mat[10]));
    }
    
    template < typename FLOAT >
    inline void orthonormal(const FLOAT v[3], FLOAT mul, FLOAT x[3], FLOAT y[3])
    {
        const FLOAT s = std::copysign((FLOAT)1, v[2]);
        /// optimized version by Marc B. Reynolds
        const FLOAT a = v[1] / ( v[2] + s );
        const FLOAT b = v[1] * a;
        const FLOAT c = v[0] * a;
        FLOAT ss = s * mul;
        x[0] = mul * ( -v[2] - b );
        x[1] = mul * c;
        x[2] = mul * v[0];
        y[0] = ss * c;
        y[1] = ss * b - mul;
        y[2] = ss * v[1];
    }

    template < typename FLOAT >
    inline void orthonormal(const FLOAT v[3], FLOAT mul, FLOAT x[3], FLOAT y[3], FLOAT z[3])
    {
        const FLOAT s = std::copysign((FLOAT)1, v[2]);
        /// optimized version by Marc B. Reynolds
        const FLOAT a = v[1] / ( v[2] + s );
        const FLOAT b = v[1] * a;
        const FLOAT c = v[0] * a;
        FLOAT ss = s * mul;
        x[0] = mul * ( -v[2] - b );
        x[1] = mul * c;
        x[2] = mul * v[0];
        y[0] = ss * c;
        y[1] = ss * b - mul;
        y[2] = ss * v[1];
        z[0] = mul * v[0];
        z[1] = mul * v[1];
        z[2] = mul * v[2];
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
    
    // set rotation to align Z with 'AB' and translate to 'A'
    void stretchAlignZ(Vector3 const& A, Vector3 const& B, float R)
    {
        float X = float(B.XX-A.XX);
        float Y = float(B.YY-A.YY);
        float Z = float(B.ZZ-A.ZZ);
        float N = invsqrt(X*X+Y*Y+Z*Z);
        float vec[3] = { N*X, N*Y, N*Z };
        float mat[16] = {
            0, 0, 0, 0,
            0, 0, 0, 0,
            X, Y, Z, 0,
            float(A.XX), float(A.YY), float(A.ZZ), 1};
        orthonormal(vec, R, mat, mat+4);
        glMultMatrixf(mat);
    }

    
    void gleRotate(Vector3 const& A, Vector3 const& B,
                   Vector3 const& C, bool inverse)
    {
        float mat[16];
        for ( int i = 0; i < 3; ++i )
        {
            if ( inverse )
            {
                mat[4*i  ] = A[i];
                mat[4*i+1] = B[i];
                mat[4*i+2] = C[i];
            }
            else
            {
                mat[i  ] = A[i];
                mat[i+4] = B[i];
                mat[i+8] = C[i];
            }
            mat[i+12]  = 0;
            mat[i*4+3] = 0;
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
        float N = R * invsqrt(X*X+Y*Y);
        X *= N;
        Y *= N;
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
        float N = invsqrt(X*X+Y*Y+Z*Z);
        float vec[3] = { N*X, N*Y, N*Z };
        float mat[16] = {
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            float(P.XX), float(P.YY), float(P.ZZ), 1};
        orthonormal(vec, R, mat, mat+4, mat+8);
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
    
    //-----------------------------------------------------------------------
    
    void circle()
    {
        glNormal3f(0, 0, 1);
        glVertexPointer(2, GL_FLOAT, 0, 4+cir_);
        glDrawArrays(GL_LINE_LOOP, 0, ncircle);
    }
    
    void disc()
    {
        glNormal3f(0, 0, 1);
        glVertexPointer(2, GL_FLOAT, 0, 2+cir_);
        glDrawArrays(GL_TRIANGLE_FAN, 0, 2+ncircle);
    }

    void discUp()
    {
        glNormal3f(0, 0, 1);
        glVertexPointer(2, GL_FLOAT, 16, cir_);
        glDrawArrays(GL_TRIANGLE_FAN, 0, 2+ncircle/2);
    }
    
    void discDown()
    {
        glNormal3f(0, 0, -1);
        glFrontFace(GL_CW);
        glVertexPointer(2, GL_FLOAT, 16, cir_);
        glDrawArrays(GL_TRIANGLE_FAN, 0, 2+ncircle/2);
        glFrontFace(GL_CCW);
    }
    
    //-----------------------------------------------------------------------
    #pragma mark - Some Platonic solids
    
    /// code to calculate and print normals
    void printNormals(GLfloat const* pts, size_t cnt)
    {
        for ( size_t i = 0; i < cnt; ++i )
        {
            GLfloat const* x = pts + 9 * i;
            Vector3 a(x[0], x[1], x[2]);
            Vector3 b(x[3], x[4], x[5]);
            Vector3 c(x[6], x[7], x[8]);
            Vector3 N = cross(b-a, c-b).normalized();
            printf("%lu %+8.5f,%+8.5f,%+8.5f,", i, N.XX, N.YY, N.ZZ);
            printf("%+8.5f,%+8.5f,%+8.5f,", N.XX, N.YY, N.ZZ);
            printf("%+8.5f,%+8.5f,%+8.5f,\n", N.XX, N.YY, N.ZZ);
        }
    }
         
    /// Tetrahedron is make of 4 triangles = 12 vertices
    void initTetrahedron(GLuint buf1, GLuint buf2, GLfloat R=1.2f)
    {
        const GLfloat S = R / M_SQRT3;
        const GLfloat Y = 2 * S;
        const GLfloat B = -M_SQRT1_2 * S;
        const GLfloat Z = -3 * B;

        // -R,-S, B
        // +R,-S, B
        //  0, Y, B
        //  0, 0, Z
        const GLfloat pts[] = {
             R,-S, B,-R,-S, B, 0, Y, B,
             R,-S, B, 0, Y, B, 0, 0, Z,
             0, Y, B,-R,-S, B, 0, 0, Z,
            -R,-S, B, R,-S, B, 0, 0, Z,
        };
        
        const GLfloat dir[] = {
            +0.00000, 0.00000,-1.00000,+0.00000, 0.00000,-1.00000,+0.00000, 0.00000,-1.00000,
            +0.81650, 0.47140, 0.33333, 0.81650, 0.47140, 0.33333, 0.81650, 0.47140, 0.33333,
            -0.81650, 0.47140, 0.33333,-0.81650, 0.47140, 0.33333,-0.81650, 0.47140, 0.33333,
            +0.00000,-0.94281, 0.33333, 0.00000,-0.94281, 0.33333, 0.00000,-0.94281, 0.33333
        };
        
        glBindBuffer(GL_ARRAY_BUFFER, buf1);
        glBufferData(GL_ARRAY_BUFFER, sizeof(pts), pts, GL_STATIC_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, buf2);
        glBufferData(GL_ARRAY_BUFFER, sizeof(dir), dir, GL_STATIC_DRAW);
    }
    
    /// The star is made of two Tetrahedrons: 8 triangles = 24 vertices
    void initStar(GLuint buf1, GLuint buf2, GLfloat R=1.2f)
    {
        const GLfloat S = R / M_SQRT3;
        const GLfloat Y = 2 * S;
        const GLfloat B = -M_SQRT1_2 * S;
        const GLfloat Z = -3.0 * B;

        const GLfloat pts[] = {
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
        
        const GLfloat dir[] = {
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
        
        glBindBuffer(GL_ARRAY_BUFFER, buf1);
        glBufferData(GL_ARRAY_BUFFER, sizeof(pts), pts, GL_STATIC_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, buf2);
        glBufferData(GL_ARRAY_BUFFER, sizeof(dir), dir, GL_STATIC_DRAW);
    }

    
    /// Cube is make of 12 triangles = 36 vertices
    void initCube(GLuint buf1, GLuint buf2, GLfloat R=0.5773502692f)
    {
        const GLfloat pts[] = {
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
        
        const GLfloat dir[] = {
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
        
        glBindBuffer(GL_ARRAY_BUFFER, buf1);
        glBufferData(GL_ARRAY_BUFFER, sizeof(pts), pts, GL_STATIC_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, buf2);
        glBufferData(GL_ARRAY_BUFFER, sizeof(dir), dir, GL_STATIC_DRAW);
    }
    
    /// Octahedron is make of 8 triangles = 24 vertices
    void initOctahedron(GLuint buf1, GLuint buf2, GLfloat R=1.46459188756f)
    {
        // Eight triangles, ordered counterclockwise
        // set size to match the volume of the unit sphere
        const GLfloat pts[] = {
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
        const GLfloat N = 1 / M_SQRT3;
        const GLfloat dir[] = {
            +N,-N, N, N,-N, N, N,-N, N,
            -N, N,-N,-N, N,-N,-N, N,-N,
            -N,-N, N,-N,-N, N,-N,-N, N,
            +N, N,-N, N, N,-N, N, N,-N,
            -N, N, N,-N, N, N,-N, N, N,
            +N,-N,-N, N,-N,-N, N,-N,-N,
            +N, N, N, N, N, N, N, N, N,
            -N,-N,-N,-N,-N,-N,-N,-N,-N
        };
        
        glBindBuffer(GL_ARRAY_BUFFER, buf1);
        glBufferData(GL_ARRAY_BUFFER, sizeof(pts), pts, GL_STATIC_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, buf2);
        glBufferData(GL_ARRAY_BUFFER, sizeof(dir), dir, GL_STATIC_DRAW);
    }

    
#if ( 0 )
    void icoFace(GLfloat* a, GLfloat* b, GLfloat* c)
    {
        glNormal3f((a[0]+b[0]+c[0])/3.0f, (a[1]+b[1]+c[1])/3.0f, (a[2]+b[2]+c[2])/3.0f);
        glVertex3fv(a);
        glVertex3fv(b);
        glVertex3fv(c);
    }
    
    void icoFaceP(GLfloat* a, GLfloat* b, GLfloat* c)
    {
        printf("%2.0f, %2.0f, %2.0f,  ", a[0], a[1], a[2]);
        printf("%2.0f, %2.0f, %2.0f,  ", b[0], b[1], b[2]);
        printf("%2.0f, %2.0f, %2.0f,\n", c[0], c[1], c[2]);
    }
    
    void icoFaceN(GLfloat* a, GLfloat* b, GLfloat* c)
    {
        Vector3 A(a[0], a[1], a[2]);
        Vector3 B(b[0], b[1], b[2]);
        Vector3 C(c[0], c[1], c[2]);
        Vector3 n = normalize(cross(B-A, C-A));
        glNormal3f(n.XX, n.YY, n.ZZ);
        glVertex3fv(a);
        glVertex3fv(b);
        glVertex3fv(c);
        printf("%+9.7f,%+9.7f,%+9.7f, ", n[0], n[1], n[2]);
        printf("%+9.7f,%+9.7f,%+9.7f, ", n[0], n[1], n[2]);
        printf("%+9.7f,%+9.7f,%+9.7f,\n", n[0], n[1], n[2]);
    }
    
    void icoFace(GLfloat* pts, size_t a, size_t b, size_t c)
    {
        icoFace(pts+3*a, pts+3*b, pts+3*c);
    }
    
    void icosahedron()
    {
        const GLfloat G = 0.5+0.5*std::sqrt(5.0);
        const GLfloat O = 1/std::sqrt(G*G+1); //0.5257311121f;
        const GLfloat T = G * O;   //0.8506508084f;
        
        // Twelve vertices of icosahedron on unit sphere
        GLfloat pts[] = {
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
        
        /* The faces are ordered with increasing Z */
        glBegin(GL_TRIANGLES);
        icoFace(pts, 5,  6, 9);
        icoFace(pts, 5, 11, 6);
        
        icoFace(pts, 6, 3,  9);
        icoFace(pts, 2, 11, 5);
        icoFace(pts, 1, 5,  9);
        icoFace(pts, 0, 6, 11);//
        
        icoFace(pts, 0, 3,  6);
        icoFace(pts, 1, 2,  5);
        
        icoFace(pts, 1, 9, 10);
        icoFace(pts, 0, 11, 8);//
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
        glEnd();
    }
#endif
    
    /// Icosahedrong with 20 triangles = 60 vertices
    void initIcosahedron(GLuint buf1, GLuint buf2, GLfloat R=1.0f)
    {
        const GLfloat T = R * 0.8506508084f;      // (1 + sqrt(5))/2
        const GLfloat O = R * 0.5257311121f;      // 1 / sqrt(1+T^2)
        
        const GLfloat pts[] = {
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
        
        const GLfloat dir[] = {
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
        
        glBindBuffer(GL_ARRAY_BUFFER, buf1);
        glBufferData(GL_ARRAY_BUFFER, sizeof(pts), pts, GL_STATIC_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, buf2);
        glBufferData(GL_ARRAY_BUFFER, sizeof(dir), dir, GL_STATIC_DRAW);
    }
    
    /// Three fins similar to the tail of a V2 rocket
    void initArrowTail(GLuint buf1, GLuint buf2, GLfloat R=0.1f,
                       GLfloat B=-0.5f, GLfloat H=-1.5f, GLfloat L=2.0f)
    {
        const GLfloat T = B + L;
        const GLfloat U = H + L;
        const GLfloat C = 0.5f;
        const GLfloat S = M_SQRT3_2;
        const GLfloat cR = R * C;
        const GLfloat sR = R * S;

        const GLfloat pts[] = {
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
        
        const GLfloat dir[] = {
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
        
        glBindBuffer(GL_ARRAY_BUFFER, buf1);
        glBufferData(GL_ARRAY_BUFFER, sizeof(pts), pts, GL_STATIC_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, buf2);
        glBufferData(GL_ARRAY_BUFFER, sizeof(dir), dir, GL_STATIC_DRAW);
    }

    //-----------------------------------------------------------------------

    void drawBuffer(GLuint buf1, GLuint buf2, unsigned cnt, GLenum mode)
    {
        glBindBuffer(GL_ARRAY_BUFFER, buf1);
        glVertexPointer(3, GL_FLOAT, 0, nullptr);
        glEnableClientState(GL_NORMAL_ARRAY);
        glBindBuffer(GL_ARRAY_BUFFER, buf2);
        glNormalPointer(GL_FLOAT, 0, nullptr);
        glDrawArrays(mode, 0, cnt);
        glDisableClientState(GL_NORMAL_ARRAY);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

    void initBuffers()
    {
        if ( !glIsBuffer(buf_[0]) )
        {
            glGenBuffers(16, buf_);
            initTetrahedron(buf_[0], buf_[1]);
            initOctahedron(buf_[2], buf_[3]);
            initIcosahedron(buf_[4], buf_[5]);
            initArrowTail(buf_[6], buf_[7]);
            initCube(buf_[8], buf_[9]);
            initStar(buf_[10], buf_[11]);
            glBindBuffer(GL_ARRAY_BUFFER, 0);
        }
    }
    
    void tetrahedron() { drawBuffer(buf_[0], buf_[1], 12, GL_TRIANGLES); }
    void octahedron() { drawBuffer(buf_[2], buf_[3], 24, GL_TRIANGLES); }
    void icosahedron() { drawBuffer(buf_[4], buf_[5], 60, GL_TRIANGLES); }

    void arrowTail() { drawBuffer(buf_[6], buf_[7], 45, GL_TRIANGLES); }
    void cube() { drawBuffer(buf_[8], buf_[9], 36, GL_TRIANGLES); }
    void star() { drawBuffer(buf_[10], buf_[11], 24, GL_TRIANGLES); }

    //-----------------------------------------------------------------------
#pragma mark - Tubes
    
    void tubeZ(GLfloat B, GLfloat T, int inc)
    {
        assert_true( B <= T );
        glBegin(GL_TRIANGLE_STRIP);
        for( size_t n = 0; n <= ncircle; n += inc )
        {
            glNormal3f(cos_(n), sin_(n), 0);
            glVertex3f(cos_(n), sin_(n), T);
            glVertex3f(cos_(n), sin_(n), B);
        }
        glEnd();
    }

    void tubeZ(GLfloat B, GLfloat rB, GLfloat T, GLfloat rT, int inc)
    {
        assert_true( B <= T );
        const GLfloat N = 1.f/sqrtf((T-B)*(T-B)+(rT-rB)*(rT-rB));
        const GLfloat C = N * (T-B);
        const GLfloat S = N * (rB-rT);
        glBegin(GL_TRIANGLE_STRIP);
        for( size_t n = 0; n <= ncircle; n += inc )
        {
            glNormal3f(C*cos_(n), C*sin_(n), S);
            glVertex3f(rT*cos_(n), rT*sin_(n), T);
            glVertex3f(rB*cos_(n), rB*sin_(n), B);
        }
        glEnd();
    }

    size_t initTube(GLuint buf1, GLuint buf2, GLfloat A, GLfloat B, size_t inc)
    {
        size_t sec = ncircle / inc;
        size_t nbf = 6 * sec + 6;     // number of coordinates
        GLfloat pos[nbf], dir[nbf];
        assert_true( A <= B );

        for( size_t n = 0, p = 0; n < nbf; n += 6, p += inc )
        {
            GLfloat C = cos_(p);
            GLfloat S = sin_(p);
            
            dir[n  ] = C;
            dir[n+1] = S;
            dir[n+2] = 0;
            
            dir[n+3] = C;
            dir[n+4] = S;
            dir[n+5] = 0;

            pos[n  ] = C;
            pos[n+1] = S;
            pos[n+2] = B;

            pos[n+3] = C;
            pos[n+4] = S;
            pos[n+5] = A;
        }
        glBindBuffer(GL_ARRAY_BUFFER, buf1);
        glBufferData(GL_ARRAY_BUFFER, nbf*sizeof(GLfloat), pos, GL_STATIC_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, buf2);
        glBufferData(GL_ARRAY_BUFFER, nbf*sizeof(GLfloat), dir, GL_STATIC_DRAW);
        return nbf;
    }
    
    /// hexagon has the same surface as a disc of radius 1.
    void initHexTube(GLuint buf1, GLuint buf2, GLfloat A, GLfloat B)
    {
        constexpr GLfloat R = 1.0996361107912678f; //std::sqrt( 2 * M_PI / ( 3 * std::sqrt(3) ));
        constexpr GLfloat C = 0.8660254037844386f; //std::sqrt(3)/2;
        constexpr GLfloat S = 0.5f;
        constexpr GLfloat H = R * C, X = R * S;
        
        const GLfloat pts[] = {
             R,  0, B,  R,  0, A,
             X,  H, B,  X,  H, A,
            -X,  H, B, -X,  H, A,
            -R,  0, B, -R,  0, A,
            -X, -H, B, -X, -H, A,
             X, -H, B,  X, -H, A,
             R,  0, B,  R,  0, A };
        
        const GLfloat dir[] = {
             1,  0, 0,  1,  0, 0,
             S,  C, 0,  S,  C, 0,
            -S,  C, 0, -S,  C, 0,
            -1,  0, 0, -1,  0, 0,
            -S, -C, 0, -S, -C, 0,
             S, -C, 0,  S, -C, 0,
             1,  0, 0,  1,  0, 0 };

        glBindBuffer(GL_ARRAY_BUFFER, buf1);
        glBufferData(GL_ARRAY_BUFFER, sizeof(pts), pts, GL_STATIC_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, buf2);
        glBufferData(GL_ARRAY_BUFFER, sizeof(dir), dir, GL_STATIC_DRAW);
    }
    
    
    void initTubeBuffers()
    {
        if ( !glIsBuffer(tub_[0]) )
        {
            glGenBuffers(24, tub_);
            const GLfloat B = -4.f, T = 256.f;
            initTube(tub_[ 0], tub_[ 1], 0, 1, 8);
            initTube(tub_[ 2], tub_[ 3], 0, 1, 4);
            initTube(tub_[ 4], tub_[ 5], 0, 1, 2);
            initTube(tub_[ 6], tub_[ 7], 0, 1, 1);
            initTube(tub_[ 8], tub_[ 9], B, T, 8);
            initTube(tub_[10], tub_[11], B, T, 4);
            initTube(tub_[12], tub_[13], B, T, 2);
            initTube(tub_[14], tub_[15], 0, T, 8);
            initTube(tub_[16], tub_[17], 0, T, 4);
            initTube(tub_[18], tub_[19], 0, T, 2);
            initHexTube(tub_[22], tub_[23], 0, 1);
            glBindBuffer(GL_ARRAY_BUFFER, 0);
        }
    }

#if ( 1 )
    // using Vertex Buffer Objects
    void tube1()     { drawBuffer(tub_[ 0], tub_[ 1], 2+ncircle/4, GL_TRIANGLE_STRIP); }
    void tube2()     { drawBuffer(tub_[ 2], tub_[ 3], 2+ncircle/2, GL_TRIANGLE_STRIP); }
    void tube4()     { drawBuffer(tub_[ 4], tub_[ 5], 2+ncircle  , GL_TRIANGLE_STRIP); }
    void tube8()     { drawBuffer(tub_[ 6], tub_[ 7], 2+ncircle*2, GL_TRIANGLE_STRIP); }
    void longTube1() { drawBuffer(tub_[ 8], tub_[ 9], 2+ncircle/4, GL_TRIANGLE_STRIP); }
    void longTube2() { drawBuffer(tub_[10], tub_[11], 2+ncircle/2, GL_TRIANGLE_STRIP); }
    void longTube4() { drawBuffer(tub_[12], tub_[13], 2+ncircle  , GL_TRIANGLE_STRIP); }
    void halfTube1() { drawBuffer(tub_[14], tub_[15], 2+ncircle/4, GL_TRIANGLE_STRIP); }
    void halfTube2() { drawBuffer(tub_[16], tub_[17], 2+ncircle/2, GL_TRIANGLE_STRIP); }
    void halfTube4() { drawBuffer(tub_[18], tub_[19], 2+ncircle  , GL_TRIANGLE_STRIP); }
    void hexTube()   { drawBuffer(tub_[22], tub_[23], 14, GL_TRIANGLE_STRIP); }
#else
    // unbuffered functions:
    void tube1()     { tubeZ(0, 1, 8); }
    void tube2()     { tubeZ(0, 1, 4); }
    void tube4()     { tubeZ(0, 1, 2); }
    void tube8()     { tubeZ(0, 1, 1); }
    void longTube1() { tubeZ(-4.f, 256.f, 8); }
    void longTube2() { tubeZ(-4.f, 256.f, 4); }
    void longTube4() { tubeZ(-4.f, 256.f, 2); }
    void halfTube1() { tubeZ(-256.f, 0.0f, 8); }
    void halfTube2() { tubeZ(-256.f, 0.0f, 4); }
    void halfTube4() { tubeZ(-256.f, 0.0f, 2); }
    void hexTube()   { hexTubeZ(0, 1); }
#endif
 
    void tubeZ(GLfloat za, GLfloat ra, gle_color ca, GLfloat zb, GLfloat rb, gle_color cb)
    {
        glBegin(GL_TRIANGLE_STRIP);
        for( size_t n = 0; n <= ncircle; ++n )
        {
            cb.load_load();
            glNormal3f(   cos_(n),    sin_(n),  0);
            glVertex3f(rb*cos_(n), rb*sin_(n), zb);
            ca.load_load();
            glNormal3f(   cos_(n),    sin_(n),  0);
            glVertex3f(ra*cos_(n), ra*sin_(n), za);
        }
        glEnd();
    }
    
    void cylinder1()
    {
        tube1();
        glTranslatef(0,0,1);
        discUp();
        glTranslatef(0,0,-1);
        discDown();
    }
    
    void cylinderZ()
    {
        glTranslatef(0,0,0.5f);
        discUp();
        glTranslatef(0,0,-1);
        tube1();
        discDown();
        glTranslatef(0,0,0.5f);
    }

    /// spherocylinder of length L, radius R, centered and aligned with axis Z
    void capsuleZ(GLfloat L, GLfloat R)
    {
        const size_t fin = ncircle >> 2;
        const size_t C = 4;
        //display strips along the side of the volume:
        for ( size_t t = 0; t < 4*fin; t += C )
        {
            //compute the transverse angles:
            GLfloat cb = cos_(t),   sb = sin_(t);
            GLfloat ca = cos_(t+C), sa = sin_(t+C);
            GLfloat cB = R * cb, sB = R * sb;
            GLfloat cA = R * ca, sA = R * sa;
            
            //draw one srip of the oval:
            glBegin(GL_TRIANGLE_STRIP);
            for ( size_t i=0; i <= fin; i += C )
            {
                GLfloat x = cos_(i), y = sin_(i);
                glNormal3f(ca*y, sa*y,     x);
                glVertex3f(cA*y, sA*y, L+R*x);
                glNormal3f(cb*y, sb*y,     x);
                glVertex3f(cB*y, sB*y, L+R*x);
            }
            for ( int i=fin; i >= 0; i -= C )
            {
                GLfloat x = -cos_(i), y = sin_(i);
                glNormal3f(ca*y, sa*y,   x);
                glVertex3f(cA*y, sA*y, R*x);
                glNormal3f(cb*y, sb*y,   x);
                glVertex3f(cB*y, sB*y, R*x);
            }
            glEnd();
        }
    }

    //-----------------------------------------------------------------------
#pragma mark - Spheres

#ifdef PLATONIC_H
    
    /// using icosahedrons to render the sphere:
    Platonic::Solid ico1(Platonic::Solid::ICOSAHEDRON, gle::finesse/6);
    Platonic::Solid ico2(Platonic::Solid::ICOSAHEDRON, gle::finesse/3);
    Platonic::Solid ico4(Platonic::Solid::ICOSAHEDRON, gle::finesse);
    Platonic::Solid ico8(Platonic::Solid::ICOSAHEDRON, gle::finesse*2);
    Platonic::Solid icoH1(Platonic::Solid::HEMISPHERE, gle::finesse/6);
    Platonic::Solid icoH2(Platonic::Solid::HEMISPHERE, gle::finesse/3);
    Platonic::Solid icoH4(Platonic::Solid::HEMISPHERE, gle::finesse);

    void drawPlatonic(Platonic::Solid & ico)
    {
        glEnableClientState(GL_NORMAL_ARRAY);
        glVertexPointer(3, GL_FLOAT, 0, ico.vertex_data());
        glNormalPointer(GL_FLOAT, 0, ico.vertex_data());
        glDrawElements(GL_TRIANGLES, 3*ico.nb_faces(), GL_UNSIGNED_INT, ico.faces_data());
        glDisableClientState(GL_NORMAL_ARRAY);
    }
    
    void sphere1U() { drawPlatonic(ico1); }
    void sphere2U() { drawPlatonic(ico2); }
    void sphere4U() { drawPlatonic(ico4); }
    void sphere8U() { drawPlatonic(ico8); }
    void hemisphere1U() { drawPlatonic(icoH1); }
    void hemisphere2U() { drawPlatonic(icoH2); }
    void hemisphere4U() { drawPlatonic(icoH4); }

    //-----------------------------------------------------------------------
    
    GLuint initIcoBuffer(GLuint buf1, GLuint buf2, Platonic::Solid & ico)
    {
        //std::clog << "initializeIco ico " << ico.nb_faces() << '\n';
        
        // upload vertex data:
        glBindBuffer(GL_ARRAY_BUFFER, buf1);
        glBufferData(GL_ARRAY_BUFFER, 3*ico.nb_vertices()*sizeof(float), ico.vertex_data(), GL_STATIC_DRAW);
        
        // upload indices:
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buf2);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, 3*ico.nb_faces()*sizeof(unsigned), ico.faces_data(), GL_STATIC_DRAW);
        return ico.nb_faces();
    }
    
    void initSphereBuffers()
    {
        if ( !glIsBuffer(ico_[0]) )
        {
            glGenBuffers(16, ico_);
            ico_nfaces[0] = initIcoBuffer(ico_[0], ico_[1], ico1);
            ico_nfaces[1] = initIcoBuffer(ico_[2], ico_[3], ico2);
            ico_nfaces[2] = initIcoBuffer(ico_[4], ico_[5], ico4);
            ico_nfaces[3] = initIcoBuffer(ico_[6], ico_[7], ico8);
            ico_nfaces[4] = initIcoBuffer(ico_[8], ico_[9], icoH1);
            ico_nfaces[5] = initIcoBuffer(ico_[10], ico_[11], icoH2);
            ico_nfaces[6] = initIcoBuffer(ico_[12], ico_[13], icoH4);
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
            glBindBuffer(GL_ARRAY_BUFFER, 0);
#if 0
            for ( int n = 0; n < 4; ++n )
                fprintf(stderr, "GLE's icosahedron %i has %u faces\n", n, ico_nfaces[n]);
#endif
        }
    }
    
    void drawIcoBuffer(GLuint nfaces, GLuint buf1, GLuint buf2)
    {
        if ( glIsBuffer(buf1) )
        {
            glEnableClientState(GL_NORMAL_ARRAY);
            // the normal in each vertex is equal to the vertex!
            glBindBuffer(GL_ARRAY_BUFFER, buf1);
            glVertexPointer(3, GL_FLOAT, 0, nullptr);
            glNormalPointer(GL_FLOAT, 0, nullptr);
            glBindBuffer(GL_ARRAY_BUFFER, 0);
            
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buf2);
            glDrawElements(GL_TRIANGLES, 3*nfaces, GL_UNSIGNED_INT, nullptr);
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
            glDisableClientState(GL_NORMAL_ARRAY);
        }
    }
    
    void sphere1() { drawIcoBuffer(ico_nfaces[0], ico_[0], ico_[1]); }
    void sphere2() { drawIcoBuffer(ico_nfaces[1], ico_[2], ico_[3]); }
    void sphere4() { drawIcoBuffer(ico_nfaces[2], ico_[4], ico_[5]); }
    void sphere8() { drawIcoBuffer(ico_nfaces[3], ico_[6], ico_[7]); }
    void hemisphere1() { drawIcoBuffer(ico_nfaces[4], ico_[8], ico_[9]); }
    void hemisphere2() { drawIcoBuffer(ico_nfaces[5], ico_[10], ico_[11]); }
    void hemisphere4() { drawIcoBuffer(ico_nfaces[4], ico_[12], ico_[13]); }
    
#else
    
    // Without using PlatonicSolid
    void initSphereBuffers() { }

    /// using trigonometric functions to draw a ball of radius 1
    void drawSphere(size_t inc)
    {
        for ( size_t n = 0; n < ncircle/2; n += inc )
        {
            real U = cos_(n), R = sin_(n);
            real L = cos_(n+inc), S = sin_(n+inc);
            glBegin(GL_TRIANGLE_STRIP);
            for ( size_t p = 0; p <= ncircle; p += inc )
            {
                glNormal3f(R*cos_(p), R*sin_(p), U);
                glVertex3f(R*cos_(p), R*sin_(p), U);
                glNormal3f(S*cos_(p), S*sin_(p), L);
                glVertex3f(S*cos_(p), S*sin_(p), L);
            }
            glEnd();
        }
    }
    
    void sphere1() { drawSphere(8); }
    void sphere2() { drawSphere(4); }
    void sphere4() { drawSphere(2); }
    void sphere8() { drawSphere(1); }

#endif
    
    void dualPassSphere1()
    {
        assert_true(glIsEnabled(GL_CULL_FACE));
        glCullFace(GL_FRONT);
        gle::sphere1();
        glCullFace(GL_BACK);
        gle::sphere1();
    }
    void dualPassSphere2()
    {
        assert_true(glIsEnabled(GL_CULL_FACE));
        glCullFace(GL_FRONT);
        gle::sphere2();
        glCullFace(GL_BACK);
        gle::sphere2();
    }
    void dualPassSphere4()
    {
        assert_true(glIsEnabled(GL_CULL_FACE));
        glCullFace(GL_FRONT);
        gle::sphere4();
        glCullFace(GL_BACK);
        gle::sphere4();
    }
    void dualPassSphere8()
    {
        assert_true(glIsEnabled(GL_CULL_FACE));
        glCullFace(GL_FRONT);
        gle::sphere8();
        glCullFace(GL_BACK);
        gle::sphere8();
    }

    /// draw a Torus of radius R and a thickness 2*T
    void torusZ(GLfloat R, GLfloat T, size_t S)
    {
        for ( size_t n = 0; n < ncircle; n += S )
        {
            GLfloat X0 = cos_(n  ), Y0 = sin_(n  );
            GLfloat X1 = cos_(n+S), Y1 = sin_(n+S);
            glBegin(GL_TRIANGLE_STRIP);
            for ( size_t p = 0; p <= ncircle; p += 2*S )
            {
                glNormal3f(X0*cos_(p), Y0*cos_(p), sin_(p));
                glVertex3f(X0*(R+T*cos_(p)), Y0*(R+T*cos_(p)), T*sin_(p));
                glNormal3f(X1*cos_(p), Y1*cos_(p), sin_(p));
                glVertex3f(X1*(R+T*cos_(p)), Y1*(R+T*cos_(p)), T*sin_(p));
            }
            glEnd();
        }
    }

    
    //-----------------------------------------------------------------------
#pragma mark - 3D primitives
    
    /**
     Draw a cylindrical band on the equator of a sphere of radius 1.
     The band is in the XY plane. The axis of the cylinder is Z.
     The band is made of triangles indicating the clockwise direction.
     */
    void drawArrowedBand(const size_t nb_triangles, float width)
    {
        GLfloat A = GLfloat(2 * M_PI / nb_triangles);
        GLfloat W = GLfloat(width * A / M_SQRT3);
        GLfloat R = 1.0f / cosf(A*0.5f);
        
        glBegin(GL_TRIANGLES);
        glNormal3f(1, 0, 0);
        glVertex3f(1, 0, W);
        glVertex3f(1, 0,-W);
        for ( size_t ii = 1; ii < nb_triangles; ++ii )
        {
            GLfloat ang = ii * A;
            GLfloat c = R * cosf(ang);
            GLfloat s = R * sinf(ang);
            
            glNormal3f(c, s, 0);
            glVertex3f(c, s, 0);
            glVertex3f(c, s, W);
            glVertex3f(c, s,-W);
        }
        glNormal3f(1, 0, 0);
        glVertex3f(1, 0, 0);
        glEnd();
    }
    
    
    void drawThreeBands(const size_t nb_triangles)
    {
        drawArrowedBand(nb_triangles, 0.25);
        glRotated(-90,1,0,0);
        drawArrowedBand(nb_triangles, 0.25);
        glRotated(90,0,1,0);
        drawArrowedBand(nb_triangles, 0.25);
    }
    
    
    void coneZ(GLfloat R, GLfloat B, GLfloat T, bool closed)
    {
        if ( closed )
        {
            glBegin(GL_TRIANGLE_FAN);
            glNormal3f( 0, 0, B );
            glVertex3f( 0, 0, B );
            for ( size_t n = 0; n <= ncircle; ++n )
                glVertex3f(R*cos_(n), -R*sin_(n), B);
            glEnd();
        }
        glBegin(GL_TRIANGLE_FAN);
        glNormal3f( 0, 0, 1 );
        glVertex3f( 0, 0, T );
        const GLfloat S = -1.f/sqrtf((T-B)*(T-B)+1);
        const GLfloat C = (B-T)*S;
        for ( size_t n = 0; n <= ncircle; ++n )
        {
            glNormal3f(C*cos_(n), C*sin_(n), S);
            glVertex3f(R*cos_(n), R*sin_(n), B);
        }
        glEnd();
    }
    
    
    GLfloat dumbbellRadius(GLfloat z)
    {
        const GLfloat PIF = GLfloat(M_PI);
        return sinf(PIF*z) * ( 1.3f + cosf(2*PIF*z) );
    }
    
    
    /**
     Draw a surface of revolution around the Z-axis.
     The surface goes from Z to Z_max, and its radius is
     given by the function `radius`(z) provided as argument.
     */
    void drawRevolution(GLfloat (*radius)(GLfloat), GLfloat Z, GLfloat Zmax, GLfloat dZ)
    {
        GLfloat R = radius(Z);
        GLfloat R0, Z0, dR, dN;
        
        while ( Z < Zmax )
        {
            Z0 = Z;
            R0 = R;
            Z += dZ;
            R = radius(Z);
            
            dR = ( R - R0 ) / dZ;
            dN = 1.0f / sqrtf( 1 + dR * dR );
            dR = dR*dN;
            
            glBegin(GL_TRIANGLE_STRIP);
            for ( size_t n = 0; n <= ncircle; ++n )
            {
                glNormal3f(dN*cos_(n), dN*sin_(n),-dR);
                glVertex3f(R *cos_(n), R *sin_(n), Z);
                glVertex3f(R0*cos_(n), R0*sin_(n), Z0);
            }
            glEnd();
        }
    }

    void dumbbell()
    {
        drawRevolution(dumbbellRadius, 0, 1, 0.0625);
    }
    
    GLfloat barrelRadius(GLfloat z)
    {
        return sinf(3.14159265359f*z);
    }
    
    void barrel()
    {
        drawRevolution(barrelRadius, 0, 1, 0.0625);
    }
    
    //-----------------------------------------------------------------------
#pragma mark - Object Placement
    
    
    /**
     draw back first, and then front of object,
     GL_CULL_FACE should be enabled
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
    void gleObject(Vector1 const& X, Vector1 const& D, const float R, void (*obj)())
    {
        glPushMatrix();
        transAlignZ(X, R, D);
        obj();
        glPopMatrix();
    }
    
    void gleObject(Vector2 const& X, Vector2 const& D, const float R, void (*obj)())
    {
        glPushMatrix();
        transAlignZ(X, R, D);
        obj();
        glPopMatrix();
    }
    
    void gleObject(Vector3 const& X, Vector3 const& D, const float R, void (*obj)())
    {
        glPushMatrix();
        transAlignZ(X, R, D);
        obj();
        glPopMatrix();
    }

    //-----------------------------------------------------------------------
#pragma mark - Tubes
    
    
    void gleTube(Vector1 const& A, Vector1 const& B, float R, void (*obj)())
    {
        glPushMatrix();
        stretchAlignZ(A, B, R);
        obj();
        glPopMatrix();
    }
    
    void gleTube(Vector2 const& A, Vector2 const& B, float R, void (*obj)())
    {
        glPushMatrix();
        stretchAlignZ(A, B, R);
        obj();
        glPopMatrix();
    }
    
    void gleTube(Vector3 const& A, Vector3 const& B, float R, void (*obj)())
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
            float2 pts[4] = { A+rad*d, A-rad*d, B+rad*d, B-rad*d };
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
        float4 col[4] = { cA, cA, cB, cB };
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
            float4 col[] = { cA, cA, cB, cB };
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
        float2 pts[6] = { b-db, b, a-da, a+da, b, b+db };
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
        float2 pts[6] = { b-db, b, a-da, a+da, b, b+db };
        float4 col[6] = { cb, cb, ca, ca, cb, cb };
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
        float2 pts[8] = { A, A-rad*dA, B, B-rad*dB,
                          A, A+rad*dA, B, B+rad*dB };
        glVertexPointer(2, GL_FLOAT, 0, pts);
        glDrawArrays(GL_TRIANGLE_FAN, 0, 4);
        glDrawArrays(GL_TRIANGLE_FAN, 4, 4);
    }
    
    void drawBar(Vector3 const& A, Vector3 const& dA,
                 Vector3 const& B, Vector3 const& dB, real rad)
    {
        Vector3 ab = normalize( A - B );
        Vector3 ea = cross(ab, dA);
        Vector3 eb = cross(ab, dB);
        float3 pts[16] = { A-rad*(dA-ea), A-rad*(dA+ea), B-rad*(dB-eb), B-rad*(dB+eb),
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
        const GLfloat S(1.0996361107912678f); //sqrt( 2 * M_PI / ( 3 * sqrt(3) ));
        const GLfloat R(diameter * S);
        const GLfloat H(R * 0.8660254037844386f); //0.5f * sqrt(3);
        const GLfloat X(R * 0.5f);
        
        Vector2 x = ( B - A ).normalized(H);
        Vector2 y = x.orthogonal(X);
        float2 pts[20] = {
            {0,0}, x+y, y+y, y-x, -x-y, -y-y, x-y, x+y,
            x+y, x-y, B-A-y-x, B-A+y-x,
            {0,0}, x+y, y+y, y-x, -x-y, -y-y, x-y, x+y };
        glVertexPointer(2, GL_FLOAT, 0, pts);
        
        glPushMatrix();
        translate(A);
        // draw hexagon around 'a':
        glDrawArrays(GL_TRIANGLE_FAN, 0, 8);
        // a band from 'a' to 'b'
        glDrawArrays(GL_TRIANGLE_FAN, 8, 4);
        translate(B-A);
        // an hexagon centered around 'b'
        glDrawArrays(GL_TRIANGLE_FAN, 12, 8);
        glPopMatrix();
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
        gle::longCone();
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
        glTranslatef(0, 0, -0.5f);
        cylinder1();
        glPopMatrix();
    }
    
    
    //-----------------------------------------------------------------------
    
    void drawArrowTail(Vector1 const& pos, Vector1 const& dir, float rad)
    {
        GLfloat dx(rad * dir.XX);
        GLfloat cx(pos.XX - dx * 0.5 );;
        GLfloat pts[12] = {cx, 0, cx-dx, -dx, cx+dx, -dx,
                           cx+2*dx, 0, cx+dx, dx, cx-dx, dx};
        glVertexPointer(2, GL_FLOAT, 0, pts);
        glDrawArrays(GL_TRIANGLE_FAN, 0, 6);
    }
    
    void drawArrowTail(Vector2 const& pos, Vector2 const& dir, float rad)
    {
        GLfloat dx(rad * dir.XX);
        GLfloat dy(rad * dir.YY);
        GLfloat cx(pos.XX - 1.5f * dx);
        GLfloat cy(pos.YY - 1.5f * dy);
        GLfloat ex(cx + 2 * dx);
        GLfloat ey(cy + 2 * dy);
        GLfloat pts[12] = {cx+dx, cy+dy, cx+dy, cy-dx, ex+dy, ey-dx,
                           ex+dx, ey+dy, ex-dy, ey+dx, cx-dy, cy+dx};
        glVertexPointer(2, GL_FLOAT, 0, pts);
        glDrawArrays(GL_TRIANGLE_FAN, 0, 6);
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
        gle::longCone();
        glPopMatrix();
    }
    
    void drawArrow(Vector2 const& A, Vector2 const& B, float R)
    {
        glPushMatrix();
        stretchAlignZ(A, B, R);
        tube1();
        glTranslatef(0, 0, 1);
        glScalef(3.0, 3.0, 3*R);
        gle::longCone();
        glPopMatrix();
    }
    
    void drawArrow(Vector3 const& A, Vector3 const& B, float R)
    {
        glPushMatrix();
        stretchAlignZ(A, B, R);
        tube1();
        glTranslatef(0, 0, 1);
        glScalef(3.0, 3.0, 3*R);
        gle::longCone();
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
    void bitmapText(const char text[], void* font, GLfloat vshift)
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
     set the current raster position to `w`
     */
    void drawText(Vector3 const& vec, const char text[], void* font, float dx)
    {
        glPushAttrib(GL_CURRENT_BIT|GL_ENABLE_BIT);
        glDisable(GL_DEPTH_TEST);
        glDisable(GL_LIGHTING);
        glDisable(GL_ALPHA_TEST);
        int L = 1;
        int H = fontHeight(font);
        int W = maxTextWidth(text, font, L);
        rasterPos(vec);
        //translate to center the bitmap:
        glBitmap(0,0,0,0,-W*dx,-H/3,nullptr);
        bitmapText(text, font, H);
        glPopAttrib();
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
    void drawText(const char text[], void* font, const gle_color bcol,
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
        glPushAttrib(GL_CURRENT_BIT|GL_ENABLE_BIT);
        glDisable(GL_DEPTH_TEST);
        glDisable(GL_LIGHTING);
        glDisable(GL_ALPHA_TEST);
        
        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();
        glLoadIdentity();
        
        glMatrixMode(GL_PROJECTION);
        glPushMatrix();
        glLoadIdentity();
        glOrtho(0, width, 0, height, 0, 1);
        
        glRasterPos2i(0, 0);
        glBitmap(0, 0, 0, 0, px, py, nullptr);
        
        if ( bcol.visible() )
        {
            glPushAttrib(GL_LIGHTING_BIT|GL_CURRENT_BIT);
            glDisable(GL_LIGHTING);
            int R = abs(lineHeight);
            int B = std::min(py, py + n_lines * lineHeight);
            int T = std::max(py, py + n_lines * lineHeight);
            
            int rec[4] = { px-R, B, px+textWidth+R, T+R+R/2+R/4 };
            
            bcol.load();
            drawNiceRectangle(rec, 3, GL_TRIANGLE_FAN);
            
            glPopAttrib();
            
            if ( position == 4 )
            {
                glLineWidth(0.5);
                drawNiceRectangle(rec, 3, GL_LINE_STRIP);
            }
        }
        
        bitmapText(text, font, lineHeight);
        
        glMatrixMode(GL_PROJECTION);
        glPopMatrix();
        glMatrixMode(GL_MODELVIEW);
        glPopMatrix();
        glPopAttrib();
    }
    
    
    //-----------------------------------------------------------------------
#pragma mark - Misc
    
    /**
     Draw an array of pixels using GL_TRIANGLE_STRIP
     
     The array rgba[] should ( nbc * width * height ) bytes,
     containing nbc-components (eg. RGBA) per pixel and
     stored by columns:
     
         load(i,j) = rgba[ nbc*(i+height*j) ]
         0 <= i < height
         0 <= j < width
     
     `pos` is the position of the top-left corner
     `dx` is the direction of the width
     `dy` is the direction of the height
     The magnitudes of `dx` and `dy` indicates the dimensions of a pixel.
     They may be of different magnitudes, and not necessarily orthogonal.
     */
    void drawPixels(int width, int height, int nbc, GLubyte rgba[], Vector2 pos, Vector2 dx, Vector2 dy)
    {
        glPushAttrib(GL_ENABLE_BIT|GL_POLYGON_BIT);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glDisable(GL_DEPTH_TEST);
        glDisable(GL_CULL_FACE);
        glDisable(GL_LIGHTING);
        
        Vector2 left, right;
        GLubyte * col = rgba;
        
        for ( int jj = 0; jj < width; ++jj )
        {
            left  = pos + dx * jj;
            right = left + dx;
            for ( int ii = 0; ii < height; ++ii )
            {
                if ( nbc == 3 )
                    glColor3ubv(col);
                else
                    glColor4ubv(col);
                glBegin(GL_TRIANGLE_STRIP);
                col += nbc;
                gleVertex(left);
                gleVertex(right);
                left  += dy;
                right += dy;
                gleVertex(left);
                gleVertex(right);
                glEnd();
            }
        }
        
        glPopAttrib();
    }
    
    //-----------------------------------------------------------------------
    
    
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

    
    void drawNiceRectangle(const int rec[4], const int rad, GLint prim)
    {
        GLfloat L(rec[0]), B(rec[1]), R(rec[2]), T(rec[3]);
        GLfloat D(rad);
        GLfloat pts[18] = {L,B+D,L+D,B,R-D,B,R,B+D,R,T-D,R-D,T,L+D,T,L,T-D,L,B+D};
        glVertexPointer(2, GL_FLOAT, 0, pts);
        glDrawArrays(prim, 0, 9);
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
        glPushAttrib(GL_ENABLE_BIT);
        glDisable(GL_LIGHTING);
        glDisable(GL_DEPTH_TEST);
        
        glBegin(GL_LINE_LOOP);
        drawRectangle(rec);
        glEnd();
        
        glPopAttrib();
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
        return a + (( std::abs(a) + b ) & 1);
    }

    void drawTiledFloor(int R, float T, float Z, gle_color col1, gle_color col2)
    {
        float H = T * 0.5;
        int Q = std::floor( double(R) * M_SQRT1_2 );
        
        if ( col1.visible() )
        {
            float U = R * T;
            col1.load_load();
            drawRectangle(-U, -U, U, U, Z);
        }
        
        col2.load_load();
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
                drawRectangle( X-H, Y-H, X+H, Y+H, Z);
                drawRectangle(-X+H,-Y+H,-X-H,-Y-H, Z);
            }
            for ( int i = copyparity(Q,y); i <= x; i+=2 )
            {
                float X = y * T;
                float Y = i * T;
                drawRectangle( X-H, Y-H, X+H, Y+H, Z);
                drawRectangle(-X+H,-Y+H,-X-H,-Y-H, Z);
            }
            for ( int i = copyparity(Q,y); i <= x; i+=2 )
            {
                float X = y * T;
                float Y = i * T;
                drawRectangle(-X-H, Y-H,-X+H, Y+H, Z);
                drawRectangle( X+H,-Y+H, X-H,-Y-H, Z);
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
            gle::longCone();
            glPopMatrix();
        }
        // display a white ball at the origin
        gle_color(1.0, 1.0, 1.0, 1.0).load_load();
        glPushMatrix();
        scale(R);
        gle::sphere4();
        glPopMatrix();
    }
    
    //-----------------------------------------------------------------------
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
        GLenum glError = glGetError();
        while ( glError != GL_NO_ERROR )
        {
            fprintf(out, "OpenGL error `%s' %s\n", errorString(glError), msg);
            glError = glGetError();
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

}

