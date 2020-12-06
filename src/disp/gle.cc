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


namespace gle
{    
    /// values of sinus over a full circle
    GLfloat si_[ncircle+1] = { 0 };

    /// values of cosinus over a full circle
    GLfloat co_[ncircle+1] = { 0 };
    
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

    // Fast method to calculate cosine and sinus over the entire circle
    /**
     C[] and S[] should be preallocated to hold 'cnt+1' values.
     Fill in with a counter-clockwise circle starting at angle `start`
    */
    void circle(size_t cnt, GLfloat C[], GLfloat S[], double rad, double start)
    {
        const double theta = 2.0 * M_PI / (double)cnt;
        const double c = std::cos(theta);
        const double s = std::sin(theta);
    
        double t;
        double x = rad * std::cos(start);
        double y = rad * std::sin(start);
        
        for( size_t n = 0; n < cnt; ++n )
        {
            C[n] = GLfloat(x);
            S[n] = GLfloat(y);
            //apply the rotation matrix
            t = x;
            x = c * x - s * y;
            y = s * t + c * y;
            //std::clog << n << " " << x << " " << y << "\n";
        }
        C[cnt] = C[0];
        S[cnt] = S[0];
    }
    
    void arc(size_t cnt, GLfloat C[], GLfloat S[], double rad,
             double start, double end, GLfloat cenX, GLfloat cenY)
    {
        const double theta = ( end - start ) / (double)cnt;
        const double c = std::cos(theta);
        const double s = std::sin(theta);
    
        double t;
        double x = rad * std::cos(start);
        double y = rad * std::sin(start);
        
        for( size_t n = 0; n <= cnt; ++n )
        {
            C[n] = GLfloat(x) + cenX;
            S[n] = GLfloat(y) + cenY;
            //apply the rotation matrix
            t = x;
            x = c * x - s * y;
            y = s * t + c * y;
            //std::clog << n << " " << x << " " << y << "\n";
        }
    }

#ifdef __AVX__
    // This works only if 'GLfloat == float'
    void circle(size_t cnt, float CS[], double rad, double start)
    {
        const double theta = 2.0 * M_PI / (double)cnt;
        const double c = std::cos(theta);
        const double s = std::sin(theta);
        const double c2 = c * c - s * s;
        const double s2 = c * s + c * s;

        vec4 cs{ c2, s2,  c2, s2};
        vec4 sc{-s2, c2, -s2, c2};
        
        const double x0 = rad * std::cos(start);
        const double y0 = rad * std::sin(start);
        vec4 pp{x0, y0, c*x0-s*y0, s*x0+c*y0};
        
        GLfloat * ptr = CS;
        GLfloat * const end = CS + 2 * cnt;
        while ( ptr < end )
        {
            store4f(ptr, pp);
            ptr += 4;
            // apply the rotation matrix
            // x = c * x - s * y;
            // y = s * y + c * y;
            pp = add4(mul4(cs, duplo4(pp)), mul4(sc, duphi4(pp)));
        }
        end[0] = x0;
        end[1] = y0;
    }
#else
    void circle(size_t cnt, GLfloat CS[], double rad, double start)
    {
        const double theta = 2.0 * M_PI / (double)cnt;
        const double c = std::cos(theta);
        const double s = std::sin(theta);
    
        double t;
        double x = rad * std::cos(start);
        double y = rad * std::sin(start);
        
        for( size_t n = 0; n < cnt; ++n )
        {
            CS[2*n  ] = GLfloat(x);
            CS[2*n+1] = GLfloat(y);
            //apply the rotation matrix
            t = x;
            x = c * x - s * y;
            y = s * t + c * y;
            //std::clog << n << " " << x << " " << y << "\n";
        }
        CS[2*cnt  ] = CS[0];
        CS[2*cnt+1] = CS[1];
    }
#endif
    
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
        circle(ncircle, co_, si_, 1);
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
        float N = invsqrt(X*X+Y*Y);
        X *= N;
        Y *= N;
        //warning! this matrix appears here transposed
        float mat[16] = {
            Y*R,  -X*R,  0,  0,
            0,       0, -R,  0,
            X,       Y,  0,  0,
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
        glBegin(GL_LINE_LOOP);
        for( size_t n = 0; n <= ncircle; ++n )
            glVertex2f(co_[n], si_[n]);
        glEnd();
    }
    
    void discUp()
    {
        glNormal3f(0, 0, 1);
        glBegin(GL_TRIANGLE_FAN);
        glVertex2f(0, 0);
        for( size_t n = 0; n <= ncircle; n+=4 )
            glVertex2f(co_[n], si_[n]);
        glEnd();
    }
    
    void discDown()
    {
        glNormal3f(0, 0, -1);
        glBegin(GL_TRIANGLE_FAN);
        glVertex2f(0, 0);
        for( size_t n = 0; n <= ncircle; n+=4 )
            glVertex2f(si_[n], co_[n]);
        glEnd();
    }
    
    void disc2()
    {
        glNormal3f(0, 0, 1);
        glBegin(GL_TRIANGLE_FAN);
        glVertex2f(0, 0);
        for( size_t n = 0; n <= ncircle; ++n )
            glVertex2f(co_[n], si_[n]);
        glEnd();
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
    void initTetrahedron(GLint buf1, GLint buf2, GLfloat R=1.0f)
    {
        const GLfloat S = R / M_SQRT3;
        const GLfloat Y = 2.0 * S;
        const GLfloat B = -M_SQRT1_2 * S;
        const GLfloat Z = -3.0 * B;

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
    
    
    /// Cube is make of 12 triangles = 36 vertices
    void initCube(GLint buf1, GLint buf2, GLfloat R=0.5773502692f)
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
    void initOctahedron(GLint buf1, GLint buf2, GLfloat R=1.46459188756f)
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
    void initIcosahedron(GLint buf1, GLint buf2, GLfloat R=1.0f)
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
    void initArrowTail(GLint buf1, GLint buf2, GLfloat R=0.1f,
                       GLfloat B=-0.5f, GLfloat H=-1.5f, GLfloat L=2.0f)
    {
        const GLfloat T = B + L;
        const GLfloat U = H + L;
        const GLfloat C = 0.5f;
        const GLfloat S = M_SQRT3 * 0.5f;
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

    void drawBuffer(GLint buf1, GLint buf2, unsigned cnt, GLenum mode)
    {
        glEnableClientState(GL_VERTEX_ARRAY);
        glBindBuffer(GL_ARRAY_BUFFER, buf1);
        glVertexPointer(3, GL_FLOAT, 0, nullptr);
        glEnableClientState(GL_NORMAL_ARRAY);
        glBindBuffer(GL_ARRAY_BUFFER, buf2);
        glNormalPointer(GL_FLOAT, 0, nullptr);
        glDrawArrays(mode, 0, cnt);
        glDisableClientState(GL_NORMAL_ARRAY);
        glDisableClientState(GL_VERTEX_ARRAY);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

    void initBuffers()
    {
        if ( !glIsBuffer(buf_[0]) )
        {
            glGenBuffers(16, buf_);
            initTetrahedron(buf_[0], buf_[1]);
            initCube(buf_[2], buf_[3]);
            initOctahedron(buf_[4], buf_[5]);
            initIcosahedron(buf_[6], buf_[7]);
            initArrowTail(buf_[8], buf_[9]);
            glBindBuffer(GL_ARRAY_BUFFER, 0);
        }
    }
    
    void tetrahedron() { drawBuffer(buf_[0], buf_[1], 12, GL_TRIANGLES); }
    void cube() { drawBuffer(buf_[2], buf_[3], 36, GL_TRIANGLES); }
    void octahedron() { drawBuffer(buf_[4], buf_[5], 24, GL_TRIANGLES); }
    void icosahedron() { drawBuffer(buf_[6], buf_[7], 60, GL_TRIANGLES); }
    void arrowTail() { drawBuffer(buf_[8], buf_[9], 45, GL_TRIANGLES); }

    //-----------------------------------------------------------------------
#pragma mark - Tubes
    
    void tubeZ(GLfloat B, GLfloat T, int inc)
    {
        assert_true( B <= T );
        glBegin(GL_TRIANGLE_STRIP);
        for( size_t n = 0; n <= ncircle; n += inc )
        {
            glNormal3f(co_[n], si_[n], 0);
            glVertex3f(co_[n], si_[n], T);
            glVertex3f(co_[n], si_[n], B);
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
            glNormal3f(C*co_[n], C*si_[n], S);
            glVertex3f(rT*co_[n], rT*si_[n], T);
            glVertex3f(rB*co_[n], rB*si_[n], B);
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
            GLfloat C = co_[p];
            GLfloat S = si_[p];
            
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
    void initHexTube(GLint buf1, GLint buf2, GLfloat A, GLfloat B)
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
            const GLfloat L = 256.f;
            initTube(tub_[ 0], tub_[ 1], 0, 1, 8);
            initTube(tub_[ 2], tub_[ 3], 0, 1, 4);
            initTube(tub_[ 4], tub_[ 5], 0, 1, 2);
            initTube(tub_[ 6], tub_[ 7], 0, 1, 1);
            initTube(tub_[ 8], tub_[ 9],-4, L, 8);
            initTube(tub_[10], tub_[11],-4, L, 4);
            initTube(tub_[12], tub_[13],-4, L, 2);
            initTube(tub_[14], tub_[15], 0, L, 8);
            initTube(tub_[16], tub_[17], 0, L, 4);
            initTube(tub_[18], tub_[19], 0, L, 2);
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
    void hexTube()  { drawBuffer(tub_[22], tub_[23], 14, GL_TRIANGLE_STRIP); }
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
    void hexTube()  { hexTubeZ(0, 1); }
#endif
 
    void tubeZ(GLfloat za, GLfloat ra, gle_color ca, GLfloat zb, GLfloat rb, gle_color cb)
    {
        glBegin(GL_TRIANGLE_STRIP);
        for( size_t n = 0; n <= ncircle; ++n )
        {
            cb.load_load();
            glNormal3f(   co_[n],    si_[n],  0);
            glVertex3f(rb*co_[n], rb*si_[n], zb);
            ca.load_load();
            glNormal3f(   co_[n],    si_[n],  0);
            glVertex3f(ra*co_[n], ra*si_[n], za);
        }
        glEnd();
    }
    
    void cylinder()
    {
        tubeZ(0, 1, 1);
        glTranslatef(0,0,1);
        discUp();
        glTranslatef(0,0,-1);
        discDown();
    }
    
    
    /// spherocylinder of length L, radius R, centered and aligned with axis Z
    void capsuleZ(GLfloat L, GLfloat R)
    {
        const size_t fin = ncircle >> 2;
        const size_t C = 4;
        //display strips along the side of the volume:
        for ( size_t t = 0; t < 4*fin; t+=C )
        {
            //compute the transverse angles:
            GLfloat cb = co_[t  ], sb = si_[t  ];
            GLfloat ca = co_[t+C], sa = si_[t+C];
            GLfloat cB = R * cb, sB = R * sb;
            GLfloat cA = R * ca, sA = R * sa;
            
            //draw one srip of the oval:
            glBegin(GL_TRIANGLE_STRIP);
            for ( size_t i=0; i <= fin; i+=C )
            {
                GLfloat x = co_[i], y = si_[i];
                glNormal3f(ca*y, sa*y,     x);
                glVertex3f(cA*y, sA*y, L+R*x);
                glNormal3f(cb*y, sb*y,     x);
                glVertex3f(cB*y, sB*y, L+R*x);
            }
            for ( int i=fin; i >= 0; i-=C )
            {
                GLfloat x = -co_[i], y = si_[i];
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
        glEnableClientState(GL_VERTEX_ARRAY);
        glEnableClientState(GL_NORMAL_ARRAY);
        glVertexPointer(3, GL_FLOAT, 0, ico.vertex_data());
        glNormalPointer(GL_FLOAT, 0, ico.vertex_data());
        glDrawElements(GL_TRIANGLES, 3*ico.nb_faces(), GL_UNSIGNED_INT, ico.faces_data());
        glDisableClientState(GL_NORMAL_ARRAY);
        glDisableClientState(GL_VERTEX_ARRAY);
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
            glEnableClientState(GL_VERTEX_ARRAY);
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
            glDisableClientState(GL_VERTEX_ARRAY);
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
            real U = co_[n], R = si_[n];
            real L = co_[n+inc], S = si_[n+inc];
            glBegin(GL_TRIANGLE_STRIP);
            for ( size_t p = 0; p <= ncircle; p += inc )
            {
                glNormal3f(R*co_[p], R*si_[p], U);
                glVertex3f(R*co_[p], R*si_[p], U);
                glNormal3f(S*co_[p], S*si_[p], L);
                glVertex3f(S*co_[p], S*si_[p], L);
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
            GLfloat X0 = co_[n  ], Y0 = si_[n  ];
            GLfloat X1 = co_[n+S], Y1 = si_[n+S];
            glBegin(GL_TRIANGLE_STRIP);
            for ( size_t p = 0; p <= ncircle; p += 2*S )
            {
                glNormal3f(X0*co_[p], Y0*co_[p], si_[p]);
                glVertex3f(X0*(R+T*co_[p]), Y0*(R+T*co_[p]), T*si_[p]);
                glNormal3f(X1*co_[p], Y1*co_[p], si_[p]);
                glVertex3f(X1*(R+T*co_[p]), Y1*(R+T*co_[p]), T*si_[p]);
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
    
    
    void cylinderZ()
    {
        const GLfloat T =  0.5;
        const GLfloat B = -0.5;
        
        glBegin(GL_TRIANGLE_FAN);
        glNormal3f( 0, 0, -1 );
        glVertex3f( 0, 0,  B );
        for ( size_t n = 0; n <= ncircle; ++n )
            glVertex3f(co_[n], -si_[n], B);
        glEnd();
        
        glBegin(GL_TRIANGLE_STRIP);
        for ( size_t n = 0; n <= ncircle; ++n )
        {
            glNormal3f(co_[n], si_[n], 0);
            glVertex3f(co_[n], si_[n], T);
            glVertex3f(co_[n], si_[n], B);
        }
        glEnd();
        
        glBegin(GL_TRIANGLE_FAN);
        glNormal3f( 0, 0, 1 );
        glVertex3f( 0, 0, T );
        for ( size_t n = 0; n <= ncircle; ++n )
            glVertex3f(co_[n], si_[n], T);
        glEnd();
    }
    
    
    void coneZ(GLfloat R, GLfloat B, GLfloat T, bool closed)
    {
        if ( closed )
        {
            glBegin(GL_TRIANGLE_FAN);
            glNormal3f( 0, 0, B );
            glVertex3f( 0, 0, B );
            for ( size_t n = 0; n <= ncircle; ++n )
                glVertex3f(R*co_[n], -R*si_[n], B);
            glEnd();
        }
        glBegin(GL_TRIANGLE_FAN);
        glNormal3f( 0, 0, 1 );
        glVertex3f( 0, 0, T );
        const GLfloat S = -1.f/sqrtf((T-B)*(T-B)+1);
        const GLfloat C = (B-T)*S;
        for ( size_t n = 0; n <= ncircle; ++n )
        {
            glNormal3f(C*co_[n], C*si_[n], S);
            glVertex3f(R*co_[n], R*si_[n], B);
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
                glNormal3f(dN*co_[n], dN*si_[n],-dR);
                glVertex3f(R *co_[n], R *si_[n], Z);
                glVertex3f(R0*co_[n], R0*si_[n], Z0);
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
            glBegin(GL_TRIANGLE_STRIP);
            gleVertex(A+rad*d);
            gleVertex(A-rad*d);
            gleVertex(B+rad*d);
            gleVertex(B-rad*d);
            glEnd();
        }
    }
    
    
    void drawBand(Vector1 const& A, real rA,
                 Vector1 const& B, real rB)
    {
        glBegin(GL_TRIANGLE_STRIP);
        gleVertex(A.XX,+rA);
        gleVertex(A.XX,-rA);
        gleVertex(B.XX,+rB);
        gleVertex(B.XX,-rB);
        glEnd();
    }
    
    void drawBand(Vector2 const& A, real rA,
                 Vector2 const& B, real rB)
    {
        Vector2 d = ( B - A ).orthogonal();
        real n = d.norm();
        if ( n > 0 )
        {
            d /= n;
            glBegin(GL_TRIANGLE_STRIP);
            gleVertex(A-rA*d);
            gleVertex(A+rA*d);
            gleVertex(B-rB*d);
            gleVertex(B+rB*d);
            glEnd();
        }
    }
    
    void drawBand(Vector1 const& a, real ra, gle_color ca,
                 Vector1 const& b, real rb, gle_color cb)
    {
        glBegin(GL_TRIANGLE_STRIP);
        ca.load();
        gleVertex(a.XX,-ra);
        gleVertex(a.XX,+ra);
        cb.load();
        gleVertex(b.XX,-rb);
        gleVertex(b.XX,+rb);
        glEnd();
    }
    
    void drawBand(Vector2 const& a, real ra, gle_color ca,
                 Vector2 const& b, real rb, gle_color cb)
    {
        Vector2 d = ( b - a ).orthogonal();
        real n = d.norm();
        if ( n > 0 )
        {
            d /= n;
            glBegin(GL_TRIANGLE_STRIP);
            ca.load();
            gleVertex(a+ra*d);
            gleVertex(a-ra*d);
            cb.load();
            gleVertex(b+rb*d);
            gleVertex(b-rb*d);
            glEnd();
        }
    }
    
    
    /**
     This will displays a rectangle if the connection is parallel,
     and a hourglass if the connection is antiparallel
     */
    void drawHourglass(Vector2 const& a, Vector2 const& da,
                Vector2 const& b, Vector2 const& db)
    {
        Vector2 pts[6] = { b-db, b, a-da, a+da, b, b+db };
        static_assert(sizeof(pts)==12*sizeof(real), "unexpected size of Vector2");
#if REAL_IS_DOUBLE
        glVertexPointer(2, GL_DOUBLE, 0, pts);
#else
        glVertexPointer(2, GL_FLOAT, 0, pts);
#endif
        glEnableClientState(GL_VERTEX_ARRAY);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 6);
        glDisableClientState(GL_VERTEX_ARRAY);
    }
    
    /**
     This will displays a rectangle if the connection is parallel,
     and a hourglass if the connection is antiparallel
     */
    void drawHourglass(Vector2 const& a, Vector2 const& da, gle_color ca,
                Vector2 const& b, Vector2 const& db, gle_color cb)
    {
        Vector2 pts[6] = { b-db, b, a-da, a+da, b, b+db };
        static_assert(sizeof(pts)==12*sizeof(real), "unexpected size of Vector2");
        GLfloat col[24];
        cb.store(col);
        cb.store(col+4);
        ca.store(col+8);
        ca.store(col+12);
        cb.store(col+16);
        cb.store(col+20);
        glEnableClientState(GL_VERTEX_ARRAY);
#if REAL_IS_DOUBLE
        glVertexPointer(2, GL_DOUBLE, 0, pts);
#else
        glVertexPointer(2, GL_FLOAT, 0, pts);
#endif
        glEnableClientState(GL_COLOR_ARRAY);
        glColorPointer(4, GL_FLOAT, 0, col);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 6);
        glDisableClientState(GL_COLOR_ARRAY);
        glDisableClientState(GL_VERTEX_ARRAY);
    }
    
    
    /**
     This will displays a rectangle if the connection is antiparallel,
     and a hourglass if the connection is parallel
     */
    void drawCross(Vector2 const& a, Vector2 const& da,
                  Vector2 const& b, Vector2 const& db, real rad)
    {
        glLineWidth(0.5);
        glBegin(GL_TRIANGLE_FAN);
        gleVertex(a);
        gleVertex(a-rad*da);
        gleVertex(b);
        gleVertex(b-rad*db);
        glEnd();
        glBegin(GL_TRIANGLE_FAN);
        gleVertex(a);
        gleVertex(a+rad*da);
        gleVertex(b);
        gleVertex(b+rad*db);
        glEnd();
    }
    
    void drawBar(Vector3 const& a, Vector3 const& da,
                Vector3 const& b, Vector3 const& db, real rad)
    {
        Vector3 ab = normalize( a - b );
        Vector3 ea = cross(ab, da);
        Vector3 eb = cross(ab, db);
        glBegin(GL_TRIANGLE_STRIP);
        gleVertex(a-rad*(da-ea));
        gleVertex(a-rad*(da+ea));
        gleVertex(b-rad*(db-eb));
        gleVertex(b-rad*(db+eb));
        glEnd();
        glBegin(GL_TRIANGLE_STRIP);
        gleVertex(a+rad*(da-ea));
        gleVertex(a+rad*(da+ea));
        gleVertex(b+rad*(db-eb));
        gleVertex(b+rad*(db+eb));
        glEnd();
        glBegin(GL_TRIANGLE_STRIP);
        gleVertex(a-rad*da);
        gleVertex(a+rad*da);
        gleVertex(b-rad*db);
        gleVertex(b+rad*db);
        glEnd();
        glBegin(GL_TRIANGLE_STRIP);
        gleVertex(a-rad*da);
        gleVertex(a+rad*da);
        gleVertex(b-rad*db);
        gleVertex(b+rad*db);
        glEnd();
    }
    
    
    /**
     Two hexagons linked by a rectangle
     hexagons have the same surface as a disc of radius 1.
     */
    void drawDumbbell(Vector2 const& A, Vector2 const& B, GLfloat diameter)
    {
        const GLfloat S = 1.0996361107912678f; //sqrt( 2 * M_PI / ( 3 * sqrt(3) ));
        const GLfloat R = diameter * S;
        const GLfloat H = R * 0.8660254037844386f; //0.5f * sqrt(3);
        const GLfloat X = R * 0.5f;
        
        Vector2 x = ( B - A ).normalized(H);
        Vector2 y = x.orthogonal(X);
        
        glPushMatrix();
        translate(A);
        
        // this is an hexagon centered around 'a':
        glBegin(GL_TRIANGLE_FAN);
        glVertex2f(0,0);
        gleVertex(x+y);
        gleVertex(2*y);
        gleVertex(-x+y);
        gleVertex(-x-y);
        gleVertex(-2*y);
        gleVertex(x-y);
        gleVertex(x+y);
        glEnd();
        
        // a band from 'a' to 'b'
        glBegin(GL_TRIANGLE_FAN);
        gleVertex(+y+x);
        gleVertex(-y+x);
        gleVertex(B-A-y-x);
        gleVertex(B-A+y-x);
        glEnd();
        
        // an hexagon centered around 'b'
        translate(B-A);
        glBegin(GL_TRIANGLE_FAN);
        glVertex2f(0,0);
        gleVertex(x+y);
        gleVertex(2*y);
        gleVertex(-x+y);
        gleVertex(-x-y);
        gleVertex(-2*y);
        gleVertex(x-y);
        gleVertex(x+y);
        glEnd();
        
        glPopMatrix();
    }
    
    //-----------------------------------------------------------------------
#pragma mark - Arrows
    
    void drawCone(Vector1 const& pos, Vector1 const& dir, const GLfloat rad)
    {
        GLfloat dx = rad*dir.XX, cx = pos.XX;
        glBegin(GL_TRIANGLE_STRIP);
        gleVertex(cx-dx   , dx);
        gleVertex(cx-dx/2 , 0 );
        gleVertex(cx+dx+dx, 0 );
        gleVertex(cx-dx   ,-dx);
        glEnd();
    }
    
    void drawCone(Vector2 const& pos, Vector2 const& dir, const GLfloat rad)
    {
        GLfloat dx = rad*dir.XX,  cx = pos.XX;
        GLfloat dy = rad*dir.YY,  cy = pos.YY;
        glBegin(GL_TRIANGLE_STRIP);
        gleVertex(cx-dx-dy, cy-dy+dx);
        gleVertex(cx-dx/2,  cy-dy/2 );
        gleVertex(cx+dx+dx, cy+dy+dy);
        gleVertex(cx-dx+dy, cy-dy-dx);
        glEnd();
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
        real cx = pos.XX;
        glBegin(GL_TRIANGLE_STRIP);
        real dx = rad * dir.XX / 2;
        gleVertex( cx-dx, -rad );
        gleVertex( cx-dx,  rad );
        gleVertex( cx+dx, -rad );
        gleVertex( cx+dx,  rad );
        glEnd();
    }
    
    void drawCylinder(Vector2 const& pos, Vector2 const& dir, float rad)
    {
        real dx = rad * dir.XX, cx = pos.XX - dx / 2;
        real dy = rad * dir.YY, cy = pos.YY - dy / 2;
        glBegin(GL_TRIANGLE_STRIP);
        gleVertex( cx+dy, cy-dx );
        gleVertex( cx-dy, cy+dx );
        gleVertex( cx+dx+dy, cy+dy-dx );
        gleVertex( cx+dx-dy, cy+dy+dx );
        glEnd();
    }
    
    void drawCylinder(Vector3 const& pos, Vector3 const& dir, float rad)
    {
        glPushMatrix();
        transAlignZ(pos, rad, dir);
        gle::cylinderZ();
        glPopMatrix();
    }
    
    
    //-----------------------------------------------------------------------
    
    void drawArrowTail(Vector1 const& pos, Vector1 const& dir, float rad)
    {
        GLfloat dx = rad * dir.XX;
        GLfloat cx = pos.XX - dx / 2;
        glBegin(GL_TRIANGLE_FAN);
        glVertex2f( cx,       0  );
        glVertex2f( cx-dx,   -dx );
        glVertex2f( cx+dx,   -dx );
        glVertex2f( cx+dx+dx, 0  );
        glVertex2f( cx+dx,    dx );
        glVertex2f( cx-dx,    dx );
        glEnd();
    }
    
    void drawArrowTail(Vector2 const& pos, Vector2 const& dir, float rad)
    {
        GLfloat dx = rad * dir.XX;
        GLfloat dy = rad * dir.YY;
        GLfloat cx = pos.XX - 1.5f * dx;
        GLfloat cy = pos.YY - 1.5f * dy;
        GLfloat ex = cx + 2 * dx;
        GLfloat ey = cy + 2 * dy;
        
        glBegin(GL_TRIANGLE_FAN);
        glVertex2f( cx+dx, cy+dy );
        glVertex2f( cx+dy, cy-dx );
        glVertex2f( ex+dy, ey-dx );
        glVertex2f( ex+dx, ey+dy );
        glVertex2f( ex-dy, ey+dx );
        glVertex2f( cx-dy, cy+dx );
        glEnd();
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
        glScaled(3.0, 3.0, 3*R);
        gle::longCone();
        glPopMatrix();
    }
    
    void drawArrow(Vector2 const& A, Vector2 const& B, float R)
    {
        glPushMatrix();
        stretchAlignZ(A, B, R);
        tube1();
        glTranslatef(0, 0, 1);
        glScaled(3.0, 3.0, 3*R);
        gle::longCone();
        glPopMatrix();
    }
    
    void drawArrow(Vector3 const& A, Vector3 const& B, float R)
    {
        glPushMatrix();
        stretchAlignZ(A, B, R);
        tube1();
        glTranslatef(0, 0, 1);
        glScaled(3.0, 3.0, 3*R);
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
            glBegin(GL_TRIANGLE_FAN);
            drawNiceRectangle(rec, 3);
            glEnd();
            
            glPopAttrib();
            
            if ( position == 4 )
            {
                glLineWidth(0.5);
                glBegin(GL_LINE_STRIP);
                drawNiceRectangle(rec, 3);
                glEnd();
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
        glVertex2i(rec[0], rec[1]);
        glVertex2i(rec[2], rec[1]);
        glVertex2i(rec[2], rec[3]);
        glVertex2i(rec[0], rec[3]);
        glVertex2i(rec[0], rec[1]);
    }
    
    
    void drawNiceRectangle(const int rec[4], const int rad)
    {
        glVertex2i(rec[0], rec[1]+rad);
        glVertex2i(rec[0]+rad, rec[1]);
        glVertex2i(rec[2]-rad, rec[1]);
        glVertex2i(rec[2], rec[1]+rad);
        glVertex2i(rec[2], rec[3]-rad);
        glVertex2i(rec[2]-rad, rec[3]);
        glVertex2i(rec[0]+rad, rec[3]);
        glVertex2i(rec[0], rec[3]-rad);
        glVertex2i(rec[0], rec[1]+rad);
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
        glBegin(GL_LINES);
        glVertex2i(16, 1);    glVertex2i(1, 16);
        glVertex2i(12, 1);    glVertex2i(1, 12);
        glVertex2i(8,  1);    glVertex2i(1,  8);
        glVertex2i(4,  1);    glVertex2i(1,  4);
        glEnd();
        
        glPopMatrix();
        glMatrixMode(GL_MODELVIEW);
        glPopMatrix();
    }
    

    //-----------------------------------------------------------------------
    
    void drawRectangle(float X1, float Y1, float X2, float Y2, float Z)
    {
        glBegin(GL_TRIANGLE_STRIP);
        glVertex3f(X1, Y1, Z);
        glVertex3f(X2, Y1, Z);
        glVertex3f(X1, Y2, Z);
        glVertex3f(X2, Y2, Z);
        glEnd();
    }
    
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
        const GLfloat S = GLfloat(size);
        const GLfloat R = S * 0.1f;
        
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

