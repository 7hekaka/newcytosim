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
#include "flute.h"

namespace gle
{
    /// function callback
    using drawCall = void (*)(GLsizei, GLfloat*, GLuint);

    /// values of cosinus, sinus over a full circle
    GLfloat cir_[2*ncircle+8] = { 0 };
    
    /// OpenGL buffers objects for streaming
    GLuint stream_[4] = { 0 };
    
    /// vertex buffer objects for tubes
    GLuint tub_[16] = { 0 };
    
    /// vertex buffer objects for hex tubes
    GLuint buf_[16] = { 0 };
    
    /// vertex buffer objects for icosahedrons
    GLuint ico_[16] = { 0 };
    
    /// number of faces in icosahedrons
    GLsizei ico_ntri[8] = { 0 };

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
        
        if ( !glIsBuffer(ico_[0]) )
        {
            initSphereBuffers();
            CHECK_GL_ERROR("gle:initSphereBuffers()");
        }
        if ( !glIsBuffer(tub_[0]) )
        {
            initTubeBuffers();
            CHECK_GL_ERROR("gle:initTubeBuffers()");
        }
        if ( !glIsBuffer(buf_[0]) )
        {
            initBuffers();
            CHECK_GL_ERROR("gle:initBuffers()");
        }
        if ( !glIsBuffer(stream_[0]) )
            glGenBuffers(4, stream_);
        CHECK_GL_ERROR("glGenBuffers(4, stream_)");
        std::atexit(release);
    }
    
    void release()
    {
        glDeleteBuffers(16, ico_);
        ico_[0] = 0;
        glDeleteBuffers(12, tub_);
        tub_[0] = 0;
        glDeleteBuffers(16, buf_);
        buf_[0] = 0;
        glDeleteBuffers(4, stream_);
        stream_[0] = 0;
    }

    //-----------------------------------------------------------------------
    #pragma mark - Compute Arc and Circle
    
    /// access to precomputed cosinus
    inline GLfloat cos_(size_t n) { return cir_[4+2*n]; }
    
    /// access to precomputed sinus
    inline GLfloat sin_(size_t n) { return cir_[5+2*n]; }

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
    
    template < typename FLOAT >
    inline void orthonormal(const FLOAT v[3], FLOAT mul, FLOAT x[3], FLOAT y[3])
    {
        const FLOAT s = std::copysign((FLOAT)1, v[2]);
        /// optimized version by Marc B. Reynolds
        const FLOAT a = v[1] / ( v[2] + s );
        const FLOAT b = v[1] * a;
        const FLOAT c = v[0] * a;
        FLOAT sul = s * mul;
        x[0] = mul * ( -v[2] - b );
        x[1] = mul * c;
        x[2] = mul * v[0];
        y[0] = sul * c;
        y[1] = sul * b - mul;
        y[2] = sul * v[1];
    }

    template < typename FLOAT >
    inline void orthonormal(const FLOAT v[3], FLOAT mul, FLOAT x[3], FLOAT y[3], FLOAT z[3])
    {
        const FLOAT s = std::copysign((FLOAT)1, v[2]);
        /// optimized version by Marc B. Reynolds
        const FLOAT a = v[1] / ( v[2] + s );
        const FLOAT b = v[1] * a;
        const FLOAT c = v[0] * a;
        FLOAT sul = s * mul;
        x[0] = mul * ( -v[2] - b );
        x[1] = mul * c;
        x[2] = mul * v[0];
        y[0] = sul * c;
        y[1] = sul * b - mul;
        y[2] = sul * v[1];
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
        float n = invsqrt(X*X+Y*Y+Z*Z);
        float vec[3] = { n*X, n*Y, n*Z };
        float mat[16] = {
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            float(P.XX), float(P.YY), float(P.ZZ), 1};
        orthonormal(vec, R, mat, mat+4, mat+8);
        glMultMatrixf(mat);
    }

    
    // rotate to align Z with 'D' and translate to center 'P'
    void transAlignZ(Vector1 const& P, float R, Vector1 const& D, float S)
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
    
    // rotate to align Z with 'D' and translate to center 'P'
    void transAlignZ(Vector2 const& P, float R, Vector2 const& D, float S)
    {
        float X = float(D.XX);
        float Y = float(D.YY);
        float n = invsqrt(X*X+Y*Y);
        X *= n;
        Y *= n;
        float mat[16] = {
            R*Y, -R*X,  0,  0,
            0,      0, -R,  0,
            S*X,  S*Y,  0,  0,
            float(P.XX), float(P.YY), 0, 1};
        glMultMatrixf(mat);
    }
    
    // rotate to align Z with 'D' and translate to center 'P'
    void transAlignZ(Vector3 const& P, float R, Vector3 const& D, float S)
    {
        float X = float(D.XX);
        float Y = float(D.YY);
        float Z = float(D.ZZ);
        float n = invsqrt(X*X+Y*Y+Z*Z);
        float vec[3] = { n*X, n*Y, n*Z };
        float mat[16] = {
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0,
            float(P.XX), float(P.YY), float(P.ZZ), 1};
        orthonormal(vec, R, mat, mat+4, mat+8);
        n = S / R;
        mat[8] *= n;
        mat[9] *= n;
        mat[10] *= n;
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
    Vector3 directionDepth()
    {
        GLfloat mat[16];
        glGetFloatv(GL_MODELVIEW_MATRIX, mat);
        return normalize(Vector3(mat[2], mat[6], mat[10]));
    }

    //-----------------------------------------------------------------------
    # pragma mark - Circle, Cone and Disc
    
    size_t setCone(flute6* flu, GLfloat B, GLfloat T, size_t inc)
    {
        const GLfloat L = T - B;
        const GLfloat Y = 1.f/sqrtf(L*L+1.f);
        const GLfloat X = Y * L;
        GLsizei i = 1;
        flu[0] = flute6{ 0, 0, T, 0, 0, 1 };
        for( size_t p = 0; p <= ncircle; p += inc )
        {
            GLfloat C = cos_(p), S = sin_(p);
            flu[i++] = flute6{ C, S, B, C*X, S*X, Y };
        }
        return i;
    }

    void circle()
    {
        glNormal3f(0, 0, 1);
        glVertexPointer(2, GL_FLOAT, 0, 4+cir_);
        glDrawArrays(GL_LINE_LOOP, 0, ncircle);
    }
    
    void disc()
    {
        glNormal3f(0, 0, 1);
        glVertexPointer(2, GL_FLOAT, 16, cir_);
        glDrawArrays(GL_TRIANGLE_FAN, 0, 2+ncircle/2);
    }

    void disc2()
    {
        glNormal3f(0, 0, 1);
        glVertexPointer(2, GL_FLOAT, 0, 2+cir_);
        glDrawArrays(GL_TRIANGLE_FAN, 0, 2+ncircle);
    }
    
    // this does not preserve the modelview transformation
    void discTop()
    {
        glNormal3f(0, 0, 1);
        glTranslatef(0, 0, 1);
        glVertexPointer(2, GL_FLOAT, 16, cir_);
        glDrawArrays(GL_TRIANGLE_FAN, 0, 2+ncircle/2);
        //glTranslatef(0, 0, -1);
    }

    void discBottom()
    {
        glNormal3f(0, 0, -1);
        glFrontFace(GL_CW);
        glVertexPointer(2, GL_FLOAT, 16, cir_);
        glDrawArrays(GL_TRIANGLE_FAN, 0, 2+ncircle/2);
        glFrontFace(GL_CCW);
    }
    
    void discBottom2()
    {
        glNormal3f(0, 0, -1);
        glFrontFace(GL_CW);
        glVertexPointer(2, GL_FLOAT, 0, 2+cir_);
        glDrawArrays(GL_TRIANGLE_FAN, 0, 2+ncircle);
        glFrontFace(GL_CCW);
    }
    
    void discZ(GLfloat R, GLfloat Z, GLfloat N)
    {
        glBegin(GL_TRIANGLE_FAN);
        glNormal3f( 0, 0, N );
        glVertex3f( 0, 0, Z );
        for ( size_t n = 0; n <= ncircle; ++n )
            glVertex3f(R*cos_(n), R*sin_(n), Z);
        glEnd();
    }
    
    void coneZ(GLfloat R, GLfloat B, GLfloat T)
    {
        glBegin(GL_TRIANGLE_FAN);
        glNormal3f( 0, 0, 1 );
        glVertex3f( 0, 0, T );
        const GLfloat L(T-B);
        const GLfloat Y = 1.f/sqrtf(L*L+1);
        const GLfloat X = Y * L;
        for ( size_t n = 0; n <= ncircle; ++n )
        {
            GLfloat S = sin_(n), C = cos_(n);
            glNormal3f(X*C, X*S, Y);
            glVertex3f(R*C, R*S, B);
        }
        glEnd();
    }

    void cone()
    {
        coneZ(1, 0, 1);
        discBottom2();
    }

    void longCone()
    {
        glTranslatef(0, 0, -1);
        glScalef(1, 1, 3);
        cone();
    }

    void shortCone()
    {
        glTranslatef(0, 0, -0.333f);
        glScalef(1, 1, 0.5);
        cone();
    }

    void truncatedCone()
    {
        tubeZ(0, 1, 1, 0.25, 4);
    }

    //-----------------------------------------------------------------------
    #pragma mark - Some Platonic solids
    
    /// code to calculate and print normals
    void printNormals(size_t cnt, GLfloat const* pts)
    {
        for ( size_t i = 0; i < cnt; ++i )
        {
            GLfloat const* x = pts + 9 * i;
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

    /// Tetrahedron is make of 4 triangles = 12 vertices
    void drawTetrahedron(drawCall func1, drawCall func2, GLuint bufN, GLuint bufP, GLfloat R=1.2f)
    {
        const GLfloat S = R / M_SQRT3;
        const GLfloat Y = 2 * S;
        const GLfloat B = -M_SQRT1_2 * S;
        const GLfloat Z = -3 * B;

        // -R,-S, B
        // +R,-S, B
        //  0, Y, B
        //  0, 0, Z
        GLfloat pts[9*4] = {
             R,-S, B,-R,-S, B, 0, Y, B,
             R,-S, B, 0, Y, B, 0, 0, Z,
             0, Y, B,-R,-S, B, 0, 0, Z,
            -R,-S, B, R,-S, B, 0, 0, Z,
        };
        //printNormals(4, pts);
        GLfloat dir[9*4] = {
            +0.00000, 0.00000,-1.00000,+0.00000, 0.00000,-1.00000,+0.00000, 0.00000,-1.00000,
            +0.81650, 0.47140, 0.33333, 0.81650, 0.47140, 0.33333, 0.81650, 0.47140, 0.33333,
            -0.81650, 0.47140, 0.33333,-0.81650, 0.47140, 0.33333,-0.81650, 0.47140, 0.33333,
            +0.00000,-0.94281, 0.33333, 0.00000,-0.94281, 0.33333, 0.00000,-0.94281, 0.33333
        };
        assert_true( sizeof(pts) == sizeof(dir) );
        func1(sizeof(dir), dir, bufN);
        func2(sizeof(pts), pts, bufP);
    }
    
    /// The star is made of two Tetrahedrons: 8 triangles = 24 vertices
    void drawStar(drawCall func1, drawCall func2, GLuint bufN, GLuint bufP, GLfloat R=1.2f)
    {
        const GLfloat S = R / M_SQRT3;
        const GLfloat Y = 2 * S;
        const GLfloat B = -M_SQRT1_2 * S;
        const GLfloat Z = -3.0 * B;

        GLfloat pts[9*8] = {
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
        
        GLfloat dir[9*8] = {
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
        func1(sizeof(dir), dir, bufN);
        func2(sizeof(pts), pts, bufP);
    }

    
    /// Cube is make of 12 triangles = 36 vertices
    void drawCube(drawCall func1, drawCall func2, GLuint bufN, GLuint bufP, GLfloat R=0.5773502692f)
    {
        GLfloat pts[9*12] = {
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
        
        GLfloat dir[9*12] = {
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
        func1(sizeof(dir), dir, bufN);
        func2(sizeof(pts), pts, bufP);
    }
    
    /// Octahedron is make of 8 triangles = 24 vertices
    void drawOctahedron(drawCall func1, drawCall func2, GLuint bufN, GLuint bufP, GLfloat R=1.46459188756f)
    {
        // Eight triangles, ordered counterclockwise
        // set size to match the volume of the unit sphere
        GLfloat pts[9*8] = {
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
        GLfloat dir[9*8] = {
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
        func1(sizeof(dir), dir, bufN);
        func2(sizeof(pts), pts, bufP);
    }

    
#if ( 0 )
    // this is used to calculate the vertices of the icosahedron
    void icoFace(GLfloat* a, GLfloat* b, GLfloat* c)
    {
        GLfloat pts[9] = {a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2]};
        GLfloat nx = (a[0]+b[0]+c[0]) / 3.0f;
        GLfloat ny = (a[1]+b[1]+c[1]) / 3.0f;
        GLfloat nz = (a[2]+b[2]+c[2]) / 3.0f;
        GLfloat n = std::sqrt( nx * nx + ny * ny + nz * nz );
        nx /= n;
        ny /= n;
        nz /= n;
        GLfloat nor[9] = {nx, ny, nz, nx, ny, nz, nx, ny, nz};
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
        GLfloat pts[3*12] = {
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
    void drawIcosahedron(drawCall func1, drawCall func2, GLuint bufN, GLuint bufP, GLfloat R=1.0f)
    {
        const GLfloat T = R * 0.8506508084f;      // (1 + sqrt(5))/2
        const GLfloat O = R * 0.5257311121f;      // 1 / sqrt(1+T^2)

        GLfloat pts[9*20] = {
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
        
        GLfloat dir[9*20] = {
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
        func1(sizeof(dir), dir, bufN);
        func2(sizeof(pts), pts, bufP);
    }
    
    /// Three fins similar to the tail of a V2 rocket
    void drawArrowTail(drawCall func1, drawCall func2, GLuint bufN, GLuint bufP, GLfloat R=0.1f,
                       GLfloat B=-0.5f, GLfloat H=-1.5f, GLfloat L=2.0f)
    {
        const GLfloat T = B + L;
        const GLfloat U = H + L;
        const GLfloat C = 0.5f;
        const GLfloat S = M_SQRT3_2;
        const GLfloat cR = R * C;
        const GLfloat sR = R * S;

        GLfloat pts[9*15] = {
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
        
        GLfloat dir[9*15] = {
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
        func1(sizeof(dir), dir, bufN);
        func2(sizeof(pts), pts, bufP);
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
    
    size_t setCube(flute3* flu, GLfloat R)
    {
        const GLfloat U = -R;
        flu[0] = flute3{R, U, U};
        flu[1] = flute3{U, U, U};
        flu[2] = flute3{R, R, U};
        flu[3] = flute3{U, R, U};
        flu[4] = flute3{U, R, R};
        flu[5] = flute3{U, U, U};
        flu[6] = flute3{U, U, R};
        flu[7] = flute3{R, U, U};
        flu[8] = flute3{R, U, R};
        flu[9] = flute3{R, R, U};
        flu[10] = flute3{R, R, R};
        flu[11] = flute3{U, R, R};
        flu[12] = flute3{R, U, R};
        flu[13] = flute3{U, U, R};
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
        const GLfloat pts[] = {
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
            a = flute3{pts[3*n  ], pts[3*n+1], pts[3*n+2]};
            b = flute3{pts[3*n+3], pts[3*n+4], pts[3*n+5]};
            c = flute3{pts[3*n+6], pts[3*n+7], pts[3*n+8]};
            flu[i++] = normalize(a+a);
            flu[i++] = normalize(a+b);
            flu[i++] = normalize(a+c);
            flu[i++] = normalize(b+c);
        }
        a = flute3{pts[3*n+3], pts[3*n+4], pts[3*n+5]};
        flu[i++] = normalize(c+c);
        flu[i++] = normalize(a+c);
        for ( n = 13; n > 2; n -= 2 )
        {
            a = flute3{pts[3*n  ], pts[3*n+1], pts[3*n+2]};
            b = flute3{pts[3*n-3], pts[3*n-2], pts[3*n-1]};
            c = flute3{pts[3*n-6], pts[3*n-5], pts[3*n-4]};
            flu[i++] = normalize(a+a);
            flu[i++] = normalize(a+b);
            flu[i++] = normalize(a+c);
            flu[i++] = normalize(b+c);
        }
        a = flute3{pts[0], pts[1], pts[2]};
        flu[i++] = normalize(c+c);
        flu[i++] = normalize(a+c);
        return i;
    }
    
    /* This moves some vertices to add an hexagonal needle to the blob */
    void modifyBlob(flute3 * flu)
    {
        const GLfloat R = 0.6f;
        const GLfloat Y = R * 0.8660254037844386f; // sqrtf(3)/2;
        const GLfloat X = R * 0.5f;
        const GLfloat Z = 0.6f;
        //{ 1, 3, 5, 7, 9, 11, 40, 41, 43, 45, 47, 49, 51 }
        // pull some vertices far away in Z
        for ( int u : { 42, 44, 46, 48, 50 } ) flu[u] = flute3{ 0, 0, 6 };
        // set the 6 vertex of an hexagon:
        for ( int u : { 3, 49 } ) flu[u] = flute3{ R, 0, Z};
        for ( int u : { 5, 47 } ) flu[u] = flute3{ X,-Y, Z};
        for ( int u : { 7, 45 } ) flu[u] = flute3{-X,-Y, Z};
        for ( int u : { 9, 43 } ) flu[u] = flute3{-R, 0, Z};
        for ( int u : {11, 41 } ) flu[u] = flute3{-X, Y, Z};
        for ( int u : {1,40,51} ) flu[u] = flute3{ X, Y, Z};
    }
    
    void initBlobBuffer(GLuint buf)
    {
        const size_t cnt = 3*64;
        glBindBuffer(GL_ARRAY_BUFFER, buf);
        glBufferData(GL_ARRAY_BUFFER, cnt*sizeof(flute3), nullptr, GL_STATIC_DRAW);
        flute3 * flu = (flute3*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
        setBlob(flu);
        setBlob(flu+64);
        modifyBlob(flu+64);
        setCube(flu+128, 0.5f);
        glUnmapBuffer(GL_ARRAY_BUFFER);
    }
    
    void drawVertexBuffer(GLenum mode, GLuint buf, GLsizei start, GLsizei cnt)
    {
        assert_true(glIsBuffer(buf));
        glBindBuffer(GL_ARRAY_BUFFER, buf);
        glVertexPointer(3, GL_FLOAT, 0, nullptr);
        glEnableClientState(GL_NORMAL_ARRAY);
        glNormalPointer(GL_FLOAT, 0, nullptr);
        glDrawArrays(mode, start, cnt);
        glDisableClientState(GL_NORMAL_ARRAY);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
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
        glDrawArrays(GL_LINE_STRIP, 0, 52);
    }

    //-----------------------------------------------------------------------

    void drawVertexNormalBuffer(GLenum mode, GLuint buf, GLsizei cnt,
                                GLsizei data_size, size_t normal_offset)
    {
        assert_true(glIsBuffer(buf));
        glBindBuffer(GL_ARRAY_BUFFER, buf);
        glVertexPointer(3, GL_FLOAT, data_size, nullptr);
        glEnableClientState(GL_NORMAL_ARRAY);
        glNormalPointer(GL_FLOAT, data_size, (void*)normal_offset);
        glDrawArrays(mode, 0, cnt);
        glDisableClientState(GL_NORMAL_ARRAY);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

    void drawTriangleBuffer(GLuint nor, GLuint pts, GLsizei cnt)
    {
        assert_true(glIsBuffer(nor));
        glBindBuffer(GL_ARRAY_BUFFER, pts);
        glVertexPointer(3, GL_FLOAT, 0, nullptr);
        glEnableClientState(GL_NORMAL_ARRAY);
        glBindBuffer(GL_ARRAY_BUFFER, nor);
        glNormalPointer(GL_FLOAT, 0, nullptr);
        glDrawArrays(GL_TRIANGLES, 0, cnt);
        glDisableClientState(GL_NORMAL_ARRAY);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }
    
    /// prepare to draw normals
    void setNormals(GLsizei, GLfloat* nor, GLuint)
    {
        glEnableClientState(GL_NORMAL_ARRAY);
        glNormalPointer(GL_FLOAT, 0, nor);
    }
    
    /// draw
    void drawTriangles(GLsizei cnt, GLfloat* pts, GLuint)
    {
        glVertexPointer(3, GL_FLOAT, 0, pts);
        glDrawArrays(GL_TRIANGLES, 0, cnt/(3*sizeof(float)));
        glDisableClientState(GL_NORMAL_ARRAY);
    }

    /// load data to GPU
    void setBuffer(GLsizei cnt, GLfloat* ptr, GLuint buf)
    {
        glBindBuffer(GL_ARRAY_BUFFER, buf);
        glBufferData(GL_ARRAY_BUFFER, cnt, ptr, GL_STATIC_DRAW);
    }
    
    void initBuffers()
    {
        glGenBuffers(16, buf_);
        drawTetrahedron(setBuffer, setBuffer, buf_[0], buf_[1]);
        drawOctahedron(setBuffer, setBuffer, buf_[2], buf_[3]);
        drawIcosahedron(setBuffer, setBuffer, buf_[4], buf_[5]);
        drawArrowTail(setBuffer, setBuffer, buf_[6], buf_[7]);
        drawCube(setBuffer, setBuffer, buf_[8], buf_[9]);
        drawStar(setBuffer, setBuffer, buf_[10], buf_[11]);
        initBlobBuffer(buf_[12]);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }
    
    void tetrahedron() { drawTriangleBuffer(buf_[0], buf_[1], 12); }
    void octahedron()  { drawTriangleBuffer(buf_[2], buf_[3], 24); }
    void icosahedron() { drawTriangleBuffer(buf_[4], buf_[5], 60); }
    //void icosahedron() { drawIcosahedron(setNormals, drawTriangles, 0, 0); }
    
    void arrowTail() { drawTriangleBuffer(buf_[6], buf_[7], 45); }
    void cube()      { drawTriangleBuffer(buf_[8], buf_[9], 36); }
    void star()      { drawTriangleBuffer(buf_[10], buf_[11], 24); }
    
    void blob()      { drawVertexBuffer(GL_TRIANGLE_STRIP, buf_[12], 0, 52); }
    void needle()    { drawVertexBuffer(GL_TRIANGLE_STRIP, buf_[12], 64, 52); }
    void smallCube() { drawVertexBuffer(GL_TRIANGLE_STRIP, buf_[12], 128, 14); }

    void loadBlobBuffer() { loadBuffer(buf_[12]); }
    void blobf()     { glDrawArrays(GL_TRIANGLE_STRIP, 0, 52); }
    void needlef()   { glDrawArrays(GL_TRIANGLE_STRIP,64, 52); }

    //-----------------------------------------------------------------------
#pragma mark - Tubes

    size_t setTube(flute6* flu, GLfloat B, GLfloat T, size_t inc)
    {
        assert_true( B <= T );
        GLsizei i = 0;
        for( size_t p = 0; p <= ncircle; p += inc )
        {
            GLfloat C = cos_(p), S = sin_(p);
            flu[i++] = flute6{ C, S, T, C, S, 0 };
            flu[i++] = flute6{ C, S, B, C, S, 0 };
        }
        return i;
    }
    
    size_t setTube(flute6* flu, GLfloat R, GLfloat B, GLfloat T, size_t inc)
    {
        assert_true( B <= T );
        GLsizei i = 0;
        for( size_t p = 0; p <= ncircle; p += inc )
        {
            GLfloat C = cos_(p), S = sin_(p);
            flu[i++] = flute6{ R*C, R*S, T, C, S, 0 };
            flu[i++] = flute6{ R*C, R*S, B, C, S, 0 };
        }
        return i;
    }

    /// hexagon has the same surface as a disc of radius rad.
    size_t setHexTube(flute6* flu, GLfloat rad, GLfloat A, GLfloat B)
    {
        constexpr GLfloat C = 0.8660254037844386f; //std::sqrt(3)/2;
        constexpr GLfloat S = 0.5f;
        const GLfloat R = rad * 1.0996361107912678f; //std::sqrt( 2 * M_PI / ( 3 * std::sqrt(3) ));
        const GLfloat Y = R * C, X = R * S;
        
        flu[0] = flute6{ R, 0, B, 1, 0, 0};
        flu[1] = flute6{ R, 0, A, 1, 0, 0};
        flu[2] = flute6{ X, Y, B, S, C, 0};
        flu[3] = flute6{ X, Y, A, S, C, 0};
        flu[4] = flute6{-X, Y, B,-S, C, 0};
        flu[5] = flute6{-X, Y, A,-S, C, 0};
        flu[6] = flute6{-R, 0, B,-1, 0, 0};
        flu[7] = flute6{-R, 0, A,-1, 0, 0};
        flu[8] = flute6{-X,-Y, B,-S, -C, 0};
        flu[9] = flute6{-X,-Y, A,-S, -C, 0};
        flu[10] = flute6{ X,-Y, B, S,-C, 0};
        flu[11] = flute6{ X,-Y, A, S,-C, 0};
        flu[12] = flute6{ R, 0, B, 1, 0, 0};
        flu[13] = flute6{ R, 0, A, 1, 0, 0};
        return 14;
    }
    
    inline size_t nbTubeTriangles(size_t inc)
    {
        return 2 + ( 2 * ncircle ) / inc;
    }
    
    constexpr size_t nbHexTubeTriangles = 14;
    
    size_t initTubeBuffer(GLuint buf, GLfloat A, GLfloat B, size_t inc)
    {
        size_t cnt = nbTubeTriangles(inc);
        glBindBuffer(GL_ARRAY_BUFFER, buf);
        glBufferData(GL_ARRAY_BUFFER, cnt*sizeof(flute6), nullptr, GL_STATIC_DRAW);
        flute6 * flu = (flute6*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
        size_t n = setTube(flu, A, B, inc);
        glUnmapBuffer(GL_ARRAY_BUFFER);
        assert_true( n <= cnt );
        return n;
    }
    
    size_t initHexTubeBuffer(GLuint buf, GLfloat R, GLfloat A, GLfloat B)
    {
        size_t cnt = nbHexTubeTriangles;
        glBindBuffer(GL_ARRAY_BUFFER, buf);
        glBufferData(GL_ARRAY_BUFFER, cnt*sizeof(flute6), nullptr, GL_STATIC_DRAW);
        flute6 * flu = (flute6*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
        size_t n = setHexTube(flu, R, A, B);
        glUnmapBuffer(GL_ARRAY_BUFFER);
        assert_true( n <= cnt );
        return n;
    }

    inline void drawTubeBuffer(GLuint buf, size_t cnt)
    {
        drawVertexNormalBuffer(GL_TRIANGLE_STRIP, buf, cnt, sizeof(flute6), 3*sizeof(GLfloat));
    }

    void initTubeBuffers()
    {
        glGenBuffers(14, tub_);
        /* The value of T limits the aspect ratio of tubes that can be drawn */
        const GLfloat B = -32.f, T = 256.f;
        initTubeBuffer(tub_[0], 0, 1, 8);
        initTubeBuffer(tub_[1],-0.03125, 1.03125, 4);
        initTubeBuffer(tub_[2], 0, 1.03125, 4);
        initTubeBuffer(tub_[3], 0, 1, 2);
        initTubeBuffer(tub_[4], 0, 1, 1);
        initTubeBuffer(tub_[5], B, T, 8);
        initTubeBuffer(tub_[6], B, T, 4);
        initTubeBuffer(tub_[7], B, T, 2);
        initTubeBuffer(tub_[8], 0, T, 8);
        initTubeBuffer(tub_[9], 0, T, 4);
        initTubeBuffer(tub_[10], 0, T, 2);
        initHexTubeBuffer(tub_[11], 1, 0, 1);
        initHexTubeBuffer(tub_[12], 0.5f, 0, 1);
        initHexTubeBuffer(tub_[13], 0.5f, 0, T);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

    // using Vertex Buffer Objects
    void tube1()     { drawTubeBuffer(tub_[0], nbTubeTriangles(8)); }
    void tube2()     { drawTubeBuffer(tub_[1], nbTubeTriangles(4)); }
    void tube3()     { drawTubeBuffer(tub_[2], nbTubeTriangles(4)); }
    void tube4()     { drawTubeBuffer(tub_[3], nbTubeTriangles(2)); }
    void tube8()     { drawTubeBuffer(tub_[4], nbTubeTriangles(1)); }
    void longTube1() { drawTubeBuffer(tub_[5], nbTubeTriangles(8)); }
    void longTube2() { drawTubeBuffer(tub_[6], nbTubeTriangles(4)); }
    void longTube4() { drawTubeBuffer(tub_[7], nbTubeTriangles(2)); }
    void halfTube1() { drawTubeBuffer(tub_[8], nbTubeTriangles(8)); }
    void halfTube2() { drawTubeBuffer(tub_[9], nbTubeTriangles(4)); }
    void halfTube4() { drawTubeBuffer(tub_[10], nbTubeTriangles(2)); }
    void hexTube()   { drawTubeBuffer(tub_[11], nbHexTubeTriangles); }
    void thinTube()  { drawTubeBuffer(tub_[12], nbHexTubeTriangles); }
    void thinLongTube() { drawTubeBuffer(tub_[13], nbHexTubeTriangles); }

    // this does not preserve the modelview transformation
    void cylinder1()
    {
        tube1();
        discBottom();
        discTop();
    }
    
    // this does not preserve the modelview transformation
    void cylinder2()
    {
        glTranslatef(0,0,-0.5f);
        discBottom();
        tube4();
        discTop();
        //glTranslatef(0,0,0.5f);
    }

    //-----------------------------------------------------------------------
    #pragma mark - Legacy 3D objects

    void tubeZ(GLfloat B, GLfloat T, int inc)
    {
        assert_true( B <= T );
        glBegin(GL_TRIANGLE_STRIP);
        for( size_t n = 0; n <= ncircle; n += inc )
        {
            GLfloat S = sin_(n), C = cos_(n);
            glNormal3f(C, S, 0);
            glVertex3f(C, S, T);
            glVertex3f(C, S, B);
        }
        glEnd();
    }

    void tubeZ(GLfloat bZ, GLfloat rB, GLfloat tZ, GLfloat rT, int inc)
    {
        assert_true( bZ <= tZ );
        const GLfloat H(tZ-bZ);
        const GLfloat N = 1.f/sqrtf(H*H+(rT-rB)*(rT-rB));
        const GLfloat tC = N * (tZ-bZ);
        const GLfloat tS = N * (rB-rT);
        glBegin(GL_TRIANGLE_STRIP);
        for( size_t n = 0; n <= ncircle; n += inc )
        {
            GLfloat S = sin_(n), C = cos_(n);
            glNormal3f(tC*C, tC*S, tS);
            glVertex3f(rT*C, rT*S, tZ);
            glVertex3f(rB*C, rB*S, bZ);
        }
        glEnd();
    }

    void tubeZ(GLfloat B, GLfloat R, gle_color col, GLfloat T, GLfloat D, gle_color loc)
    {
        glBegin(GL_TRIANGLE_STRIP);
        for( size_t n = 0; n <= ncircle; ++n )
        {
            GLfloat S = sin_(n), C = cos_(n);
            col.load_load();
            glNormal3f(  C,   S, 0);
            glVertex3f(D*C, D*S, B);
            loc.load_load();
            glNormal3f(  C,   S, 0);
            glVertex3f(R*C, R*S, T);
        }
        glEnd();
    }
    
    /// draw spherocylinder of radius R, of axis Z with Z in [B, T]
    void capsuleZ(GLfloat B, GLfloat T, GLfloat R)
    {
        const size_t fin = ncircle >> 2;
        const size_t inc = 4;
        //display strips along the side of the volume:
        for ( size_t t = 0; t < 4*fin; t += inc )
        {
            //compute the transverse angles:
            GLfloat cb = cos_(t),     sb = sin_(t);
            GLfloat ca = cos_(t+inc), sa = sin_(t+inc);
            GLfloat cB = R * cb, sB = R * sb;
            GLfloat cA = R * ca, sA = R * sa;
            
            //draw one srip of the oval:
            glBegin(GL_TRIANGLE_STRIP);
            for ( size_t i=0; i <= fin; i += inc )
            {
                GLfloat x = cos_(i), y = sin_(i);
                glNormal3f(ca*y, sa*y, x);
                glVertex3f(cA*y, sA*y, T+R*x);
                glNormal3f(cb*y, sb*y, x);
                glVertex3f(cB*y, sB*y, T+R*x);
            }
            for ( int i=fin; i >= 0; i -= inc )
            {
                GLfloat x = -cos_(i), y = sin_(i);
                glNormal3f(ca*y, sa*y, x);
                glVertex3f(cA*y, sA*y, B+R*x);
                glNormal3f(cb*y, sb*y, x);
                glVertex3f(cB*y, sB*y, B+R*x);
            }
            glEnd();
        }
    }
    
    /// draw a Torus of radius R and a thickness 2*T
    void torusZ(GLfloat R, GLfloat T, size_t inc)
    {
        for ( size_t n = 0; n < ncircle; n += inc )
        {
            GLfloat X0 = cos_(n    ), Y0 = sin_(n    );
            GLfloat X1 = cos_(n+inc), Y1 = sin_(n+inc);
            glBegin(GL_TRIANGLE_STRIP);
            for ( size_t p = 0; p <= ncircle; p += 2*inc )
            {
                GLfloat S = sin_(p), C = cos_(p);
                glNormal3f(X0*C, Y0*C, S);
                glVertex3f(X0*(R+T*C), Y0*(R+T*C), T*S);
                glNormal3f(X1*C, Y1*C, S);
                glVertex3f(X1*(R+T*C), Y1*(R+T*C), T*S);
            }
            glEnd();
        }
    }
    
    void ellipseZ(GLfloat rX, GLfloat rY, GLfloat rZ)
    {
        GLfloat iX(1.f/rX), iY(1.f/rY), iZ(1.f/rZ);
        /*
         A vector orthogonal to the ellipse surface at position (X, Y, Z) is
         ( X / rX^2, Y / rY^2, Z / rZ^2 )
          */
        for ( size_t n = 0; n < ncircle/2; ++n )
        {
            GLfloat uC = cos_(n  ), uS = sin_(n  );
            GLfloat lC = cos_(n+1), lS = sin_(n+1);
            GLfloat uX = uS * rX, uY = uS * rY, uZ = uC * rZ;
            GLfloat lX = lS * rX, lY = lS * rY, lZ = lC * rZ;
            GLfloat xu = uS * iX, yu = uS * iY, zu = uC * iZ;
            GLfloat xl = lS * iX, yl = lS * iY, zl = lC * iZ;
            glBegin(GL_TRIANGLE_STRIP);
            for ( size_t p = 0; p <= ncircle; ++p )
            {
                GLfloat S = sin_(p), C = cos_(p);
                glNormal3f(C*xu, S*yu, zu);
                glVertex3f(C*uX, S*uY, uZ);
                glNormal3f(C*xl, S*yl, zl);
                glVertex3f(C*lX, S*lY, lZ);
            }
            glEnd();
        }
    }

    void ellipse_circleZ(GLfloat rX, GLfloat rY, GLfloat rZ, GLfloat u)
    {
        GLfloat R(std::sqrt((1-u)*(1+u)));
        GLfloat iX(1.f/rX), iY(1.f/rY), iZ(u/rZ);
        GLfloat vX(R*rX), vY(R*rY), vZ(u*rZ);
        glBegin(GL_LINE_LOOP);
        for ( size_t n = 0; n <= ncircle; ++n )
        {
            GLfloat S = sin_(n), C = cos_(n);
            glNormal3f(C*iX, S*iY, iZ);
            glVertex3f(C*vX, S*vY, vZ);
        }
        glEnd();
    }

    //-----------------------------------------------------------------------
#pragma mark - Spheres
    
    GLuint initIcoBuffer(GLuint buf1, GLuint buf2, Platonic::Solid & ico)
    {
        constexpr size_t DOUZE = 3*sizeof(float);
        //std::clog << "initializeIco ico " << ico.nb_faces() << '\n';
        glBindBuffer(GL_ARRAY_BUFFER, buf1);
        glBufferData(GL_ARRAY_BUFFER, DOUZE*ico.nb_vertices(), nullptr, GL_STATIC_DRAW);
        void * glb = glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
        ico.store_vertices((float*)glb);
        glUnmapBuffer(GL_ARRAY_BUFFER);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        // upload indices:
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buf2);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, DOUZE*ico.nb_faces(), ico.faces_data(), GL_STATIC_DRAW);
        //fprintf(stderr, "icosahedron buffer %i has %u faces\n", buf1, ico.nb_faces());
        return ico.nb_faces();
    }
    
    void initSphereBuffers()
    {
        /// using icosahedrons to render the sphere:
        Platonic::Solid ico1(Platonic::Solid::ICOSAHEDRON, gle::finesse/6);
        Platonic::Solid ico2(Platonic::Solid::ICOSAHEDRON, gle::finesse/3);
        Platonic::Solid ico4(Platonic::Solid::ICOSAHEDRON, gle::finesse);
        Platonic::Solid ico8(Platonic::Solid::ICOSAHEDRON, gle::finesse*2);
        Platonic::Solid icoH1(Platonic::Solid::HEMISPHERE, gle::finesse/6);
        Platonic::Solid icoH2(Platonic::Solid::HEMISPHERE, gle::finesse/3);
        Platonic::Solid icoH4(Platonic::Solid::HEMISPHERE, gle::finesse);
        glGenBuffers(16, ico_);
        ico_ntri[0] = initIcoBuffer(ico_[0], ico_[1], ico1);
        ico_ntri[1] = initIcoBuffer(ico_[2], ico_[3], ico2);
        ico_ntri[2] = initIcoBuffer(ico_[4], ico_[5], ico4);
        ico_ntri[3] = initIcoBuffer(ico_[6], ico_[7], ico8);
        ico_ntri[4] = initIcoBuffer(ico_[8], ico_[9], icoH1);
        ico_ntri[5] = initIcoBuffer(ico_[10], ico_[11], icoH2);
        ico_ntri[6] = initIcoBuffer(ico_[12], ico_[13], icoH4);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }
    
    void drawIcoBuffer(GLsizei nfaces, GLuint buf1, GLuint buf2)
    {
        assert_true(glIsBuffer(buf1));
        glEnableClientState(GL_NORMAL_ARRAY);
        glBindBuffer(GL_ARRAY_BUFFER, buf1);
        glVertexPointer(3, GL_FLOAT, 0, nullptr);
        // the normal in each vertex is equal to the vertex!
        glNormalPointer(GL_FLOAT, 0, nullptr);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buf2);
        glDrawElements(GL_TRIANGLES, 3*nfaces, GL_UNSIGNED_INT, nullptr);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glDisableClientState(GL_NORMAL_ARRAY);
    }
    
    void sphere1() { drawIcoBuffer(ico_ntri[0], ico_[0], ico_[1]); }
    void sphere2() { drawIcoBuffer(ico_ntri[1], ico_[2], ico_[3]); }
    void sphere4() { drawIcoBuffer(ico_ntri[2], ico_[4], ico_[5]); }
    void sphere8() { drawIcoBuffer(ico_ntri[3], ico_[6], ico_[7]); }
    void hemisphere1() { drawIcoBuffer(ico_ntri[4], ico_[8], ico_[9]); }
    void hemisphere2() { drawIcoBuffer(ico_ntri[5], ico_[10], ico_[11]); }
    void hemisphere4() { drawIcoBuffer(ico_ntri[6], ico_[12], ico_[13]); }
    
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
    
    /**
     Draw a cylindrical band on the equator of a sphere of radius 1.
     The band is in the XY plane. The axis of the cylinder is Z.
     The band is made of triangles indicating the clockwise direction.
     */
    void arrowedBand(const size_t nb_triangles, float width)
    {
        GLfloat A(2 * M_PI / nb_triangles);
        GLfloat W(width * A / M_SQRT3);
        GLfloat R(1.0f / cosf(A*0.5f));
        
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
    
    
    void threeBands(const size_t nb_triangles)
    {
        arrowedBand(nb_triangles, 0.25);
        glRotated(-90,1,0,0);
        arrowedBand(nb_triangles, 0.25);
        glRotated(90,0,1,0);
        arrowedBand(nb_triangles, 0.25);
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
                GLfloat S = sin_(n), C = cos_(n);
                glNormal3f(dN*C, dN*S,-dR);
                glVertex3f(R *C, R *S, Z);
                glVertex3f(R0*C, R0*S, Z0);
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
        flute2 pts[8] = { A, A-rad*dA, B, B-rad*dB,
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
        const GLfloat S(1.0996361107912678f); //sqrt( 2 * M_PI / ( 3 * sqrt(3) ));
        const GLfloat R(diameter * S);
        const GLfloat H(R * 0.8660254037844386f); //0.5f * sqrt(3);
        const GLfloat X(R * 0.5f);
        
        Vector2 x = ( B - A ).normalized(H);
        Vector2 y = x.orthogonal(X);
        flute2 pts[20] = {
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
        
        drawRectangle(rec);
        
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

