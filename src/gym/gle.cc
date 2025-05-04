// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University

#include <cmath>
#include <cctype>
#include <cstdlib>
#include "assert_macro.h"
#include "gle.h"
#include "tesselator.h"
#include "simd.h"
#include "simd_float.h"
#include "simd_math.h"
#include "gym_flute.h"
#include "gym_check.h"
#include "gym_view.h"
#include "gym_vect.h"
#include "gym_draw.h"

#include "vector1.h"
#include "vector2.h"
#include "vector3.h"

namespace gle
{
    /// values of cosine and sine over two full revolutions
    float circle_[4*pi_4half+8] = { 0 };
        
    /// vertex buffer objects for static draw
    GLuint buf_[4] = { 0 };

    /// offset for object's data stored in buffers
    unsigned tubes_[32] = { 0 };

    /// offset for object's data with 6 floats / vertex
    unsigned cubes_[10] = { 0 };
    
    /// offset for object's data with 3 floats / vertex
    unsigned blobs_[4] = { 0 };
    
    /// offset for object's data stored in buffers
    unsigned discs_[8] = { 0 };

    /// offset of first vertex coordinate in buffer object
    unsigned ico_pts_[14] = { 0 };
    /// offset, in bytes, of first index in buffer object
    unsigned ico_idx_[14] = { 0 };
    /// number of triangles making up the faces in icosahedrons
    unsigned ico_cnt_[14] = { 0 };
    /// another number of faces
    unsigned ico_mid_[14] = { 0 };

    const unsigned FOOT = 9;
    /// index used for the dome:
    const unsigned DOME = 10;
    const unsigned ROOF = 11;
    /// index used for the icoid object:
    const unsigned ICOID = 12;
    
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
    void set_arc(size_t cnt, float ptr[], double radius,
                 double start, double delta, double cX, double cY)
    {
#ifdef __SSE3__
        return set_arc_SSE(cnt, ptr, radius, start, delta, cX, cY);
#else
        const double C = std::cos(delta);
        const double S = std::sin(delta);

        double x = radius * std::cos(start);
        double y = radius * std::sin(start);
        
        for( size_t n = 0; n < cnt; ++n )
        {
            ptr[  2*n] = float(x + cX);
            ptr[1+2*n] = float(y + cY);
            //apply the rotation matrix
            double t = x;
            x = C * x - S * y;
            y = S * t + C * y;
            //std::clog << n << " " << x << " " << y << "\n";
        }
        ptr[  2*cnt] = float(x + cX);
        ptr[1+2*cnt] = float(y + cY);
#endif
    }
    
    /** Compute cosine and sine over [0, 2*PI],  setting 2 + 4 * cnt pairs of [cosine, sine] values */
    void compute_circle(size_t cnt, float ptr0[], double radius, double start)
    {
        const double delta = M_PI_2 / cnt;
        const double C = std::cos(delta);
        const double S = std::sin(delta);

        double x = radius * std::cos(start);
        double y = radius * std::sin(start);

        float * ptr1 = ptr0 + 2 * cnt;
        float * ptr2 = ptr0 + 4 * cnt;
        float * ptr3 = ptr0 + 6 * cnt;
        
        for( size_t n = 0; n < cnt; ++n )
        {
            float fx = float(x);
            float fy = float(y);
            ptr0[0] = +fx;
            ptr0[1] = +fy;
            ptr0 += 2;
            ptr1[0] = -fy;
            ptr1[1] = +fx;
            ptr1 += 2;
            ptr2[0] = -fx;
            ptr2[1] = -fy;
            ptr2 += 2;
            ptr3[0] = +fy;
            ptr3[1] = -fx;
            ptr3 += 2;
            // apply rotation by angle 'delta':
            double t = x;
            x = C * x - S * y;
            y = S * t + C * y;
        }
        ptr3[0] = float(x);
        ptr3[1] = float(y);
    }
    
    void compute_arc(size_t cnt, float ptr[], double radius, double start,
                     double angle, double cX, double cY)
    {
        set_arc(cnt, ptr, radius, start, angle/(cnt-1), cX, cY);
    }

    //-----------------------------------------------------------------------
    #pragma mark - Some 3D objects: Cubes, Octahedron, Icosahedron...

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
        // top hexagon
        flu[0] = {-R, 0, T, 0, 0, 1};
        flu[1] = {-X,-Y, T, 0, 0, 1};
        flu[2] = {-X, Y, T, 0, 0, 1};
        flu[3] = { X,-Y, T, 0, 0, 1};
        flu[4] = { X, Y, T, S, C, 0};
        // 6 sides
        flu[5] = { R, 0, T, 1, 0, 0};
        flu[6] = { R, 0, B, 1, 0, 0};
        flu[7] = { X,-Y, T, S,-C, 0};
        flu[8] = { X,-Y, B, S,-C, 0};
        flu[9] = {-X,-Y, T,-S,-C, 0};
        flu[10] = {-X,-Y, B,-S,-C, 0};
        flu[11] = {-R, 0, T,-1, 0, 0};
        flu[12] = {-R, 0, B,-1, 0, 0};
        flu[13] = {-X, Y, T,-S, C, 0};
        flu[14] = {-X, Y, B,-S, C, 0};
        flu[15] = { X, Y, T, S, C, 0};
        flu[16] = { X, Y, B, S, C, 0};
        // bottom hexagon
        flu[17] = { R, 0, B, 1, 0, 0};
        flu[18] = {-X, Y, B, 0, 0,-1};
        flu[19] = { X,-Y, B, 0, 0,-1};
        flu[20] = {-R, 0, B, 0, 0,-1};
        flu[21] = {-X,-Y, B, 0, 0,-1};
        return 22;
    }

    /// Tetrahedron is made of 4 triangles and 12 vertices
    size_t setTetrahedron(flute6* flt, float R=1.2f)
    {
        const float S = R / M_SQRT3;
        const float Y = 2 * S;
        const float B = -M_SQRT1_2 * S;
        const float Z = -3 * B;

        const float T = 1 / 3.0;
        const float D = M_SQRT2 / M_SQRT3;
        const float H = 2 * M_SQRT2 / 3.0;
        const float U = M_SQRT2 / 3.0;

        // -R,-S, B
        // +R,-S, B
        //  0, Y, B
        //  0, 0, Z
        flt[0] = { R,-S, B, 0, 0,-1};
        flt[1] = {-R,-S, B, 0, 0,-1};
        flt[2] = { 0, Y, B, 0, 0,-1};
        flt[3] = { R,-S, B, D, U, T};
        flt[4] = { 0, Y, B, D, U, T};
        flt[5] = { 0, 0, Z, D, U, T};
        flt[6] = { 0, Y, B,-D, U, T};
        flt[7] = {-R,-S, B,-D, U, T};
        flt[8] = { 0, 0, Z,-D, U, T};
        flt[9] = {-R,-S, B, 0,-H, T};
        flt[10] = { R,-S, B, 0,-H, T};
        flt[11] = { 0, 0, Z, 0,-H, T};
        return 12;
    }
    
    /// inversed Tetrahedrons by central symmetry
    size_t invTetrahedron(flute6* flt, float R=1.2f)
    {
        const float S = R / M_SQRT3;
        const float Y = 2 * S;
        const float B = -M_SQRT1_2 * S;
        const float Z = -3.0 * B;

        const float T = 1 / 3.0;
        const float D = M_SQRT2 / M_SQRT3;
        const float H = 2 * M_SQRT2 / 3.0;
        const float U = M_SQRT2 / 3.0;

        // reversed tetrahedron by central symmetry
        flt[0] = { R, S,-B, 0, 0, 1};
        flt[1] = {-R, S,-B, 0, 0, 1};
        flt[2] = { 0,-Y,-B, 0, 0, 1};
        flt[3] = { 0,-Y,-B,-D,-U,-T};
        flt[4] = {-R, S,-B,-D,-U,-T};
        flt[5] = { 0, 0,-Z,-D,-U,-T};
        flt[6] = { R, S,-B, D,-U,-T};
        flt[7] = { 0,-Y,-B, D,-U,-T};
        flt[8] = { 0, 0,-Z, D,-U,-T};
        flt[9] = {-R, S,-B, 0, H,-T};
        flt[10] = { R, S,-B, 0, H,-T};
        flt[11] = { 0, 0,-Z, 0, H,-T};
        return 12;
    }

    
    /// A cube made of 12 triangles with normals = 36 vertices
    size_t setCubeFaces(flute6* flt, float X, float Y, float Z)
    {
        // bottom plate
        flt[0] = { X, Y,-Z, 0, 0,-1};
        flt[1] = {-X,-Y,-Z, 0, 0,-1};
        flt[2] = {-X, Y,-Z, 0, 0,-1};
        flt[3] = { X, Y,-Z, 0, 0,-1};
        flt[4] = { X,-Y,-Z, 0, 0,-1};
        flt[5] = {-X,-Y,-Z, 0, 0,-1};
        // top plate
        flt[6] = {-X, Y, Z, 0, 0, 1};
        flt[7] = {-X,-Y, Z, 0, 0, 1};
        flt[8] = { X,-Y, Z, 0, 0, 1};
        flt[9] = { X, Y, Z, 0, 0, 1};
        flt[10] = {-X, Y, Z, 0, 0, 1};
        flt[11] = { X,-Y, Z, 0, 0, 1};
        // sides
        flt[12] = { X, Y, Z, 1, 0, 0};
        flt[13] = { X,-Y,-Z, 1, 0, 0};
        flt[14] = { X, Y,-Z, 1, 0, 0};
        flt[15] = { X,-Y,-Z, 1, 0, 0};
        flt[16] = { X, Y, Z, 1, 0, 0};
        flt[17] = { X,-Y, Z, 1, 0, 0};
        flt[18] = { X, Y, Z, 0, 1, 0};
        flt[19] = { X, Y,-Z, 0, 1, 0};
        flt[20] = {-X, Y,-Z, 0, 1, 0};
        flt[21] = { X, Y, Z, 0, 1, 0};
        flt[22] = {-X, Y,-Z, 0, 1, 0};
        flt[23] = {-X, Y, Z, 0, 1, 0};
        flt[24] = {-X,-Y,-Z,-1, 0, 0};
        flt[25] = {-X,-Y, Z,-1, 0, 0};
        flt[26] = {-X, Y, Z,-1, 0, 0};
        flt[27] = {-X,-Y,-Z,-1, 0, 0};
        flt[28] = {-X, Y, Z,-1, 0, 0};
        flt[29] = {-X, Y,-Z,-1, 0, 0};
        flt[30] = { X,-Y, Z, 0,-1, 0};
        flt[31] = {-X,-Y,-Z, 0,-1, 0};
        flt[32] = { X,-Y,-Z, 0,-1, 0};
        flt[33] = { X,-Y, Z, 0,-1, 0};
        flt[34] = {-X,-Y, Z, 0,-1, 0};
        flt[35] = {-X,-Y,-Z, 0,-1, 0};
        return 36;
    }
    
    /// set triangles for exploded cube where faces are translated outward
    size_t setExplodedCube(flute6* flt, float X, float Y, float Z, float XR, float YR, float ZR)
    {
        // bottom plate
        flt[0] = { X, Y,-ZR, 0, 0,-1};
        flt[1] = {-X,-Y,-ZR, 0, 0,-1};
        flt[2] = {-X, Y,-ZR, 0, 0,-1};
        flt[3] = { X, Y,-ZR, 0, 0,-1};
        flt[4] = { X,-Y,-ZR, 0, 0,-1};
        flt[5] = {-X,-Y,-ZR, 0, 0,-1};
        // top plate
        flt[6] = {-X, Y, ZR, 0, 0, 1};
        flt[7] = {-X,-Y, ZR, 0, 0, 1};
        flt[8] = { X,-Y, ZR, 0, 0, 1};
        flt[9] = { X, Y, ZR, 0, 0, 1};
        flt[10] = {-X, Y, ZR, 0, 0, 1};
        flt[11] = { X,-Y, ZR, 0, 0, 1};
        // sides
        flt[12] = { XR, Y, Z, 1, 0, 0};
        flt[13] = { XR,-Y,-Z, 1, 0, 0};
        flt[14] = { XR, Y,-Z, 1, 0, 0};
        flt[15] = { XR,-Y,-Z, 1, 0, 0};
        flt[16] = { XR, Y, Z, 1, 0, 0};
        flt[17] = { XR,-Y, Z, 1, 0, 0};
        flt[18] = { X, YR, Z, 0, 1, 0};
        flt[19] = { X, YR,-Z, 0, 1, 0};
        flt[20] = {-X, YR,-Z, 0, 1, 0};
        flt[21] = { X, YR, Z, 0, 1, 0};
        flt[22] = {-X, YR,-Z, 0, 1, 0};
        flt[23] = {-X, YR, Z, 0, 1, 0};
        flt[24] = {-XR,-Y,-Z,-1, 0, 0};
        flt[25] = {-XR,-Y, Z,-1, 0, 0};
        flt[26] = {-XR, Y, Z,-1, 0, 0};
        flt[27] = {-XR,-Y,-Z,-1, 0, 0};
        flt[28] = {-XR, Y, Z,-1, 0, 0};
        flt[29] = {-XR, Y,-Z,-1, 0, 0};
        flt[30] = { X,-YR, Z, 0,-1, 0};
        flt[31] = {-X,-YR,-Z, 0,-1, 0};
        flt[32] = { X,-YR,-Z, 0,-1, 0};
        flt[33] = { X,-YR, Z, 0,-1, 0};
        flt[34] = {-X,-YR, Z, 0,-1, 0};
        flt[35] = {-X,-YR,-Z, 0,-1, 0};
        return 36;
    }

    /// Cube is made of 4 linestrips of 4 points each
    size_t setCubeEdges(flute3* flt, float X, float Y, float Z)
    {
        // vertical edges
        flt[0] = { X, Y, -Z};
        flt[1] = { X, Y,  Z};
        flt[2] = { X,-Y, -Z};
        flt[3] = { X,-Y,  Z};
        flt[4] = {-X,-Y, -Z};
        flt[5] = {-X,-Y,  Z};
        flt[6] = {-X, Y, -Z};
        flt[7] = {-X, Y,  Z};
        // top plate
        flt[8] = { X, Y, Z};
        flt[9] = {-X, Y, Z};
        flt[10] = { X,-Y, Z};
        flt[11] = {-X,-Y, Z};
        flt[12] = { X, Y, Z};
        flt[13] = { X,-Y, Z};
        flt[14] = {-X, Y, Z};
        flt[15] = {-X,-Y, Z};
        // bottom plate
        flt[16] = { X, Y, -Z};
        flt[17] = {-X, Y, -Z};
        flt[18] = { X,-Y, -Z};
        flt[19] = {-X,-Y, -Z};
        flt[20] = { X, Y, -Z};
        flt[21] = { X,-Y, -Z};
        flt[22] = {-X, Y, -Z};
        flt[23] = {-X,-Y, -Z};
        return 24;
    }
    
    /// Octahedron is made of 8 triangles and 24 vertices
    size_t setOctahedron(flute6* flt, float R=1.46459188756f)
    {
        // the default size is set to match the volume of the unit sphere
        // triangles ordered counterclockwise
        const float N = 1 / M_SQRT3;
        size_t i = 0;
        // base at Z = 0, faces orientated up:
        flt[i++] = {-R, 0, 0,  0, 0,-1};
        flt[i++] = { 0,-R, 0,  0, 0,-1};
        flt[i++] = { 0, R, 0,  0, 0,-1};
        flt[i++] = { 0, R, 0,  0, 0,-1};
        flt[i++] = { 0,-R, 0,  0, 0,-1};
        flt[i++] = { R, 0, 0,  0, 0,-1};
        // lower size Z < 0
        flt[i++] = { 0, 0,-R, -N, N,-N};
        flt[i++] = {-R, 0, 0, -N, N,-N};
        flt[i++] = { 0, R, 0, -N, N,-N};
        flt[i++] = { 0, 0,-R, -N,-N,-N};
        flt[i++] = { 0,-R, 0, -N,-N,-N};
        flt[i++] = {-R, 0, 0, -N,-N,-N};
        flt[i++] = { 0, 0,-R,  N, N,-N};
        flt[i++] = { 0, R, 0,  N, N,-N};
        flt[i++] = { R, 0, 0,  N, N,-N};
        flt[i++] = { 0, 0,-R,  N,-N,-N};
        flt[i++] = { R, 0, 0,  N,-N,-N};
        flt[i++] = { 0,-R, 0,  N,-N,-N};
        // upper size Z > 0
        flt[i++] = { 0, 0, R, -N,-N, N};
        flt[i++] = {-R, 0, 0, -N,-N, N};
        flt[i++] = { 0,-R, 0, -N,-N, N};
        flt[i++] = {-R, 0, 0, -N, N, N};
        flt[i++] = { 0, 0, R, -N, N, N};
        flt[i++] = { 0, R, 0, -N, N, N};
        flt[i++] = { R, 0, 0,  N,-N, N};
        flt[i++] = { 0, 0, R,  N,-N, N};
        flt[i++] = { 0,-R, 0,  N,-N, N};
        flt[i++] = { 0, 0, R,  N, N, N};
        flt[i++] = { R, 0, 0,  N, N, N};
        flt[i++] = { 0, R, 0,  N, N, N};
        // base at Z = 0, faces orientated down:
        flt[i++] = {-R, 0, 0,  0, 0,-1};
        flt[i++] = { 0, R, 0,  0, 0,-1};
        flt[i++] = { 0,-R, 0,  0, 0,-1};
        flt[i++] = { 0, R, 0,  0, 0,-1};
        flt[i++] = { R, 0, 0,  0, 0,-1};
        flt[i++] = { 0,-R, 0,  0, 0,-1};
        return i; // 12+12+6+6 = 36 vertices
    }

    
#if ( 0 )
    // this is used to calculate the vertices of the icosahedron
    void icoFace(float* a, float* b, float* c)
    {
        float nx = (a[0]+b[0]+c[0]) / 3.0f;
        float ny = (a[1]+b[1]+c[1]) / 3.0f;
        float nz = (a[2]+b[2]+c[2]) / 3.0f;
        float n = std::sqrt( nx * nx + ny * ny + nz * nz );
        nx /= n;
        ny /= n;
        nz /= n;
        flute6* pts = gym::mapBufferV3N3(3);
        pts[0] = {a[0], a[1], a[2], nx, ny, nz};
        pts[1] = {b[0], b[1], b[2], nx, ny, nz};
        pts[2] = {c[0], c[1], c[2], nx, ny, nz};
        gym::unmapBufferV3N3();
        gym::ref_view();
        gym::drawTriangleStrip(0, 3);
        if ( 1 ) {
            printf("%2.0f, %2.0f, %2.0f, %+9.7f, %+9.7f, %+9.7f\n", a[0], a[1], a[2], nx, ny, nz);
            printf("%2.0f, %2.0f, %2.0f, %+9.7f, %+9.7f, %+9.7f\n", b[0], b[1], b[2], nx, ny, nz);
            printf("%2.0f, %2.0f, %2.0f, %+9.7f, %+9.7f, %+9.7f\n", c[0], c[1], c[2], nx, ny, nz);
        }
    }
    
    void icoFace(float* pts, size_t a, size_t b, size_t c)
    {
        icoFace(pts+3*a, pts+3*b, pts+3*c);
    }
    
    void ICOSAHEDRON()
    {
        const float G = 0.5+0.5*std::sqrt(5.0);
        const float H = 1/std::sqrt(G*G+1); //0.5257311121f;
        const float T = G * H;   //0.8506508084f;
        
        // Twelve vertices of icosahedron on unit sphere
        float pts[3*12] = {
            +T,  H,  0, // 0
            -T, -H,  0, // 1
            +T, -H,  0, // 3
            -T,  H,  0, // 2
            +H,  0,  T, // 4
            -H,  0, -T, // 5
            +H,  0, -T, // 6
            -H,  0,  T, // 7
            +0,  T,  H, // 8
            +0, -T, -H, // 9
            +0, -T,  H, // 10
            +0,  T, -H  // 11
        };
        
        glEnableClientState(GL_NORMAL_ARRAY);
        /* The faces are ordered with increasing Z */
        icoFace(pts, 5,  6, 9);
        icoFace(pts, 5, 11, 6);
        
        icoFace(pts, 6, 2,  9);
        icoFace(pts, 3, 11, 5);
        icoFace(pts, 1, 5,  9);
        icoFace(pts, 0, 6, 11);
        
        icoFace(pts, 0, 2,  6);
        icoFace(pts, 1, 3,  5);
        
        icoFace(pts, 1, 9, 10);
        icoFace(pts, 0, 11, 8);
        icoFace(pts, 8, 11, 3);
        icoFace(pts, 9, 2, 10);
        
        icoFace(pts, 0, 4,  2);
        icoFace(pts, 1, 7,  3);
        
        icoFace(pts, 0, 8,  4);
        icoFace(pts, 1, 10, 7);
        icoFace(pts, 2, 4, 10);
        icoFace(pts, 7, 8,  3);
        
        icoFace(pts, 4, 8,  7);
        icoFace(pts, 4, 7, 10);
        glDisableClientState(GL_NORMAL_ARRAY);
    }
#endif
    
    /// Icosahedron with 20 triangles = 60 vertices, using true face-normals
    size_t setIcosahedron(flute6* flt, float R=1.0f)
    {
        const float T = R * 0.8506508084f;      // (1 + sqrt(5))/2
        const float H = R * 0.5257311121f;      // 1 / sqrt(1+T^2)
        const float N = 1 / M_SQRT3; // 0.5773503
        const float X = 0.3568221;
        const float Y = 0.9341724;

        size_t i = 0;
        flt[i++] = {-H,  0, -T, +0,-X,-Y};
        flt[i++] = { H,  0, -T, +0,-X,-Y};
        flt[i++] = { 0, -T, -H, +0,-X,-Y};
        
        flt[i++] = {-H,  0, -T, +0,+X,-Y};
        flt[i++] = { 0,  T, -H, +0,+X,-Y};
        flt[i++] = { H,  0, -T, +0,+X,-Y};
        
        flt[i++] = { H,  0, -T, +N,-N,-N};
        flt[i++] = { T, -H,  0, +N,-N,-N};
        flt[i++] = { 0, -T, -H, +N,-N,-N};
        
        flt[i++] = {-T,  H,  0, -N,+N,-N};
        flt[i++] = { 0,  T, -H, -N,+N,-N};
        flt[i++] = {-H,  0, -T, -N,+N,-N};
        
        flt[i++] = {-T, -H,  0, -N,-N,-N};
        flt[i++] = {-H,  0, -T, -N,-N,-N};
        flt[i++] = { 0, -T, -H, -N,-N,-N};
        
        flt[i++] = { T,  H,  0, +N,+N,-N};
        flt[i++] = { H,  0, -T, +N,+N,-N};
        flt[i++] = { 0,  T, -H, +N,+N,-N};
        
        flt[i++] = { T,  H,  0, +Y,+0,-X};
        flt[i++] = { T, -H,  0, +Y,+0,-X};
        flt[i++] = { H,  0, -T, +Y,+0,-X};
        
        flt[i++] = {-T, -H,  0, -Y,+0,-X};
        flt[i++] = {-T,  H,  0, -Y,+0,-X};
        flt[i++] = {-H,  0, -T, -Y,+0,-X};
        
        flt[i++] = {-T, -H,  0, -X,-Y,+0};
        flt[i++] = { 0, -T, -H, -X,-Y,+0};
        flt[i++] = { 0, -T,  H, -X,-Y,+0};
        
        flt[i++] = { T,  H,  0, +X,+Y,+0};
        flt[i++] = { 0,  T, -H, +X,+Y,+0};
        flt[i++] = { 0,  T,  H, +X,+Y,+0};
        
        flt[i++] = { 0,  T,  H, -X,+Y,+0};
        flt[i++] = { 0,  T, -H, -X,+Y,+0};
        flt[i++] = {-T,  H,  0, -X,+Y,+0};
        
        flt[i++] = { 0, -T, -H, +X,-Y,+0};
        flt[i++] = { T, -H,  0, +X,-Y,+0};
        flt[i++] = { 0, -T,  H, +X,-Y,+0};
        
        flt[i++] = { T,  H,  0, +Y,+0,+X};
        flt[i++] = { H,  0,  T, +Y,+0,+X};
        flt[i++] = { T, -H,  0, +Y,+0,+X};
        
        flt[i++] = {-T, -H,  0, -Y,+0,+X};
        flt[i++] = {-H,  0,  T, -Y,+0,+X};
        flt[i++] = {-T,  H,  0, -Y,+0,+X};
        
        flt[i++] = { T,  H,  0, +N,+N,+N};
        flt[i++] = { 0,  T,  H, +N,+N,+N};
        flt[i++] = { H,  0,  T, +N,+N,+N};
        
        flt[i++] = {-T, -H,  0, -N,-N,+N};
        flt[i++] = { 0, -T,  H, -N,-N,+N};
        flt[i++] = {-H,  0,  T, -N,-N,+N};
        
        flt[i++] = { T, -H,  0, +N,-N,+N};
        flt[i++] = { H,  0,  T, +N,-N,+N};
        flt[i++] = { 0, -T,  H, +N,-N,+N};
        
        flt[i++] = {-H,  0,  T, -N,+N,+N};
        flt[i++] = { 0,  T,  H, -N,+N,+N};
        flt[i++] = {-T,  H,  0, -N,+N,+N};
        
        flt[i++] = { H,  0,  T, +0,+X,+Y};
        flt[i++] = { 0,  T,  H, +0,+X,+Y};
        flt[i++] = {-H,  0,  T, +0,+X,+Y};
        
        flt[i++] = { H,  0,  T, +0,-X,+Y};
        flt[i++] = {-H,  0,  T, +0,-X,+Y};
        flt[i++] = { 0, -T,  H, +0,-X,+Y};
        assert_true( i == 60 );
        return i;
    }
    
    /// this only sets vertices, skipping normals
    size_t setIcoidBuffer(flute3* flu, Tesselator::INDEX* idx)
    {
        const float U = 1.0;
        const float Z = std::sqrt(0.2);
        const float C = std::cos(M_PI*0.4);
        const float S = std::sin(M_PI*0.4);
        const float D = C*C - S*S;
        const float T = C*S + S*C;
        const float H = std::sqrt(1 + Z*Z);

        // Twelve vertices of icosahedron on unit sphere
        flu[0] = { 0,  0, -H};
        flu[1] = { U,  0, -Z};
        flu[2] = { C, -S, -Z};
        flu[3] = { D, -T, -Z};
        flu[4] = { D,  T, -Z};
        flu[5] = { C,  S, -Z};
        flu[6] = {-D, -T,  Z};
        flu[7] = {-C, -S,  Z};
        flu[8] = {-U,  0,  Z};
        flu[9] = {-C,  S,  Z};
        flu[10] = {-D,  T,  Z};
        flu[11] = { 0,  0,  H};

        const size_t F = 30; // number of points in the triangle strip
        Tesselator::INDEX strip[F] = {0,1,2,6,7,11,
            11,8,7,3,2,0, 0,3,4,8,9,11, 11,10,9,5,4,0, 0,5,1,10,6,11};
        memcpy(idx, strip, sizeof(strip));
        return F;
    }
    
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
        flt[i++] = {cR,-sR, B,  0, -1, 0};
        flt[i++] = { 1,  0, H,  0, -1, 0};
        flt[i++] = { 1,  0, U,  0, -1, 0};
        
        flt[i++] = {cR,-sR, B,  0, -1, 0};
        flt[i++] = { 1,  0, U,  0, -1, 0};
        flt[i++] = { 0,  0, T,  0, -1, 0};
        
        flt[i++] = {cR, sR, B,  0, +1, 0};
        flt[i++] = { 0,  0, T,  0, +1, 0};
        flt[i++] = { 1,  0, U,  0, +1, 0};
        
        flt[i++] = {cR, sR, B,  0, +1, 0};
        flt[i++] = { 1,  0, U,  0, +1, 0};
        flt[i++] = { 1,  0, H,  0, +1, 0};
        
        flt[i++] = {cR,-sR, B,  S, -C, 0};
        flt[i++] = { 0,  0, T,  S, -C, 0};
        flt[i++] = {-C, -S, U,  S, -C, 0};
        
        flt[i++] = {cR,-sR, B,  S, -C, 0};
        flt[i++] = {-C, -S, U,  S, -C, 0};
        flt[i++] = {-C, -S, H,  S, -C, 0};
        
        flt[i++] = {-R,  0, B, -S,  C, 0};
        flt[i++] = {-C, -S, H, -S,  C, 0};
        flt[i++] = {-C, -S, U, -S,  C, 0};
        
        flt[i++] = {-R,  0, B, -S,  C, 0};
        flt[i++] = {-C, -S, U, -S,  C, 0};
        flt[i++] = { 0,  0, T, -S,  C, 0};
        
        flt[i++] = {cR, sR, B,  S,  C, 0};
        flt[i++] = {-C,  S, H,  S,  C, 0};
        flt[i++] = {-C,  S, U,  S,  C, 0};
        
        flt[i++] = {cR, sR, B,  S,  C, 0};
        flt[i++] = {-C,  S, U,  S,  C, 0};
        flt[i++] = { 0,  0, T,  S,  C, 0};
        
        flt[i++] = {-R,  0, B, -S, -C, 0};
        flt[i++] = { 0,  0, T, -S, -C, 0};
        flt[i++] = {-C,  S, U, -S, -C, 0};
        
        flt[i++] = {-R,  0, B, -S, -C, 0};
        flt[i++] = {-C,  S, U, -S, -C, 0};
        flt[i++] = {-C,  S, H, -S, -C, 0};
        
        flt[i++] = {cR, sR, B,  C, -S,-1};
        flt[i++] = {-R,  0, B,  C, -S,-1};
        flt[i++] = {-C,  S, H,  C, -S,-1};
        
        flt[i++] = {-R,  0, B,  C,  S,-1};
        flt[i++] = {cR,-sR, B,  C,  S,-1};
        flt[i++] = {-C, -S, H,  C,  S,-1};
        
        flt[i++] = {cR,-sR, B, -1,  0,-1};
        flt[i++] = {cR, sR, B, -1,  0,-1};
        flt[i++] = { 1,  0, H, -1,  0,-1};
        assert_true( i == 45 );
        return i;
    }
    
    /// set 12 pentagon as on a football constructed as a truncated icosahedron
    size_t setFootballPentagons(flute6* flu, float R, float P)
    {
        const float K = std::sqrt(0.2);
        float C = std::cos(M_PI * 0.4); // 0.31
        float S = std::sin(M_PI * 0.4); // 0.95
        float D = C*C - S*S; // double angle -0.80
        float T = C*S + C*S; // double angle +0.59
        float H = std::sqrt(1+K*K);

        /* Twelve vertices of icosahedron on unit sphere,
        ordered to follow a helical path */
        flute3 vex[13] = {
            { 0,  0,  H },
            {-1,  0,  K },
            {-C,  S,  K },
            {-D,  T,  K },
            {-D, -T,  K },
            {-C, -S,  K },
            { D, -T, -K },
            { D,  T, -K },
            { C,  S, -K },
            { 1,  0, -K },
            { C, -S, -K },
            { 0,  0, -H },
            { 1,  0, -K },
        };

        C *= P;
        S *= P;
        D *= P;
        T *= P;
        size_t i = 0;
        for ( int u = 0; u < 12; ++u )
        {
            flute3 Z = normalize(vex[u]);
            flute3 X = normalize(vex[u+1] - dot(vex[u+1], Z) * Z);
            flute3 Y = cross(Z, X);

            //float pentagon[] = { R,0, C,S, D,T, D,-T, C,-S };
            flute3 M = normalize(Z+P*X);
            flute3 N = normalize(Z+C*X-S*Y);
            flute3 O = normalize(Z+D*X-T*Y);
            flute3 P = normalize(Z+D*X+T*Y);
            flute3 Q = normalize(Z+C*X+S*Y);
            
            flu[i+0] = { R * M, M };
            flu[i+1] = { R * M, M };
            flu[i+2] = { R * N, N };
            flu[i+3] = { R * Z, Z };
            flu[i+4] = { R * O, O };
            flu[i+5] = { R * P, P };
            flu[i+6] = { R * P, P };
            flu[i+7] = { R * Z, Z };
            flu[i+8] = { R * Q, Q };
            flu[i+9] = { R * M, M };
            flu[i+10] = { R * M, M };
            flu[i+11] = { R * M, M };
            i += 12;
        }
        return i;
    }
    
    //-----------------------------------------------------------------------
    
    static size_t sizeCubeBuffers()
    {
        return ( 12*2 + 36 + 60 + 45 + 36*2 + 22*3 );
    }
    
    size_t setCubeBuffers(flute6* ptr, flute6* const ori)
    {
        unsigned i = 0, s = ptr - ori;
        cubes_[0] = i+s; i += setTetrahedron(ptr);
        cubes_[1] = i+s; i += invTetrahedron(ptr+i);
        cubes_[2] = i+s; i += setOctahedron(ptr+i);
        cubes_[3] = i+s; i += setIcosahedron(ptr+i);
        cubes_[4] = i+s; i += setArrowTail(ptr+i);
        cubes_[5] = i+s; i += setCubeFaces(ptr+i, 1.0f, 1.0f, 1.0f);
        cubes_[6] = i+s; i += setCubeFaces(ptr+i, 0.57735f, 0.57735f, 0.57735f);
        cubes_[7] = i+s; i += setHexTube(ptr+i, 0, 1, 1.0f);  // hexTube
        cubes_[8] = i+s; i += setHexTube(ptr+i, 0, 1, 0.5f);  // thinTube
        cubes_[9] = i+s; i += setHexTube(ptr+i, 0, 256.f, 0.5f); // thinLongTube
        assert_true( i <= sizeCubeBuffers() );
        return i;
    }
    
    void cubeEdges(float width)
    {
        gym::bindBufferV3(buf_[0], 1, 3*blobs_[3]);
        gym::drawLines(width, 0, 24);
        gym::cleanupV();
    }
    
    void cubeVerticalEdges(float width)
    {
        gym::bindBufferV3(buf_[0], 1, 3*blobs_[3]);
        gym::drawLines(width, 0, 8);
        gym::cleanupV();
    }

    void doCubeTriangles(size_t inx, GLsizei cnt)
    {
        gym::bindBufferV3N3(buf_[0]);
        gym::drawTriangles(cubes_[inx], cnt);
        gym::cleanupVN();
    }

    void doCubeStrip(size_t inx, GLsizei cnt)
    {
        gym::bindBufferV3N3(buf_[0]);
        gym::drawTriangleStrip(cubes_[inx], cnt);
        gym::cleanupVN();
    }

    void tetrahedron() { doCubeTriangles(0, 12); }
    void upsideTetra() { doCubeTriangles(1, 12); }
    void star()        { doCubeTriangles(0, 24); } // union of 2 tetrahedrons

    void octahedron()  { doCubeTriangles(2+6, 24); }
    void invPyramid()  { doCubeTriangles(2, 18); }
    void pyramid()     { doCubeTriangles(2+18, 18); }
    void icosahedron() { doCubeTriangles(3, 60); }
    
    void arrowTail() { doCubeTriangles(4, 45); }
    void cube()      { doCubeTriangles(5, 36); }
    void cubeFaces() { doCubeTriangles(5, 12); }
    void smallCube() { doCubeTriangles(6, 36); }

    void hexTube()      { doCubeStrip(7, 22); }
    void thinTube()     { doCubeStrip(8, 22); }
    void thinLongTube() { doCubeStrip(9, 22); }
    
    void ICOSAHEDRON()
    {
        flute6* tmp = gym::mapBufferV3N3(60);
        setIcosahedron(tmp);
        gym::unmapBufferV3N3();
        gym::drawTriangles(0, 60);
        gym::cleanupVN();
    }

    //-----------------------------------------------------------------------
    #pragma mark - Blobs = 3D objects without normals

    size_t makeBlob(flute3* flu)
    {
        constexpr float R = 1.f, U = -1.f, H(M_SQRT2);
        /* start from a centered cube, rotated appropriately
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
            flu[i++] = (a+a);
            flu[i++] = (a+b);
            flu[i++] = (a+c);
            flu[i++] = (b+c);
        }
        a = {pts[3*n+3], pts[3*n+4], pts[3*n+5]};
        flu[i++] = (c+c);
        flu[i++] = (a+c);
        flu[i++] = (a+c);
        a = {pts[3*13], pts[3*13+1], pts[3*13+2]};
        flu[i++] = (a+a);
        for ( n = 13; n > 2; n -= 2 )
        {
            a = {pts[3*n  ], pts[3*n+1], pts[3*n+2]};
            b = {pts[3*n-3], pts[3*n-2], pts[3*n-1]};
            c = {pts[3*n-6], pts[3*n-5], pts[3*n-4]};
            flu[i++] = (a+a);
            flu[i++] = (a+b);
            flu[i++] = (a+c);
            flu[i++] = (b+c);
        }
        a = {pts[0], pts[1], pts[2]};
        flu[i++] = (c+c);
        flu[i++] = (a+c);
        assert_true(i==54);
        return 54;
    }
    
    size_t setBlob(flute3* flu)
    {
        const size_t FF = 54;
        float U = 1.0;
        float H = 0.5;
        float S = M_SQRT1_2;
        float O = std::sqrt(1/3.0); //0.57735
        float T = std::sqrt(2/3.0); //0.816497
        
        flute3 pts[26] = {
            {T, O, 0}, {H, S, H}, {U, 0, 0}, {S, 0, S}, {T, -O, 0}, {H, -S, H},
            {0, -U, 0}, {-H, -S, H}, {-T, -O, 0}, {-S, 0, S}, {-U, 0, 0}, {-H, S, H},
            {-T, O, 0}, {0, U, 0}, {-H, S, -H}, {H, S, -H}, {0, O, -T}, {S, 0, -S},
            {0, 0, -U}, {H, -S, -H}, {0, -O, -T}, {-H, -S, -H}, {-S, 0, -S}, {0, O, T},
            {0, 0, U}, {0, -O, T}
        };

        int inx[FF] = {
            0, 1, 2, 3, 4, 5,
            6, 7, 8, 9, 10, 11,
            12, 13, 14, 15, 16, 17,
            18, 19, 20, 21, 18, 22,
            16, 14, 14, 12, 12, 14,
            10, 22, 8, 21, 6, 19,
            4, 17, 2, 15, 0, 13,
            1, 11, 23, 9, 24, 7,
            25, 5, 24, 3, 23, 1
        };

        for ( size_t i = 0; i < FF; ++i )
            flu[i] = normalize(pts[inx[i]]);

        return FF;
    }

    /* This moves some vertices to add an hexagonal needle to the blob */
    void modifyBlob(flute3 * flu, float F)
    {
        const float R = 0.6f;
        const float Y = R * 0.8660254037844386f; // sqrtf(3)/2;
        const float X = R * 0.5f;
        const float Z = 0.6f;
        //{ 1, 3, 5, 7, 9, 11, 40, 41, 43, 45, 47, 49, 51 }
        // pull some vertices far away in Z
        for ( int u : { 44, 46, 48, 50, 52 } ) flu[u] = { 0, 0, F };
        // set the 6 vertex of an hexagon:
        for ( int u : { 3, 51 } ) flu[u] = { R, 0, Z};
        for ( int u : { 5, 49 } ) flu[u] = { X,-Y, Z};
        for ( int u : { 7, 47 } ) flu[u] = {-X,-Y, Z};
        for ( int u : { 9, 45 } ) flu[u] = {-R, 0, Z};
        for ( int u : {11, 43 } ) flu[u] = {-X, Y, Z};
        for ( int u : {1,42,53} ) flu[u] = { X, Y, Z};
    }
    
    size_t setPin(flute3* flu)
    {
        size_t i = setBlob(flu);
        modifyBlob(flu, 10);
        return i;
    }

    void thing()
    {
        flute3* flu = gym::mapBufferV3(64);
        setBlob(flu);
        modifyBlob(flu, 10);
        gym::unmapBufferV3N0();

        gym::color_both(1,1,1,1);
        for ( unsigned u : { 1, 3, 5, 7, 9, 11 } )
            gym::drawPoints(8, u, 1);
        gym::color_both(1,0,1,1);
        for ( unsigned u : { 44, 46, 48, 50, 52 } )
            gym::drawPoints(8, u, 1);
        gym::drawLineStrip(1, 0, 54);
        
        gym::enableCullFace(GL_BACK);
        gym::color_both(0,1,1,0.75);
        gym::drawTriangleStrip(0, 26);
        gym::color_both(1,1,0,0.75);
        gym::drawTriangleStrip(28, 26);
        gym::cleanupVN();
    }
    
    void nothing()
    {
    }

    static size_t sizeBlobBuffers()
    {
        return ( 54 + 54 + 16 + 24 );
    }

    size_t setBlobBuffers(flute3* ptr, flute3* const ori)
    {
        unsigned i = 0, s = ptr - ori;
        blobs_[0] = i+s; i += setBlob(ptr+i);
        blobs_[1] = i+s; i += setPin(ptr+i);
        blobs_[2] = i+s; i += setCuboid(ptr+i, 1.0);
        blobs_[3] = i+s; i += setCubeEdges(ptr+i, 1.0f, 1.0f, 1.0f);
        assert_true( i <= sizeBlobBuffers() );
        return i;
    }

    void doBlobStrip(size_t inx, GLsizei cnt)
    {
        gym::bindBufferV3N0(buf_[0], 0);
        gym::drawTriangleStrip(blobs_[inx], cnt);
        gym::cleanupVN();
    }

    void blob()   { doBlobStrip(0, 54); }
    void needle() { doBlobStrip(1, 54); }
    void cuboid() { doBlobStrip(2, 14); }
    //void icoidS() { doBlobStrip(3, 32); }

    //-----------------------------------------------------------------------
    #pragma mark - 2D Circle
    
    size_t setSquare(flute2* flu, float X, float Y)
    {
        flu[0] = { X, Y};
        flu[1] = {-X, Y};
        flu[2] = {-X,-Y};
        flu[3] = { X,-Y};
        flu[4] = { X, Y};
        return 5;
    }

    size_t setCircle(flute2* flu, unsigned inc, float R)
    {
        size_t i = 0;
        for ( unsigned n = 0; n <= pi_4half; n += inc )
            flu[i++] = {R*cos_(n), R*sin_(n)};
        return i;
    }
    
    size_t setCircleZ(flute3* flu, unsigned inc, float R, float Z, float N)
    {
        size_t i = 0;
        for ( unsigned n = 0; n <= pi_4half; n += inc )
        {
            float C = cos_(n), S = sin_(n);
            flu[i++] = {R*C, R*S, Z};
        }
        return i;
    }
    
    /// set buffer for strokeCapsule
    size_t setCapsuleStroke(flute2* flu, float R)
    {
        float Y = R * 0.5 * M_SQRT3;
        float H = R * 0.5;
        flu[0] = { 0, 1+R };
        flu[1] = { H, 1+Y };
        flu[2] = { Y, 1+H };
        flu[3] = { R, 1 };
        flu[4] = { R, -1 };
        flu[5] = { Y, -1-H };
        flu[6] = { H, -1-Y };
        flu[7] = { 0, -1-R };
        flu[8] = { -H, -1-Y };
        flu[9] = { -Y, -1-H };
        flu[10] = { -R, -1 };
        flu[11] = { -R, 1 };
        flu[12] = { -Y, 1+H };
        flu[13] = { -H, 1+Y };
        flu[14] = { 0, 1+R };
        return 15;
    }
    
    /// set buffer for paintCapsule
    size_t setCapsulePaint(flute2* flu, float R)
    {
        float Y = R * 0.5 * M_SQRT3;
        float H = R * 0.5;
        flu[0] = { 0, 1+R };
        flu[1] = { -H, 1+Y };
        flu[2] = { H, 1+Y };
        flu[3] = { -Y, 1+H };
        flu[4] = { Y, 1+H };
        flu[5] = { -R, 1 };
        flu[6] = { R, 1 };
        flu[7] = { -R, -1 };
        flu[8] = { R, -1 };
        flu[9] = { -Y, -1-H };
        flu[10] = { Y, -1-H };
        flu[11] = { -H, -1-Y };
        flu[12] = { H, -1-Y };
        flu[13] = { 0, -1-R };
        flu[14] = { 0, 1+R };
        return 15;
    }
    
    /// set buffer for strokeCross
    size_t setCrossStroke(flute2* flu)
    {
        float W = 0.2;
        float L = 0.75;
        float R = 0.75;
        flu[0] = { -W, -L };
        flu[1] = { W, -L };
        flu[2] = { W, R-W };
        flu[3] = { L, R-W };
        flu[4] = { L, R+W };
        flu[5] = { W, R+W };
        flu[6] = { W, R+L };
        flu[7] = { -W, R+L };
        flu[8] = { -W, R+W };
        flu[9] = { -L, R+W };
        flu[10] = { -L, R-W };
        flu[11] = { -W, R-W };
        flu[12] = { -W, -L };
        return 13;
    }

    /// set buffer for paintCross
    size_t setCrossPaint(flute2* flu)
    {
        float W = 0.2;
        float L = 0.75;
        float R = 0.75;
        flu[0] = { -L, R+W };
        flu[1] = { -L, R-W };
        flu[2] = { -W, R+W };
        flu[3] = { -W, R-W };
        flu[4] = { -W, R+L };
        flu[5] = { -W, -L };
        flu[6] = { W, R+L };
        flu[7] = { W, -L };
        flu[8] = { W, R+W };
        flu[9] = { W, R-W };
        flu[10] = { L, R+W };
        flu[11] = { L, R-W };
        return 12;
    }

    static size_t sizeCircBuffers()
    {
        return 56 + 4 * pi_4half;
    }
    
    size_t setCircBuffers(flute2* ptr, flute2* const ori)
    {
        unsigned i = 0, s = ptr - ori;
        discs_[0] = i+s; i += setCircle(ptr+i, 1, 1);
        discs_[1] = i+s; i += setCircle(ptr+i, 2, 1);
        discs_[2] = i+s; i += setSquare(ptr+i, 1, 1);
        discs_[4] = i+s; i += setCapsuleStroke(ptr+i, 0.5);
        discs_[5] = i+s; i += setCapsulePaint(ptr+i, 0.5);
        discs_[6] = i+s; i += setCrossStroke(ptr+i);
        discs_[7] = i+s; i += setCrossPaint(ptr+i);
        assert_true( i <= sizeCircBuffers() );
        return i;
    }
    
    void paintHalo(float R0, float R1)
    {
        flute2 *buf = gym::mapBufferV2(2*pi_4half+2);
        flute2 *ptr = buf;
        for ( unsigned j = 0; j <= pi_4half; ++j )
        {
            float C = cos_(j), S = sin_(j);
            ptr[0] = { R0*C, R0*S };
            ptr[1] = { R1*C, R1*S };
            ptr += 2;
        }
        gym::unmapBufferV2();
        assert_true( ptr <= buf+2*pi_4half+2 );
        gym::drawTriangleStrip(0, ptr-buf);
        //gym::drawPoints(width, 0, ptr-buf);
    }

    void paintCapsule(float L, float R, float rad, unsigned inc)
    {
        flute2 *buf = gym::mapBufferV2(pi_4half+4);
        flute2 *ptr = buf;
        *ptr++ = { R+rad, 0 };
        for ( unsigned j = inc; j <= pi_1half; j += inc )
        {
            float C = cos_(j), S = sin_(j);
            ptr[0] = { rad*C + R,  rad*S };
            ptr[1] = { rad*C + R, -rad*S };
            ptr += 2;
        }
        for ( unsigned j = pi_1half; j < pi_2half; j += inc )
        {
            float C = cos_(j), S = sin_(j);
            ptr[0] = { rad*C + L,  rad*S };
            ptr[1] = { rad*C + L, -rad*S };
            ptr += 2;
        }
        *ptr++ = { L-rad, 0 };
        gym::unmapBufferV2();
        assert_true( ptr <= buf+pi_4half+4 );
        gym::drawTriangleStrip(0, ptr-buf);
        //gym::drawPoints(width, 0, ptr-buf);
    }

    void strokeCapsule(float L, float R, float rad, float width, unsigned inc)
    {
        flute2 *buf = gym::mapBufferV2(pi_4half+4);
        flute2 *ptr = buf;
        for ( unsigned j = pi_1half; j <= pi_3half; j += inc )
            *ptr++ = { rad*cos_(j) + L, rad*sin_(j) };
        for ( unsigned j = pi_3half; j <= pi_5half; j += inc )
            *ptr++ = { rad*cos_(j) + R, rad*sin_(j) };
        *ptr++ = { L, rad };
        gym::unmapBufferV2();
        assert_true( ptr <= buf+pi_4half+4 );
        gym::drawLineStrip(width, 0, ptr-buf);
        //gym::drawPoints(width, 0, ptr-buf);
    }
        
    void paintBicapsule(float L, float R, float rad, float G, float H, unsigned inc)
    {
        const float F = std::max(rad - G, 0.f);
        const float X = F + H;
        flute2 *buf = gym::mapBufferV2(2*pi_4half+4);
        flute2 *ptr = buf;
        *ptr++ = { R+rad, 0 };
        for ( unsigned j = inc; j <= pi_1half; j += inc )
        {
            float C = cos_(j), S = sin_(j);
            ptr[0] = { rad*C + R,  rad*S };
            ptr[1] = { rad*C + R, -rad*S };
            ptr += 2;
        }
        for ( unsigned j = pi_1half; j <= pi_2half; j += inc )
        {
            float C = cos_(j), S = sin_(j);
            ptr[0] = { F*C + X, +F*S + G };
            ptr[1] = { F*C + X, -F*S - G };
            ptr += 2;
        }
        for ( unsigned j = 0; j < pi_1half; j += inc )
        {
            float C = cos_(j), S = sin_(j);
            ptr[0] = { F*C - X, +F*S + G };
            ptr[1] = { F*C - X, -F*S - G };
            ptr += 2;
        }
        for ( unsigned j = pi_1half; j < pi_2half; j += inc )
        {
            float C = cos_(j), S = sin_(j);
            ptr[0] = { rad*C + L,  rad*S };
            ptr[1] = { rad*C + L, -rad*S };
            ptr += 2;
        }
        *ptr++ = { L-rad, 0 };
        gym::unmapBufferV2();
        assert_true( ptr <= buf+2*pi_4half+4 );
        gym::drawTriangleStrip(0, ptr-buf);
        //gym::drawPoints(width, 0, ptr-buf);
    }

    void strokeBicapsule(float L, float R, float rad, float G, float H, float width, unsigned inc)
    {
        const float F = std::max(rad - G, 0.f);
        const float X = F + H;
        flute2 *buf = gym::mapBufferV2(2*pi_4half+8);
        flute2 *ptr = buf;
        for ( unsigned j = pi_1half; j <= pi_3half; j += inc )
            *ptr++ = { rad*cos_(j) + L, rad*sin_(j) };
        for ( unsigned j = pi_3half; j <= pi_4half; j += inc )
            *ptr++ = { F*cos_(j) - X, F*sin_(j) - G };
        for ( unsigned j = pi_2half; j <= pi_3half; j += inc )
            *ptr++ = { F*cos_(j) + X, F*sin_(j) - G };
        for ( unsigned j = pi_3half; j <= pi_5half; j += inc )
            *ptr++ = { rad*cos_(j) + R, rad*sin_(j) };
        for ( unsigned j = pi_1half; j <= pi_2half; j += inc )
            *ptr++ = { F*cos_(j) + X, F*sin_(j) + G };
        for ( unsigned j = 0; j <= pi_1half; j += inc )
            *ptr++ = { F*cos_(j) - X, F*sin_(j) + G };
        *ptr++ = { L, rad };
        gym::unmapBufferV2();
        assert_true( ptr <= buf+2*pi_4half+8 );
        gym::drawLineStrip(width, 0, ptr-buf);
        //gym::drawPoints(width, 0, ptr-buf);
    }

    //-----------------------------------------------------------------------
    #pragma mark - Tubes primitives

    inline unsigned nbTrianglesCylinder(unsigned inc, int top = 0)
    {
        return 2 + 2 * (( pi_6half + top*pi_2half ) / inc );
    }

    /// set triangle strip for a tube of radius 1 at Z=B, and R at Z=T, closed at Z=B
    /** Note that if T < B, the triangles will be given clockwise */
    size_t setCylinder(flute6* flu, unsigned inc, float B, float T, float R, int top = 0)
    {
        const float H(T-B), W(R-1);
        const float Z = std::copysign(1.f, H);
        const float tg(1.f/sqrtf(H*H+W*W));
        const float X(tg*H);
        const float Y(tg*W);
        unsigned i = 0;
        unsigned p = pi_2half;
        // bottom disc, from PI to 0:
        flu[i++] = { 1, 0, B, 0, 0, -Z };
        while ( p > inc )
        {
            p -= inc;
            float C = -cos_(p), S = sin_(p);
            flu[i++] = { C,-S, B, 0, 0, -Z };
            flu[i++] = { C, S, B, 0, 0, -Z };
        }
        assert_true( p == inc );
        flu[i++] = { -1, 0, B, 0, 0, -Z };
        // repeat point to adjust normal:
        flu[i++] = { -1, 0, B, X, 0, Y };
        flu[i++] = { -R, 0, T, X, 0, Y };
        // sides, from 0 to 2*PI:
        while ( p < pi_4half )
        {
            float C = -cos_(p), S = sin_(p);
            flu[i++] = { C,   S,   B, X*C, X*S, Y };
            flu[i++] = { C*R, S*R, T, X*C, X*S, Y };
            p += inc;
        }
        flu[i++] = { -1, 0, B, X, 0, Y };
        flu[i++] = { -R, 0, T, X, 0, Y };
        assert_true( p == pi_4half );
        if ( top )
        {
            // repeat point to adjust normal:
            flu[i++] = { -R, 0, T, 0, 0, Z };
            // top disc, from 2*PI to PI:
            while ( p > pi_2half + inc )
            {
                p -= inc;
                float C = -cos_(p), S = sin_(p);
                flu[i++] = { C*R, S*R, T, 0, 0, Z };
                flu[i++] = { C*R,-S*R, T, 0, 0, Z };
            }
            assert_true( p == pi_2half + inc );
            flu[i++] = { R, 0, T, 0, 0, Z };
        }
        size_t j = nbTrianglesCylinder(inc, top);
        //std::clog << top << "   " << i << " " << j << "\n";
        assert_true( i == j );
        return i;
    }

    /// return essentially finesse * 24 / inc
    inline int nbTrianglesTube(unsigned inc)
    {
        return 2 + 2 * ( pi_4half / inc );
    }

    /// set triangle strip for a tube of constant radius 1 with Z in [B, T]
    size_t setTube(flute6* flu, unsigned inc, float B, float T)
    {
        assert_true( B <= T );
        size_t i = 0;
        for( unsigned p = 0; p <= pi_4half; p += inc )
        {
            float C = cos_(p), S = sin_(p);
            flu[i++] = { C, S, T, C, S, 0 };
            flu[i++] = { C, S, B, C, S, 0 };
        }
        size_t j = nbTrianglesTube(inc);
        assert_true( i <= j );
        return i;
    }
    
    /// set triangle strip for a tube of constant radius R with Z in [B, T]
    size_t setTube(flute6* flu, unsigned inc, float B, float T, float R)
    {
        assert_true( B <= T );
        size_t i = 0;
        for( unsigned p = 0; p <= pi_4half; p += inc )
        {
            float C = cos_(p), S = sin_(p);
            flu[i++] = { R*C, R*S, T, C, S, 0 };
            flu[i++] = { R*C, R*S, B, C, S, 0 };
        }
        return i;
    }
    
    /// set triangle strip for a cone of radius rB at Z=B, and rT at Z=T
    size_t setCone(flute6* flu, unsigned inc, float B, float rB, float T, float rT)
    {
        assert_true( B <= T );
        const float H(T-B), W(rT-rB);
        const float t(std::copysign(1.f,H)/sqrtf(H*H+W*W));
        const float x(t*H);
        const float y(t*W);
        size_t i = 0;
        for( unsigned n = 0; n <= pi_4half; n += inc )
        {
            float S = sin_(n), C = cos_(n);
            flu[i++] = { rT*C, rT*S, T, x*C, x*S, y };
            flu[i++] = { rB*C, rB*S, B, x*C, x*S, y };
        }
        return i;
    }
    
    /// return essentially finesse * 24 / inc
    inline unsigned nbTrianglesHelix(int turns)
    {
        return 2 + pi_4half * ( 2 * turns + 1 );
    }

    /// set triangle strip for a spiral, with Z in [Z, T]
    size_t setHelix(flute6* flu, float Z, float T, int turns)
    {
        assert_true( Z <= T );
        float W = ( T - Z ) / ( 1 + turns * 2 );
        float dZ = 2.0 * W / pi_4half;
        unsigned i = 0;
        // add quarter-turn with no increment in Z
        for ( unsigned p = pi_3half; p < pi_4half; ++p )
        {
            float C = cos_(p), S = sin_(p);
            flu[i++] = { C, S, Z+W, C, S, 0 };
            flu[i++] = { C, S, Z, C, S, 0 };
        }
        for ( int t = 0; t < turns; ++t )
        for ( unsigned p = 0; p < pi_4half; ++p )
        {
            float C = cos_(p), S = sin_(p);
            flu[i++] = { C, S, Z+W, C, S, 0 };
            flu[i++] = { C, S, Z, C, S, 0 };
            Z += dZ;
        }
        // add quarter-turn with no increment in Z
        for ( unsigned p = 0; p <= pi_1half; ++p )
        {
            float C = cos_(p), S = sin_(p);
            flu[i++] = { C, S, Z+W, C, S, 0 };
            flu[i++] = { C, S, Z, C, S, 0 };
        }
        unsigned j = nbTrianglesHelix(turns);
        assert_true( i <= j );
        return i;
    }

    /// set triangle strip for a disc at Z, with given normal in Z
    size_t setDisc(flute6* flu, unsigned inc, float Z, float N)
    {
        unsigned i = 0;
        flu[i++] = { 1, 0, Z, 0, 0, N };
        for( unsigned n = inc; n < pi_2half; n += inc )
        {
            float S = N * sin_(n), C = cos_(n);
            flu[i++] = { C,  S, Z, 0, 0, N };
            flu[i++] = { C, -S, Z, 0, 0, N };
        }
        flu[i++] = {-1, 0, Z, 0, 0, N };
        return i;
    }
    
    /// set triangle strip for a hollow disc at Z, with radius [ 1, R ] assuming R > 1
    size_t setRing(flute6* flu, unsigned inc, float R, float Z, float N)
    {
        unsigned i = 0;
        flu[i++] = { 1, 0, Z, 0, 0, N };
        flu[i++] = { R, 0, Z, 0, 0, N };
        for( unsigned n = inc; n < pi_4half; n += inc )
        {
            float S = sin_(n), C = cos_(n);
            flu[i++] = { C, S, Z, 0, 0, N };
            flu[i++] = { R*C, R*S, Z, 0, 0, N };
        }
        flu[i++] = { 1, 0, Z, 0, 0, N };
        flu[i++] = { R, 0, Z, 0, 0, N };
        return i;
    }

    /**
     Set triangles for a surface of revolution around the Z-axis.
     The surface goes from Z = [ B, T ] with increments dZ, and its radius is
     given by the function `radius(Z)` provided as argument.
     */
    unsigned setRevolution(flute6 * flu, float (*radius)(float), float B, const float T, float dZ)
    {
        unsigned i = 0;
        float Z = B, R, Q = radius(B);
        while ( B < T )
        {
            Z += dZ;
            R = radius(Z);
            
            float dR = ( R - Q ) / dZ;
            float dN = 1.0f / sqrtf( 1 + dR * dR );
            dR = -dR * dN;
            
            for ( unsigned n = 0; n <= pi_4half; ++n )
            {
                float S = sin_(n), C = cos_(n);
                flu[i++] = {R*C, R*S, Z, dN*C, dN*S, dR};
                flu[i++] = {Q*C, Q*S, B, dN*C, dN*S, dR};
            }
            flu[i++] = {Q, 0, B, dN, dN, dR};
            flu[i++] = {R, 0, Z, dN, dN, dR};
            B = Z;
            Q = R;
        }
        assert_true( i == 2*(2+pi_4half)*(1+std::ceil((T-B)/dZ)) );
        return i;
    }
    
    //-----------------------------------------------------------------------
    #pragma mark - Tubes buffers
    
    static size_t sizeTubeBuffers()
    {
        return 4 + 52 * pi_4half;  // this is empirical!
    }
    
    size_t setTubeBuffers(flute6* ptr, flute6* const ori)
    {
        /* The value of T limits the aspect ratio of tubes that can be drawn */
        const float B = -32.f, T = 256.f, E = 0.03125;
        unsigned i = 0, s = ptr - ori;
        tubes_[0] = i+s; i += setTube(ptr+i, 1, 0, 1);
        tubes_[1] = i+s; i += setTube(ptr+i, 2, 0, 1);
        tubes_[2] = i+s; i += setTube(ptr+i, 4, 0, 1);
        tubes_[3] = i+s; i += setTube(ptr+i, 8, 0, 1);

        tubes_[4] = i+s; i += setTube(ptr+i, 4, 0, 1+E); // tubeS
        tubes_[5] = i+s; i += setTube(ptr+i, 4,-E, 1+E); // tubeM
        tubes_[6] = i+s; i += setTube(ptr+i, 4,-E, 1);   // tubeE
        
        tubes_[7] = i+s; i += setTube(ptr+i, 1, B, T); // longTube1
        tubes_[8] = i+s; i += setTube(ptr+i, 2, B, T); // longTube2
        tubes_[9] = i+s; i += setTube(ptr+i, 4, B, T); // longTube4
        
        tubes_[10] = i+s; i += setTube(ptr+i, 1, 0, T); // halfTube1
        tubes_[11] = i+s; i += setTube(ptr+i, 2, 0, T); // halfTube2
        tubes_[12] = i+s; i += setTube(ptr+i, 4, 0, T); // halfTube4
        
        tubes_[13] = i+s; i += setCylinder(ptr+i, 1, 0, 1, 1, 1); //cylinder1
        tubes_[14] = i+s; i += setCylinder(ptr+i, 1, -1, 1, 1, 1); //cylinderC
        tubes_[15] = i+s; i += setCylinder(ptr+i, 2, 0, 1, 1, 0); //shutTube2
        tubes_[16] = i+s; i += setCylinder(ptr+i, 2, 0, T, 1, 0); //shutLongTube2

        tubes_[17] = i+s; i += setCylinder(ptr+i, 1, 0, 1, 0); // cone1
        tubes_[18] = i+s; i += setCylinder(ptr+i, 2, 0, 1, 0); // cone2
        tubes_[19] = i+s; i += setCylinder(ptr+i, 2,-1, 1, 0); // coneC
        tubes_[20] = i+s; i += setCylinder(ptr+i, 2, -1, 2, 0); // longCone
        tubes_[21] = i+s; i += setCone(ptr+i, 2, 0, 1, 1, 0.5); // cutCone
        tubes_[22] = i+s; i += setHelix(ptr+i, 0, 1, 5);

        tubes_[24] = i+s; i += setDisc(ptr+i, 1, 0, 1);
        tubes_[25] = i+s; i += setDisc(ptr+i, 2, 0, 1);
        tubes_[26] = i+s; i += setDisc(ptr+i, 1, 1, 1);
        tubes_[27] = i+s; i += setDisc(ptr+i, 2, 1, 1);
        tubes_[28] = i+s; i += setDisc(ptr+i, 1, 0, -1);
        tubes_[29] = i+s; i += setDisc(ptr+i, 2, 0, -1);
        tubes_[30] = i+s; i += setRing(ptr+i, 2, M_SQRT2, 0, 1);
        tubes_[31] = i+s; i += setRing(ptr+i, 2, 1.19, 0, 1);
        size_t j = sizeTubeBuffers();
        assert_true( i <= j );
        return i;
    }

    inline void doTubeStrip(size_t inx, unsigned cnt)
    {
        gym::bindBufferV3N3(buf_[0]);
        gym::drawTriangleStrip(tubes_[inx], cnt);
        gym::cleanupVN();
    }
    
    inline void doStripedTubeStrip(float line_width, size_t inx, unsigned cnt, gym_color col)
    {
        gym::bindBufferV3N3(buf_[0]);
        gym::drawTriangleStrip(tubes_[inx], cnt);
        gym::rebind();
        gym::disableLighting();
        gym::color(col);
        gym::drawLines(line_width, tubes_[inx], cnt);
        gym::restoreLighting();
        gym::cleanupVN();
    }

    // using Vertex Buffer Objects
    void tube1()         { doTubeStrip(0, nbTrianglesTube(1)); }
    void tube2()         { doTubeStrip(1, nbTrianglesTube(2)); }
    void tube4()         { doTubeStrip(2, nbTrianglesTube(4)); }
    void tube8()         { doTubeStrip(3, nbTrianglesTube(8)); }

    void tubeS()         { doTubeStrip(4, nbTrianglesTube(4)); }
    void tubeM()         { doTubeStrip(5, nbTrianglesTube(4)); }
    void tubeE()         { doTubeStrip(6, nbTrianglesTube(4)); }
    
    void longTube1()     { doTubeStrip(7, nbTrianglesTube(1)); }
    void longTube2()     { doTubeStrip(8, nbTrianglesTube(2)); }
    void longTube4()     { doTubeStrip(9, nbTrianglesTube(4)); }
    
    void halfTube1()     { doTubeStrip(10, nbTrianglesTube(1)); }
    void halfTube2()     { doTubeStrip(11, nbTrianglesTube(2)); }
    void halfTube4()     { doTubeStrip(12, nbTrianglesTube(4)); }
    
    void cylinder1()     { doTubeStrip(13, nbTrianglesCylinder(1, 1)); }
    void cylinderC()     { doTubeStrip(14, nbTrianglesCylinder(1, 1)); }
    void shutTube2()     { doTubeStrip(15, nbTrianglesCylinder(2, 0)); }
    void shutLongTube2() { doTubeStrip(16, nbTrianglesCylinder(2, 0)); }

    void cone1()         { doTubeStrip(17, nbTrianglesCylinder(1)); }
    void cone2()         { doTubeStrip(18, nbTrianglesCylinder(2)); }
    void coneC()         { doTubeStrip(19, nbTrianglesCylinder(2)); }
    void longCone()      { doTubeStrip(20, nbTrianglesCylinder(2)); }
    void cutCone()       { doTubeStrip(21, nbTrianglesTube(2)); }
    void helix()         { doTubeStrip(22, nbTrianglesHelix(5)); }

    void disc1()         { doTubeStrip(24, pi_4half); }
    void disc2()         { doTubeStrip(25, pi_4half/2); }
    void discTop1()      { doTubeStrip(26, pi_4half); }
    void discTop2()      { doTubeStrip(27, pi_4half/2); }
    void discBottom1()   { doTubeStrip(28, pi_4half); }
    void discBottom2()   { doTubeStrip(29, pi_4half/2); }
    void ring()          { doTubeStrip(30, 2+pi_4half); }
    void thinRing()      { doTubeStrip(31, 2+pi_4half); }

    void stripedTube(float w, gym_color col) { doStripedTubeStrip(w, 3, nbTrianglesTube(8), col); }

    inline void bindBufferV2() { gym::bindBufferV2(buf_[0]); }
    
    void stripes(float w) { gym::bindBufferV3(buf_[0], 2); gym::drawLines(w, tubes_[2], nbTrianglesTube(4)); }
    void circle1(float w) { bindBufferV2(); gym::drawLineStrip(w, discs_[0], 1+pi_4half); }
    void circle2(float w) { bindBufferV2(); gym::drawLineStrip(w, discs_[1], 1+pi_4half/2); }
    void dottedCircle(float w) { bindBufferV2(); gym::drawPoints(w, discs_[1], 1+pi_4half/2); }
    void square1(float w) { bindBufferV2(); gym::drawLineStrip(w, discs_[2], 5); }

    void circle(float R, float w) { gym::scale(R); circle1(w); gym::pull_ref(); }

    void strokeCapsule(float w) { bindBufferV2(); gym::drawLineStrip(w, discs_[4], 15); }
    void paintCapsule() { bindBufferV2(); gym::drawTriangleStrip(discs_[5], 15); }
    void strokeCross(float w) { bindBufferV2(); gym::drawLineStrip(w, discs_[6], 12); }
    void paintCross() { bindBufferV2(); gym::drawTriangleStrip(discs_[7], 12); }

    //-----------------------------------------------------------------------
#pragma mark - Spheres made from refined Icosahedrons
    
    /// using icosahedrons to render the sphere:
    static unsigned setIcoBuffer(Tesselator& ico, float*& ptr, Tesselator::INDEX*& idx)
    {
        //fprintf(stderr, "setIcoBuffer %i: %u %u\n", i, ico.max_vertices(), ico.num_vertices());
        assert_true(ico.num_vertices() <= 65535);
        ico.store_vertices(ptr);
        ptr += 3 * ico.num_vertices();
        
        unsigned cnt = 3 * ico.num_faces();
        memcpy(idx, ico.face_data(), ico.face_data_size());
        idx += cnt;
        return cnt;
    }
    
    void drawIcoBuffer(GLsizei pts, size_t inx, GLsizei cnt)
    {
        assert_true(buf_[0]); assert_true(buf_[1]);
        // the normal in each vertex is equal to the vertex!
        gym::bindBufferV3N0(buf_[0], pts);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buf_[1]);
        glDrawElements(GL_TRIANGLES, cnt, GL_UNSIGNED_SHORT, (void*)(sizeof(GLushort)*inx));
        gym::unbind2();
        gym::cleanupVN();
    }
    
    /** Assuming faces of the icosahedron have been properly sorted */
    void drawFootball(GLsizei pts, size_t inx, GLsizei cnt, gym_color col, GLsizei mid)
    {
        assert_true(buf_[0]); assert_true(buf_[1]);
        // the normal in each vertex is equal to the vertex!
        gym::bindBufferV3N0(buf_[0], pts);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buf_[1]);
        glDrawElements(GL_TRIANGLES, mid, GL_UNSIGNED_SHORT, (void*)(sizeof(GLushort)*inx));
        gym::color_both(col);
        glDrawElements(GL_TRIANGLES, cnt-mid, GL_UNSIGNED_SHORT, (void*)(sizeof(GLushort)*(inx+mid)));
        gym::unbind2();
        gym::cleanupVN();
    }
    
    /** Assuming faces of the icosahedron have been properly sorted */
    void drawFootballT(GLsizei pts, size_t inx, GLsizei cnt, gym_color col, GLsizei mid)
    {
        assert_true(buf_[0]); assert_true(buf_[1]);
        // the normal in each vertex is equal to the vertex!
        gym::bindBufferV3N0(buf_[0], pts);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buf_[1]);
        glDrawElements(GL_TRIANGLES, mid, GL_UNSIGNED_SHORT, (void*)(sizeof(GLushort)*inx));
        size_t edg = std::ceil(0.33*std::sqrt((cnt-mid)/180));
        edg = 180 * edg * edg;
        glDrawElements(GL_TRIANGLES, edg, GL_UNSIGNED_SHORT, (void*)(sizeof(GLushort)*(inx+cnt-edg)));
        gym::color_both(col);
        glDrawElements(GL_TRIANGLES, cnt-mid-edg, GL_UNSIGNED_SHORT, (void*)(sizeof(GLushort)*(inx+mid)));
        gym::unbind2();
        gym::cleanupVN();
    }

    void dualPassIcoBuffer(GLsizei pts, size_t inx, GLsizei cnt)
    {
        assertLighting();
        assertCullFace();
        assert_true(buf_[0]); assert_true(buf_[1]);
        // the normal in each vertex is equal to the vertex!
        gym::bindBufferV3N0(buf_[0], pts);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buf_[1]);
        gym::switchCullFace(GL_FRONT);
        glDrawElements(GL_TRIANGLES, cnt, GL_UNSIGNED_SHORT, (void*)(sizeof(GLushort)*inx));
        gym::switchCullFace(GL_BACK);
        glDrawElements(GL_TRIANGLES, cnt, GL_UNSIGNED_SHORT, (void*)(sizeof(GLushort)*inx));
        gym::unbind2();
        gym::cleanupVN();
    }

    void drawStripBuffer(GLsizei pts, size_t inx, GLsizei cnt)
    {
        assert_true(buf_[0]); assert_true(buf_[1]);
        gym::bindBufferV3N0(buf_[0], pts);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buf_[1]);
        glDrawElements(GL_TRIANGLE_STRIP, cnt, GL_UNSIGNED_SHORT, (void*)inx);
        gym::drawPoints(4, 0, 12);
        gym::unbind2();
        gym::cleanupVN();
    }

    void sphere1() { drawIcoBuffer(ico_pts_[0], ico_idx_[0], ico_cnt_[0]); }
    void sphere2() { drawIcoBuffer(ico_pts_[1], ico_idx_[1], ico_cnt_[1]); }
    void sphere4() { drawIcoBuffer(ico_pts_[2], ico_idx_[2], ico_cnt_[2]); }
    void sphere8() { drawIcoBuffer(ico_pts_[3], ico_idx_[3], ico_cnt_[3]); }
    
    void hemisphere1() { drawIcoBuffer(ico_pts_[4], ico_idx_[4], ico_cnt_[4]); }
    void hemisphere2() { drawIcoBuffer(ico_pts_[5], ico_idx_[5], ico_cnt_[5]); }

    void pin() { drawIcoBuffer(ico_pts_[6], ico_idx_[6], ico_cnt_[6]); }
    void droplet() { drawIcoBuffer(ico_pts_[7], ico_idx_[7], ico_cnt_[7]); }
    void droplet2() { drawIcoBuffer(ico_pts_[8], ico_idx_[8], ico_cnt_[8]); }
    
    void dome() { drawIcoBuffer(ico_pts_[DOME], ico_idx_[DOME], ico_cnt_[DOME]); }
    void roof() { drawIcoBuffer(ico_pts_[ROOF], ico_idx_[ROOF], ico_cnt_[ROOF]); }

    // dual pass routines can be used to draw transparent spheres:
    void dualPassSphere1() { dualPassIcoBuffer(ico_pts_[0], ico_idx_[0], ico_cnt_[0]); }
    void dualPassSphere2() { dualPassIcoBuffer(ico_pts_[1], ico_idx_[1], ico_cnt_[1]); }
    void dualPassSphere4() { dualPassIcoBuffer(ico_pts_[2], ico_idx_[2], ico_cnt_[2]); }
    void dualPassSphere8() { dualPassIcoBuffer(ico_pts_[3], ico_idx_[3], ico_cnt_[3]); }

    void football() { drawFootball(ico_pts_[FOOT], ico_idx_[FOOT], ico_cnt_[FOOT], gym_color(0,0,0), ico_mid_[FOOT]); }

    void football1(gym_color col) { drawFootball(ico_pts_[FOOT], ico_idx_[FOOT], ico_cnt_[FOOT], col, ico_mid_[FOOT]); }
    void footballT(gym_color col) { drawFootballT(ico_pts_[FOOT], ico_idx_[FOOT], ico_cnt_[FOOT], col, ico_mid_[FOOT]); }
    void football2(gym_color col) { drawFootball(ico_pts_[1], ico_idx_[1], ico_cnt_[1], col, ico_mid_[1]); }
    void football4(gym_color col) { drawFootball(ico_pts_[2], ico_idx_[2], ico_cnt_[2], col, ico_mid_[2]); }
    void football8(gym_color col) { drawFootball(ico_pts_[3], ico_idx_[3], ico_cnt_[3], col, ico_mid_[3]); }

    void icoid() { drawStripBuffer(ico_pts_[ICOID], ico_idx_[ICOID], ico_cnt_[ICOID]); }
    // this draws without loading/initializing the buffer
    void icoidF() { glDrawElements(GL_TRIANGLE_STRIP, ico_cnt_[ICOID], GL_UNSIGNED_SHORT, (void*)ico_idx_[ICOID]); }
    
    void createBuffers()
    {
        Tesselator ico[10];
        ico[0].buildIcosahedron(finesse*9);
        ico[1].buildIcosahedron(finesse*6);
        ico[2].buildIcosahedron(finesse*3);
        ico[3].buildIcosahedron(finesse);

        ico[4].buildHemisphere(finesse*3);
        ico[5].buildHemisphere(finesse);
        ico[6].buildPin(finesse*3);
        ico[7].buildDroplet(finesse*3);
        ico[8].buildDroplet(finesse);
        ico[FOOT].buildFootball(finesse*9);

        size_t F = 32; // reserve for ICOID (30)
        size_t S = 12;
        for ( int i = 0; i < 10; ++i )
        {
            F += 3 * ico[i].max_faces();
            S += 3 * ico[i].max_vertices();
        }
        const int DOME_ICO = 5;
        // reserve for DOME
        F += 3 * ico[DOME_ICO].max_faces();
        S += 3 * ico[DOME_ICO].max_vertices();
        // reserve for ROOF
        F += 3 * ico[DOME_ICO].max_faces();
        S += 3 * ico[DOME_ICO].max_vertices();

        // required buffer size:
        size_t T = 6 * sizeTubeBuffers();
        size_t C = 6 * sizeCubeBuffers();
        size_t B = 3 * sizeBlobBuffers();
        size_t O = 2 * sizeCircBuffers();

        // request buffer for vertex data:
        glBufferData(GL_ARRAY_BUFFER, (S+T+C+B+O)*sizeof(float), nullptr, GL_STATIC_DRAW);
        float* ptr = (float*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
        
        // request buffer for index data:
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, F*sizeof(GLushort), nullptr, GL_STATIC_DRAW);
        GLushort* idx = (GLushort*)glMapBuffer(GL_ELEMENT_ARRAY_BUFFER, GL_WRITE_ONLY);
        static_assert(std::is_same<Tesselator::INDEX, GLushort>::value, "Index type mismatch");

        float*const ptr0 = ptr;
        GLushort*const idx0 = idx;
        
        ptr += 6 * setTubeBuffers((flute6*)ptr, (flute6*)ptr0);
        //fprintf(stderr, "setTubeBuffers : %li %li\n", ptr-ptr0, t); float* sub=ptr;
        assert_true( ptr <= ptr0 + T );

        ptr += 6 * setCubeBuffers((flute6*)ptr, (flute6*)ptr0);
        //fprintf(stderr, "setCubeBuffer : %li %li\n", ptr-sub, c); sub=ptr;
        assert_true( ptr <= ptr0 + C + T );

        for ( int i = 0; i < 10; ++i )
        {
            ico_pts_[i] = ptr - ptr0;
            ico_idx_[i] = idx - idx0;
            unsigned cnt = setIcoBuffer(ico[i], ptr, idx);
            ico_cnt_[i] = cnt;
            ico_mid_[i] = cnt - 3 * ico[i].num_foot_faces();
        }

        assert_true(DOME >= 10);
        float * dome = ptr;
        ico_pts_[DOME] = ptr - ptr0;
        ico_idx_[DOME] = idx - idx0;
        ico_cnt_[DOME] = setIcoBuffer(ico[DOME_ICO], ptr, idx);
        Tesselator::scale3D(ico_cnt_[DOME], dome, 1.f, 1.f, 0.5f);
        
        assert_true(ROOF >= DOME);
        float * roof = ptr;
        ico_pts_[ROOF] = ptr - ptr0;
        ico_idx_[ROOF] = idx - idx0;
        ico_cnt_[ROOF] = setIcoBuffer(ico[DOME_ICO], ptr, idx);
        Tesselator::scale3D(ico_cnt_[ROOF], roof, 1.f, 1.f, 0.33f);

        assert_true(ICOID >= ROOF);
        ico_pts_[ICOID] = ptr - ptr0;
        ico_idx_[ICOID] = idx - idx0;
        ico_cnt_[ICOID] = setIcoidBuffer((flute3*)ptr, idx);
        ptr += 3*12;
        idx += ico_cnt_[ICOID];
        
        //fprintf(stderr, "setIcoBuffer : %li %li -- %li\n", ptr-sub, s, idx-idx0); sub=ptr;
        assert_true( ptr <= ptr0 + S + C + T );
        assert_true( idx <= idx0 + F );

        ptr += 3 * setBlobBuffers((flute3*)ptr, (flute3*)ptr0);
        //fprintf(stderr, "setBlobBuffers : %li %li\n", ptr-sub, b); sub=ptr;
        assert_true( ptr <= ptr0 + B + S + C + T );

        // align pointer:
        ptr += (ptr-ptr0) & 1;
        
        ptr += 2 * setCircBuffers((flute2*)ptr, (flute2*)ptr0);
        //fprintf(stderr, "setCircBuffers : %li %li\n", ptr-sub, o);

        assert_true( ptr <= ptr0 + O + B + S + C + T );
        glUnmapBuffer(GL_ELEMENT_ARRAY_BUFFER);
        glUnmapBuffer(GL_ARRAY_BUFFER);
    }

    void initBuffers()
    {
        glGenBuffers(2, buf_);
        glBindBuffer(GL_ARRAY_BUFFER, buf_[0]);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buf_[1]);
        gle::createBuffers();
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
            _exit(1);
        }
#endif
        // circle_[] covers 2 revolutions = 4 * PI
        //set_arc(2*pi_4half, circle_, 1, 0, 2*M_PI/pi_4half, 0, 0);
        compute_circle(pi_1half, circle_, 1, 0);
        compute_circle(pi_1half, circle_+4*pi_2half, 1, 0);

        if ( !glIsBuffer(buf_[0]) )
        {
            initBuffers();
            gym::initStreams();
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
    #pragma mark - Sphere decorations

    /**
     Draw a cylindrical band on the equator of a sphere of radius 1.
     The band is in the XY plane. The axis of the cylinder is Z.
     The band is made of triangles indicating the clockwise direction.
     */
    void arrowStrip(float width, const unsigned inc)
    {
        float A(M_PI * inc / pi_1half);
        float W(width * A / M_SQRT3);
        float R(1.0f / cosf(A*0.5f));
        
        flute6 * flu = gym::mapBufferV3N3(3*(1+pi_4half/(2*inc)));
        unsigned i = 0;
        flu[i++] = {R, 0, W, 1, 0, 0};
        flu[i++] = {R, 0,-W, 1, 0, 0};
        for ( unsigned n = 0; n < pi_4half; n += 2*inc )
        {
            float c = R * cos_(n);
            float s = R * sin_(n);
            flu[i++] = {c, s, 0, c, s, 0};
            flu[i++] = {c, s, W, c, s, 0};
            flu[i++] = {c, s,-W, c, s, 0};
        }
        flu[i++] = {R, 0, 0, 1, 0, 0};
        gym::unmapBufferV3N3();
        gym::drawTriangles(0, i);
        gym::cleanupVN();
    }
    
    
    void threeArrowStrip(float w, const unsigned inc)
    {
        arrowStrip(w, inc);
        gym::rotateX(0, -1);
        arrowStrip(w, inc);
        gym::rotateY(0, 1);
        arrowStrip(w, inc);
    }
    
    /** The path followed by the seam of a tennis ball */
    unsigned setSeamCurve(flute3 * flu, float R, float B, float W, unsigned inc)
    {
        float A = R - B;
        float D = 2 * std::sqrt(A*B);
        unsigned cnt = 0;
        for ( unsigned i = 0; i <= pi_4half; i += inc )
        {
            float C = cos_(i), S = sin_(i);
            float C2 = C*C-S*S, S2 = S*C+C*S;
            float C3 = C*C2-S*S2, S3 = S*C2+C*S2;
            flute3 P = { A*S+B*S3, A*C-B*C3, D*C2 };
            flute3 T = { A*C+3*B*C3, -A*S+3*B*S3, -2*D*S2 };
            flute3 N = W * cross(P, T);
            flu[cnt++] = P + N;
            flu[cnt++] = P - N;
        }
        return cnt;
    }

    /** The path followed by the seam of a tennis ball */
    void baseballSeamCurve(float R, float W)
    {
        flute3 * flu = gym::mapBufferV3(2*pi_4half);
        unsigned cnt = setSeamCurve(flu, R, R*0.3, W/28, 1);
        gym::unmapBufferV3();
        gym::drawLineStrip(4, 0, cnt);
    }

    /** The path followed by the seam of a tennis ball */
    void tennisballSeamCurve(float R, float W)
    {
        flute3 * flu = gym::mapBufferV3(2*pi_4half);
        unsigned cnt = setSeamCurve(flu, R, R*0.3, W/28, 1);
        gym::unmapBufferV3();
        gym::drawTriangleStrip(0, cnt);
    }
    
    void baseball()
    {
        gym::color_front(1,1,1);
        sphere1();
        gym::color_front(1,0,0);
        baseballSeamCurve(1, 2);
    }

    void tennisball()
    {
        gym::color_front(1,1,0);
        sphere1();
        gym::color_front(1,1,1);
        tennisballSeamCurve(1.02, 1);
    }
    
    //-----------------------------------------------------------------------
    #pragma mark - 3D volume primitives
    
    /// draw a Torus of radius R and a thickness 2*T
    void torusZ(float R, float T, unsigned inc)
    {
        for ( unsigned n = 0; n < pi_4half; n += inc )
        {
            flute6 * flu = gym::mapBufferV3N3(2+pi_4half);
            float X0 = cos_(n), X1 = cos_(n+inc);
            float Y0 = sin_(n), Y1 = sin_(n+inc);
            size_t i = 0;
            for ( unsigned p = 0; p <= pi_4half; p += 2*inc )
            {
                float S = sin_(p), C = cos_(p);
                flu[i++] = {X0*(R+T*C), Y0*(R+T*C), T*S, X0*C, Y0*C, S};
                flu[i++] = {X1*(R+T*C), Y1*(R+T*C), T*S, X1*C, Y1*C, S};
            }
            gym::unmapBufferV3N3();
            gym::drawTriangleStrip(0, i);
        }
        gym::cleanupVN();
    }
    
    /**
     Draw a surface of revolution around the Z-axis.
     The surface goes from Z = [ B, T ] with increments dZ, and its radius is
     given by the function `radius`(Z) provided as argument.
     */
    void drawRevolution(float (*radius)(float), float B, const float T, float dZ)
    {
        unsigned cnt = 2*(2+pi_4half)*(1+std::ceil((T-B)/dZ));
        flute6 * flu = gym::mapBufferV3N3(cnt+4);
        unsigned i = setRevolution(flu, radius, B, T, dZ);
        assert_true( i <= cnt );
        gym::unmapBufferV3N3();
        gym::drawTriangleStrip(0, cnt);
        gym::cleanupVN();
    }

    /// some volume of revolution with axis along Z
    float barrelRadius(float z) { return sinf(M_PI*z); }
    void barrel() { drawRevolution(barrelRadius, 0, 1, 0.0625); }

    float dumbbellRadius(float z) { return sinf(M_PI*z) * (1.3f+cosf(2*M_PI*z)); }
        void dumbbell() { drawRevolution(dumbbellRadius, 0, 1, 0.0625); }

    void dualPassBarrel() { gym::dualPass(barrel); }

    //-----------------------------------------------------------------------
    
    void drawBand(Vector2 const& A, Vector2 const& B, real rad)
    {
        Vector2 d = ( B - A ).orthogonal();
        real n = d.norm();
        if ( n > 0 )
        {
            d *= rad / n;
            flute2 * pts = gym::mapBufferV2(4);
            pts[0] = A + d;
            pts[1] = A - d;
            pts[2] = B + d;
            pts[3] = B - d;
            gym::unmapBufferV2();
            gym::ref_view();
            gym::drawTriangleStrip(0, 4);
        }
    }
    
    void drawSpikyBand(Vector2 const& A, Vector2 const& B, real rad)
    {
        Vector2 t = B - A;
        Vector2 d = t.orthogonal();
        real n = t.norm();
        if ( n > 0 )
        {
            d *= rad / n;
            t *= rad / n;
            flute2 * pts = gym::mapBufferV2(6);
            pts[0] = A - t;
            pts[1] = A - d;
            pts[2] = A + d;
            pts[3] = B - d;
            pts[4] = B + d;
            pts[5] = B + t;
            gym::unmapBufferV2();
            gym::ref_view();
            gym::drawTriangleStrip(0, 6);
        }
    }

    
    void drawBand(Vector1 const& A, float rA,
                  Vector1 const& B, float rB)
    {
        float AX(A.XX);
        float BX(B.XX);
        flute2 * pts = gym::mapBufferV2(4);
        pts[0] = { AX, rA };
        pts[1] = { AX,-rA };
        pts[2] = { BX, rB };
        pts[3] = { BX,-rB };
        gym::unmapBufferV2();
        gym::ref_view();
        gym::drawTriangleStrip(0, 4);
    }
    
    void drawBand(Vector2 const& A, float rA,
                  Vector2 const& B, float rB)
    {
        Vector2 d = ( B - A ).orthogonal();
        real n = d.norm();
        if ( n > 0 )
        {
            float dX(d.XX/n), dY(d.YY/n);
            float AX(A.XX), AY(A.YY);
            float BX(B.XX), BY(B.YY);
            flute2 * pts = gym::mapBufferV2(4);
            pts[0] = { AX+rA*dX, AY+rA*dY };
            pts[2] = { AX-rA*dX, AY-rA*dY };
            pts[4] = { BX+rB*dX, BY+rB*dY };
            pts[7] = { BX-rB*dX, BY-rB*dY };
            gym::unmapBufferV2();
            gym::ref_view();
            gym::drawTriangleStrip(0, 4);
        }
    }
    
    void drawBand(Vector1 const& A, float rA, gym_color cA,
                  Vector1 const& B, float rB, gym_color cB)
    {
        float AX(A.XX);
        float BX(B.XX);
        flute6 * flu = gym::mapBufferC4V2(4);
        flu[0] = { cA, AX, -rA };
        flu[1] = { cA, AX,  rA };
        flu[2] = { cB, BX, -rB };
        flu[3] = { cB, BX,  rB };
        gym::unmapBufferC4V2();
        gym::ref_view();
        gym::drawTriangleStrip(0, 4);
        gym::cleanupCV();
    }
    
    void drawBand(Vector2 const& A, float rA, gym_color cA,
                  Vector2 const& B, float rB, gym_color cB)
    {
        Vector2 d = ( B - A ).orthogonal();
        real n = d.norm();
        if ( n > 0 )
        {
            float dX(d.XX/n), dY(d.YY/n);
            float AX(A.XX), AY(A.YY);
            float BX(B.XX), BY(B.YY);
            flute6 * flu = gym::mapBufferC4V2(6);
            flu[0] = { cA, AX+rA*dX, AY+rA*dY };
            flu[1] = { cA, AX-rA*dX, AY-rA*dY };
            flu[2] = { cB, BX+rB*dX, BY+rB*dY };
            flu[3] = { cB, BX-rB*dX, BY-rB*dY };
            gym::unmapBufferC4V2();
            gym::ref_view();
            gym::drawTriangleStrip(0, 4);
            gym::cleanupCV();
        }
    }
    
    
    /**
     This will displays a rectangle if the connection is parallel,
     and a hourglass if the connection is antiparallel
     */
    void drawHourglass(Vector2 const& A, Vector2 const& dA,
                       Vector2 const& B, Vector2 const& dB)
    {
        flute2 * pts = gym::mapBufferV2(6);
        pts[0] = B-dB;
        pts[1] = B;
        pts[2] = A-dA;
        pts[3] = A+dA;
        pts[4] = B;
        pts[5] = B+dB;
        gym::unmapBufferV2();
        gym::ref_view();
        gym::drawTriangleStrip(0, 6);
    }
    
    /**
     This will displays a rectangle if the connection is parallel,
     and a hourglass if the connection is antiparallel
     */
    void drawHourglass(Vector2 const& a, Vector2 const& da, gym_color cA,
                       Vector2 const& b, Vector2 const& db, gym_color cB)
    {
        flute6 * flu = gym::mapBufferC4V2(6);
        flu[0] = { cB, b-db };
        flu[1] = { cB, b };
        flu[2] = { cA, a-da };
        flu[3] = { cA, a+da };
        flu[4] = { cB, b };
        flu[5] = { cB, b+db };
        gym::unmapBufferC4V2();
        gym::drawTriangleStrip(0, 6);
        gym::cleanupCV();
    }
    
    
    /**
     This will displays a rectangle if the connection is antiparallel,
     and a hourglass if the connection is parallel
     */
    void drawCross(Vector2 const& A, Vector2 const& dA,
                   Vector2 const& B, Vector2 const& dB, real rad)
    {
        flute2 * pts = gym::mapBufferV2(8);
        pts[0] = A-rad*dA;
        pts[1] = A;
        pts[2] = B;
        pts[3] = B-rad*dB;
        pts[4] = A+rad*dA;
        pts[5] = A;
        pts[6] = B;
        pts[7] = B+rad*dB;
        gym::unmapBufferV2();
        gym::ref_view();
        gym::drawTriangleStrip(0, 4);
        gym::drawTriangleStrip(4, 4);
        gym::cleanupV();
    }
    
    void drawBar(Vector3 const& A, Vector3 const& dA,
                 Vector3 const& B, Vector3 const& dB, real rad)
    {
        Vector3 ab = normalize( A - B );
        Vector3 ea = cross(ab, dA);
        Vector3 eb = cross(ab, dB);
        flute3 * pts = gym::mapBufferV3(16);
        pts[0] = A-rad*(dA-ea);
        pts[1] = A-rad*(dA+ea);
        pts[2] = B-rad*(dB-eb);
        pts[3] = B-rad*(dB+eb);
        pts[4] = A+rad*(dA-ea);
        pts[5] = A+rad*(dA+ea);
        pts[6] = B+rad*(dB-eb);
        pts[7] = B+rad*(dB+eb);
        pts[8] = A-rad*dA;
        pts[9] = A+rad*dA;
        pts[10] = B-rad*dB;
        pts[11] = B+rad*dB;
        pts[12] = A-rad*dA;
        pts[13] = A+rad*dA;
        pts[14] = B-rad*dB;
        pts[15] = B+rad*dB;
        gym::unmapBufferV3();
        gym::ref_view();
        gym::drawTriangleStrip(0, 4);
        gym::drawTriangleStrip(4, 4);
        gym::drawTriangleStrip(8, 4);
        gym::drawTriangleStrip(12, 4);
        gym::cleanupV();
    }
    
    
    /**
     Two hexagons linked by a rectangle
     hexagons have the same surface as a disc of radius 1.
     */
    void drawDumbbell(Vector2 const& A, Vector2 const& B, float diameter)
    {
        const float S(1.0996361107912678f); //sqrt( 2 * M_PI / ( 3 * sqrt(3) ));
        const float R(diameter * S);
        Vector2 x = ( B - A ).normalized(R*0.8660254037844386f);
        Vector2 y = x.orthogonal(R*0.5f);
        flute2 * pts = gym::mapBufferV2(12);
        pts[0] = A-x-y;
        pts[1] = A-x+y;
        pts[2] = A-y-y;
        pts[3] = A+y+y;
        pts[4] = A+x-y;
        pts[5] = A+x+y;
        pts[6] = B-x-y;
        pts[7] = B-x+y;
        pts[8] = B-y-y;
        pts[9] = B+y+y;
        pts[10] = B+x-y;
        pts[11] = B+x+y;
        gym::unmapBufferV2();
        gym::ref_view();
        gym::drawTriangleStrip(0, 12);
        gym::cleanupV();
    }
    
    void drawTipi(real* ref, int inx, real rad)
    {
        flute3 S(ref+3*inx);
        flute3 A(ref+3*inx+1);
        flute3 B(ref+3*inx+2);
        flute3 C(ref+3*inx+3);
        flute6 * pts = gym::mapBufferV3N3(18);
        pts[0] = {S, B};
        gym::unmapBufferV3N3();
        gym::ref_view();
        gym::drawTriangles(0, 4);
        gym::cleanupVN();
    }

    //-----------------------------------------------------------------------
#pragma mark - Arrows
    
    void drawCone(Vector1 const& pos, Vector1 const& dir, const float rad)
    {
        float dx = rad * dir.XX, cx = pos.XX;
        flute2 * pts = gym::mapBufferV2(4);
        pts[0] = {cx-dx, dx};
        pts[1] = {cx-dx/2, 0};
        pts[2] = {cx+dx+dx, 0};
        pts[3] = {cx-dx,-dx };
        gym::unmapBufferV2();
        gym::ref_view();
        gym::drawTriangleStrip(0, 4);
    }
    
    void drawCone(Vector2 const& pos, Vector2 const& dir, const float rad)
    {
        float dx(rad*dir.XX), cx(pos.XX);
        float dy(rad*dir.YY), cy(pos.YY);
        float dxy = dx + dy, dyx = dy - dx;
        flute2 * pts = gym::mapBufferV2(4);
        pts[0] = {cx-dxy, cy-dyx};
        pts[1] = {cx-dx/2, cy-dy/2};
        pts[2] = {cx+2*dx, cy+2*dy};
        pts[3] = {cx+dyx, cy-dxy};
        gym::unmapBufferV2();
        gym::ref_view();
        gym::drawTriangleStrip(0, 4);
    }
    
    void drawCone(Vector3 const& pos, Vector3 const& dir, const float rad)
    {
        gym::transAlignZ(pos, rad, dir);
        longCone();
    }
    
    //-----------------------------------------------------------------------
    
    void drawCylinder(Vector1 const& pos, Vector1 const& dir, float rad)
    {
        float cx(pos.XX);
        float dx(rad * dir.XX * 0.5);
        flute2 * pts = gym::mapBufferV2(4);
        pts[0] = {cx-dx,-rad};
        pts[1] = {cx-dx, rad};
        pts[2] = {cx+dx,-rad};
        pts[3] = {cx+dx, rad};
        gym::unmapBufferV2();
        gym::ref_view();
        gym::drawTriangleStrip(0, 4);
    }
    
    void drawCylinder(Vector2 const& pos, Vector2 const& dir, float rad)
    {
        float dx(rad * dir.XX), cx(pos.XX - dx * 0.5);
        float dy(rad * dir.YY), cy(pos.YY - dy * 0.5);
        flute2 * pts = gym::mapBufferV2(4);
        pts[0] = {cx+dy, cy-dx};
        pts[1] = {cx-dy, cy+dx};
        pts[2] = {cx+dx+dy, cy+dy-dx};
        pts[3] = {cx+dx-dy, cy+dy+dx};
        gym::unmapBufferV2();
        gym::ref_view();
        gym::drawTriangleStrip(0, 4);
    }
    
    void drawCylinder(Vector3 const& pos, Vector3 const& dir, float rad)
    {
        gym::transAlignZ(pos, rad, dir);
        cylinderC();
    }
    
    
    //-----------------------------------------------------------------------
    
    void drawArrowTail(Vector1 const& pos, Vector1 const& dir, float rad)
    {
        float dx(rad * dir.XX);
        float cx(pos.XX - dx * 0.5 );
        flute2 * pts = gym::mapBufferV2(6);
        pts[0] = {cx-dx, -dx};
        pts[1] = {cx+dx, -dx};
        pts[2] = {cx, 0};
        pts[3] = {cx+2*dx, 0};
        pts[4] = {cx-dx, dx};
        pts[5] = {cx+dx, dx};
        gym::unmapBufferV2();
        gym::ref_view();
        gym::drawTriangleStrip(0, 6);
    }
    
    void drawArrowTail(Vector2 const& pos, Vector2 const& dir, float rad)
    {
        float dx(rad * dir.XX);
        float dy(rad * dir.YY);
        float cx(pos.XX - 1.5f * dx);
        float cy(pos.YY - 1.5f * dy);
        float ex(cx + 2 * dx);
        float ey(cy + 2 * dy);
        flute2 * pts = gym::mapBufferV2(6);
        pts[0] = {cx+dy, cy-dx};
        pts[1] = {ex+dy, ey-dx};
        pts[2] = {cx+dx, cy+dy};
        pts[3] = {ex+dx, ey+dy};
        pts[4] = {cx-dy, cy+dx};
        pts[5] = {ex-dy, ey+dx};
        gym::unmapBufferV2();
        gym::ref_view();
        gym::drawTriangleStrip(0, 6);
    }
    
    void drawArrowTail(Vector3 const& pos, Vector3 const& dir, float rad)
    {
        gym::transAlignZ(pos, rad, dir);
        arrowTail();
    }
    
    //-----------------------------------------------------------------------
    void drawArrow(Vector1 const& A, Vector1 const& B, float R)
    {
        gym::stretchAlignZ(A, B, R);
        tube1();
        gym::shift(0, 0, 1-3*R);
        gym::scale(3.0, 3.0, 6*R);
        cone2();
    }
    
    void drawArrow(Vector2 const& A, Vector2 const& B, float R)
    {
        gym::stretchAlignZ(A, B, R);
        tube1();
        gym::shift(0, 0, 1-3*R);
        gym::scale(3.0, 3.0, 6*R);
        cone2();
    }
    
    void drawArrow(Vector3 const& A, Vector3 const& B, float R)
    {
        gym::stretchAlignZ(A, B, R);
        tube1();
        gym::transAlignZ(B, 2*R, B-A);
        gym::scale(1.0, 1.0, 2.0);
        cone2();
    }
    
    //-----------------------------------------------------------------------
    void drawAxes(const float L, int dim)
    {
        const float R(L*0.1f);
        // display a white ball at the origin
        gym::color_front(1, 1, 1);
        gym::pull_ref();
        gym::scale(R);
        sphere2();
        for (int d = 0; d < dim; ++d)
        {
            gym::pull_ref();
            switch(d)
            {
                case 0:
                    gym::color_front(1, 0, 0);
                    gym::rotateY(0, 1);
                    break;
                case 1:
                    gym::color_front(0, 1, 0);
                    gym::rotateX(0, -1);
                    gym::rotateZ(-1, 0);
                    break;
                case 2:
                    gym::color_front(0, 0, 1);
                    gym::rotateZ(0, -1);
                    break;
            }
            gym::shift(0, 0, 8*R);
            gym::scale(1.5*R, 1.5*R, 3*R);
            cone2();
            gym::shift(0, 0, -2.4);
            gym::scale(0.26, 0.26, 2.4);
            tube1();
        }
        gym::load_ref();
    }
    
    /// draw two rectangles parallel with the XY plane
    void strokeCuboid(Vector3 const& A, Vector3 const& B, float w)
    {
        float AX = A.XX, AY = A.YY, AZ = A.ZZ;
        float BX = B.XX, BY = B.YY, BZ = B.ZZ;
        flute3 * flu = gym::mapBufferV3(10);
        flu[0] = { AX, AY, AZ };
        flu[1] = { BX, AY, AZ };
        flu[2] = { BX, BY, AZ };
        flu[3] = { AX, BY, AZ };
        flu[4] = { AX, AY, AZ };
        flu[5] = { AX, AY, BZ };
        flu[6] = { BX, AY, BZ };
        flu[7] = { BX, BY, BZ };
        flu[8] = { AX, BY, BZ };
        flu[9] = { AX, AY, BZ };
        gym::unmapBufferV3();
        gym::drawLineStrip(w, 0, 5);
        gym::drawLineStrip(w, 5, 5);
    }

    /// draw faces of cuboid of axis [A,B] and size `rad`
    void paintCuboid(Vector3 const& A3, Vector3 const& B3, float rad)
    {
        flute3 A(A3), B(B3);
        flute3 X, Y, AB(B-A);
        float n = norm(AB);
        gym::orthonormal(AB, n, X, Y);
        X = X * ( rad / n );
        Y = Y * ( rad / n );
        flute6 * flu = gym::mapBufferV3N3(16);
        flu[0] = { A + X - Y, X };
        flu[1] = { A + X + Y, X };
        flu[2] = { B + X - Y, X };
        flu[3] = { B + X + Y, X };
        flu[4] = { B - X - Y,-X };
        flu[5] = { B - X + Y,-X };
        flu[6] = { A - X - Y,-X };
        flu[7] = { A - X + Y,-X };
        
        flu[ 8] = { B - X + Y, Y };
        flu[ 9] = { B + X + Y, Y };
        flu[10] = { A - X + Y, Y };
        flu[11] = { A + X + Y, Y };
        flu[12] = { A - X - Y,-Y };
        flu[13] = { A + X - Y,-Y };
        flu[14] = { B - X - Y,-Y };
        flu[15] = { B + X - Y,-Y };
        gym::unmapBufferV3N3();
        gym::drawTriangleStrip(0, 8);
        gym::drawTriangleStrip(8, 8);
        gym::cleanupVN();
    }

    /// draw faces of cuboid of axis [A,B] and size `rad`
    void paintTetrahedron(Vector3 const& X3, Vector3 const& Y3, Vector3 const& Z3,
                          Vector3 const& P3)
    {
        flute3 X(X3), Y(Y3), Z(Z3), P(P3);
        // translate to place center of base at P
        P = P - ( X + Y + Z ) * 0.33333333f;
        flute3 * flu = gym::mapBufferV3(8);
        flu[0] = { P + X };
        flu[1] = { P + Y };
        flu[2] = { P + Z };
        flu[3] = { P };
        flu[4] = { P + X };
        flu[5] = { P + Y };
        gym::unmapBufferV3();
        gym::drawTriangleStrip(0, 6);
        gym::cleanupV();
    }

    /// draw triangular prism extending over [P3,Q3] with directions XYZ and ABC
    void paintPrism(Vector3 const& X3, Vector3 const& Y3, Vector3 const& Z3, Vector3 const& P3,
                    Vector3 const& A3, Vector3 const& B3, Vector3 const& C3, Vector3 const& Q3)
    {
        flute3 X(X3), Y(Y3), Z(Z3), P(P3);
        flute3 A(A3), B(B3), C(C3), Q(Q3);
        // translate to place center of base at P
        P = P - ( X + Y + Z ) * 0.33333333f;
        Q = Q - ( A + B + C ) * 0.33333333f;
        flute3 * flu = gym::mapBufferV3(10);
        flu[0] = { P + X };
        flu[1] = { P + Y };
        flu[2] = { Q + A };
        flu[3] = { Q + B };
        flu[4] = { Q + C };
        flu[5] = { P + Y };
        flu[6] = { P + Z };
        flu[7] = { P + X };
        flu[8] = { Q + C };
        flu[9] = { Q + A };
        gym::unmapBufferV3();
        gym::drawTriangleStrip(0, 10);
        gym::cleanupV();
    }

    /// draw triangular prism extending over [P3,Q3] with directions XYZ and ABC
    void paintSpikyPrism(Vector3 const& X3, Vector3 const& Y3, Vector3 const& Z3, Vector3 const& P3,
                         Vector3 const& A3, Vector3 const& B3, Vector3 const& C3, Vector3 const& Q3)
    {
        flute3 X(X3), Y(Y3), Z(Z3), P(P3);
        flute3 A(A3), B(B3), C(C3), Q(Q3);
        // translate to place center of base at P
        P = P - ( X + Y + Z ) * 0.33333333f;
        Q = Q - ( A + B + C ) * 0.33333333f;
        flute3 * flu = gym::mapBufferV3(14);
        flu[0] = { P + X };
        flu[1] = { P + Y };
        flu[2] = { Q + A };
        flu[3] = { Q + B };
        flu[4] = { Q };
        flu[5] = { Q + C };
        flu[6] = { Q + A };
        flu[7] = { P + Z };
        flu[8] = { P + X };
        flu[9] = { P };
        flu[10] = { P + Y };
        flu[11] = { P + Z };
        flu[12] = { Q + B };
        flu[13] = { Q + C };
        gym::unmapBufferV3();
        gym::drawTriangleStrip(0, 14);
        gym::cleanupV();
    }

}

