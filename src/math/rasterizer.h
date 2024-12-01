// Cytosim was created by Francois Nedelec. Copyright Cambridge University 2020
#ifndef RASTERIZER_H
#define RASTERIZER_H

#include "real.h"
#include "vector1.h"
#include "vector2.h"
#include "vector3.h"

/// 2D and 3D rasterizer
/**
 The various rasterizer methods call a given function `func(P)` for every
 point `P` with INTEGER coordinates that are inside a certain volume.
 
 The volume can be specified:
 - as a polygon described by a list of points, using paintPolygon?D(polygon) 
 - as a cylinder specified by two points [P,Q] and a scalar 'radius'
 .
 [P,Q] or the points defining the polygon do not need to be integers.
 
 Note that in 3D, the functions do not rasterize the cylinder exactly, but a
 generalized cylinder of axis [PQ], with a rectangular or hexagonal crosssection.
 This slightly larger volume contains all the points located at distance `radius`
 or less from [PQ].

 F.Nedelec, EMBL 2002-2017, Cambridge 2019-- nedelec@slcu.cam.ac.uk

*/
namespace Rasterizer 
{
    /// type used by the Rasterizer, for which single precision is sufficient
    typedef float FLOAT;

    /// a point in 2D
    struct Vertex2
    {
        /// coordinates of the point
        FLOAT XX, YY;
        
        Vertex2()
        {
        }
        
        Vertex2(Vector2 const& vec)
        {
            XX = (FLOAT)vec.XX;
            YY = (FLOAT)vec.YY;
        }
        
        void print(std::ostream& os) const
        {
            os << "( " << std::setw(9) << XX;
            os << "  " << std::setw(9) << YY << " )";
        }
    };

    /// a point in 2D, with slopes in Z
    struct Vertex2dZ
    {
        /// coordinates of the point
        FLOAT XX, YY;
        
        /// derivatives of the line with respect to Z
        FLOAT dX, dY;
        
        Vertex2dZ()
        {
        }

        void fix(FLOAT x, FLOAT dx, FLOAT y, FLOAT dy)
        {
            XX = x;
            YY = y;
            dX = dx;
            dY = dy;
        }
        
        void move()
        {
            XX += dX;
            YY += dY;
        }
        
        void print(std::ostream& os) const
        {
            os << "( " << std::setw(9) << XX;
            os << "  " << std::setw(9) << YY << " )";
        }
    };

    
    /// a point in 3D, with a bitfield for connectivity
    struct Vertex3
    {
        /// coordinates of the point
        FLOAT XX, YY, ZZ;
        
        /// bit-field used to describe the connectivity between the points.
        /**
         Two points A and B are connected if ( A.UU & B.UU ) using the bit-wise AND.
         With an integer, this limits the number of edges to 32.
         With a long integer, this limits the number of edges to 64.
         A bigger integer could be used if needed.
        */
        unsigned UU;
        
        void set(Vector3 const& vec, unsigned u)
        {
            XX = (FLOAT)vec.XX;
            YY = (FLOAT)vec.YY;
            ZZ = (FLOAT)vec.ZZ;
            UU = u;
        }
        
        void print(std::ostream& os) const
        {
            os << "( " << std::setw(9) << XX;
            os << "  " << std::setw(9) << YY;
            os << "  " << std::setw(9) << ZZ << " )";
        }
    };

    //------------------------------------ 1D --------------------------------------
#pragma mark - 1D

    /// Rasterizer function in 1D
    void paintThickLine1D(void (*paint)(int, int, int, int, void*),
                          void* arg,             ///< last argument to paint()
                          const Vector1& P,      ///< segment end point
                          const Vector1& Q,      ///< other segment end point
                          real  radius,          ///< half-width of painted area
                          const Vector1& offset, ///< phase of the grid
                          const Vector1& delta); ///< period for the grid
    
    
    //------------------------------------ 2D --------------------------------------
#pragma mark - 2D
    
    /// compute convex hull of points in the plane
    int convexHull2D(int n_pts, Vertex2 pts[]);

    /// Paint a polygon in 2D
    /**
     paintPolygon2D() calls paintPoint(x,y,zz) for every point (x,y) of
     integral coordinates, which are inside the polygon given in xy[]. 
     The polygon should be convex, and ordered anti-clockwise.
     */
    void paintPolygon2D(void (*paint)(int, int, int, int, void*),
                        void * arg,           ///< last argument to paint
                        size_t n_pts,         ///< number of points
                        const Vertex2[],      ///< the 2D points ( x0 y0 x1 y1 ...)
                        int zz);              ///< third coordinate, passed as argument to paint()

    
    /// Paint the inside of a rectangle with edges parallel to the segment PQ
    void paintRectangle(void (*paint)(int, int, int, int, void*),
                        void* arg,            ///< last argument to paint
                        Vector2 P,            ///< segment end point [dim=3]
                        Vector2 Q,            ///< other segment end point [dim=3]
                        real  iPQ,            ///< 1 / (length of PQ)
                        real  radius);        ///< half-width of painted area
    
    
    /// Paint the inside of a rectangle with edges parallel to the segment PQ
    void paintRectangle(void (*paint)(int, int, int, int, void*),
                        void* arg,             ///< last argument to paint
                        Vector2 P,             ///< segment end point
                        Vector2 Q,             ///< other segment end point
                        real  iPQ,             ///< 1 / (length of PQ)
                        const real radius,     ///< radius of cylinder
                        const Vector2& offset, ///< phase of the grid
                        const Vector2& delta); ///< period for the grid
    
    /// Paint a 2D rectangular volume with edges parallel to the main axes
    /**
     The painted volume is square and aligned with the principal axes (X, Y, Z)
     It contains all the points at a distance 'radius' or less from the segment [P,Q].
     This is the fastest rasterizer, but the volume can be much greater than that of the cylinder.
     However, the volume is nearly optimal if PQ is aligned with one of the main axis, 
     and paintBox3D is then the best choice.
     */
    void paintBox2D(void (*paint)(int, int, int, int, void*),
                    void* arg,             ///< last argument to paint
                    const Vector2& P,      ///< segment end point
                    const Vector2& Q,      ///< other segment end point
                    real radius,           ///< radius of cylinder
                    const Vector2& offset, ///< phase of the grid
                    const Vector2& delta); ///< period for the grid
    
    
    //------------------------------------ 3D --------------------------------------
#pragma mark - 3D
         
    /// Paint a 3D polygon for which the edges of the convex hull are known
    /**
     The polygon is the convex hull of the 'nbpts' vertices given in pts[].
     Each Vertex contains coordinates and information on the connectivity to other points.
     The connections between Vertices are the edge of the 3D polygon.
     */
    void paintPolygon3D(void (*paint)(int, int, int, int, void*),
                        void * arg,        ///< last argument to paint
                        size_t n_pts,      ///< number of points
                        Vertex3 pts[]);    ///< coordinates + connectivity
    
    
    /// Paint a 3D cylinder with square section, aligned with the segment [P,Q]
    /**
     A volume is painted around the segment [P,Q], containing the cylinder of
     all the points located at a distance 'radius' or less from [P,Q].
     The volume is a right cylinder with a square section.
     */
    void paintCuboid(void (*paint)(int, int, int, int, void*),
                     void* arg,             ///< last argument to paint
                     Vector3 P,             ///< segment end point
                     Vector3 Q,             ///< other segment end point
                     real  iPQ,             ///< 1 / (length of PQ)
                     real  radius,          ///< radius of cylinder
                     const Vector3& offset, ///< phase of the grid
                     const Vector3& delta); ///< period for the grid
    
    
    /// Paint a 3D cylinder with hexagonal section, aligned with the segment [P,Q]
    /**
     A volume is painted around points [P,Q], which contains the cylinder of
     all the points at a distance 'radius' or less from the segment [P,Q].
     The volume is a right cylinder with hexagonal section.
     This is a tighter approximation of the cylinder than the square cylinder of paintCuboid.
     */
    void paintHexagonalPrism(void (*paint)(int, int, int, int, void*),
                             void* arg,             ///< last argument to paint
                             Vector3 P,             ///< segment end point
                             Vector3 Q,             ///< other segment end point
                             real  iPQ,             ///< 1 / (length of PQ)
                             real  radius,          ///< radius of cylinder
                             const Vector3& offset, ///< phase of the grid
                             const Vector3& delta); ///< period for the grid

    
    /// Paint a 3D rectangular volume with edges parallel to the main axes
    /**
     The painted volume is square and its edges are parallel to the principal axes (X, Y, Z)
     It contains all the points at a distance 'radius' or less from the segment [P,Q].
     This is the fastest rasterizer, but the volume can be much greater than that of the cylinder,
     in particular in the case where PQ >> radius, and PQ is oriented along a diagonal.
     However, the volume is nearly optimal if PQ is almost aligned with one of the main axis, 
     and paintBox3D is then a better choice than the rasterizers that paint a cylinder,
     because it is much faster.
     */
    void paintBox3D(void (*paint)(int, int, int, int, void*),
                    void * arg,            ///< last argument to paint
                    const Vector3& P,      ///< segment end point
                    const Vector3& Q,      ///< other segment end point
                    real radius,           ///< radius of cylinder
                    const Vector3& offset, ///< phase of the grid
                    const Vector3& delta); ///< period for the grid

}

#endif

