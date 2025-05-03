// Cytosim was created by Francois Nedelec. Copyright Cambridge University 2020
#include "assert_macro.h"
#include "rasterizer.h"
#include "vector2.h"
#include "vector3.h"
#include <cmath>

/**
 DISPLAY enables some code useful for visual debugging,
 and should be defined for compiling test_rasterizer.cc
 */
#ifdef DISPLAY
#  include "gym_flute.h"
#  include "gym_view.h"
#  include "gym_draw.h"
bool rasterizer_draws = false;
#endif

//==============================================================================
//                                TEMPLATED
//==============================================================================

/// accessory data swap function
template<typename VEC>
void swap(VEC& A, VEC& B) { VEC T = A; A = B; B = T; }

/// Compute the Convex Hull of a set of points
/**
 on entry, pts[] contains 'n_pts' points with coordinates members XX and YY.
 on exit, pts[] is the anti-clockwise polygon hull, starting with the point of lowest Y.
 
 @returns The number of points in the convex hull (at most n_pts).
 */
template<typename VEC>
static int convexHull(int n_pts, VEC pts[])
{
    //---------- find bottom and top points:
    int inx = 0, top = 0;
    Rasterizer::FLOAT y_bot = pts[0].YY;
    Rasterizer::FLOAT y_top = pts[0].YY;
    
    for ( int n = 1; n < n_pts; ++n )
    {
        if ( pts[n].YY < y_bot || ( pts[n].YY == y_bot  &&  pts[n].XX > pts[inx].XX ) )
        {
            inx = n;
            y_bot = pts[n].YY;
        }
        if ( pts[n].YY > y_top || ( pts[n].YY == y_top  &&  pts[n].XX < pts[top].XX ) )
        {
            top = n;
            y_top = pts[n].YY;
        }
    }
    
    if ( inx == top )  //all points are equal ?
        return 1;
    
    // put the bottom point at index zero:
    if ( inx )
    {
        swap(pts[0], pts[inx]);
        if ( top == 0 )
            top = inx;
    }
    // reset
    inx = 0;
    
    // wrap upward on the right side of the hull
    int nxt;
    while ( 1 )
    {
        Rasterizer::FLOAT pX = pts[inx].XX;
        Rasterizer::FLOAT pY = pts[inx].YY;
        ++inx;
        
        nxt = top;
        Rasterizer::FLOAT dx = pts[top].XX - pX;
        Rasterizer::FLOAT dy = pts[top].YY - pY;
        
        for ( int n = inx; n < n_pts; ++n )
        {
            Rasterizer::FLOAT dxt = pts[n].XX - pX;
            Rasterizer::FLOAT dyt = pts[n].YY - pY;
            // keep if slope is lower:
            if ( dxt * dy > dyt * dx )
            {
                nxt = n;
                dx = dxt;
                dy = dyt;
            }
        }
        
        // if we reached bottom point, sweep is complete
        if ( nxt == 0 )
            break;
        
        swap(pts[inx], pts[nxt]);
        
        // if we reached topmost point, change reference for downward sweep:
        if ( nxt == top )
            top = 0;
        else if ( inx == top )
            top = nxt;  // this compensates the swap
    }
    
    return inx;
}


int Rasterizer::convexHull2D(int n_pts, Rasterizer::Vertex2 pts[])
{
    return convexHull(n_pts, pts);
}


/**
 pts[] stores a counter-clockwise polygon, with pts[0] the point of lowest Y.
 */
template < typename VEC >
static void paintPolygon(void (*paint)(int, int, int, int, void*), void * arg,
                              const size_t n_pts, const VEC pts[],
                              const int zz)
{
#ifdef DISPLAY
    if ( rasterizer_draws )
    {
        gym::ref_view();
        gym::color(0, 0, 1);
        flute3 * flu = gym::mapBufferV3(n_pts+1);
        for ( size_t i = 0; i < n_pts; ++i )
            flu[i] = { pts[i].XX, pts[i].YY, float(zz) };
        flu[n_pts] = { pts[0].XX, pts[0].YY, float(zz) };
        gym::unmapBufferV3();
        gym::drawLineStrip(1, 0, n_pts+1);
        gym::drawPoints(7, 0, 1);
    }
#endif

#if ( 0 )
    // print polygon:
    for ( size_t n = 0; n < n_pts; ++n )
        pts[n].print(std::clog);
    std::clog << " " << zz << '\n';
#endif
  
    // start with bottom point:
    Rasterizer::FLOAT xxR = pts[0].XX, yyR = pts[0].YY;
    Rasterizer::FLOAT xxL = pts[0].XX, yyL = pts[0].YY;

    // start on the Y-line just above the bottom point
    int yy = (int)std::ceil(yyR);
        
    size_t iR = 1;
    size_t iL = n_pts-1;
    
    // finish Y-lines for right and left sides:
    int yR = (int)std::floor(pts[iR].YY);
    int yL = (int)std::floor(pts[iL].YY);
    
    // left and right sided slopes dX/dY:
    Rasterizer::FLOAT dxR = ( pts[iR].XX - xxR ) / ( pts[iR].YY - yyR );
    Rasterizer::FLOAT dxL = ( pts[iL].XX - xxL ) / ( pts[iL].YY - yyL );

    xxR += dxR * ( yy - yyR );
    xxL += dxL * ( yy - yyL );
    
    //printf("\nY%i ", yy);
    while ( iR <= iL )
    {
        // index of the last line without changing edges:
        int top = std::min(yL, yR);
        
        for ( ; yy <= top; ++yy )
        {
            //printf("[%.2f %.2f]Y%i ", xxL, xxR, yy);
            int inf = (int) std::ceil(xxL);
            int sup = (int)std::floor(xxR);
            paint(inf, sup, yy, zz, arg);
            xxL += dxL;
            xxR += dxR;
        }
        
        //printf("y%i yL%i yR%i ", yy, yL, yR);
        // update point on right edge:
        if ( yR < yy )
        {
            xxR = pts[iR].XX;
            yyR = pts[iR].YY;
            do {
                if ( ++iR > iL ) return;
                yR = (int)std::floor(pts[iR].YY);
            } while ( yR < yy );
            //printf("R%iL%i ", yR, yL);
            dxR = ( pts[iR].XX - xxR ) / ( pts[iR].YY - yyR );
            xxR += dxR * ( yy - yyR );
        }
        
        // update point on left edge:
        if ( yL < yy )
        {
            xxL = pts[iL].XX;
            yyL = pts[iL].YY;
            do {
                if ( --iL < iR ) return;
                yL = (int)std::floor(pts[iL].YY);
            } while ( yL < yy );
            //printf("L%iR%i ", yL, yR);
            dxL = ( pts[iL].XX - xxL ) / ( pts[iL].YY - yyL );
            xxL += dxL * ( yy - yyL );
        }
    }
}

void Rasterizer::paintPolygon2D(void (*paint)(int, int, int, int, void*), void * arg,
                                const size_t n_pts, const Vertex2 pts[],
                                const int zz)
{
    paintPolygon(paint, arg, n_pts, pts, zz);
}

//==============================================================================
//                                   1D
//==============================================================================
#pragma mark - 1D

void Rasterizer::paintThickLine1D(void (*paint)(int, int, int, int, void*), void * arg,
                                  const Vector1& P, const Vector1& Q,
                                  const real radius, const Vector1& offset, const Vector1& delta)
{
    real L = std::min(P.XX, Q.XX);
    real R = std::max(P.XX, Q.XX);

    int inf = (int) ceil( ( L - radius - offset.XX ) * delta.XX );
    int sup = (int)floor( ( R + radius - offset.XX ) * delta.XX );

    paint(inf, sup, 0, 0, arg);
}


//==============================================================================
//                                     2D
//==============================================================================
#pragma mark - 2D

void Rasterizer::paintRectangle(void (*paint)(int, int, int, int, void*), void * arg,
                                Vector2 P, Vector2 Q, const real iPQ,
                                const real radius)
{
    Vector2 A, B;
    // swap to place P lower than Q
    if ( P.YY < Q.YY )
    {
        B = ( radius * iPQ ) * ( Q - P );
    }
    else
    {
        B = ( radius * iPQ ) * ( P - Q );
        A = P;
        P = Q;
        Q = A;
    }
    
    A = B + Vector2(B.YY, -B.XX);
    B = B - Vector2(B.YY, -B.XX);

    Vertex2 pts[5];

    // compose the 4 corners of the rectangle:
    pts[0] = ( P - A );
    pts[1] = ( P - B );
    pts[2] = ( Q + A );
    pts[3] = ( Q + B );
    pts[4] = pts[0];
    
    // always start from index of the lowest corner
    paintPolygon(paint, arg, 4, pts+(P.XX<Q.XX), 0);
}


void Rasterizer::paintRectangle(void (*paint)(int, int, int, int, void*), void * arg,
                                Vector2 P, Vector2 Q, const real iPQ,
                                const real radius, const Vector2& offset, const Vector2& delta)
{
    Vector2 A, B;
    // swap to place P lower than Q
    if ( P.YY < Q.YY )
    {
        B = ( radius * iPQ ) * ( Q - P );
        P = P - offset;
        Q = Q - offset;
    }
    else
    {
        A = P;
        B = ( radius * iPQ ) * ( P - Q );
        P = Q - offset;
        Q = A - offset;
    }
    
    A = B + Vector2(B.YY, -B.XX);
    B = B - Vector2(B.YY, -B.XX);

    Vertex2 pts[5];

    // compose the 4 corners of the rectangle:
    pts[0] = ( P - A ).e_mul(delta);
    pts[1] = ( P - B ).e_mul(delta);
    pts[2] = ( Q + A ).e_mul(delta);
    pts[3] = ( Q + B ).e_mul(delta);
    pts[4] = pts[0];
    
    // always start from index of the lowest corner
    paintPolygon(paint, arg, 4, pts+(P.XX<Q.XX), 0);
}


void Rasterizer::paintBox2D(void (*paint)(int, int, int, int, void*), void * arg,
                            const Vector2& P, const Vector2& Q, const real radius,
                            const Vector2& offset, const Vector2& delta )
{
    int inf[2], sup[2];
    
    for ( int d = 0; d < 2; ++d )
    {
        real i = std::min(P[d], Q[d]);
        real s = std::max(P[d], Q[d]);
        inf[d] = (int) std::ceil( ( i - radius - offset[d] ) * delta[d] );
        sup[d] = (int)std::floor( ( s + radius - offset[d] ) * delta[d] );
    }
    
    for ( int yy = inf[1]; yy <= sup[1]; ++yy )
        paint(inf[0], sup[0], yy, 0, arg);
}


//==============================================================================
//                                     3D
//==============================================================================
#pragma mark - 3D


/// qsort function comparing the Z component of two points
static int compareVertex3(const void * a, const void * b)
{
    Rasterizer::FLOAT az = ((Rasterizer::Vertex3 const*)(a))->ZZ;
    Rasterizer::FLOAT bz = ((Rasterizer::Vertex3 const*)(b))->ZZ;
    
    return ( az > bz ) - ( bz > az );
}


void Rasterizer::paintPolygon3D(void (*paint)(int, int, int, int, void*), void * arg,
                                const size_t n_pts, Vertex3 pts[])
{
    assert_true( n_pts > 1 );
    
#ifdef DISPLAY
    if ( rasterizer_draws )
    {
        gym::ref_view();
        gym::color(0, 1, 1);
        const size_t sup = n_pts * (n_pts-1);
        flute3* flu = gym::mapBufferV3(sup);
        flute3* ptr = flu;
        for ( size_t n = 0; n < n_pts; ++n )
        {
            for ( size_t u = n+1; u < n_pts; ++u )
            {
                if ( pts[n].UU  &  pts[u].UU )
                {
                    ptr[0] = { pts[n].XX, pts[n].YY, pts[n].ZZ };
                    ptr[1] = { pts[u].XX, pts[u].YY, pts[u].ZZ };
                    ptr += 2;
                }
            }
        }
        assert_true( ptr <= flu+sup );
        gym::unmapBufferV3();
        gym::drawLines(0.5, 0, ptr-flu);
    }
#endif
    
    //order the points by increasing Z values:
    qsort(pts, n_pts, sizeof(Vertex3), &compareVertex3);
    
    //we can normally only cross four sides of a parallelogram in 3D
    //but in some degenerate cases, it can be more
    constexpr unsigned TOP = 8;
    Vertex2dZ xydz[TOP];
    
    int above = 0;
    int zzz = (int) std::ceil( pts[0].ZZ );
    
    while ( ++above < n_pts )
    {
        //printf("restart at zz %4i\n", zz );
        
        //find the first point strictly above the plane Z = zz:
        //the index of this point is (above-1)
        while ( pts[above].ZZ <= zzz )
        {
            if ( ++above >= n_pts )
                return;
        }
        
        //the next time we have to recalculate the lines
        //is when pts[above] will be below the plane Z = ZZtop:
        int ZZtop = (int) std::ceil( pts[above].ZZ );
        
        //number of edges crossing the plane at Z=zz;
        int edg = 0;
        //set-up all the edges of the solid 3D polygon, where one vertex
        //is below the plane, and the other point above:
        for ( int ii = 0;     ii < above; ++ii )
        for ( int jj = above; jj < n_pts; ++jj )
        {
            //test if [ii, jj] constitute an edge of the solid 3D polygon:
            if ( pts[ii].UU  &  pts[jj].UU )
            {
                double dZ = pts[jj].ZZ - pts[ii].ZZ;
                
                if ( dZ > FLT_EPSILON )
                {
                    double dXdZ = ( pts[jj].XX - pts[ii].XX ) / dZ;
                    double dYdZ = ( pts[jj].YY - pts[ii].YY ) / dZ;
                    double offZ = zzz - pts[ii].ZZ;
                    xydz[edg].fix(pts[ii].XX + dXdZ * offZ, dXdZ,
                                  pts[ii].YY + dYdZ * offZ, dYdZ);
                    ++edg;
                    assert_true( edg < TOP );
                }
            }
        }
        //printf("zz %3i : nbp = %i\n", zz, nbl);

        // the edges of the convex solid polygon should not intersect,
        // so we can take the convex hull only once in most cases:
        int nbp = convexHull(edg, xydz);
        
        // nbp = number of points in the hull. If this is not equal to the
        // the number of edges in the intersection, we have to redo the Hull:
        for ( ; ( nbp != edg ) && ( zzz < ZZtop ); ++zzz )
        {
            paintPolygon(paint, arg, nbp, xydz, zzz);
            
            //update X, Y coordinates according to the slopes, to reach Z+1:
            for ( int i = 0; i < edg; ++i )
                xydz[i].move();

            //rebuild the convex hull around the points xydz[]:
            nbp = convexHull(edg, xydz);
            //printf("zz %3i : hull nbp = %i\n", zzz, nbp);
        }

        // from here on, we do not need to remake the convex hull:
        for ( ; zzz < ZZtop; ++zzz )
        {
            paintPolygon(paint, arg, nbp, xydz, zzz);
            
            //update X, Y coordinates according to the slopes, to reach Z+1:
            for ( int i = 0; i < edg; ++i )
                xydz[i].move();
        }
    }
}


void Rasterizer::paintCuboid(void (*paint)(int, int, int, int, void*), void * arg,
                             Vector3 P, Vector3 Q, const real iPQ,
                             const real radius, const Vector3& offset, const Vector3& delta )
{
    Vector3 A, B;
    Vector3 PQ = ( Q - P ) * iPQ;
    
    // make an orthogonal basis with norm = radius * sqrt(2):
#if ( 1 )
    //std::clog << std::scientific << PQ.normSqr() << '\n';
    PQ.orthonormal(A, B, radius*M_SQRT2);
#else
    A = PQ.orthogonal(radius*M_SQRT2);
    B = cross(PQ, A);
#endif
    
    PQ *= radius;
    A = A.e_mul(delta);
    B = B.e_mul(delta);
    
    /*
     Extend segment PQ by `radius` on each side,
     and convert to the grid's coordinates:
     grid coordinates = ( coordinates - min ) * delta
     */
    P = ( P - PQ - offset ).e_mul(delta);
    Q = ( Q + PQ - offset ).e_mul(delta);
    
    // set the vertex of a generalized cylinder aligned along PQ
    Vertex3 pts[8];
    
    /*
     Below, the last arguments is a bitfield defining which point
     are connected to form the edges of the polygonal volume.
     */
    pts[0].set(P + A, 0b000000010011);
    pts[1].set(P + B, 0b000000100110);
    pts[2].set(P - A, 0b000001001100);
    pts[3].set(P - B, 0b000010001001);
    pts[4].set(Q + A, 0b001100010000);
    pts[5].set(Q + B, 0b011000100000);
    pts[6].set(Q - A, 0b110001000000);
    pts[7].set(Q - B, 0b100110000000);
    
    //paint the volume:
    paintPolygon3D(paint, arg, 8, pts);
}


/**
 Paint a cylinder of Hexagonal base.
 The hexagon covering the unit disc has vertices:
     A (  0, -2*a )
     B (  b,   -a )
     C (  b,    a )
     D (  0,  2*a ) = -A
     E ( -b,    a ) = -B
     F ( -b,   -a ) = -C
 with b = sqrt(3) * a:
     b = 1
     a = 1 / sqrt(3)
 */
void Rasterizer::paintHexagonalPrism(void (*paint)(int, int, int, int, void*), void * arg,
                                     Vector3 P, Vector3 Q, const real iPQ,
                                     const real radius, const Vector3& offset, const Vector3& delta)
{
    Vector3 A, B, C;
    Vector3 PQ = ( Q - P ) * iPQ;
    PQ.orthonormal(A, C, radius);
    PQ *= radius;

    // build the vertices of the Hexagon
    A *= 2.0 / M_SQRT3;
    B = C + A * 0.5;
    C = B - A;
    
    A = A.e_mul(delta);
    B = B.e_mul(delta);
    C = C.e_mul(delta);
    
    /*
     Extend segment PQ by `radius` on each side,
     and convert to the grid's coordinates:
     grid coordinates = ( coordinates - min ) * delta
    */
    P = ( P - PQ - offset ).e_mul(delta);
    Q = ( Q + PQ - offset ).e_mul(delta);
    
    // set the vertex of a generalized cylinder aligned along PQ
    Vertex3 pts[12];
    
    /*
     Below, the last arguments is a bitfield defining which point
     are connected to form the edges of the polygonal volume.
    */
    pts[ 0].set(P + A, 0b000000000001000011);
    pts[ 1].set(P + B, 0b000000000010000110);
    pts[ 2].set(P + C, 0b000000000100001100);
    pts[ 3].set(P - A, 0b000000001000011000);
    pts[ 4].set(P - B, 0b000000010000110000);
    pts[ 5].set(P - C, 0b000000100000100001);
    pts[ 6].set(Q + A, 0b000011000001000000);
    pts[ 7].set(Q + B, 0b000110000010000000);
    pts[ 8].set(Q + C, 0b001100000100000000);
    pts[ 9].set(Q - A, 0b011000001000000000);
    pts[10].set(Q - B, 0b110000010000000000);
    pts[11].set(Q - C, 0b100001100000000000);
    
    paintPolygon3D(paint, arg, 12, pts);
}


void Rasterizer::paintBox3D(void (*paint)(int, int, int, int, void*), void * arg,
                            const Vector3& P, const Vector3& Q, const real radius,
                            const Vector3& offset, const Vector3& delta )
{
    int inf[3], sup[3];
    
    for ( int d = 0; d < 3; ++d )
    {
        real i = std::min(P[d], Q[d]);
        real s = std::max(P[d], Q[d]);
        inf[d] = (int) std::ceil( ( i - radius - offset[d] ) * delta[d] );
        sup[d] = (int)std::floor( ( s + radius - offset[d] ) * delta[d] );
    }
    
    for ( int zz = inf[2]; zz <= sup[2]; ++zz )
    for ( int yy = inf[1]; yy <= sup[1]; ++yy )
        paint(inf[0], sup[0], yy, zz, arg);
}

