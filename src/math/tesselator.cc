// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "assert_macro.h"
#include "tesselator.h"
#include <algorithm>
#include <cmath>


// return 'a < b'
bool Tesselator::Vertex::smaller(unsigned a, unsigned b)
{
    if ( weight_[a] < weight_[b] )
        return true;
    return ( weight_[a] > 0 && weight_[a] == weight_[b]
            && index_[a] > index_[b] );
}

void Tesselator::Vertex::swap(unsigned a, unsigned b)
{
    unsigned p = index_[b];
    unsigned w = weight_[b];
    index_[b] = index_[a];
    weight_[b] = weight_[a];
    index_[a] = p;
    weight_[a] = w;
}

void Tesselator::Vertex::set(unsigned A, unsigned wA,
                             unsigned B, unsigned wB,
                             unsigned C, unsigned wC )
{
    index_[0] = A;
    weight_[0] = wA;
    assert_true( wA >= 0 );
    
    index_[1] = B;
    weight_[1] = wB;
    assert_true( wB >= 0 );
    
    index_[2] = C;
    weight_[2] = wC;
    assert_true( wC >= 0 );
    
    if ( smaller(0,1) ) swap(0,1);
    if ( smaller(1,2) ) swap(1,2);
    if ( smaller(0,1) ) swap(0,1);
}

/**
 The indices should be sorted for this to work
 */
bool Tesselator::Vertex::equivalent(unsigned A, unsigned wA, unsigned B, unsigned wB) const
{
    assert_true( weight_[0] >= weight_[1] );
    assert_true( weight_[1] >= weight_[2] );
    // order by larger weigth and smaller point index
    if ( wB > wA || ( wA == wB && A > B ))
    {
        unsigned T = A;
        A = B;
        B = T;
        unsigned wT = wA;
        wA = wB;
        wB = wT;
    }
    unsigned S = sum_weights();
    unsigned T = wA + wB;
    return (   T*weight_[0]==S*wA && ( weight_[0]==0 || index_[0]==A )
            && T*weight_[1]==S*wB && ( weight_[1]==0 || index_[1]==B ));
}


void Tesselator::Vertex::print(unsigned inx, FILE* f) const
{
    fprintf(f, "P%i = ( ", inx);
    if ( weight_[2] == 0 && weight_[1] == 0 )
    {
        fprintf(f, "%i %u", index_[0], weight_[0]);
    }
    else if ( weight_[2] == 0 )
    {
        fprintf(f, "%i %u", index_[0], weight_[0]);
        fprintf(f, "  %i %u", index_[1], weight_[1]);
    }
    else
    {
        fprintf(f, "%i %u", index_[0], weight_[0]);
        fprintf(f, "  %i %u", index_[1], weight_[1]);
        fprintf(f, "  %i %u", index_[2], weight_[2]);
    }
    fprintf(f, " )");
}


//------------------------------------------------------------------------------
#pragma mark - Solid

void Tesselator::build()
{
    num_corners_  = 0;
    corners_      = nullptr;
    
    max_vertices_ = 0;
    num_vertices_ = 0;
    vertices_     = nullptr;

    max_faces_    = 0;
    num_faces_    = 0;
    faces_        = nullptr;
    
    max_edges_    = 0;
    num_edges_    = 0;
    edges_        = nullptr;
    
    vex_ = nullptr;
    
    kind_ = 0;
    halfZ_ = 0;
    
    for ( int i = 0; i < 3; ++i )
        length_[i] = 0;
    length_[3] = 1;
}


void Tesselator::allocate(unsigned V, unsigned E, unsigned F, unsigned N)
{
    unsigned nv = V + E * (N-1) + F * ((N-1)*(N-2)/2);
    unsigned ne = E * N + F * 3 * ((N-1)*N/2);
    unsigned nf = F * N * N;
    max_vertices_ = nv;
    vertices_     = new Vertex[nv];
    max_edges_    = ne;
    max_faces_    = nf;
    faces_        = new unsigned[3*nf];
}


void Tesselator::build(Polyhedra kind, unsigned div, int make)
{
    build();
    if ( div > 0 )
    {
        switch( kind )
        {
            case UNSET: break;
            case TETRAHEDRON: initTetrahedron(div); break;
            case OCTAHEDRON: initOctahedron(div); break;
            case ICOSAHEDRON: initIcosahedron(div); break;
            case HEMISPHERE: initHemisphere(div); break;
            case DICE: initDice(0.7, 0.5, 0.5, 0.3, div, div); break;
        }
        
        assert_true( num_vertices_ <= max_vertices_ );
        assert_true( num_faces_ <= max_faces_ );
        if ( make )
            setVertices();
        if ( make & 2 )
            setEdges();
    }
}

void Tesselator::setVertices()
{
    delete[] vex_;
    vex_ = new float[3*max_vertices_];
    store_vertices(vex_);
}

Tesselator::~Tesselator()
{
    delete[] corners_;
    delete[] vertices_;
    delete[] edges_;
    delete[] faces_;
}


//------------------------------------------------------------------------------
#pragma mark -

/**
 register a new derived-vertex:
 */
unsigned Tesselator::addVertex(unsigned A, unsigned wA, unsigned B, unsigned wB, unsigned C, unsigned wC)
{
    unsigned i = num_vertices_++;
    assert_true(i < max_vertices_);
    vertices_[i].set(A, wA, B, wB, C, wC);
    //vertices_[i].print(i, std::cerr);
    return i;
}


/**
 return index of the Vertex corresponding to given interpolation
 of max_vertices_ if not found
 */
unsigned Tesselator::findEdgeVertex(unsigned A, unsigned wA, unsigned B, unsigned wB) const
{
    /**
     We limit the search to the vertices that are on the edges of
     the original Platonic solid, since they are the only one to be duplicated
     Yet, we have a linear search that may be limiting for very fine divisions
     */
    for ( unsigned i = 0; i < num_edge_vertices_; ++i )
    {
        if ( vertices_[i].equivalent(A, wA, B, wB) )
            return i;
    }
    return max_vertices_;
}

/**
 return index of the Vertex corresponding to given interpolation
 */
unsigned Tesselator::getEdgeVertex(unsigned A, unsigned wA, unsigned B, unsigned wB) const
{
    unsigned i = findEdgeVertex(A, wA, B, wB);
    if ( i >= max_vertices_ )
    {
        fprintf(stderr, "Tesselator internal error: non-existent edge vertex\n");
        return 0;
    }
    return i;
}


/**
 return the Vertex identical to `X` if it exists, or store it otherwise.
 
 We only check existence if X is on an edge, since only then may it
 has been created previously while processing another face.
 */
unsigned Tesselator::makeVertex(unsigned A, unsigned wA, unsigned B, unsigned wB, unsigned C, unsigned wC)
{
    assert_true( wA > 0 && wB > 0 && wC >= 0 );
    if ( wC == 0 )
        return getEdgeVertex(A, wA, B, wB);
    return addVertex(A, wA, B, wB, C, wC);
}


//------------------------------------------------------------------------------
#pragma mark -

/// add intermediate vertices between `a` and `b`
void Tesselator::refineEdge(unsigned a, unsigned b, unsigned div)
{
    // check if this edge has been processes already:
    unsigned i = findEdgeVertex(a, div-1, b, 1);
    if ( i >= max_vertices_ )
    {
        for ( unsigned u = 1; u < div; ++u )
        {
            addVertex(a, div-u, b, u, 0, 0);
            num_edge_vertices_ = num_vertices_;
        }
    }
}


void Tesselator::addFace(unsigned a, unsigned b, unsigned c)
{
    //printf("face %i %i %i\n", a, b, c);
    assert_true( a!=b && b!=c && a!=c );
    if ( halfZ_ )
    {
        // reject face if 2 corners are outside
        FLOAT pos[3];
        int out = 0;
        interpolate(vertex(a), pos, 0); out += ( pos[2] * halfZ_ > 0 );
        interpolate(vertex(b), pos, 0); out += ( pos[2] * halfZ_ > 0 );
        interpolate(vertex(c), pos, 0); out += ( pos[2] * halfZ_ > 0 );
        if ( out > 1 )
            return;
    }
    faces_[3*num_faces_  ] = a;
    faces_[3*num_faces_+1] = b;
    faces_[3*num_faces_+2] = c;
    ++num_faces_;
}


void Tesselator::refineFace(unsigned* line, unsigned a, unsigned b, unsigned c, unsigned div)
{
    for ( unsigned u = 0; u <= div; ++u )
        line[u] = getEdgeVertex(a, div-u, b, u);
    
    for ( unsigned ii = 1; ii <= div; ++ii )
    {
        unsigned x = getEdgeVertex(a, div-ii, c, ii);
        addFace(line[0], line[1], x);
        line[0] = x;
        
        for ( unsigned jj = 1; ii+jj <= div; ++jj )
        {
            x = makeVertex(b, jj, c, ii, a, div-ii-jj);
            addFace(line[jj-1], line[jj], x);
            addFace(line[jj], line[jj+1], x);
            line[jj] = x;
        }
    }
}

void Tesselator::refineQuad(unsigned* line, unsigned quad[4], unsigned div)
{
    unsigned A = quad[0];
    unsigned B = quad[1];
    unsigned C = quad[2];
    unsigned D = quad[3];

    for ( unsigned u = 0; u <= div; ++u )
        line[u] = getEdgeVertex(A, div-u, B, u);
    
    for ( unsigned ii = 1; ii <= div; ++ii )
    {
        unsigned x = getEdgeVertex(A, div-ii, C, ii);
        addFace(line[0], line[1], x);
        line[0] = x;
        unsigned stop = div-ii;
        for ( unsigned jj = 1; jj <= stop; ++jj )
        {
            x = makeVertex(B, jj, C, ii, A, div-ii-jj);
            addFace(line[jj-1], line[jj], x);
            addFace(line[jj], line[jj+1], x);
            line[jj] = x;
        }
        for ( unsigned jj = stop+1; jj < div; ++jj )
        {
            x = makeVertex(D, jj-stop, C, stop+ii-jj, B, div-ii);
            addFace(line[jj-1], line[jj], x);
            addFace(line[jj], line[jj+1], x);
            line[jj] = x;
        }
        x = getEdgeVertex(D, ii, B, div-ii);
        addFace(line[div-1], line[div], x);
        line[div] = x;
    }
}


/// unfinished
void Tesselator::refineStrip(unsigned cnt, unsigned inx[], unsigned div)
{
    unsigned* line = new unsigned[cnt*(div+1)];

    for ( unsigned c = 0; c < cnt; ++c )
    {
        unsigned A = inx[0];
        unsigned B = inx[1];
        unsigned C = inx[2];

        for ( unsigned u = 0; u <= div; ++u )
            line[div*c+u] = getEdgeVertex(A, div-u, C, u);
    }
    
    for ( unsigned ii = 1; ii <= div; ++ii )
    {
    }
    delete[] line;
}


//------------------------------------------------------------------------------
#pragma mark -

void Tesselator::setCorners(unsigned n_vex, FLOAT vex[][3], unsigned div)
{
    delete[] corners_;
    num_corners_ = n_vex;
    corners_ = new Corner[num_corners_];
    
    for ( unsigned c = 0; c < n_vex; ++c )
    {
        corners_[c].init(c, vex[c][0], vex[c][1], vex[c][2]);
        // also create the corresponding 'derived' vertex:
        addVertex(c, div, 0, 0, 0, 0);
    }
    num_edge_vertices_ = n_vex;
}


void Tesselator::refineTriangles(unsigned n_fac, unsigned fac[][3], unsigned div)
{
    for ( unsigned f = 0; f < n_fac; ++f )
    {
        unsigned a = fac[f][0];
        unsigned b = fac[f][1];
        unsigned c = fac[f][2];
        refineEdge(a, b, div);
        refineEdge(b, c, div);
        refineEdge(c, a, div);
    }
    assert_true(num_vertices_ <= max_vertices_);
    
    unsigned* line = new unsigned[div+1];
    for ( unsigned f = 0; f < n_fac; ++f )
        refineFace(line, fac[f][0], fac[f][1], fac[f][2], div);
    delete[] line;
}


void Tesselator::initTetrahedron(unsigned div)
{
    kind_ = TETRAHEDRON;
    constexpr FLOAT F_SQRT3 = 1.7320508075688772935274463415059f;
    FLOAT a = 1.0/3.0;
    FLOAT b = M_SQRT2/3.0;
    FLOAT c = M_SQRT2/F_SQRT3;
    
    // Four vertices on unit sphere
    FLOAT vex[4][3] = {
        { 0, 2*b, -a},
        {-c,  -b, -a},
        { c,  -b, -a},
        { 0,  0,   1},
    };
    
    // Faces are ordered for OpenGL's default rule: Counter-Clockwise = facing out
    unsigned fac[4][3] = {
        {0, 2, 1},
        {1, 3, 0},
        {0, 3, 2},
        {1, 2, 3}
    };
    
    allocate(4, 6, 4, div);
    setCorners(4, vex, div);
    refineTriangles(4, fac, div);
}


void Tesselator::initOctahedron(unsigned div)
{
    kind_ = OCTAHEDRON;
    // Eight vertices on unit sphere
    FLOAT vex[6][3] = {
        { 0,  0,  1},
        { 0,  0, -1},
        { 1,  0,  0},
        {-1,  0,  0},
        { 0, -1,  0},
        { 0,  1,  0},
    };
    
    // Faces are ordered for OpenGL's default rule: Counter-Clockwise = facing out
    unsigned fac[8][3] = {
        {2, 0, 4},
        {1, 3, 5},
        {0, 3, 4},
        {1, 5, 2},
        {3, 0, 5},
        {1, 2, 4},
        {0, 2, 5},
        {1, 4, 3}
    };
    
    allocate(6, 12, 8, div);
    setCorners(6, vex, div);
    refineTriangles(8, fac, div);
}


void Tesselator::initIcosahedron(unsigned div)
{
    kind_ = ICOSAHEDRON;
    const FLOAT G = 0.5+0.5*std::sqrt(5.0);
    const FLOAT Z = 1 / std::sqrt(G*G+1.0);
    const FLOAT T = G * Z;
    
    // Twelve vertices of icosahedron on unit sphere
    FLOAT vex[12][3] = {
        { T,  Z,  0},
        {-T, -Z,  0},
        {-T,  Z,  0},
        { T, -Z,  0},
        { Z,  0,  T},
        {-Z,  0, -T},
        { Z,  0, -T},
        {-Z,  0,  T},
        { 0,  T,  Z},
        { 0, -T, -Z},
        { 0, -T,  Z},
        { 0,  T, -Z}
    };
    
    // Faces are ordered for OpenGL's default rule: Counter-Clockwise = facing out
    unsigned fac[20][3] = {
        {0,  3,  6},
        {1,  7,  2},
        {0,  4,  3},
        {1,  2,  5},
        {0,  8,  4},
        {1,  5,  9},
        {0, 11,  8},
        {1,  9, 10},
        {0,  6, 11},
        {1, 10,  7},
        {4,  8,  7},
        {5,  6,  9},
        {2,  7,  8},
        {6,  3,  9},
        {2,  8, 11},
        {9,  3, 10},
        {5,  2, 11},
        {3,  4, 10},
        {6,  5, 11},
        {4,  7, 10},
    };
    
    allocate(12, 30, 20, div);
    setCorners(12, vex, div);
    refineTriangles(20, fac, div);
}

void Tesselator::initIcosahedronRotated(unsigned div)
{
    kind_ = ICOSAHEDRON;
    const FLOAT Z = std::sqrt(0.2);
    const FLOAT C = std::cos(5*M_PI_2);
    const FLOAT S = std::sin(5*M_PI_2);
    const FLOAT D = C*C - S*S;
    const FLOAT T = C*S*2;
    
    // Twelve vertices of icosahedron on unit sphere
    FLOAT vex[12][3] = {
        { 0,  0,  1},
        { 1,  0,  Z},
        { C,  S,  Z},
        { D,  T,  Z},
        { D, -T,  Z},
        { C, -S,  Z},
        {-D,  T, -Z},
        {-C,  S, -Z},
        {-1,  0, -Z},
        {-C, -S, -Z},
        {-D, -S, -Z},
        { 0,  0, -1}
    };
    
    // Faces are ordered for OpenGL's default rule: Counter-Clockwise = facing out
    unsigned fac[20][3] = {
        {0,  1,  2},
        {0,  2,  3},
        {0,  3,  4},
        {0,  4,  5},
        {0,  5,  1} ,
        {1,  6,  2},
        {2,  7,  3},
        {3,  8,  4},
        {4,  9,  5},
        {5, 10,  1} ,
        {6,  7,  2},
        {7,  8,  3},
        {8,  9,  4},
        {9, 10,  5},
        {10, 6,  1} ,
        {11, 7,  6},
        {11, 8,  7},
        {11, 9,  8},
        {11, 10, 9},
        {11, 6, 10}
    };
    
    allocate(12, 30, 20, div);
    setCorners(12, vex, div);
    refineTriangles(20, fac, div);
}

void Tesselator::initHemisphere(unsigned div)
{
    kind_ = HEMISPHERE;
    const FLOAT G = 0.5+0.5*std::sqrt(5.0);
    const FLOAT Z = 1 / std::sqrt(G*G+1.0);
    const FLOAT T = G * Z;
    
    // Twelve vertices of icosahedron on unit sphere
    FLOAT vex[12][3] = {
        {-Z,  T,  0},
        { Z, -T,  0},
        {-Z, -T,  0},
        { Z,  T,  0},
        { 0,  Z,  T}, // 4
        { 0, -Z, -T},
        { 0,  Z, -T},
        { 0, -Z,  T}, // 7
        {-T,  0,  Z},
        { T,  0, -Z},
        { T,  0,  Z},
        {-T,  0, -Z}
    };
    
    /* Remove any face involving vertex 4 or 7 */
    // Faces are ordered for OpenGL's default rule: Counter-Clockwise = facing out
    unsigned fac[12][3] = {
        {0,  3,  6},
      //{1,  7,  2},
      //{0,  4,  3},
        {1,  2,  5},
      //{0,  8,  4},
        {1,  5,  9},
        {0, 11,  8},
        {1,  9, 10},
        {0,  6, 11},
      //{1, 10,  7},
      //{4,  8,  7},
        {5,  6,  9},
      //{2,  7,  8},
        {6,  3,  9},
        {2,  8, 11},
        {9,  3, 10},
        {5,  2, 11},
      //{3,  4, 10},
        {6,  5, 11},
      //{4,  7, 10}
    };
    
    halfZ_ = 1;
    // we can skip 5 triangles which are entirely in Z > 0
    allocate(12, 21, 12, div);
    setCorners(12, vex, div);
    refineTriangles(12, fac, div);
}


/**
 This divides the edges by 'div' and the 6 square faces bi 'vid'
 */
void Tesselator::initDice(FLOAT X, FLOAT Y, FLOAT Z, FLOAT R, unsigned div, unsigned vid)
{
    kind_ = DICE;
    length_[0] = X;
    length_[1] = Y;
    length_[2] = Z;
    length_[3] = R;
    
    const FLOAT Xr = X + R;
    const FLOAT Yr = Y + R;
    const FLOAT Zr = Z + R;

    FLOAT vex[24][3] = {
        {+Xr, Y,-Z}, { Xr, Y, Z}, { Xr,-Y,-Z}, { Xr,-Y, Z},
        {-Xr,-Y,-Z}, {-Xr,-Y, Z}, {-Xr, Y,-Z}, {-Xr, Y, Z},
        {+X, Yr,-Z}, {-X, Yr,-Z}, { X, Yr, Z}, {-X, Yr, Z},
        {+X,-Yr, Z}, {-X,-Yr, Z}, { X,-Yr,-Z}, {-X,-Yr,-Z},
        {+X, Y, Zr}, {-X, Y, Zr}, { X,-Y, Zr}, {-X,-Y, Zr},
        {+X, Y,-Zr}, { X,-Y,-Zr}, {-X, Y,-Zr}, {-X,-Y,-Zr}
    };
    
    unsigned quad[18][4] = {
        { 0, 1, 2, 3 }, { 4, 5, 6, 7 },         // faces at +X and -X
        { 8, 9, 10, 11 }, { 12, 13, 14, 15 },   // faces at +Y and -Y
        { 16, 17, 18, 19 }, { 20, 21, 22, 23 }, // faces at +Z and -Z
        // parallel to the Z axis
        { 8, 10, 0, 1 }, { 2, 3, 14, 12 },
        { 15, 13, 4, 5 }, { 6, 7, 9, 11 },
        // parallel to the Y axis
        { 0, 2, 20, 21 }, { 22, 23, 6, 4 },
        { 7, 5, 17, 19 }, { 16, 18, 1, 3 },
        // parallel to the X axis
        { 10, 11, 16, 17 }, { 18, 19, 12, 13 },
        { 14, 15, 21, 23 }, { 20, 22, 8, 9 }
    };
    
    unsigned fac[12][3] = {
        { 0, 20, 8 }, { 1, 10, 16 },
        { 2, 14, 21 }, { 3, 18, 12 },
        { 4, 23, 15 }, { 5, 13, 19 },
        { 7, 17, 11 }, { 6, 9, 22 }
    };

    allocate(24, 90, 44, div);
    setCorners(24, vex, div);
    for ( unsigned q = 0; q < 18; ++q )
    {
        unsigned a = quad[q][0];
        unsigned b = quad[q][1];
        unsigned c = quad[q][2];
        unsigned d = quad[q][3];
        refineEdge(a, b, div);
        refineEdge(b, c, div);
        refineEdge(c, a, div);
        refineEdge(b, d, div);
        refineEdge(c, d, div);
    }
    refineTriangles(8, fac, div);
    
    unsigned* line = new unsigned[div+1];
    // 4 edges parallel to the Z axis
    for ( int n = 6; n < 10; ++n )
        refineQuad(line, quad[n], div);
    
    // 4 edges parallel to the Y axis
    for ( int n = 10; n < 14; ++n )
        refineQuad(line, quad[n], div);

    // 4 edges parallel to the X axis
    for ( int n = 14; n < 18; ++n )
        refineQuad(line, quad[n], div);

    // faces
    for ( int n = 0; n < 6; ++n )
        refineQuad(line, quad[n], vid);
    
    delete[] line;
}


//------------------------------------------------------------------------------
#pragma mark - Solid vertices

template < typename REAL >
void scale(REAL& X, REAL& Y, REAL& Z, REAL n)
{
    X *= n;
    Y *= n;
    Z *= n;
}


template < typename REAL >
void project(REAL* X)
{
    REAL n = 1.0 / std::sqrt(X[0]*X[0] + X[1]*X[1] + X[2]*X[2]);
    X[0] *= n;
    X[1] *= n;
    X[2] *= n;
}


template < typename REAL >
void projectDice(REAL* X, const REAL len[4])
{
    REAL px = std::min(len[0], std::max(-len[0], X[0]));
    REAL py = std::min(len[1], std::max(-len[1], X[1]));
    REAL pz = std::min(len[2], std::max(-len[2], X[2]));

    REAL n = 0;
    n += ( X[0] - px ) * ( X[0] - px );
    n += ( X[1] - py ) * ( X[1] - py );
    n += ( X[2] - pz ) * ( X[2] - pz );
    n = len[3] / sqrt(n);
    X[0] = px + n * ( X[0] - px );
    X[1] = py + n * ( X[1] - py );
    X[2] = pz + n * ( X[2] - pz );
}


void Tesselator::interpolate(Vertex const& vex, double ptr[3], int half) const
{
    double X = 0, Y = 0, Z = 0, S = 0;
    for ( int i = 0; i < 3; ++i )
    {
        FLOAT W = (FLOAT)vex.weight_[i];
        if ( W )
        {
            Corner const& v = corners_[vex.index_[i]];
            X += W * v.pos_[0];
            Y += W * v.pos_[1];
            Z += W * v.pos_[2];
            S += W;
        }
    }
    assert_true( S > 0 );

    if ( half * Z > 0 )
        Z = 0;

    scale(X, Y, Z, 1.0/S);
    ptr[0] = X;
    ptr[1] = Y;
    ptr[2] = Z;
}


void Tesselator::interpolate(Vertex const& vex, float ptr[3], int half) const
{
    float X = 0, Y = 0, Z = 0, S = 0;
    for ( int i = 0; i < 3; ++i )
    {
        float W = (float)vex.weight_[i];
        if ( W )
        {
            Corner const& v = corners_[vex.index_[i]];
            X += W * v.pos_[0];
            Y += W * v.pos_[1];
            Z += W * v.pos_[2];
            S += W;
        }
    }
    assert_true( S > 0 );

    if ( half * Z > 0 )
        Z = 0;
    
    scale(X, Y, Z, 1.f/S);
    ptr[0] = X;
    ptr[1] = Y;
    ptr[2] = Z;
}


void Tesselator::store_vertices(float * vec) const
{
    for ( unsigned n = 0; n < num_vertices_; ++n )
        interpolate(vertices_[n], vec+3*n, halfZ_);
    
    if ( kind_ == DICE )
    {
        for ( unsigned n = 0; n < num_vertices_; ++n )
            projectDice(vec+3*n, length_);
    }
    else
    {
        for ( unsigned n = 0; n < num_vertices_; ++n )
            project(vec+3*n);
    }
}


void Tesselator::store_vertices(double * vec) const
{
    for ( unsigned n = 0; n < num_vertices_; ++n )
        interpolate(vertices_[n], vec+3*n, halfZ_);
    
    if ( kind_ == DICE )
    {
        double len[4] = { length_[0], length_[1], length_[2], length_[3] };
        for ( unsigned n = 0; n < num_vertices_; ++n )
            projectDice(vec+3*n, len);
    }
    else
    {
        for ( unsigned n = 0; n < num_vertices_; ++n )
            project(vec+3*n);
    }
}


//------------------------------------------------------------------------------
#pragma mark - edges

void Tesselator::addEdge(unsigned a, unsigned b)
{
    edges_[2*num_edges_  ] = a;
    edges_[2*num_edges_+1] = b;
    ++num_edges_;
}


void Tesselator::setEdges()
{
    delete[] edges_;
    edges_ = new unsigned[2*max_edges_];
    num_edges_ = 0;
    
    // build edges from the faces:
    for ( unsigned i = 0; i < num_faces_; ++i )
    {
        unsigned a = faces_[3*i  ];
        unsigned b = faces_[3*i+1];
        unsigned c = faces_[3*i+2];
        if ( a < b ) addEdge(a, b);
        if ( b < c ) addEdge(b, c);
        if ( c < a ) addEdge(c, a);
    }
    
    assert_true( num_edges_ <= max_edges_ );
}

/**
 http://paulbourke.net/dataformats/ply
 */
void Tesselator::exportPLY(FILE* f) const
{
    fprintf(f, "ply\n");
    fprintf(f, "format ascii 1.0\n");
    fprintf(f, "comment made by Cytosim\n");
    fprintf(f, "element vertex %u\n", num_vertices_);
    fprintf(f, "property float x\n");
    fprintf(f, "property float y\n");
    fprintf(f, "property float z\n");
    fprintf(f, "element face %u\n", num_faces_);
    fprintf(f, "property list uchar int vertex_index\n");
//    fprintf(f, "property uint vertex0\n");
//    fprintf(f, "property uint vertex1\n");
//    fprintf(f, "property uint vertex2\n");
    fprintf(f, "end_header\n");
    
    for ( unsigned u = 0; u < num_vertices_; ++u )
        fprintf(f, "%5.3f %5.3f %5.3f\n", vex_[3*u], vex_[3*u+1], vex_[3*u+2]);
    
    for ( unsigned u = 0; u < num_faces_; ++u )
        fprintf(f, "3 %u %u %u\n", faces_[3*u], faces_[3*u+1], faces_[3*u+2]);
}


/**
 https://en.wikipedia.org/wiki/STL_(file_format)
 */
void Tesselator::exportSTL(FILE* f) const
{
    float n[3] = { 0 };
    char buf[128] = { 0 };
    fwrite(buf, 80, 1, f); // 80 bytes header
    fwrite(&num_faces_, 1, 4, f); // 4 bytes
   
    for ( unsigned u = 0; u < num_faces_; ++u )
    {
        float const* a = face_vertex0(u);
        float const* b = face_vertex1(u);
        float const* c = face_vertex2(u);
        // calculate normal of triangle:
        n[0] = ( b[1] - a[1] ) * ( b[2] - a[2] ) - ( c[2] - a[2] ) * ( b[1] - a[1] );
        n[1] = ( b[2] - a[2] ) * ( b[0] - a[0] ) - ( c[0] - a[0] ) * ( b[2] - a[2] );
        n[2] = ( b[0] - a[0] ) * ( b[1] - a[1] ) - ( c[1] - a[1] ) * ( b[0] - a[0] );

        fwrite(n, 3, 4, f); // normal
        fwrite(a, 3, 4, f);
        fwrite(b, 3, 4, f);
        fwrite(c, 3, 4, f);
        fwrite(buf, 1, 2, f);
    }
}
