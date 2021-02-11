// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "assert_macro.h"
#include "platonic.h"
#include <cmath>

//------------------------------------------------------------------------------
#pragma mark - Solid static members

unsigned Platonic::Solid::nb_vertices(Polyhedra k, unsigned N)
{
    static const unsigned V[] = { 4, 6, 12, 12 };
    static const unsigned F[] = { 4, 8, 20, 12 };
    static const unsigned E[] = { 6, 12, 30, 21 };
    return V[k] + E[k] * (N-1) + F[k] * ((N-1)*(N-2)/2);
}

unsigned Platonic::Solid::nb_faces(Polyhedra k, unsigned N)
{
    static const unsigned F[] = { 4, 8, 20, 12 };
    return F[k] * N * N;
}

unsigned Platonic::Solid::nb_edges(Polyhedra k, unsigned N)
{
    static const unsigned F[] = { 4, 8, 20, 12 };
    static const unsigned E[] = { 6, 12, 30, 21 };
    return E[k] * N + F[k] * 3 * ((N-1)*N/2);
}

/**
 We estimate L from:
 - the area of the sphere is 4*PI
 - the number of triangles on it is nb_faces()
 - the area of an equilateral triangle of side L is std::sqrt(3)*(L/2)^2
 */
Platonic::FLOAT Platonic::Solid::length_edge(Polyhedra k, unsigned N)
{
    return 4 * std::sqrt( M_PI / ( std::sqrt(3) * nb_faces(k,N) ) );
}

//------------------------------------------------------------------------------
#pragma mark - Vertex

// return 'a < b'
bool Platonic::Vertex::smaller(unsigned a, unsigned b)
{
    if ( weight_[a] < weight_[b] )
        return true;
    return ( weight_[a] > 0 && weight_[a] == weight_[b]
            && vertex_[a]->inx_ > vertex_[b]->inx_ );
}

void Platonic::Vertex::swap(unsigned a, unsigned b)
{
    Corner  * p = vertex_[b];
    unsigned w = weight_[b];
    vertex_[b] = vertex_[a];
    weight_[b] = weight_[a];
    vertex_[a] = p;
    weight_[a] = w;
}

void Platonic::Vertex::set(Corner* A, unsigned wA,
                           Corner* B, unsigned wB,
                           Corner* C, unsigned wC )
{
    vertex_[0] = A;
    weight_[0] = wA;
    assert_true( wA == 0 || A );
    
    vertex_[1] = B;
    weight_[1] = wB;
    assert_true( wB == 0 || B );
    
    vertex_[2] = C;
    weight_[2] = wC;
    assert_true( wC == 0 || C );
    
    if ( smaller(0,1) ) swap(0,1);
    if ( smaller(1,2) ) swap(1,2);
    if ( smaller(0,1) ) swap(0,1);
}

/**
 The indices should be sorted for this to work
 */
bool Platonic::Vertex::equivalent(Corner* A, unsigned wA, Corner* B, unsigned wB) const
{
    assert_true( weight_[0] >= weight_[1] );
    assert_true( weight_[1] >= weight_[2] );
    // order by larger weigth and smaller point index
    if ( wB > wA || ( wA == wB && A->inx_ > B->inx_ ))
    {
        Corner* T = A;
        A = B;
        B = T;
        unsigned wT = wA;
        wA = wB;
        wB = wT;
    }
    unsigned S = sum_weights();
    unsigned T = wA + wB;
    return (   T*weight_[0]==S*wA && ( weight_[0]==0 || vertex_[0]==A )
            && T*weight_[1]==S*wB && ( weight_[1]==0 || vertex_[1]==B ));
}


void Platonic::Vertex::store_pos(double C[3]) const
{
    double X = 0, Y = 0, Z = 0;
    assert_true( sum_weights() > 0 );
    
    for ( int i = 0; i < 3; ++i )
    {
        Corner const* v = vertex_[i];
        if ( v )
        {
            FLOAT W = (FLOAT)weight_[i];
            X += W * v->pos_[0];
            Y += W * v->pos_[1];
            Z += W * v->pos_[2];
        }
    }
    
    // normalize:
    double n = X*X + Y*Y + Z*Z;
    if ( n > 0 )
        n = 1.0 / std::sqrt(n);
    else
        n = 1.0 / sum_weights();
    C[0] = X * n;
    C[1] = Y * n;
    C[2] = Z * n;
}

void Platonic::Vertex::store_pos(float ptr[3], int half) const
{
    float X = 0, Y = 0, Z = 0;
    assert_true( sum_weights() > 0 );
    
    for ( int i = 0; i < 3; ++i )
    {
        Corner const* v = vertex_[i];
        if ( v )
        {
            float W = (float)weight_[i];
            X += W * v->pos_[0];
            Y += W * v->pos_[1];
            Z += W * v->pos_[2];
        }
    }
    
    if ( half * Z > 0 )
        Z = 0;
    
    //normalize:
    float n = X*X + Y*Y + Z*Z;
    if ( n > 0 )
        n = 1.f / std::sqrt(n);
    else
        n = 1.f / sum_weights();
    ptr[0] = X * n;
    ptr[1] = Y * n;
    ptr[2] = Z * n;
}

void Platonic::Vertex::print(unsigned inx, std::ostream& out) const
{
    out << "P" << inx << " = ( ";
    if ( weight_[2] == 0 && weight_[1] == 0 )
    {
        out << vertex_[0]->inx_ << " " << weight_[0] ;
    }
    else if ( weight_[2] == 0 )
    {
        out << vertex_[0]->inx_ << " " << weight_[0] << ", ";
        out << vertex_[1]->inx_ << " " << weight_[1];
    }
    else
    {
        out << vertex_[0]->inx_ << " " << weight_[0] << ", ";
        out << vertex_[1]->inx_ << " " << weight_[1] << ", ";
        out << vertex_[2]->inx_ << " " << weight_[2];
    }
    out << " )" << std::endl;
}

//------------------------------------------------------------------------------
#pragma mark - Solid

void Platonic::Solid::build()
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
    
    coordinates_  = nullptr;
}


Platonic::Solid::Solid(Polyhedra K, unsigned div, bool make_edges)
{
    build();
    
    if ( div > 0 )
    {
        max_vertices_ = nb_vertices(K, div);
        vertices_     = new Vertex[max_vertices_];
        max_faces_    = nb_faces(K, div);
        faces_        = new unsigned[3*max_faces_];
        max_edges_    = nb_edges(K, div);
        setVertices(K, div);
        if ( make_edges )
            setEdges();
    }
}


Platonic::Solid::~Solid()
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
unsigned Platonic::Solid::addVertex(Corner* A, unsigned wA, Corner* B, unsigned wB, Corner* C, unsigned wC)
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
unsigned Platonic::Solid::findEdgeVertex(Corner* A, unsigned wA, Corner* B, unsigned wB) const
{
    /**
     We limit the search to the vertices that are on the edges of
     the original PlatonicSolid, since they are the only one to be duplicated
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
unsigned Platonic::Solid::getEdgeVertex(Corner* A, unsigned wA, Corner* B, unsigned wB) const
{
    unsigned i = findEdgeVertex(A, wA, B, wB);
    if ( i >= max_vertices_ )
    {
        fprintf(stderr, "Platonic::Solid internal error: non-existent edge vertex\n");
        return 0;
    }
    return i;
}


/**
 return the Vertex identical to `X` if it exists, or store it otherwise.
 
 We only check existence if X is on an edge, since only then may it
 has been created previously while processing another face.
 */
unsigned Platonic::Solid::makeVertex(Corner* A, unsigned wA, Corner* B, unsigned wB, Corner* C, unsigned wC)
{
    if ( wC == 0 )
        return getEdgeVertex(A, wA, B, wB);
    return addVertex(A, wA, B, wB, C, wC);
}


//------------------------------------------------------------------------------
#pragma mark -

/// add intermediate vertices between `a` and `b`
void Platonic::Solid::refineEdge(unsigned a, unsigned b, unsigned div)
{
    Corner *A = &corners_[a];
    Corner *B = &corners_[b];
    
    // check if this edge has been processes already:
    unsigned i = findEdgeVertex(A, div-1, B, 1);
    if ( i >= max_vertices_ )
    {
        for ( unsigned u = 1; u < div; ++u )
        {
            addVertex(A, div-u, B, u, nullptr, 0);
            num_edge_vertices_ = num_vertices_;
        }
    }
}


void Platonic::Solid::addFace(unsigned a, unsigned b, unsigned c, int half)
{
    //printf("face %i %i %i\n", a, b, c);
    assert_true( a!=b && b!=c && a!=c );
    if ( half )
    {
        // reject face if 2 corners are outside
        FLOAT pos[3];
        int out = 0;
        vertex(a).store_pos(pos, 0); out += ( pos[2] * half > 0 );
        vertex(b).store_pos(pos, 0); out += ( pos[2] * half > 0 );
        vertex(c).store_pos(pos, 0); out += ( pos[2] * half > 0 );
        if ( out > 1 )
            return;
    }
    faces_[3*num_faces_  ] = a;
    faces_[3*num_faces_+1] = b;
    faces_[3*num_faces_+2] = c;
    ++num_faces_;
}


void Platonic::Solid::refineFace(unsigned a, unsigned b, unsigned c, unsigned div, int half)
{
    Corner *A = &corners_[a];
    Corner *B = &corners_[b];
    Corner *C = &corners_[c];
    
    unsigned* line = new unsigned[div+1];
    
    for ( unsigned u = 0; u <= div;  ++u )
        line[u] = getEdgeVertex(A, div-u, B, u);
    
    for ( unsigned ii = 1; ii <= div; ++ii )
    {
        unsigned X = getEdgeVertex(A, div-ii, C, ii);
        
        addFace(line[0], line[1], X, half);
        line[0] = X;
        
        for ( unsigned jj = 1; ii+jj <= div; ++jj )
        {
            X = makeVertex(B, jj, C, ii, A, div-ii-jj);
            addFace(line[jj-1], line[jj], X, half);
            addFace(line[jj], line[jj+1], X, half);
            line[jj] = X;
        }
    }
    delete[] line;
}


//------------------------------------------------------------------------------
#pragma mark -

void Platonic::Solid::refineTriangles(unsigned n_vex, FLOAT vex[][3],
                                      unsigned n_fac, unsigned fac[][3],
                                      unsigned div, int half = 0)
{
    delete[] corners_;
    num_corners_ = n_vex;
    corners_ = new Corner[num_corners_];
    
    for ( unsigned c = 0; c < n_vex; ++c )
    {
        corners_[c].init(c, vex[c][0], vex[c][1], vex[c][2]);
        // also create the vertex 'derived' from this point
        addVertex(corners_+c, div, nullptr, 0, nullptr, 0);
    }
    num_edge_vertices_ = n_vex;

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
    
    for ( unsigned f = 0; f < n_fac; ++f )
        refineFace(fac[f][0], fac[f][1], fac[f][2], div, half);
}


void Platonic::Solid::initTetrahedron(unsigned div)
{
    constexpr FLOAT M_SQRT3 = 1.7320508075688772935274463415059;
    FLOAT a = 1.0/3.0;
    FLOAT b = M_SQRT2/3.0;
    FLOAT c = M_SQRT2/M_SQRT3;
    
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
    
    refineTriangles(4, vex, 4, fac, div);
}


void Platonic::Solid::initOctahedron(unsigned div)
{
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
    
    refineTriangles(6, vex, 8, fac, div);
}


void Platonic::Solid::initIcosahedron(unsigned div)
{
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
    
    refineTriangles(12, vex, 20, fac, div);
}

void Platonic::Solid::initIcosahedronRotated(unsigned div)
{
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
        {11, 6, 10},
    };
    
    refineTriangles(12, vex, 20, fac, div);
}

void Platonic::Solid::initHemisphere(unsigned div)
{
    const FLOAT G = 0.5+0.5*std::sqrt(5.0);
    const FLOAT Z = 1 / std::sqrt(G*G+1.0);
    const FLOAT T = G * Z;
    
    // Twelve vertices of icosahedron on unit sphere
    FLOAT vex[12][3] = {
        { T,  Z,  0},
        {-T, -Z,  0},
        {-T,  Z,  0},
        { T, -Z,  0},
        { Z,  0,  T}, // 4
        {-Z,  0, -T},
        { Z,  0, -T},
        {-Z,  0,  T}, // 7
        { 0,  T,  Z},
        { 0, -T, -Z},
        { 0, -T,  Z},
        { 0,  T, -Z}
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
      //{4,  7, 10},
    };
    
    // we can skip 5 triangles which are entirely in Z > 0
    refineTriangles(12, vex, 12, fac, div, 1);
}
//------------------------------------------------------------------------------
#pragma mark - vertices

void Platonic::Solid::setVertices(Polyhedra kind, unsigned div)
{
    assert_true(div > 0);
    
    switch( kind )
    {
        case TETRAHEDRON: initTetrahedron(div); break;
        case OCTAHEDRON: initOctahedron(div); break;
        case ICOSAHEDRON: initIcosahedron(div); break;
        case HEMISPHERE: initHemisphere(div); break;
    }
    
    assert_true( num_vertices_ <= max_vertices_ );
    assert_true( num_faces_ <= max_faces_ );
    
    delete[] coordinates_;
    coordinates_ = new float[3*max_vertices_];
        
    int half = ( kind == HEMISPHERE );
    for ( unsigned n = 0; n < num_vertices_; ++n )
        vertices_[n].store_pos(coordinates_+3*n, half);
}


void Platonic::Solid::store_vertices(float * vec) const
{
    for ( unsigned n = 0; n < num_vertices_; ++n )
        vertices_[n].store_pos(vec+3*n, 0);
}


void Platonic::Solid::store_vertices(double * vec) const
{
    for ( unsigned n = 0; n < num_vertices_; ++n )
        vertices_[n].store_pos(vec+3*n);
}


//------------------------------------------------------------------------------
#pragma mark - edges

void Platonic::Solid::addEdge(unsigned a, unsigned b)
{
    edges_[2*num_edges_  ] = a;
    edges_[2*num_edges_+1] = b;
    ++num_edges_;
}


void Platonic::Solid::setEdges()
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


