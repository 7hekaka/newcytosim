// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "assert_macro.h"
#include "tesselator.h"
#include <algorithm>
#include <cmath>


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
    unsigned S = weight_[0] + weight_[1] + weight_[2];
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

Tesselator::Tesselator()
{
    apices_     = nullptr;
    num_apices_ = 0;
    
    vertices_     = nullptr;
    max_vertices_ = 0;
    num_vertices_ = 0;
    
    faces_     = nullptr;
    max_faces_ = 0;
    num_faces_ = 0;
    
    edges_     = nullptr;
    max_edges_ = 0;
    num_edges_ = 0;
    
    vex_ = nullptr;
    
    kind_ = 0;
    
    dim_[0] = 0;
    dim_[1] = 0;
    dim_[2] = 0;
    dim_[3] = 1;
}


void Tesselator::setGeometry(int K, unsigned V, unsigned E, unsigned F, unsigned N)
{
    assert_true(N > 0);
    kind_ = K;
    num_apices_ = V;
    max_vertices_ = V + E * (N-1) + F * ((N-1)*(N-2))/2;
    max_edges_ = E * N + F * 3 * ((N-1)*N)/2;
    max_faces_ = F * N * N;
}


void Tesselator::allocate()
{
    assert_true(apices_ == nullptr);
    assert_true(vertices_ == nullptr);
    assert_true(faces_ == nullptr);
    apices_ = new Apex[num_apices_];
    vertices_ = new Vertex[max_vertices_];
    faces_ = new INDEX[3*max_faces_];
}


void Tesselator::destroy()
{
    delete[] apices_;
    delete[] vertices_;
    delete(edges_);
    delete(faces_);
    delete(vex_);
    vertices_ = nullptr;
    apices_ = nullptr;
    faces_ = nullptr;
    edges_ = nullptr;
    vex_ = nullptr;
}


void Tesselator::setVertexCoordinates()
{
    assert_true( num_vertices_ <= max_vertices_ );
    assert_true( num_faces_ <= max_faces_ );
    delete(vex_);
    vex_ = new float[3*max_vertices_];
    store_vertices(vex_);
}


void scale_(size_t num, float* ptr, float X, float Y, float Z)
{
    for ( unsigned n = 0; n < num; ++n )
    {
        ptr[3*n  ] *= X;
        ptr[3*n+1] *= Y;
        ptr[3*n+2] *= Z;
    }
}


/** This transforms the sphere into a 'pin'-like smooth surface */
template < typename FLOAT >
static void dropletify_(size_t num, FLOAT* ptr)
{
    const FLOAT m = 0.75f;
    for ( unsigned n = 0; n < num; ++n )
    {
        FLOAT Z = ptr[3*n+2];
        FLOAT W = 0.5f * ( 1.f + std::tanh(-m*Z) );
        ptr[3*n  ] *= W;
        ptr[3*n+1] *= W;
    }
}

struct VertexIZ
{
    unsigned vex;
    float val;
};


/// qsort function comparing the Z component of two points
static int compareVertexZ(const void * A, const void * B)
{
    float a = ((VertexIZ const*)(A))->val;
    float b = ((VertexIZ const*)(B))->val;
    return ( a > b ) - ( a < b );
}


/** Sort vertices in Z order and reassign the face indices */
void Tesselator::sortVertices()
{
    VertexIZ * map = new VertexIZ[num_vertices_];
    unsigned * per = new unsigned[num_vertices_];

    for ( unsigned i = 0; i < num_vertices_; ++i )
    {
        Vertex const& vex = vertices_[i];
        FLOAT S = 1.0 / ( vex.weight_[0] + vex.weight_[1] + vex.weight_[2] );
        FLOAT a = vex.weight_[0] * S;
        FLOAT b = vex.weight_[1] * S;
        FLOAT c = vex.weight_[2] * S;
        FLOAT* pA = apices_[vex.index_[0]].pos_;
        FLOAT* pB = apices_[vex.index_[1]].pos_;
        FLOAT* pC = apices_[vex.index_[2]].pos_;
        FLOAT Z = a * pA[2] + b * pB[2] + c * pC[2];
        map[i] = { i, Z };
    }

    qsort(map, num_vertices_, sizeof(VertexIZ), &compareVertexZ);
    
    // copy vertices according to the new order
    Vertex * ver = new Vertex[max_vertices_];
    for ( unsigned i = 0; i < num_vertices_; ++i )
    {
        ver[i] = vertices_[map[i].vex];
        per[map[i].vex] = i;
    }
    delete[] vertices_;
    vertices_ = ver;
    
    // update vertex indices in faces:
    for ( unsigned i = 0; i < 3*num_faces_; ++i )
        faces_[i] = per[faces_[i]];
    delete[] map;
    delete[] per;
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
    //printf("vertex %p %u\n", this, wA+wB+wC);
    assert_true( wA > 0 && wB > 0 && wC >= 0 );
    if ( wC == 0 )
        return getEdgeVertex(A, wA, B, wB);
    return addVertex(A, wA, B, wB, C, wC);
}


//------------------------------------------------------------------------------
#pragma mark -

/// add intermediate vertices between `a` and `b`
void Tesselator::cutEdge(unsigned a, unsigned b, unsigned div)
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
    unsigned f = 3 * num_faces_;
    faces_[f  ] = (INDEX)a;
    faces_[f+1] = (INDEX)b;
    faces_[f+2] = (INDEX)c;
    assert_true(faces_[f  ] == a);
    assert_true(faces_[f+1] == b);
    assert_true(faces_[f+2] == c);
    ++num_faces_;
}


void Tesselator::cutFace(unsigned* line, unsigned a, unsigned b, unsigned c, unsigned div)
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

void Tesselator::cutQuad(unsigned* line, unsigned quad[4], unsigned div)
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


/*
 /// unfinished
 void Tesselator::cutStrip(unsigned cnt, unsigned inx[], unsigned div)
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
 */

void Tesselator::setApices(FLOAT vex[][3], unsigned div)
{
    for ( unsigned c = 0; c < num_apices_; ++c )
    {
        apices_[c].init(c, vex[c][0], vex[c][1], vex[c][2]);
        // also create the corresponding 'derived' vertex:
        addVertex(c, div, 0, 0, 0, 0);
    }
    num_edge_vertices_ = num_apices_;
}



void Tesselator::refineTriangles(unsigned n_fac, unsigned fac[][3], unsigned div)
{
    for ( unsigned f = 0; f < n_fac; ++f )
    {
        unsigned a = fac[f][0];
        unsigned b = fac[f][1];
        unsigned c = fac[f][2];
        cutEdge(a, b, div);
        cutEdge(b, c, div);
        cutEdge(c, a, div);
    }
    assert_true(num_vertices_ <= max_vertices_);
    
    unsigned* line = new unsigned[div+1];
    for ( unsigned f = 0; f < n_fac; ++f )
        cutFace(line, fac[f][0], fac[f][1], fac[f][2], div);
    delete[] line;
}


//------------------------------------------------------------------------------
#pragma mark -


void Tesselator::construct(Tesselator::Polyhedra kind, unsigned div, int make)
{
    switch( kind )
    {
        case UNSET: break;
        case TETRAHEDRON: buildTetrahedron(div, make); break;
        case OCTAHEDRON: buildOctahedron(div, make); break;
        case ICOSAHEDRON: buildIcosahedron(div, make); break;
        case ICOSAHEDRONX: buildIcosahedronX(div, make); break;
        case HEMISPHERE: buildHemisphere(div, make); break;
        case CYLINDER: buildCylinder(div, make); break;
        case DICE: buildDice(0.7, 0.5, 0.5, 0.3, div, div, make); break;
        case DROPLET: buildDroplet(div, make); break;
    }
}


void Tesselator::buildTetrahedron(unsigned div, int make)
{
    div = std::max(div, 1U);
    setGeometry(TETRAHEDRON, 4, 6, 4, div);
    
    constexpr FLOAT a = 1.0/3.0;
    constexpr FLOAT b = M_SQRT2/3.0;
    constexpr FLOAT c = M_SQRT2/1.7320508075688772935274463415059; //that is SQRT3
    
    // Four vertices on unit sphere
    FLOAT vex[4][3] = {
        { 0, 2*b, -a},
        {-c,  -b, -a},
        { c,  -b, -a},
        { 0,  0,   1},
    };
    
    // ordered faces: Counter-Clockwise = facing out
    unsigned fac[4][3] = {
        {0, 2, 1},
        {1, 3, 0},
        {0, 3, 2},
        {1, 2, 3}
    };
    
    allocate();
    setApices(vex, div);
    refineTriangles(4, fac, div);
    if ( make & 2 ) setVertexCoordinates();
    if ( make & 4 ) setEdges();
}


void Tesselator::buildOctahedron(unsigned div, int make)
{
    div = std::max(div, 1U);
    setGeometry(OCTAHEDRON, 6, 12, 8, div);
    
    // Eight vertices on unit sphere
    FLOAT vex[6][3] = {
        { 0,  0,  1},
        { 0,  0, -1},
        { 1,  0,  0},
        {-1,  0,  0},
        { 0, -1,  0},
        { 0,  1,  0},
    };
    
    // ordered faces: Counter-Clockwise = facing out
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
    
    allocate();
    setApices(vex, div);
    refineTriangles(8, fac, div);
    if ( make & 2 ) setVertexCoordinates();
    if ( make & 4 ) setEdges();
}


void Tesselator::buildIcosahedronX(unsigned div, int make)
{
    div = std::max(div, 1U);
    setGeometry(ICOSAHEDRON, 12, 30, 20, div);
    
    const FLOAT G = 0.5+0.5*std::sqrt(5.0);
    const FLOAT Z = 1 / std::sqrt(G*G+1.0);
    const FLOAT T = G * Z;
    
    // Twelve vertices of icosahedron on unit sphere
    FLOAT vex[12][3] = {
        { T,  Z,  0},
        {-T, -Z,  0},
        { T, -Z,  0},
        {-T,  Z,  0},
        { Z,  0,  T},
        {-Z,  0, -T},
        { Z,  0, -T},
        {-Z,  0,  T},
        { 0,  T,  Z},
        { 0, -T, -Z},
        { 0, -T,  Z},
        { 0,  T, -Z}
    };
    
    // ordered faces: Counter-Clockwise = facing out
    unsigned fac[20][3] = {
        {0,  2,  6},
        {1,  7,  3},
        {0,  4,  2},
        {1,  3,  5},
        {0,  8,  4},
        {1,  5,  9},
        {0, 11,  8},
        {1,  9, 10},
        {0,  6, 11},
        {1, 10,  7},
        {4,  8,  7},
        {5,  6,  9},
        {3,  7,  8},
        {6,  2,  9},
        {3,  8, 11},
        {9,  2, 10},
        {5,  3, 11},
        {2,  4, 10},
        {6,  5, 11},
        {4,  7, 10},
    };
    
    allocate();
    setApices(vex, div);
    refineTriangles(20, fac, div);
    if ( make & 2 ) setVertexCoordinates();
    if ( make & 4 ) setEdges();
}

/** The faces are draw in order of increasing Z */
void Tesselator::buildIcosahedron(unsigned div, int make)
{
    div = std::max(div, 1U);
    setGeometry(ICOSAHEDRON, 12, 30, 20, div);
    
    const FLOAT Z = std::sqrt(0.2);
    const FLOAT C = std::cos(M_PI/2.5);
    const FLOAT S = std::sin(M_PI/2.5);
    const FLOAT D = C*C - S*S;
    const FLOAT T = C*S + C*S;
    const FLOAT H = std::sqrt(1 + Z*Z);

    // Twelve vertices of icosahedron
    FLOAT vex[12][3] = {
        { 0,  0, -H},
        { 1,  0, -Z},
        { C, -S, -Z},
        { D, -T, -Z},
        { D,  T, -Z},
        { C,  S, -Z},
        {-D, -T,  Z},
        {-C, -S,  Z},
        {-1,  0,  Z},
        {-C,  S,  Z},
        {-D,  T,  Z},
        { 0,  0,  H}
    };
    
    // ordered faces: Counter-Clockwise = facing out
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
    
    allocate();
    setApices(vex, div);
    refineTriangles(20, fac, div);
    if ( make & 2 ) setVertexCoordinates();
    if ( make & 4 ) setEdges();
}

static void triangle(unsigned i[3], unsigned a, unsigned b, unsigned c)
{
    i[0] = a;
    i[1] = b;
    i[2] = c;
}

/**
 This extends a half icosahedron by adding rows of equilateral triangles
 The faces are draw in order of increasing Z */
void Tesselator::buildCylinder(unsigned div, int make)
{
    div = std::max(div, 1U);
    setGeometry(CYLINDER, 21, 55, 35, div);

    const FLOAT Z = std::sqrt(0.2);
    const FLOAT C = std::cos(M_PI/2.5);
    const FLOAT S = std::sin(M_PI/2.5);
    const FLOAT D = C*C - S*S;
    const FLOAT T = C*S + C*S;
    const FLOAT H = 3 * Z; // first row added
    const FLOAT W = 5 * Z; // second row added

    // Vertices for half-icosahedron & tubular extension
    FLOAT vex[21][3] = {
        { 0,  0, -1},
        { 1,  0, -Z},
        { C, -S, -Z},
        { D, -T, -Z},
        { D,  T, -Z},
        { C,  S, -Z},
        {-D, -T,  Z},
        {-C, -S,  Z},
        {-1,  0,  Z},
        {-C,  S,  Z},
        {-D,  T,  Z},
        { C, -S,  H},
        { D, -T,  H},
        { D,  T,  H},
        { C,  S,  H},
        { 1,  0,  H},
        {-C, -S,  W},
        {-1,  0,  W},
        {-C,  S,  W},
        {-D,  T,  W},
        {-D, -T,  W},
    };
    
    // ordered faces: Counter-Clockwise = facing out
    unsigned fac[40][3] = {
        {0,  1,  2},
        {0,  2,  3},
        {0,  3,  4},
        {0,  4,  5},
        {0,  5,  1},
    };
    int f = 5;
    // every row adds 10 triangles
    for ( unsigned x : { 0, 5, 10 } )
    {
        for ( unsigned i = 1; i <= 5; ++i )
        {
            unsigned ix = i + x;
            unsigned nx = (i==5) ? 1+x : i+x+1;
            unsigned px = (i==1) ? 5+x : i+x-1;
            triangle(fac[f++], nx, ix, ix+5);
            triangle(fac[f++], ix, px+5, ix+5);
        }
    }
    allocate();
    setApices(vex, div);
    refineTriangles(f, fac, div);
    sortVertices();
    if ( make & 2 ) setVertexCoordinates();
    if ( make & 4 ) setEdges();
}

void Tesselator::buildDroplet(unsigned div, int make)
{
    div = std::max(div, 1U);
    buildCylinder(div, make);
    kind_ = DROPLET;
}

static void copyme(Tesselator::FLOAT x[3], Tesselator::FLOAT a[3])
{
    x[0] = a[0];
    x[1] = a[1];
    x[2] = a[2];
}

static void middle(Tesselator::FLOAT x[3], Tesselator::FLOAT a[3], Tesselator::FLOAT b[3])
{
    x[0] = 0.5 * ( a[0] + b[0] );
    x[1] = 0.5 * ( a[1] + b[1] );
    x[2] = 0.5 * ( a[2] + b[2] );
}

void Tesselator::buildHemisphere(unsigned div, int make)
{
    div = std::max(div, 1U);
    setGeometry(HEMISPHERE, 26, 75, 40, div);
    
    const FLOAT P = M_PI/10.;
    const FLOAT D = M_PI/2.5;
    FLOAT C[5], S[5];
    for ( int i = 0; i < 5; ++i )
    {
        FLOAT a = P + D*i;
        C[i] = std::cos(a);
        S[i] = std::sin(a);
    }

    const FLOAT Z = std::sqrt(0.2);
    const FLOAT H = std::sqrt(1 + Z*Z);

    // Twelve vertices of icosahedron
    FLOAT ico[12][3] = {
        { 0,  0, -H},
        { C[0],-S[0], -Z},
        { C[1],-S[1], -Z},
        { C[2],-S[2], -Z},
        { C[3],-S[3], -Z},
        { C[4],-S[4], -Z},
        {-C[3], S[3],  Z},
        {-C[4], S[4],  Z},
        {-C[0], S[0],  Z},
        {-C[1], S[1],  Z},
        {-C[2], S[2],  Z},
        { 0,  0,  H}
    };
    
    int j = 0, f = 0;
    FLOAT vex[26][3];
    unsigned fac[40][3];
    
    copyme(vex[j++], ico[0]);
    for ( int i = 1; i <= 5; ++i )
    {
        middle(vex[j++], ico[0], ico[i]);
        triangle(fac[f++], 0, i, i<5?i+1:1);
    }
    for ( int i = 1; i <= 5; ++i )
    {
        unsigned n = i<5 ? i+1 : 1;
        middle(vex[j++], ico[i], ico[n]);
        triangle(fac[f++], n, i, i+5);
    }
    for ( int i = 1; i <= 5; ++i )
    {
        unsigned p = i<2 ? 5 : i-1;
        copyme(vex[j++], ico[i]);
        triangle(fac[f++], i, p+5, i+10);
        triangle(fac[f++], i, i+10, i+5);
    }
    for ( int i = 1; i <= 5; ++i )
    {
        unsigned n = i<5 ? i+1 : 1;
        middle(vex[j++], ico[n], ico[i+5]);
        triangle(fac[f++], n+10, i+5, i+15);
    }
    for ( int i = 1; i <= 5; ++i )
    {
        unsigned p = i<2 ? 5 : i-1;
        middle(vex[j++], ico[i], ico[i+5]);
        triangle(fac[f++], i+10, p+15, i+20);
        triangle(fac[f++], i+5, i+20, i+15);
        triangle(fac[f++], i+5, i+10, i+20);
    }
    allocate();
    setApices(vex, div);
    refineTriangles(f, fac, div);
    if ( make & 2 ) setVertexCoordinates();
    if ( make & 4 ) setEdges();
}


/**
 This divides the edges by 'div' and the 6 square faces bi 'vid'
 */
void Tesselator::buildDice(FLOAT X, FLOAT Y, FLOAT Z, FLOAT R, unsigned div, unsigned vid, int make)
{
    setGeometry(DICE, 24, 90, 44, div);
    
    dim_[0] = X;
    dim_[1] = Y;
    dim_[2] = Z;
    dim_[3] = R;
    
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
    
    allocate();
    setApices(vex, div);
    for ( unsigned q = 0; q < 18; ++q )
    {
        unsigned a = quad[q][0];
        unsigned b = quad[q][1];
        unsigned c = quad[q][2];
        unsigned d = quad[q][3];
        cutEdge(a, b, div);
        cutEdge(b, c, div);
        cutEdge(c, a, div);
        cutEdge(b, d, div);
        cutEdge(c, d, div);
    }
    refineTriangles(8, fac, div);
    
    unsigned* line = new unsigned[div+1];
    // 4 edges parallel to the Z axis
    for ( int n = 6; n < 10; ++n )
        cutQuad(line, quad[n], div);
    
    // 4 edges parallel to the Y axis
    for ( int n = 10; n < 14; ++n )
        cutQuad(line, quad[n], div);
    
    // 4 edges parallel to the X axis
    for ( int n = 14; n < 18; ++n )
        cutQuad(line, quad[n], div);
    
    // faces
    for ( int n = 0; n < 6; ++n )
        cutQuad(line, quad[n], vid);
    
    delete[] line;
    if ( make & 2 ) setVertexCoordinates();
    if ( make & 4 ) setEdges();
}


//------------------------------------------------------------------------------
#pragma mark - Solid vertices

template < typename REAL >
static void scale_(REAL& X, REAL& Y, REAL& Z, REAL n)
{
    X *= n;
    Y *= n;
    Z *= n;
}


template < typename REAL >
static void projectSphere(REAL* X)
{
    REAL n = 1.0 / std::sqrt(X[0]*X[0] + X[1]*X[1] + X[2]*X[2]);
    X[0] *= n;
    X[1] *= n;
    X[2] *= n;
}


template < typename REAL >
static void projectCylinder(REAL* X)
{
    if ( X[2] < 0 )
    {
        REAL n = 1.0 / std::sqrt(X[0]*X[0] + X[1]*X[1] + X[2]*X[2]);
        X[0] *= n;
        X[1] *= n;
        X[2] *= n;
    }
    else
    {
        REAL n = 1.0 / std::sqrt(X[0]*X[0] + X[1]*X[1]);
        X[0] *= n;
        X[1] *= n;
    }
}


template < typename REAL >
static void projectDice(REAL* X, const REAL len[4])
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


void Tesselator::interpolate(Vertex const& vex, double ptr[3]) const
{
    double S = 1.0 / ( vex.weight_[0] + vex.weight_[1] + vex.weight_[2] );
    double a = (double)vex.weight_[0] * S;
    double b = (double)vex.weight_[1] * S;
    double c = (double)vex.weight_[2] * S;
    auto pA = apices_[vex.index_[0]].pos_;
    auto pB = apices_[vex.index_[1]].pos_;
    auto pC = apices_[vex.index_[2]].pos_;
    ptr[0] = a * pA[0] + b * pB[0] + c * pC[0];
    ptr[1] = a * pA[1] + b * pB[1] + c * pC[1];
    ptr[2] = a * pA[2] + b * pB[2] + c * pC[2];
}


void Tesselator::interpolate(Vertex const& vex, float ptr[3]) const
{
    float S = 1.0 / ( vex.weight_[0] + vex.weight_[1] + vex.weight_[2] );
    float a = (float)vex.weight_[0] * S;
    float b = (float)vex.weight_[1] * S;
    float c = (float)vex.weight_[2] * S;
    auto pA = apices_[vex.index_[0]].pos_;
    auto pB = apices_[vex.index_[1]].pos_;
    auto pC = apices_[vex.index_[2]].pos_;
    ptr[0] = a * pA[0] + b * pB[0] + c * pC[0];
    ptr[1] = a * pA[1] + b * pB[1] + c * pC[1];
    ptr[2] = a * pA[2] + b * pB[2] + c * pC[2];
}


void Tesselator::store_vertices(float * vec) const
{
    for ( unsigned n = 0; n < num_vertices_; ++n )
        interpolate(vertices_[n], vec+3*n);
    
    if ( kind_ == DICE )
    {
        float len[4] = { (float)dim_[0], (float)dim_[1], (float)dim_[2], (float)dim_[3] };
        for ( unsigned n = 0; n < num_vertices_; ++n )
            projectDice(vec+3*n, len);
    }
    else if ( kind_ == CYLINDER || kind_ == DROPLET )
    {
        for ( unsigned n = 0; n < num_vertices_; ++n )
            projectCylinder(vec+3*n);
        if ( kind_ == DROPLET )
            dropletify_(num_vertices_, vec);
    }
    else
    {
        for ( unsigned n = 0; n < num_vertices_; ++n )
            projectSphere(vec+3*n);
    }
}


void Tesselator::store_vertices(double * vec) const
{
    for ( unsigned n = 0; n < num_vertices_; ++n )
        interpolate(vertices_[n], vec+3*n);
    
    if ( kind_ == DICE )
    {
        double len[4] = { dim_[0], dim_[1], dim_[2], dim_[3] };
        for ( unsigned n = 0; n < num_vertices_; ++n )
            projectDice(vec+3*n, len);
    }
    else if ( kind_ == CYLINDER || kind_ == DROPLET )
    {
        for ( unsigned n = 0; n < num_vertices_; ++n )
            projectCylinder(vec+3*n);
        if ( kind_ == DROPLET )
            dropletify_(num_vertices_, vec);
    }
    else
    {
        for ( unsigned n = 0; n < num_vertices_; ++n )
            projectSphere(vec+3*n);
    }
}


//------------------------------------------------------------------------------
#pragma mark - edges

void Tesselator::addEdge(unsigned a, unsigned b)
{
    if ( num_edges_ < max_edges_ )
    {
        //printf("edge: %i %i\n", a, b);
        size_t e = 2 * num_edges_++;
        edges_[e  ] = (INDEX)a;
        edges_[e+1] = (INDEX)b;
        assert_true((INDEX)a == a);
        assert_true((INDEX)b == b);
    }
    else
        printf("Tesselator::addEdge overflow (%i): %i %i\n", max_edges_, a, b);

}


void Tesselator::setEdges()
{
    delete[] edges_;
    edges_ = new INDEX[2*max_edges_];
    num_edges_ = 0;
    
    // build edges from the faces:
    for ( unsigned i = 0; i < num_faces_; ++i )
    {
        unsigned a = faces_[3*i  ];
        unsigned b = faces_[3*i+1];
        unsigned c = faces_[3*i+2];
        //printf("face %i: %i %i %i\n", i, a, b, c);
        /*
         Here we rely on the fact that every edge will be covered by two faces,
         and thus we add this test to keep only half of them.
         This may fail if the volume is not closed.
         */
        if ( a < b ) addEdge(a, b);
        if ( b < c ) addEdge(b, c);
        if ( c < a ) addEdge(c, a);
    }
    //printf("ico: %i vertices, %i edges\n", num_faces_, num_edges_);
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
