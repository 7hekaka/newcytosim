// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

#ifndef TESSELATOR_H
#define TESSELATOR_H

#include <cstdio>
#include <math.h>

/// Provides a triangulation of a surface made by refining a Platonic solid
/**
 Platonic solids made of triangles are refined by subdividing the faces
 into smaller triangles. The faces are re-assembled to cover the surface
 without duplicating vertices.
 
 The level of refinement is set by an integer N > 0, corresponding to the
 number of section in which each edge of the original Platonic solid is divided.
 */
class Tesselator
{
    /// Disabled copy constructor
    Tesselator(Tesselator const&);
    
    /// disabled assignment operator
    Tesselator& operator = (const Tesselator&);

public:
    
    /// floating type used for calculations
    typedef float FLOAT;

    /// integer type used for vertex indexing
    typedef unsigned short INDEX;

    /// starting shapes
    enum Polyhedra { UNSET=0, TETRAHEDRON, OCTAHEDRON, ICOSAHEDRON, 
        ICOSAHEDRONX, HEMISPHERE, FOOTBALL, CYLINDER, DICE, DROPLET, PIN };
    
    /// One of the vertex of the unrefined template model
    struct Apex
    {
        /// Coordinates in space
        FLOAT pos_[3];
        
        /// an index to identify this vertex
        unsigned inx_;
        
        Apex()
        { inx_=0; pos_[0]=0; pos_[1]=0; pos_[2]=0; }
        
        void init(unsigned n, FLOAT x, FLOAT y, FLOAT z)
        { inx_=n; pos_[0]=x; pos_[1]=y; pos_[2]=z; }
    };
    
    /// A vertex is interpolated from 3 Apex
    class Vertex
    {
    public:
        
        /// pointers to the apices being interpolated
        unsigned index_[3];
        
        /// Coefficients of the interpolation, before normalization
        unsigned weight_[3];
        
        /// check if weights are equal
        bool equivalent(unsigned, unsigned, unsigned, unsigned) const;
        
        ///
        Vertex() { index_[0]=-1; index_[1]=-1; index_[2]=-1; }
                
        void set(unsigned, unsigned);
        void set(unsigned, unsigned, unsigned, unsigned, unsigned, unsigned);

        unsigned weight(int x) const { return weight_[x]; }
        
        void print(unsigned, FILE*) const;
    };
    
    
private:
    /// dimensions
    FLOAT dim_[4];
    
    /// Array of coordinates of all vertices
    float * vex_;
    
    /// Array of primary vertices of the geometry
    Apex * apices_;
    
    /// Array of derived vertices
    Vertex * vertices_;
    
    /// Array of indices of the points making the edges
    INDEX * edges_;
    
    /// Array of indices of the points making the faces
    INDEX * faces_;
    
    /// number of primary vertices
    unsigned num_apices_;
    
    unsigned num_vertices_, max_vertices_;
    
    /// number of vertices on the edges between primary corners
    unsigned num_edge_vertices_;
    
    /// number of faces
    unsigned num_faces_, max_faces_;
    
    /// number of edges
    unsigned num_edges_, max_edges_;
    
    /// type
    int kind_;
    
    /// number of subdivision in an edge
    int rank_;

    unsigned findEdgeVertex(unsigned, unsigned, unsigned, unsigned) const;
    unsigned getEdgeVertex(unsigned, unsigned, unsigned, unsigned) const;
    unsigned addVertex(unsigned, unsigned, unsigned, unsigned, unsigned, unsigned);
    unsigned makeVertex(unsigned, unsigned, unsigned, unsigned, unsigned, unsigned);
    
    void setApices(FLOAT vex[][3], unsigned div);
    void refineTriangles(unsigned, unsigned fac[][3], unsigned div);
    
    void addFace(unsigned, unsigned, unsigned);
    void addEdge(unsigned, unsigned);
    void cutEdge(unsigned a, unsigned b, unsigned div);
    void cutEdges(unsigned, unsigned[][2], unsigned div);
    void cutFace(unsigned a, unsigned b, unsigned c, unsigned div, unsigned*);
    void cutQuad(unsigned quad[4], unsigned div, unsigned*);
    void cutStrip(unsigned inx[6], unsigned div, unsigned*);
    
    void allocate();
    void destroy();
    void setGeometry(int K, unsigned V, unsigned E, unsigned F, unsigned div);

    template <typename REAL> void interpolateVertex(REAL vec[3], Vertex const&) const;
    template <typename REAL> void projectFootball(REAL vec[3], Vertex const&) const;
    template <typename REAL> void projectVertices(REAL*) const;

public:

    /// build as empty structure
    Tesselator();
    
    /// destructor
    ~Tesselator() { destroy(); }
    
    /// build as polyhedra refined by order `div`
    void construct(Polyhedra, unsigned div, int make = 0);

    void buildTetrahedron(unsigned div);
    void buildOctahedron(unsigned div);
    void buildIcosahedron(unsigned div);
    void buildIcosahedronS(unsigned div);
    void buildIcosahedronX(unsigned div);
    void buildCylinder(unsigned div);
    void buildHemisphere(unsigned div);
    void buildDice(FLOAT X, FLOAT Y, FLOAT Z, FLOAT R, unsigned div, unsigned vid);
    void buildFootball(unsigned div);
    void buildDroplet(unsigned div);
    void buildPin(unsigned div);

    /// set array of indices that define the edges
    void setEdges();
    /// sort faces
    void sortFaces(unsigned);
    /// sort vertices in Z and update faces
    void sortVertices();

    /// calculate coordinates of vertices used in vertex_data()
    void setVertexCoordinates();

    /// reference to derived vertex `ii`
    Vertex& vertex(int i) const { return vertices_[i]; }
    
    /// copy coordinates of points to given array
    void store_vertices(float* vec) const;
    
    /// copy coordinates of points to given array
    void store_vertices(double* vec) const;

    /** Scale all 3D points by (X, Y, Z) */
    static void scale3D(size_t num, float* ptr, float X, float Y, float Z);

    /// number of derived vertices
    unsigned max_vertices() const { return max_vertices_; }
       
    /// number of faces (each face is a triangle of 3 vertices)
    unsigned max_faces() const { return max_faces_; }

    /// number of subdivisions
    unsigned rank() const { return rank_; }
    
    /// number of derived vertices
    unsigned num_vertices() const { return num_vertices_; }
    
    /// return pointer to array of coordinates of vertices, initialized in setVertices()
    const float* vertex_data() const { return vex_; }
    
    /// address of coordinates for vertex `v` ( `v < num_vertices()` )
    const float* vertex_data(int v) const { return vex_ + 3 * v; }
    
    /// number of points in the edges = 2 * nb-of-edges
    unsigned num_edges() const { return num_edges_; }
    
    /// array of indices to the vertices in each edge (2 per edge)
    INDEX * edge_data() const { return edges_; }
    
    /// number of faces (each face is a triangle of 3 vertices)
    unsigned num_faces() const { return num_faces_; }
    
    /// return threshold for triangles belonging to the football's pentagons
    unsigned foot_rank() const { return rint(rank_*0.3333); }
    
    /// number of faces making up the football's pentagons
    unsigned num_foot_faces() const { unsigned R = foot_rank(); return 60*(R*R); }

    /// array of indices to the vertices in each face (3 vertices per face)
    INDEX * face_data() const { return faces_; }
    
    /// number of faces (each face is a triangle of 3 vertices)
    size_t face_data_size() const { return num_faces_ * 3 * sizeof(INDEX); }

    /// export ascii PLY format
    void exportPLY(FILE *) const;
    
    /// export binary STL format
    void exportSTL(FILE *) const;

    
    /// return address of first vertex of edge `e`
    const float* edge_vertex0(unsigned e) const { return vex_ + 3 * edges_[2*e]; }
    /// return address of second vertex of edge `e`
    const float* edge_vertex1(unsigned e) const { return vex_ + 3 * edges_[2*e+1]; }
    
    /// return address of first vertex of face `f`
    const float* face_vertex0(unsigned f) const { return vex_ + 3 * faces_[3*f]; }
    /// return address of second vertex of face `f`
    const float* face_vertex1(unsigned f) const { return vex_ + 3 * faces_[3*f+1]; }
    /// return address of third vertex of face `f`
    const float* face_vertex2(unsigned f) const { return vex_ + 3 * faces_[3*f+2]; }
};

#endif
