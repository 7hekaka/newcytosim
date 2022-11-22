// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University

#ifndef TESSELATOR_H
#define TESSELATOR_H

#include <cstdio>

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
    Tesselator& operator =(const Tesselator&);

public:
    
    /// floating type used for calculations
    typedef float FLOAT;

    /// floating type used for calculations
    typedef unsigned short INDEX;

    /// starting shapes
    enum Polyhedra { UNSET=0, TETRAHEDRON=1, OCTAHEDRON=2, ICOSAHEDRON=3,
        ICOSAHEDRONX=4, HEMISPHERE=5, CYLINDER=6, DICE=7, DROPLET=8 };
    
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
    private:
        
        void swap(unsigned, unsigned);
        
        bool smaller(unsigned, unsigned);
        
    public:
        
        /// pointers to the apices being interpolated
        unsigned index_[3];
        
        /// Coefficients of the interpolation, before normalization
        unsigned weight_[3];
        
        /// check if weights are equal
        bool equivalent(unsigned, unsigned, unsigned, unsigned) const;
        
        ///
        Vertex() { index_[0]=-1; index_[1]=-1; index_[2]=-1; }
                
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
    
    ///
    int kind_;
    
    unsigned findEdgeVertex(unsigned, unsigned, unsigned, unsigned) const;
    unsigned getEdgeVertex(unsigned, unsigned, unsigned, unsigned) const;
    unsigned addVertex(unsigned, unsigned, unsigned, unsigned, unsigned, unsigned);
    unsigned makeVertex(unsigned, unsigned, unsigned, unsigned, unsigned, unsigned);
    
    void setApices(FLOAT vex[][3], unsigned div);
    void refineTriangles(unsigned, unsigned fac[][3], unsigned div);
    
    void addFace(unsigned, unsigned, unsigned);
    void addEdge(unsigned, unsigned);
    void cutEdge(unsigned a, unsigned b, unsigned div);
    void cutFace(unsigned*, unsigned a, unsigned b, unsigned c, unsigned div);
    void cutQuad(unsigned*, unsigned quad[4], unsigned div);
    void cutStrip(unsigned cnt, unsigned inx[], unsigned div);
    
    void allocate();
    void destroy();
    void setGeometry(int K, unsigned V, unsigned E, unsigned F, unsigned div);

    void interpolate(Vertex const&, float vec[3]) const;
    void interpolate(Vertex const&, double vec[3]) const;
    
public:

    /// build as empty structure
    Tesselator();
    
    /// destructor
    ~Tesselator() { destroy(); }
    
    /// build as polyhedra refined by order `div`
    void construct(Polyhedra, unsigned div, int make = 0);

    void buildTetrahedron(unsigned div, int make = 1);
    void buildOctahedron(unsigned div, int make = 1);
    void buildIcosahedron(unsigned div, int make = 1);
    void buildIcosahedronX(unsigned div, int make = 1);
    void buildCylinder(unsigned div, int make = 1);
    void buildHemisphere(unsigned div, int make = 1);
    void buildDice(FLOAT X, FLOAT Y, FLOAT Z, FLOAT R, unsigned div, unsigned vid, int make);
    void buildDroplet(unsigned div, int make = 1);

    /// set array of indices that define the edges
    void setEdges();
    /// calculate coordinates of vertices used in vertex_data()
    void sortVertices();
    /// calculate coordinates of vertices used in vertex_data()
    void setVertexCoordinates();

    /// reference to derived vertex `ii`
    Vertex& vertex(int i) const { return vertices_[i]; }
    
    /// copy coordinates of points to given array
    void store_vertices(float* vec) const;
    
    /// copy coordinates of points to given array
    void store_vertices(double* vec) const;
    
    /// transform sphere into a droplet-like surface
    void dropletify(float*, float Z);
    void dropletify(double*, double Z);

    
    /// number of derived vertices
    unsigned max_vertices() const { return max_vertices_; }
       
    /// number of faces (each face is a triangle of 3 vertices)
    unsigned int max_faces() const { return max_faces_; }

    
    /// number of derived vertices
    unsigned num_vertices() const { return num_vertices_; }
    
    /// return pointer to array of coordinates of vertices, initialized in setVertices()
    const float* vertex_data() const { return vex_; }
    
    /// address of coordinates for vertex `v` ( `v < num_vertices()` )
    const float* vertex_data(int v) const { return vex_ + 3 * v; }
    
    /// number of points in the edges = 2 * nb-of-edges
    unsigned int num_edges() const { return num_edges_; }
    
    /// array of indices to the vertices in each edge (2 per edge)
    INDEX * edge_data() const { return edges_; }
    
    /// number of faces (each face is a triangle of 3 vertices)
    unsigned int num_faces() const { return num_faces_; }
    
    /// array of indices to the vertices in each face (3 vertices per face)
    INDEX * face_data() const { return faces_; }
    
    /// export ascii PLY format
    void exportPLY(FILE *) const;
    
    /// export binary STL format
    void exportSTL(FILE *) const;

    
    /// return address of first vertex of edge `e`
    const float* edge_vertex0(int e) const { return vex_ + 3 * edges_[2*e]; }
    /// return address of second vertex of edge `e`
    const float* edge_vertex1(int e) const { return vex_ + 3 * edges_[2*e+1]; }
    /// return address of first vertex of face `f`
    const float* face_vertex0(int f) const { return vex_ + 3 * faces_[3*f]; }
    /// return address of second vertex of face `f`
    const float* face_vertex1(int f) const { return vex_ + 3 * faces_[3*f+1]; }
    /// return address of third vertex of face `f`
    const float* face_vertex2(int f) const { return vex_ + 3 * faces_[3*f+2]; }
};

#endif
