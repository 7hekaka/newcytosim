// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef PLATONIC_H
#define PLATONIC_H

#include <iostream>

namespace Platonic
{
    /// floating type used by Platonic
    typedef float FLOAT;

    /// One of the corner of a Platonic solid
    struct Corner
    {
        /// Coordinates in space
        FLOAT pos_[3];
        
        /// an index to identify this vertex
        unsigned inx_;
        
        Corner()
        { inx_=0; pos_[0]=0; pos_[1]=0; pos_[2]=0; }
        
        ~Corner() {}
        
        void init(unsigned n, FLOAT x, FLOAT y, FLOAT z)
        { inx_=n; pos_[0]=x; pos_[1]=y; pos_[2]=z; }
    };
    
    /// A vertex is interpolated from 3 Corners
    class Vertex
    {
    private:
        
        void swap(unsigned, unsigned);
        
        bool smaller(unsigned, unsigned);

    public:
        
        /// pointers to the corners being interpolated
        unsigned index_[3];
        
        /// Weights of the interpolation
        unsigned weight_[3];
        
        /// check if weights are equal
        bool equivalent(unsigned, unsigned, unsigned, unsigned) const;
        
        ///
        Vertex() { index_[0]=-1; index_[1]=-1; index_[2]=-1; }
        
        //Vertex(Corner *, unsigned, Corner *, unsigned, Corner *, unsigned) { set(); }
        
        ~Vertex() {}
        
        void set(unsigned, unsigned, unsigned, unsigned, unsigned, unsigned);

        unsigned weight(int x) const { return weight_[x]; }
        
        unsigned sum_weights() const { return weight_[0]+weight_[1]+weight_[2]; }
        
        void     print(unsigned, std::ostream&) const;
    };
    
    
    /// A refined polyhedra to serve as a tesselation
    /**
     Platonic solids made of triangles can be refined by subdividing the faces
     into smaller triangles. The faces are re-assembled together into the Solid
     without duplicating points.
     
     The level of refinement is set by an integer N > 0, corresponding to the 
     number of section in which each edge of the original Platonic solid is divided.
     
     */
    class Solid
    {
    public:
        /// regular polyhedra made of triangles
        enum Polyhedra { UNSET=0, TETRAHEDRON=1, OCTAHEDRON=2, ICOSAHEDRON=3, HEMISPHERE=4, DICE=5 };

    private:
        /// dimensions
        FLOAT    length_[4];
        
        /// Array of coordinates of all vertices
        float  * coordinates_;
        
        /// Array of primary vertices of the geometry
        Corner * corners_;
        
        /// Array of derived vertices
        Vertex * vertices_;
        
        /// Array of indices of the points making the edges
        unsigned *edges_;
        
        /// Array of indices of the points making the faces
        unsigned *faces_;
        
        /// number of primary vertices
        unsigned num_corners_;
        
        unsigned num_vertices_, max_vertices_;
        
        /// number of vertices on the edges between primary corners
        unsigned num_edge_vertices_;
        
        /// number of faces
        unsigned num_faces_, max_faces_;
        
        /// number of edges
        unsigned num_edges_, max_edges_;
        
        /// defining which Z half space for hemisphere
        FLOAT halfZ_;
        
        /// 
        int kind_;
        
        unsigned findEdgeVertex(unsigned, unsigned, unsigned, unsigned) const;
        unsigned getEdgeVertex(unsigned, unsigned, unsigned, unsigned) const;
        unsigned addVertex(unsigned, unsigned, unsigned, unsigned, unsigned, unsigned);
        unsigned makeVertex(unsigned, unsigned, unsigned, unsigned, unsigned, unsigned);
        
        void setCorners(unsigned, FLOAT vex[][3], unsigned div);
        void refineTriangles(unsigned, unsigned fac[][3], unsigned div);
        
        void addFace(unsigned, unsigned, unsigned);
        void addEdge(unsigned, unsigned);
        void refineEdge(unsigned a, unsigned b, unsigned div);
        void refineFace(unsigned*, unsigned a, unsigned b, unsigned c, unsigned div);
        void refineQuad(unsigned*, unsigned quad[4], unsigned div);
        void refineStrip(unsigned cnt, unsigned inx[], unsigned div);
        
        void build();
        void build(Polyhedra K, unsigned div);
        void allocate(unsigned nv, unsigned ne, unsigned nf, unsigned div);

        void interpolate(Vertex const&, float vec[3], int half) const;
        void interpolate(Vertex const&, double vec[3], int half) const;

    public:
        
        void initTetrahedron(unsigned div);
        void initOctahedron(unsigned div);
        void initIcosahedron(unsigned div);
        void initIcosahedronRotated(unsigned div);
        void initHemisphere(unsigned div);
        void initDice(FLOAT X, FLOAT Y, FLOAT Z, FLOAT R, unsigned div);

        /// build as polyhedra refined by order `div`
        Solid(Polyhedra, unsigned div, int make = 0);
        
        /// build as empty structure
        Solid() { build(); }

        /// destructor
        ~Solid();
        
        /// set array of indices that define the edges
        void     setEdges();
        
        /// number of derived vertices
        unsigned nb_vertices() const { return num_vertices_; }
        
        /// reference to derived vertex `ii`
        Vertex&  vertex(int i) const { return vertices_[i]; }
        
        /// copy coordinates of points to given array
        void     store_vertices(float* vec) const;
        
        /// copy coordinates of points to given array
        void     store_vertices(double* vec) const;

        /// return pointer to array of coordinates of vertices
        const float* vertex_data()      const { return coordinates_; }
        
        /// address of coordinates for vertex `v` ( `v < nb_vertices()` )
        const float* vertex_data(int v) const { return coordinates_ + 3 * v; }
        
        
        /// number of points in the edges = 2 * nb-of-edges
        unsigned int nb_edges()        const { return num_edges_; }
        
        /// array of indices to the vertices in each edge (2 per edge)
        unsigned int* edge_indices()   const { return edges_; }
        
        /// return address of first vertex of edge `e`
        unsigned edge_index0(int e)    const { return edges_[2*e]; }
        
        /// return address of second vertex of edge `e`
        unsigned edge_index1(int e)    const { return edges_[2*e+1]; }

        /// return address of first vertex of edge `e`
        const float* edge_vertex0(int e) const { return coordinates_ + 3 * edges_[2*e]; }
        
        /// return address of second vertex of edge `e`
        const float* edge_vertex1(int e) const { return coordinates_ + 3 * edges_[2*e+1]; }
        
        
        /// number of faces (each face is a triangle of 3 vertices)
        unsigned int nb_faces()        const { return num_faces_; }
        
        /// array of indices to the vertices in each face (3 vertices per face)
        unsigned int* face_data()     const { return faces_; }
        
        /// return address of first vertex of face `f`
        const float* face_data0(int f) const { return coordinates_ + 3 * faces_[3*f]; }
        /// return address of second vertex of face `f`
        const float* face_data1(int f) const { return coordinates_ + 3 * faces_[3*f+1]; }
        /// return address of third vertex of face `f`
        const float* face_data2(int f) const { return coordinates_ + 3 * faces_[3*f+2]; }
        
        /// return index of first vertex in face `f`
        unsigned face_indx0(int f)     const { return faces_[3*f]; }
        /// return index of second vertex in face `f`
        unsigned face_indx1(int f)     const { return faces_[3*f+1]; }
        /// return index of third vertex in face `f`
        unsigned face_indx2(int f)     const { return faces_[3*f+2]; }
    };
}

#endif
