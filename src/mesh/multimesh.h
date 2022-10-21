// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University

#include <cstdio>
#include "real.h"


/// MultiMesh structure from Christopher Betty's Multitrack
class MultiMesh
{    
private:

    /// number of points
    size_t n_points;
    
    /// number of faces
    size_t n_faces;
    
    /// coordinates of points
    real * points;
    
    /// indices of points in the faces
    unsigned * faces;
    
    /// info on the faces
    int * labels;
    
public:
    
    /// constructor
    MultiMesh();
    
    /// desctructor
    ~MultiMesh();
    
    /// number of points
    size_t nbPoints() { return n_points; }
    
    /// number of points
    size_t nbFaces() { return n_faces; }

    /// release memory
    void release();
    
    /// read from file
    int read(char const* filename);
    
    /// read from file
    int read_ascii(FILE * file);
    
    /// read from file
    int read_binary(FILE * file);
    
    /// OpenGL picking function
    unsigned pick() const;
    
    /// draw vertex points using OpenGL
    void drawPoints(float size, const float color[4]) const;
    
    /// display triangles using OpenGL
    void drawFaces(const float color[4], int selected) const;
    
    /// display triangles using OpenGL
    void drawFaces(const float dir[3], const float color[4], int selected) const;

};
