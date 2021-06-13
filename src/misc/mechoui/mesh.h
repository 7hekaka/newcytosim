// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University

#include "vector3.h"
class MechouiParam;

/// something in 3D
class Mesh
{
    
    /// name of OpenGL buffer
    mutable unsigned buffer;
    
private:

    /// number of points
    size_t n_points;
    
    /// number of faces
    size_t n_faces;
    
    /// coordinates of points
    float * points;
    
    /// indices of points in the faces
    unsigned * faces;
    
    /// info on the faces
    int * labels;
    
public:
    
    /// constructor
    Mesh();
    
    /// desctructor
    ~Mesh();
    
    /// number of points
    size_t nbPoints() { return n_points; }
    
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
    
    /// display using OpenGL
    void draw(MechouiParam const&) const;

};
