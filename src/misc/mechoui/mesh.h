// Cytosim 3.0 -  Copyright Francois Nedelec et al.  EMBL 2015

#include "vector3.h"
class MechouiParam;

/// something in 3D
class Mesh
{
    
    /// name of OpenGL buffer
    mutable unsigned buffer;
    
private:

    /// number of points
    unsigned   n_points;
    
    /// number of faces
    unsigned   n_faces;
    
    /// coordinates of points
    float    * points;
    
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
    unsigned nbPoints() { return n_points; }
    
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
    void display(MechouiParam const&) const;

};
