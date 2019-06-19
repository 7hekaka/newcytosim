// Cytosim 3.0 -  Copyright Francois Nedelec et al.  EMBL 2015

#include "mesh.h"
#include <stdlib.h>
#include "opengl.h"
#include "mechoui_param.h"
#include "gle_color_list.h"




Mesh::Mesh()
{
    n_points = 0;
    n_faces = 0;
    points = 0;
    faces = 0;
    labels = 0;
    buffer = 0;
}


Mesh::~Mesh()
{
    release();
}


void Mesh::release()
{
    if ( points ) delete(points);
    if ( faces ) delete(faces);
    if ( labels ) delete(labels);
    points = 0;
    faces = 0;
    labels = 0;
    if ( buffer )
        glDeleteBuffers(1, &buffer);
}


int Mesh::read_ascii(FILE * file)
{
    char str[1024], * ptr;
    
    release();
    
    // read vertices:
    if ( 0 == fgets(str, sizeof(str), file) )
        return 1;
    n_points = strtol(str, 0, 10);
    points = new float[3*n_points];
    
    for ( unsigned n = 0; n < n_points; ++n )
    {
        if ( 0 == fgets(str, sizeof(str), file) )
            return 1;
        points[3*n  ] = strtof(str, &ptr);
        points[3*n+1] = strtof(ptr, &ptr);
        points[3*n+2] = strtof(ptr, &ptr);
    }
    
    // read faces:
    if ( 0 == fgets(str, sizeof(str), file) )
        return 1;
    n_faces = strtol(str, 0, 10);
    faces = new unsigned[3*n_faces];
    labels = new int[2*n_faces];

    for ( unsigned n = 0; n < n_faces; ++n )
    {
        if ( 0 == fgets(str, sizeof(str), file) )
            return 1;
        faces[3*n  ] = strtol(str, &ptr, 10);
        faces[3*n+1] = strtol(ptr, &ptr, 10);
        faces[3*n+2] = strtol(ptr, &ptr, 10);
        labels[2*n  ] = strtol(ptr, &ptr, 10);
        labels[2*n+1] = strtol(ptr, &ptr, 10);
    }

    return 0;
}


int Mesh::read_binary(FILE * file)
{
    double d[3];
    size_t s[3];
    int i[2];
    
    release();
    
    // read vertices:
    if ( 1 != fread(s, sizeof(size_t), 1, file) )
        return 1;
    
    n_points = s[0];
    points = new float[3*n_points];
    
    for ( unsigned n = 0; n < n_points; ++n )
    {
        if ( 3 != fread(d, sizeof(double), 3, file) )
            return 2;
        points[3*n  ] = d[0];
        points[3*n+1] = d[1];
        points[3*n+2] = d[2];
    }
    
    // read faces:
    if ( 1 != fread(s, sizeof(size_t), 1, file) )
        return 3;
    
    n_faces = s[0];
    faces = new unsigned[3*n_faces];
    labels = new int[2*n_faces];
    
    for ( unsigned n = 0; n < n_faces; ++n )
    {
        fread(s, sizeof(size_t), 3, file);
        if ( 2 != fread(i, sizeof(int), 2, file) )
            return 4;
        faces[3*n  ]  = s[0];
        faces[3*n+1]  = s[1];
        faces[3*n+2]  = s[2];
        labels[2*n  ] = i[0];
        labels[2*n+1] = i[1];
    }
    
    return 0;
}




int Mesh::read(char const* filename)
{
    FILE * f = fopen(filename, "r");
    if ( f )
    {
        if ( ferror(f) )
        {
            fclose(f);
            return 2;
        }
/*
        std::clog << " Reading `" << filename << "':";
*/
        bool binary = false;
        std::string str(filename);
        size_t pos = str.find(".");
        if ( pos != std::string::npos )
            binary = ( "rec" == str.substr(pos+1) );

        int res = binary?read_binary(f):read_ascii(f);
        fclose(f);
/*
        if ( res )
            std::clog << "   error " << res << std::endl;
        else
            std::clog << "   " << n_points << " points" << " and " << n_faces << " faces" << std::endl;
*/
        return res;
    }
    std::clog << " Cannot read `" << filename << std::endl;
    return 1;
}



/// structure to depth-sort triangles
struct Triangle
{
public:
    unsigned a, b, c;
    float z;
};


/// function to sort Triangles according to their 'z'
static int closer(const void * ap, const void * bp)
{
    Triangle const* a = static_cast<const Triangle*>(ap);
    Triangle const* b = static_cast<const Triangle*>(bp);
    
    if ( a->z > b->z ) return  1;
    if ( a->z < b->z ) return -1;
    return 0;
}


/// OpenGL display function
void Mesh::display(MechouiParam const& pam) const
{
    glEnable(GL_NORMALIZE);
    glEnableClientState(GL_VERTEX_ARRAY);
    
    //Create one buffer
    if ( buffer == 0 )
        glGenBuffers(1, &buffer);
    
    //Make the new VBO active
    glBindBuffer(GL_ARRAY_BUFFER, buffer);
    
    //Upload vertex data to the video device
    glBufferData(GL_ARRAY_BUFFER, 3*n_points*sizeof(float), points, GL_STREAM_DRAW);

    glVertexPointer(3, GL_FLOAT, 0, 0);

    if ( pam.point_style )
    {
        glPointSize(pam.point_size);
        pam.point_color.load();
        glDrawArrays(GL_POINTS, 0, n_points);
    }
    
    // get current modelview transformation:
    GLfloat mat[16];
    glGetFloatv(GL_MODELVIEW_MATRIX, mat);
    
    // extract axis corresponding to vertical direction:
    GLfloat ver[] = { mat[2], mat[6], mat[10] };

    glEnable(GL_LIGHTING);
    glDisable(GL_CULL_FACE);
    glDepthMask(GL_TRUE);
    glMateriali (GL_FRONT_AND_BACK, GL_SHININESS, 64);

    Triangle * tris = new Triangle[n_faces];
    unsigned n_tris = 0;
    
    for ( unsigned n = 0; n < n_faces; ++n )
    {
        float * a = points+3*faces[3*n  ];
        float * b = points+3*faces[3*n+1];
        float * c = points+3*faces[3*n+2];
        
        int cell = labels[2*n];
        int boundary = labels[2*n+1];
        
        bool transparent = boundary > 0;
        
        // if one cell is selected, all the other ones are transparent:
        if ( pam.selected )
        {
            transparent = ( cell != pam.selected && boundary != pam.selected );
            cell = pam.selected;
        }
       
        if ( transparent )
        {
            // this is an internal triangle, store in 'tris'
            GLfloat g[] = { a[0]+b[0]+c[0], a[1]+b[1]+c[1], a[2]+b[2]+c[2] };
            tris[n_tris].a = faces[3*n  ];
            tris[n_tris].b = faces[3*n+1];
            tris[n_tris].c = faces[3*n+2];
            tris[n_tris].z = g[0] * ver[0] + g[1] * ver[1] + g[2] * ver[2];
            ++n_tris;
        }
        else
        {
            // set identical back and front material properties
            gle_color col = gle::bright_color(cell);
            glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, col.data());
            glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, col.data());

            // calculate the normal to the triangle:
            GLfloat ab[] = { b[0]-a[0], b[1]-a[1], b[2]-a[2] };
            GLfloat ac[] = { c[0]-a[0], c[1]-a[1], c[2]-a[2] };
            GLfloat ns[] = { ab[1]*ac[2]-ab[2]*ac[1], ab[2]*ac[0]-ab[0]*ac[2], ab[0]*ac[1]-ab[1]*ac[0] };
            
            glNormal3fv(ns);
            glDrawElements(GL_TRIANGLES, 3, GL_UNSIGNED_INT, faces+3*n);
        }
    }
    
    // prepare for transparency
    if ( pam.face_color.transparent() )
        glDepthMask(GL_FALSE);

    // set identical back and front material properties
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, pam.face_color.data());
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, pam.face_color.data());
    glMateriali (GL_FRONT_AND_BACK, GL_SHININESS, 64);
    
    // depth-sort triangles:
    qsort(tris, n_tris, sizeof(Triangle), &closer);
    
    // render transparent triangles
    for ( unsigned i = 0; i < n_tris; ++i )
    {
        float * a = points + 3*tris[i].a;
        float * b = points + 3*tris[i].b;
        float * c = points + 3*tris[i].c;
        
        GLfloat ab[] = { b[0]-a[0], b[1]-a[1], b[2]-a[2] };
        GLfloat ac[] = { c[0]-a[0], c[1]-a[1], c[2]-a[2] };
        GLfloat ns[] = { ab[1]*ac[2]-ab[2]*ac[1], ab[2]*ac[0]-ab[0]*ac[2], ab[0]*ac[1]-ab[1]*ac[0] };
        
        glNormal3fv(ns);
        glDrawElements(GL_TRIANGLES, 3, GL_UNSIGNED_INT, tris+i);
    }
    
    delete[] tris;
    glDepthMask(GL_TRUE);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glDisableClientState(GL_VERTEX_ARRAY);
}


/// return index of cell on which the mouse-click occurred
unsigned Mesh::pick() const
{
    const GLsizei buf_size = 1024;
    GLuint buf[buf_size] = { 0 };

    glEnableClientState(GL_VERTEX_ARRAY);
    glVertexPointer(3, GL_FLOAT, 0, points);

    glDisable(GL_CULL_FACE);

    glSelectBuffer(buf_size, buf);
    glRenderMode(GL_SELECT);
    glInitNames();
    glPushName(0);
    
    for ( unsigned n = 0; n < n_faces; ++n )
    {
        glLoadName(labels[2*n]);
        glDrawElements(GL_TRIANGLES, 3, GL_UNSIGNED_INT, faces+3*n);
    }
    
    glDisableClientState(GL_VERTEX_ARRAY);
    glPopName();
    
    GLint n_hits = glRenderMode(GL_RENDER);
    GLuint z_min = buf[1];
    GLuint hit = buf[3];
    for ( unsigned i=0; i < n_hits && 4*i < buf_size; ++i )
    {
        GLuint z = buf[4*i+1];
        if ( z < z_min )
        {
            z_min = z;
            hit = buf[4*i+3];
        }
    }
    return hit;
}

