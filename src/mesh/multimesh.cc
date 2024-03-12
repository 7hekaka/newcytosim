// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University

#include "assert_macro.h"
#include "multimesh.h"

#include "opengl.h"
#include "gym_flute.h"
#include "gym_color_list.h"
#include "gym_draw.h"
#include "gym_view.h"
#include "gym_cap.h"


MultiMesh::MultiMesh()
{
    n_points = 0;
    n_faces = 0;
    points = nullptr;
    faces = nullptr;
    labels = nullptr;
}


MultiMesh::~MultiMesh()
{
    release();
}


void MultiMesh::release()
{
    delete(points);
    delete(faces);
    delete(labels);
    points = nullptr;
    faces = nullptr;
    labels = nullptr;
}


int MultiMesh::read_ascii(FILE * file)
{
    char str[1024], * ptr;
    
    release();
    
    // read vertices:
    if ( !fgets(str, sizeof(str), file) )
        return 1;
    n_points = strtol(str, nullptr, 10);
    points = new real[3*n_points];
    
    for ( size_t n = 0; n < n_points; ++n )
    {
        if ( !fgets(str, sizeof(str), file) )
            return 1;
        points[3*n  ] = strtof(str, &ptr);
        points[3*n+1] = strtof(ptr, &ptr);
        points[3*n+2] = strtof(ptr, &ptr);
    }
    
    // read faces:
    if ( !fgets(str, sizeof(str), file) )
        return 1;
    n_faces = strtol(str, nullptr, 10);
    faces = new unsigned[3*n_faces];
    labels = new int[2*n_faces];

    for ( size_t n = 0; n < n_faces; ++n )
    {
        if ( !fgets(str, sizeof(str), file) )
            return 1;
        faces[3*n  ] = strtol(str, &ptr, 10);
        faces[3*n+1] = strtol(ptr, &ptr, 10);
        faces[3*n+2] = strtol(ptr, &ptr, 10);
        labels[2*n  ] = strtol(ptr, &ptr, 10);
        labels[2*n+1] = strtol(ptr, &ptr, 10);
    }

    return 0;
}


int MultiMesh::read_binary(FILE * file)
{
    double d[3];
    size_t s[3];
    int i[2];
    
    release();
    
    // read vertices:
    if ( 1 != fread(s, sizeof(size_t), 1, file) )
        return 1;
    
    n_points = s[0];
    points = new real[3*n_points];
    
    for ( size_t n = 0; n < n_points; ++n )
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
    
    for ( size_t n = 0; n < n_faces; ++n )
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


int MultiMesh::read(char const* filename)
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
            std::clog << "   error " << res << '\n';
        else
            std::clog << "   " << n_points << " points" << " and " << n_faces << " faces" << '\n';
*/
        return res;
    }
    std::clog << " Cannot read `" << filename << '\n';
    return 1;
}


/// OpenGL display function
void MultiMesh::drawPoints(float point_size, const float color[4]) const
{
    float * flu = gym::mapFloatBuffer(3*n_points);
    // Upload vertex data to the video device
    for ( size_t i = 0; i < 3*n_points; ++i )
        flu[i] = points[i];
    gym::unmapBufferV3();
    gym::color(color);
    gym::disableLighting();
    gym::drawPoints(point_size, 0, n_points);
    gym::enableLighting();
    gym::cleanupV();
}


/// simple OpenGL display function
void MultiMesh::drawFaces(const float color[4], int selected) const
{
    flute6 * flu = gym::mapBufferV3N3(3*n_faces);
    flute6 * ptr = flu;
    // Upload vertex data to the video device, and calculate normals
    for ( size_t i = 0; i < n_faces; ++i )
    {
        int cell = labels[2*i];
        int boundary = labels[2*i+1];
        bool transparent = ( boundary > 0 );
        // if one cell is selected, all the other ones are transparent:
        if ( selected )
        {
            transparent = ( cell != selected && boundary != selected );
            cell = selected;
        }
        real * ap = points + 3 * faces[3*i  ];
        real * bp = points + 3 * faces[3*i+1];
        real * cp = points + 3 * faces[3*i+2];
        flute3 a(ap), b(bp), c(cp);
        // calculate the normal to the triangle:
        flute3 ab = b - a, ac = c - a;
        flute3 nor = cross(ab, ac);
        *ptr++ = { a, nor };
        *ptr++ = { b, nor };
        *ptr++ = { c, nor };
    }
    gym::unmapBufferV3N3();
    gym::enableLighting();
    gym::disableCullFace();
    gym::openDepthMask();
    float black[4] = { 0, 0, 0, 1 };
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, color);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, black);
    glMateriali(GL_FRONT_AND_BACK, GL_SHININESS, 64);
    gym::drawTriangles(0, ptr-flu);
    gym::cleanupVN();
}



/// function to sort Triangles according to their 'z'
static int compareFluteTriangles(const void * A, const void * B)
{
    float a = static_cast<const float*>(A)[3];
    float b = static_cast<const float*>(B)[3];
    
    return ( a > b ) - ( b > a );
}

/// OpenGL display function with depth-sorted transparent faces
void MultiMesh::drawFaces(const float dir[3], const float color[4], int selected) const
{
    flute6 * flu = gym::mapBufferV3N3(3*n_faces);
    flute6 * ptr = flu;
    flute6 * end = flu + 3 * n_faces;
    
    // process all faces making separate triangles
    for ( size_t i = 0; i < n_faces; ++i )
    {
        int cell = labels[2*i];
        int boundary = labels[2*i+1];
        bool transparent = ( boundary > 0 );
        // if one cell is selected, all the other ones are transparent:
        if ( selected )
        {
            transparent = ( cell != selected && boundary != selected );
            cell = selected;
        }
        real * ap = points + 3 * faces[3*i  ];
        real * bp = points + 3 * faces[3*i+1];
        real * cp = points + 3 * faces[3*i+2];

        flute3 a(ap), b(bp), c(cp);
        if ( transparent )
        {
            end -= 3;
            float Z = dot(a + b + c, dir);
            end[0] = { a, Z, 0, 0 };
            end[1] = { b, 0, 0, 0 };
            end[2] = { c, 0, 0, 0 };
        }
        else
        {
            ptr[0] = { a, 0, 0, 0 };
            ptr[1] = { b, 0, 0, 0 };
            ptr[2] = { c, 0, 0, 0 };
            ptr += 3;
        }
    }
    assert_true(end == ptr);
    
    // number of opaque triangles:
    size_t n_tri = ( ptr - flu ) / 3;
    // depth-sort transparent triangles:
    qsort(ptr, n_faces-n_tri, 3*sizeof(flute6), &compareFluteTriangles);
    //printf("n_faces %lu  n_tri %lu\n", n_faces, n_tri);
    
    // calculate normals for all triangles
    end = flu + 3 * n_faces;
    for ( ptr = flu; ptr < end; ptr += 3 )
    {
        flute3& a(reinterpret_cast<flute3&>(ptr[0]));
        flute3& b(reinterpret_cast<flute3&>(ptr[1]));
        flute3& c(reinterpret_cast<flute3&>(ptr[2]));
        // calculate a normal vector to the triangle:
        flute3 ab = b - a, ac = c - a;
        flute3 nor = cross(ab, ac);
        ptr[0] = { a, nor };
        ptr[1] = { b, nor };
        ptr[2] = { c, nor };
    }

    gym::unmapBufferV3N3();
    gym::enableLighting();
    gym::enableCullFace(GL_BACK);
    gym::openDepthMask();

    // set identical back and front material properties
    //gym_color col = gym::bright_color(cell);
    float black[4] = { 0, 0, 0, 1 };
    float glass[4] = { 1, 1, 1, 0.3 };
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, color);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, black);
    glMateriali(GL_FRONT_AND_BACK, GL_SHININESS, 64);

    gym::drawTriangles(0, 3*n_tri);
    
    // render transparent triangles
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, glass);
    gym::disableCullFace();
    gym::closeDepthMask();
    gym::drawTriangles(3*n_tri, 3*(n_faces-n_tri));
    gym::openDepthMask();
    gym::cleanupVN();
}



/// return index of cell on which the mouse-click occurred
unsigned MultiMesh::pick() const
{
    float * flu = gym::mapFloatBuffer(3*n_points);
    // Upload vertex data to the video device
    for ( size_t i = 0; i < 3*n_points; ++i )
        flu[i] = points[i];
    gym::unmapBufferV3();
#if 0
    gym::clearPixels(0,0,0.25,1);
    gym::color(1,1,0);
    gym::enableCullFace(GL_BACK);
    for ( size_t n = 0; n < n_faces; ++n )
        glDrawElements(GL_TRIANGLES, 3, GL_UNSIGNED_INT, faces+3*n);
    gym::disableLighting();
    glFlush();
    return 0;
#endif
    const GLsizei buf_size = 1024;
    GLuint buf[buf_size] = { 0 };
    glSelectBuffer(buf_size, buf);
    glRenderMode(GL_SELECT);
    glInitNames();
    glPushName(0);
    for ( size_t n = 0; n < n_faces; ++n )
    {
        glLoadName(labels[2*n]);
        glDrawElements(GL_TRIANGLES, 3, GL_UNSIGNED_INT, faces+3*n);
    }
    glPopName();
    
    GLint n_hits = glRenderMode(GL_RENDER);
    GLuint z_min = buf[1];
    GLuint hit = buf[3];
    for ( GLint i=0; i < n_hits && 4*i < buf_size; ++i )
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

