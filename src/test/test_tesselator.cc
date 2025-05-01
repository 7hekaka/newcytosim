// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include "gle.h"
#include "gym_color.h"
#include "glut.h"
#include "glapp.h"
#include "gym_cap.h"
#include "gym_check.h"
#include "gym_flute.h"
#include "gym_draw.h"
#include "gym_view.h"
#include "fg_stroke.h"
#include "tesselator.h"
#include <cstdio>

int style = 0;
int kind = 2;
int rank = 1;

bool showPlane = false;
bool showNames = false;
bool showPoints = false;
bool showEdges = false;
bool showFaces = true;

Tesselator * ico = nullptr;

GLint cull_test = true;

GLuint buffers[2] = { 0 };

void flip_cap(GLenum cap)
{
    GLint i = glIsEnabled(cap);
    if ( i )
        glDisable(cap);
    else
        glEnable(cap);
    gym::printCaps("flip");
}

//------------------------------------------------------------------------------
void initVBO();

void reset(int K, int R)
{
    kind = abs(K);
    rank = R;
    if ( ico )
        delete ico;
    ico = new Tesselator();
    ico->construct((Tesselator::Polyhedra)K, R, 2);

    char tmp[128];
    snprintf(tmp, sizeof(tmp), "%i div, %i points, %i faces",
             R, ico->num_vertices(), ico->num_faces());
    glApp::setMessage(tmp);
    initVBO();
    
}

FILE * openFile(const char name[])
{
    FILE * f = fopen(name, "w");

    if ( !f || ferror(f) )
    {
        glApp::flashText("input file could not be opened");
        return nullptr;
    }
    if ( ferror(f) )
    {
        fclose(f);
        glApp::flashText("input file opened with error");
        return nullptr;
    }
    return f;
}

void exportPLY()
{
    FILE * f = openFile("mesh.ply");
    if ( f ) {
        ico->exportPLY(f);
        fclose(f);
        glApp::flashText("exported `mesh.ply'");
    }
}

void exportSTL()
{
    FILE * f = openFile("mesh.stl");
    if ( f ) {
        ico->exportSTL(f);
        fclose(f);
        glApp::flashText("exported `mesh.stl'");
    }
}

//------------------------------------------------------------------------------

void drawPlane()
{
    glColor4f(0.5f, 0.5f, 0.5f, 1.0f);
    GLfloat pts[8] = {1, 1,-1, 1, 1,-1,-1,-1};
    glVertexPointer(2, GL_FLOAT, 0, pts);
    glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
}


void drawFacesArray()
{
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    glVertexPointer(3, GL_FLOAT, 0, ico->vertex_data());
    glNormalPointer(GL_FLOAT, 0, ico->vertex_data());
    static_assert(std::is_same<Tesselator::INDEX, GLushort>::value, "Index type mismatch");
    glDrawElements(GL_TRIANGLES, 3*ico->num_faces(), GL_UNSIGNED_SHORT, ico->face_data());
    glDisableClientState(GL_NORMAL_ARRAY);
}


void initVBO()
{
    glGenBuffers(2, buffers);
#if 0
    // copy vertex data to device memory
    glBindBuffer(GL_ARRAY_BUFFER, buffers[0]);
    glBufferData(GL_ARRAY_BUFFER, 3*ico->num_vertices()*sizeof(float), ico->vertex_data(), GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
#else
    // calculate vertex data directly into device memory
    glBindBuffer(GL_ARRAY_BUFFER, buffers[0]);
    glBufferData(GL_ARRAY_BUFFER, 3*ico->num_vertices()*sizeof(float), nullptr, GL_STATIC_DRAW);
    void * glb = glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
    ico->store_vertices((float*)glb);
    glUnmapBuffer(GL_ARRAY_BUFFER);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
#endif
    // create a new VBO for vertex indices defining the triangles
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffers[1]);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, 3*ico->num_faces()*sizeof(Tesselator::INDEX), ico->face_data(), GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}


void drawFacesVBO()
{
    glEnableClientState(GL_NORMAL_ARRAY);
    glEnableClientState(GL_VERTEX_ARRAY);

    glBindBuffer(GL_ARRAY_BUFFER, buffers[0]);
    glVertexPointer(3, GL_FLOAT, 0, nullptr);
    glNormalPointer(GL_FLOAT, 0, nullptr);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffers[1]);
    static_assert(std::is_same<Tesselator::INDEX, GLushort>::value, "Index type mismatch");
    glDrawElements(GL_TRIANGLES, 3*ico->num_faces(), GL_UNSIGNED_SHORT, nullptr);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    
    glDisableClientState(GL_NORMAL_ARRAY);
}


void drawEdges()
{
#if 0
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    drawFacesArray();
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
#else
    if ( ico->num_edges() == 0 )
        ico->setEdges();
    glVertexPointer(3, GL_FLOAT, 0, ico->vertex_data());
    glDrawElements(GL_LINES, 2*ico->num_edges(), GL_UNSIGNED_SHORT, ico->edge_data());
#endif
}

void namePoints(View& view)
{
    const float S = 1.03;
    gym::disableLighting();
    gym::disableAlphaTest();
    gym::cancelRotation();
    char tmp[128];
    for ( unsigned i=0; i < ico->num_vertices(); ++i )
    {
        float scale = 1.f;
        Tesselator::Vertex & dv = ico->vertex(i);
        float col[4] = {1.f, 0.f, 0.f, 1.f};
        if ( dv.weight(2) == 0 )
        { col[1] = 1; scale = 1.4142f; }
        if ( dv.weight(1) == 0 )
        { col[2] = 1; scale = 2.f; }
        
        gym::color(col);
        const float* ptr = ico->vertex_data(i);
        snprintf(tmp, sizeof(tmp), "%u", i);
        view.strokeString(S*ptr[0], S*ptr[1], S*ptr[2], tmp, scale);
    }
    gym::restoreAlphaTest();
    gym::restoreLighting();
    gym::load_ref();
}

void drawPoints()
{
    glColor3f(1, 1, 1);
    if ( style == 0 )
    {
        glBindBuffer(GL_ARRAY_BUFFER, buffers[0]);
        glVertexPointer(3, GL_FLOAT, 0, nullptr);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }
    else
    {
        glEnableClientState(GL_VERTEX_ARRAY);
        glVertexPointer(3, GL_FLOAT, 0, ico->vertex_data());
    }
    gym::drawPoints(10, 0, ico->num_vertices());
    glDisableClientState(GL_NORMAL_ARRAY);
}

void drawSkeleton()
{
    glLineWidth(1);
    glPointSize(10);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glColor4f(0, 1, 1, 0.5f);
    glEnable(GL_LIGHTING);
    //gle::sphere1();
    //gle::needle();
    gle::football();
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glDisable(GL_LIGHTING);
}

int display(View& view)
{
    view.openDisplay();
    glShadeModel(GL_FLAT);
    GLfloat blue[4] = { 0, 0, 1, 1 };
    GLfloat pink[4] = { 1.0, 0.0, 0.7, 1 };
    glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, blue);
    glMaterialfv(GL_BACK, GL_AMBIENT_AND_DIFFUSE, pink);

    //drawSkeleton();
    glDepthMask(GL_TRUE);
    if ( showPlane )
    {
        glDisable(GL_LIGHTING);
        drawPlane();
    }
    if ( showFaces )
    {
        glEnable(GL_LIGHTING);
        if ( cull_test )
            glEnable(GL_CULL_FACE);
        glCullFace(GL_BACK);
        if ( style == 1 )
            drawFacesArray();
        else
            drawFacesVBO();
        glDisable(GL_CULL_FACE);
    }
    if ( showEdges )
    {
        glDisable(GL_LIGHTING);
        glLineWidth(0.5);
        glColor3f(1, 1, 1);
        drawEdges();
    }
    if ( showPoints )
    {
        glDisable(GL_LIGHTING);
        drawPoints();
    }
    if ( showNames )
    {
        glDisable(GL_LIGHTING);
        gym::color(1, 1, 1);
        namePoints(view);
    }
    view.closeDisplay();
    return 0;
}

//------------------------------------------------------------------------------

void processNormalKey(unsigned char c, int x, int y)
{
    switch (c)
    {
        case ' ': reset(kind, rank); break;
        case 'i': reset(Tesselator::ICOSAHEDRON, rank); break;
        case 'I': reset(Tesselator::ICOSAHEDRONX, rank); break;
        case 'o': reset(Tesselator::OCTAHEDRON, rank); break;
        case 'd': reset(Tesselator::DICE, rank); break;
        case 'h': reset(Tesselator::HEMISPHERE, rank); break;
        case 'a': reset(Tesselator::DROPLET, rank); break;
        case 'c': reset(Tesselator::CYLINDER, rank); break;
        case 't': reset(Tesselator::TETRAHEDRON, rank); break;
        case ']': reset(kind, rank+1); break;
        case '}': reset(kind, rank+16); break;
        case '[': reset(kind, std::max(rank-1, 1)); break;
        case '{': reset(kind, std::max(rank-16, 1)); break;
        
        case 'y': exportPLY(); return;
        case 'Y': exportSTL(); return;
            
        case 'e': showEdges = !showEdges; break;
        case 'f': showFaces = !showFaces; break;
        case 'n': showNames = !showNames; break;
        case 'p': showPoints = !showPoints; break;
        case 'v': showPlane = !showPlane; break;
        case 'C': cull_test = !cull_test; break;
        case 'D': flip_cap(GL_DEPTH_TEST); break;

        case 's': style = (style+1) % 2; glApp::flashText("style = %i", style); break;
        default: glApp::processNormalKey(c,x,y); return;
    }
    glApp::postRedisplay();
}

//------------------------------------------------------------------------------

int main(int argc, char* argv[])
{
    glutInit(&argc, argv);
    glApp::setDimensionality(3);
    glApp::normalKeyFunc(processNormalKey);
    glApp::newWindow(display);
    glApp::attachMenu();
    glApp::setScale(3);
    gle::initialize();
    reset(kind, rank);
    glutMainLoop();
}
