// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include "gle.h"
#include "gym_text.h"
#include "gym_flute.h"
#include "glut.h"
#include "glapp.h"
#include "tesselator.h"
#include <cstdio>

int style = 0;
int kind = 2;
int rank = 1;

bool showPlane = true;
bool showNames = false;
bool showPoints = false;
bool showEdges = false;
bool showFaces = true;

Tesselator * ico = nullptr;

GLuint buffers[2] = { 0 };

//------------------------------------------------------------------------------
void initVBO();

void reset()
{
    if ( ico )
        delete ico;
    ico = new Tesselator();
    ico->build((Tesselator::Polyhedra)kind, rank, 1);
    ico->setVertices();
    
    char tmp[128];
    snprintf(tmp, sizeof(tmp), "%i div, %i points, %i faces",
             rank, ico->num_vertices(), ico->num_faces());
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
    glColor3f(0.25f, 0.25f, 0.25f);
    GLfloat pts[8] = {1, 1,-1, 1, 1,-1,-1,-1};
    glVertexPointer(2, GL_FLOAT, 0, pts);
    glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
}


void drawFacesArray()
{
    glEnableClientState(GL_NORMAL_ARRAY);
    glVertexPointer(3, GL_FLOAT, 0, ico->vertex_data());
    glNormalPointer(GL_FLOAT, 0, ico->vertex_data());
    glDrawElements(GL_TRIANGLES, 3*ico->num_faces(), GL_UNSIGNED_INT, ico->face_data());
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
    // calculate vertex data into device memory
    glBindBuffer(GL_ARRAY_BUFFER, buffers[0]);
    glBufferData(GL_ARRAY_BUFFER, 3*ico->num_vertices()*sizeof(float), nullptr, GL_STATIC_DRAW);
    void * glb = glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
    ico->store_vertices((float*)glb);
    glUnmapBuffer(GL_ARRAY_BUFFER);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
#endif
    // create a new VBO for vertex indices defining the triangles
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffers[1]);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, 3*ico->num_faces()*sizeof(unsigned), ico->face_data(), GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}


void drawFacesVBO()
{
    glEnableClientState(GL_NORMAL_ARRAY);

    glBindBuffer(GL_ARRAY_BUFFER, buffers[0]);
    glVertexPointer(3, GL_FLOAT, 0, nullptr);
    glNormalPointer(GL_FLOAT, 0, nullptr);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buffers[1]);
    glDrawElements(GL_TRIANGLES, 3*ico->num_faces(), GL_UNSIGNED_INT, nullptr);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    
    glDisableClientState(GL_NORMAL_ARRAY);
}


void drawFaces()
{
    if ( style == 1 )
        drawFacesArray();
    else
        drawFacesVBO();
}

void drawEdges()
{
    glVertexPointer(3, GL_FLOAT, 0, ico->vertex_data());
    glDrawElements(GL_LINES, 2*ico->num_edges(), GL_UNSIGNED_INT, ico->edge_data());
}

void drawNames()
{
    char tmp[128];
    for ( unsigned i=0; i < ico->num_vertices(); ++i )
    {
        Tesselator::Vertex & dv = ico->vertex(i);
        if ( dv.weight(2) == 0  &&  dv.weight(1) == 0 )
            glColor3f(1.f, 1.f, 1.f);
        else if ( dv.weight(2) == 0 )
            glColor3f(0.f, 1.f, 0.f);
        else
            glColor3f(.5f, .5f, .5f);
        
        const float* v = ico->vertex_data(i);
        Vector3 pos(v[0], v[1], v[2]);
        snprintf(tmp, sizeof(tmp), "%u", i);
        gym::drawText(pos, tmp, GLUT_BITMAP_8_BY_13, 0.5);
    }
}

void drawPoints()
{
    glPointSize(10);
    glColor3f(1, 1, 1);
    glBindBuffer(GL_ARRAY_BUFFER, buffers[0]);
    glVertexPointer(3, GL_FLOAT, 0, nullptr);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glDrawArrays(GL_POINTS, 0, ico->num_vertices());
}

void display(View& view, int)
{
    view.openDisplay();
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    
    if ( 0 )
    {
        glLineWidth(1);
        glPointSize(10);
        //glPolygonMode(GL_FRONT, GL_LINE);
        glColor4f(0, 1, 1, 0.5f);
        glEnable(GL_LIGHTING);
        //glEnable(GL_CULL_FACE);
        glCullFace(GL_BACK);
        //gle::sphere1();
        gle::needle();
        glDisable(GL_CULL_FACE);
        glDisable(GL_LIGHTING);
        return;
    }

    if ( showPlane )
    {
        glEnable(GL_LIGHTING);
        drawPlane();
    }
    if ( showFaces )
    {
        glColor3f(0, 0, 0.75f);
        glEnable(GL_LIGHTING);
        glEnable(GL_CULL_FACE);
        glCullFace(GL_BACK);
        drawFaces();
        glDisable(GL_CULL_FACE);
    }
    if ( showEdges )
    {
        glDisable(GL_LIGHTING);
        glLineWidth(0.5);
        glColor3f(1, 1, 1);
#if 1
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        drawFaces();
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
#else
        if ( ico->num_edges() == 0 )
            ico->setEdges();
        drawEdges();
#endif
    }
    if ( showPoints )
    {
        glDisable(GL_LIGHTING);
        drawPoints();
    }
    if ( showNames )
    {
        glDisable(GL_LIGHTING);
        drawNames();
    }
    view.closeDisplay();
}

//------------------------------------------------------------------------------

void processNormalKey(unsigned char c, int x, int y)
{
    switch (c)
    {
        case ' ': break; // update the Platonic
        case 'k': kind = ( kind + 1 ) % 6; reset(); break;
        case 'i': kind = Tesselator::ICOSAHEDRON; reset(); break;
        case 'o': kind = Tesselator::OCTAHEDRON; reset(); break;
        case 'd': kind = Tesselator::DICE; reset(); break;
        case ']': rank += 1; reset(); break;
        case '}': rank += 16; reset(); break;
        case '[': rank = std::max(rank-1, 1); reset(); break;
        case '{': rank = std::max(rank-16, 1); reset(); break;
        case 'y': exportPLY(); return;
        case 'Y': exportSTL(); return;
        case 'e': showEdges = !showEdges; break;
        case 'f': showFaces = !showFaces; break;
        case 'n': showNames = !showNames; break;
        case 'p': showPoints = !showPoints; break;
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
    glApp::attachMenu(GLUT_RIGHT_BUTTON);
    glApp::normalKeyFunc(processNormalKey);
    glApp::newWindow(display);
    glApp::setScale(3);
    gle::initialize();
    reset();
    glEnableClientState(GL_VERTEX_ARRAY);
    glutMainLoop();
}
