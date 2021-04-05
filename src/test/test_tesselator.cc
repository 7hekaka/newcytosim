// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include "gle.h"
#include "glut.h"
#include "glapp.h"
#include "tesselator.h"
#include <cstdio>

using namespace gle;

int style = 0;
int kind = 2;
int rank = 1;

bool showPlane = true;
bool showIndices = false;
bool showVertices = false;
bool showEdges = false;
bool showFaces = true;

Tesselator * ico = nullptr;

GLuint glpts = 0;
GLuint gldir = 0;

//------------------------------------------------------------------------------
void initVBO();

void reset()
{
    if ( ico )
        delete ico;
    ico = new Tesselator((Tesselator::Polyhedra)kind, rank, 1);

    char tmp[128];
    snprintf(tmp, sizeof(tmp), "%i div, %i points, %i faces",
             rank, ico->nb_vertices(), ico->nb_faces());
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

void processNormalKey(unsigned char c, int x, int y)
{
    switch (c) 
    {
        case ' ': break; // update the Platonic
        case 'k': kind = ( kind + 1 ) % 6; reset(); break;
        case ']': rank += 1; reset(); break;
        case '}': rank += 16; reset(); break;
        case '[': rank = std::max(rank-1, 1); reset(); break;
        case '{': rank = std::max(rank-16, 1); reset(); break;
        case 'y': exportPLY(); return;
        case 'Y': exportSTL(); return;
        case 'e': showEdges = !showEdges; break;
        case 't': showFaces = !showFaces; break;
        case 'i': showIndices = !showIndices; break;
        case 'p': showVertices = !showVertices; break;
        case 's': style = (style+1) % 2; glApp::flashText("style = %i", style); break;
        default: glApp::processNormalKey(c,x,y); return;
    }
    glApp::postRedisplay();
}

//------------------------------------------------------------------------------
void drawPlane()
{
    glColor3f(0.25f, 0.25f, 0.25f);
    glEnableClientState(GL_VERTEX_ARRAY);
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
    glDrawElements(GL_TRIANGLES, 3*ico->nb_faces(), GL_UNSIGNED_INT, ico->face_data());
    glDisableClientState(GL_NORMAL_ARRAY);
    //glDisableClientState(GL_VERTEX_ARRAY);
}


void initVBO()
{
    glGenBuffers(1, &glpts);
    glGenBuffers(1, &gldir);
    //Create a new VBO for the vertex information
#if 0
    glBindBuffer(GL_ARRAY_BUFFER, glpts);
    glBufferData(GL_ARRAY_BUFFER, 3*ico->nb_vertices()*sizeof(float), ico->vertex_data(), GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
#else
    glBindBuffer(GL_ARRAY_BUFFER, glpts);
    glBufferData(GL_ARRAY_BUFFER, 3*ico->nb_vertices()*sizeof(float), nullptr, GL_STATIC_DRAW);
    void * glb = glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
    ico->store_vertices((float*)glb);
    glUnmapBuffer(GL_ARRAY_BUFFER);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
#endif
    //Create a new VBO for the indices
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gldir);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, 3*ico->nb_faces()*sizeof(unsigned), ico->face_data(), GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}


void drawFacesVBO()
{
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);

    glBindBuffer(GL_ARRAY_BUFFER, glpts);
    glVertexPointer(3, GL_FLOAT, 0, nullptr);
    glNormalPointer(GL_FLOAT, 0, nullptr);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, gldir);
    glDrawElements(GL_TRIANGLES, 3*ico->nb_faces(), GL_UNSIGNED_INT, nullptr);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    
    glDisableClientState(GL_NORMAL_ARRAY);
    //glDisableClientState(GL_VERTEX_ARRAY);
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
    glColor3f(1,1,1);
    glLineWidth(0.5);
    glBegin(GL_LINES);
    for ( unsigned i = 0; i < ico->nb_edges(); ++i )
    {
        glVertex3fv(ico->edge_vertex0(i));
        glVertex3fv(ico->edge_vertex1(i));
    }
    glEnd();
}

void nameVertices()
{
    char tmp[128];
    for ( unsigned i=0; i < ico->nb_vertices(); ++i )
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
        gle::drawText(pos, tmp, GLUT_BITMAP_8_BY_13, 0.5);
    }
}

void drawVertices()
{
    glPointSize(10);
    glColor3f(1, 1, 1);
    glEnableClientState(GL_VERTEX_ARRAY);
    glBindBuffer(GL_ARRAY_BUFFER, glpts);
    glVertexPointer(3, GL_FLOAT, 0, nullptr);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glDrawArrays(GL_POINTS, 0, ico->nb_vertices());
    glDisableClientState(GL_VERTEX_ARRAY);
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
        drawPlane();
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
#if 1
        glDisable(GL_LIGHTING);
        glLineWidth(0.25);
        glColor3f(1, 1, 1);
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        drawFaces();
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
#else
        if ( ico->nb_edges() == 0 )
            ico->setEdges();
        drawEdges();
#endif
    }
    if ( showVertices )
    {
        glDisable(GL_LIGHTING);
        drawVertices();
    }
    if ( showIndices )
    {
        glDisable(GL_LIGHTING);
        nameVertices();
    }
    view.closeDisplay();
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
    glutMainLoop();
}
