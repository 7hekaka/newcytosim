// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
// Test class Platonic::Solid

#include "gle.h"
#include "glut.h"
#include "glapp.h"
#include "platonic.h"
#include <cstdio>

using namespace gle;

int style = 2;
int kind = 2;
int rank = 1;

bool showVertices = false;
bool showEdges = false;
bool showFaces = true;

Platonic::Solid * ico = nullptr;

GLuint glpts = 0;
GLuint gldir = 0;

//------------------------------------------------------------------------------
void initVBO();

void setPlatonic()
{
    if ( ico )
        delete ico;
    
    ico = new Platonic::Solid((Platonic::Solid::Polyhedra)kind, rank, 1);

    char tmp[128];
    snprintf(tmp, sizeof(tmp), "%i div, %i points", rank, ico->nb_vertices());
    glApp::setMessage(tmp);
    
    initVBO();
}

void processNormalKey(unsigned char c, int x, int y)
{
    switch (c) 
    {
        case 'k': kind = ( kind + 1 ) % 4;      break;
        case ']': rank += 1;                    break;
        case '[': rank = std::max(rank-1, 1);   break;
        case '}': rank += 16;                   break;
        case '{': rank = std::max(rank-16, 1);  break;
        case 's': style = ( style + 1 ) % 4;    break;
        case 'e': showEdges = !showEdges;       break;
        case 't': showFaces = !showFaces;       break;
        case 'p': showVertices = !showVertices; break;
        case ' ': break; // update the Platonic
        default: glApp::processNormalKey(c,x,y); return;
    }
    
    setPlatonic();
    
    glApp::flashText("style = %i", style);
    glutPostRedisplay();
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

void drawFaces1()
{
    glBegin(GL_TRIANGLES);
    for ( unsigned i = 0; i < ico->nb_faces(); ++i )
    {
        glNormal3fv(ico->face_data0(i));
        glVertex3fv(ico->face_data0(i));
        glNormal3fv(ico->face_data1(i));
        glVertex3fv(ico->face_data1(i));
        glNormal3fv(ico->face_data2(i));
        glVertex3fv(ico->face_data2(i));
    }
    glEnd();
}


void drawFaces2()
{
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_NORMAL_ARRAY);
    glVertexPointer(3, GL_FLOAT, 0, ico->vertex_data());
    glNormalPointer(GL_FLOAT, 0, ico->vertex_data());
    glDrawElements(GL_TRIANGLES, 3*ico->nb_faces(), GL_UNSIGNED_INT, ico->faces_data());
    glDisableClientState(GL_NORMAL_ARRAY);
    glDisableClientState(GL_VERTEX_ARRAY);
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
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, 3*ico->nb_faces()*sizeof(unsigned), ico->faces_data(), GL_STATIC_DRAW);
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
    glDisableClientState(GL_VERTEX_ARRAY);
}


void drawFaces()
{
    if ( style == 1 )
        drawFaces1();
    else if ( style == 2 )
        drawFaces2();
    else if ( style == 3 )
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

void drawVertices()
{
    char tmp[128];
    for ( unsigned i=0; i < ico->nb_vertices(); ++i )
    {
        Platonic::Vertex & dv = ico->vertex(i);
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

void display(View&, int)
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    
    drawPlane();
    if ( showFaces )
    {
        glColor3f(0, 0, 0.75f);
        glEnable(GL_LIGHTING);
        glEnable(GL_CULL_FACE);
        glCullFace(GL_BACK);
        drawFaces();
        glDisable(GL_CULL_FACE);
        glDisable(GL_LIGHTING);
    }
    if ( showEdges )
    {
        if ( ico->nb_edges() == 0 )
            ico->setEdges();
        drawEdges();
    }
    if ( showVertices )
    {
        drawVertices();
    }
    glutReportErrors();
}


//------------------------------------------------------------------------------

int main(int argc, char* argv[])
{
    glutInit(&argc, argv);
    glApp::setDimensionality(3);
    glApp::attachMenu(GLUT_RIGHT_BUTTON);
    glApp::normalKeyFunc(processNormalKey);
    glApp::createWindow(display);
    glApp::setScale(3);
  
    setPlatonic();
    glutMainLoop();
}

/*
int main(int argc, char* argv[])
 {
     Platonic::Solid::Polyhedra kind = Platonic::Solid::ICOSAHEDRON;
     printf("-------------- order 2:\n");
     Platonic::Solid T0(kind, 2);
     printf("-------------- order 3:\n");
     Platonic::Solid T1(kind, 2);
     printf("-------------- order 4:\n");
     Platonic::Solid T2(kind, 4);
     printf("-------------- order 8:\n");
     Platonic::Solid T3(kind, 8);
     printf("-------------- order 16:\n");
     Platonic::Solid T4(kind, 16);
     printf("-------------- order 32:\n");
     Platonic::Solid T5(kind, 32);
     printf("-------------- order 64:\n");
     Platonic::Solid T6(kind, 64);
}
*/

