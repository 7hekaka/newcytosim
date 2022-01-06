// Cytosim was created by Francois Nedelec. Copyright 2007-2022 Cambridge University.

/*
 This is a test for glApp
 mouse driven zoom and rotation with Quaternions and GLU unproject
 Francois Nedelec nedelec@embl.de,  Oct. 2002, modified Jan. 2006
*/

#include "glossary.h"
#include "glapp.h"
#include "glut.h"
#include "gle.h"

Vector3 origin(0,0,0), position(0,0,0);


void display(View& view, int)
{
    view.openDisplay();
    glEnable(GL_LIGHTING);
    glDepthMask(GL_TRUE);

    // icosahedron:
    gle_color(1.0, 0.0, 1.0).load_front();
    gle_color(0.0, 0.0, 0.1).load_back();
    gle::icosahedron();
    //gle::ICOSAHEDRON();
    //gle::icoid();
    if ( 0 )
    {
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        gle_color(1, 1, 1).load_both();
        glDisable(GL_DEPTH_TEST);
        glLineWidth(0.5);
        gle::icoid();
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glEnable(GL_DEPTH_TEST);
    }
    
    if ( 0 )
    {
        const float rad = 0.1f;
        gle_color(1.0, 1.0, 1.0).load_both();
        glPushMatrix();
        gle::transScale(origin, rad);
        gle::sphere2();
        glPopMatrix();

        gle_color(0.0, 1.0, 0.0).load_both();
        glPushMatrix();
        gle::transScale(position, rad);
        gle::sphere2();
        glPopMatrix();
        
        gle_color(1.0, 0.0, 1.0).load_both();
        glPushMatrix();
        gle::scale(rad);
        gle::sphere2();
        glPopMatrix();
    }
    if ( 1 )
    {
        // transparent cube
        glDisable(GL_LIGHTING);
        glDepthMask(GL_FALSE);
        glLineWidth(1);
        glPushMatrix();
        gle::scale(1.732f);
        gle_color(1.0, 1.0, 1.0, 0.1).load();
        gle::cube();
        gle_color(1.0, 1.0, 1.0).load();
        //gle::wireCube();
        glPopMatrix();
    }
    view.closeDisplay();
}


///set callback for shift-click, with unprojected click position
void processMouseClick(int, int, const Vector3 & a, int)
{
    origin = a;
    glApp::postRedisplay();
}

///set callback for shift-drag, with unprojected mouse and click positions 
void processMouseDrag(int, int, Vector3 & a, const Vector3 & b, int)
{
    origin = a;
    position = b;
    glApp::postRedisplay();
}


int main(int argc, char* argv[])
{
    glutInit(&argc, argv);
    glApp::setDimensionality(3);
    glApp::actionFunc(processMouseClick);
    glApp::actionFunc(processMouseDrag);
    glApp::attachMenu(GLUT_RIGHT_BUTTON);
    glApp::newWindow(display);
    glApp::setScale(4);
    gle::initialize();
    glutMainLoop();
}
