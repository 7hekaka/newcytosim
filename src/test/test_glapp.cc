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
#include "gym_view.h"
#include "gym_draw.h"
#include "gym_vect.h"


Vector3 origin(0,0,0), position(0,0,0);

void drawObject()
{
    glEnable(GL_LIGHTING);
    gym::color_front(1, 0, 1, 1);
    gym::color_back(0, 0, 0.1, 1);
    gle::droplet();
    glLineWidth(0.5);
    glColor4f(0, 1, 1, 0.5f);
    glDisable(GL_LIGHTING);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    gle::droplet();
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

int display(View& view)
{
    view.openDisplay();
    gym::enableDepthTest();
    gym::enableLighting();
    gym::openDepthMask();

    drawObject();
    if ( 0 )
    {
        // icosahedron:
        gym::color_front(1, 0, 1, 0.5);
        gym::color_back(0, 0, 0.1, 0.5);
        gle::icosahedron();
        //gle::ICOSAHEDRON();
    }
    if ( 0 )
    {
        const float rad = 0.1f;
        gym::color_both(1, 1, 1, 1);
        gym::transScale(origin, rad);
        gle::sphere2();

        gym::color_both(0, 1, 0, 1);
        gym::transScale(position, rad);
        gle::sphere2();
        
        gym::color_both(1, 0, 1, 1);
        gym::transScale(0, 0, 0, rad);
        gle::sphere2();
    }
    if ( 0 )
    {
        gym::pull_ref();
        gym::scale(1.732f);
        // transparent cube
        gym::color_both(1, 1, 1, 0.5);
        gle::cubeEdges(1);
        gym::closeDepthMask();
        gym::enableLighting();
        gym::color_both(1, 1, 1, 0.1);
        gle::cube();
        gym::restoreLighting();
        gym::openDepthMask();
    }
    view.closeDisplay();
    return 0;
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
    glApp::newWindow(display);
    glApp::attachMenu();
    glApp::setScale(4);
    gle::initialize();
    glutMainLoop();
}
