// Cytosim 3.0 -  Copyright Francois Nedelec et al.  EMBL 2015
//------------------------------------------------------------------------------
//      Basic Mesh display - October 28 2015
//

#include "gle.h"
#include "glut.h"
#include "glapp.h"
#include "glossary.h"
#include "filepath.h"
#include "mechoui_param.h"
#include "gym_draw.h"
#include "gym_view.h"
#include "multimesh.h"

MultiMesh mesh;
MechouiParam pam;

bool animate = false;
unsigned file_index = 0;
std::vector<std::string> file_list;

//------------------------------------------------------------------------------

void load(size_t index)
{
    //printf("%lu / %lu:\n", index, file_list.size());
    if ( index < file_list.size() )
    {
        file_index = index;
        std::string file = file_list[index];
        mesh.read(file.c_str());
        glApp::flashText("%s: %lu points", file.c_str(), mesh.nbPoints());
        glApp::postRedisplay();
    }
    else
        animate = false;
}


/// callback for keyboard
void processNormalKey(unsigned char c, int x, int y)
{
    switch (c)
    {
        case 9: //TAB
            pam.selected = ( pam.selected + 1 ) & 7;
            glApp::flashText("selected %i", pam.selected);
            break;
        case ' ':
            animate = ! animate;
            break;
        case 'c':
            pam.selected = ( pam.selected + 1 ) % 17;
            break;
        case 'p':
            pam.point_style = ! pam.point_style;
            glApp::flashText("point_style %i", pam.point_style);
            break;
        case 'f':
            pam.point_style = ! pam.point_style;
            glApp::flashText("point_style %i", pam.point_style);
            break;
        case ',':
        case '<':
            load(file_index-1);
            break;
        case '.':
        case '>':
            load(file_index+1);
            break;
        case 'z':
            glApp::resetView();
            glApp::flashText("Reset view");
            break;
        case 'Z':
            load(0);
            break;
        case 'R':
            pam.write(std::cout);
            break;
        default:
            glApp::processNormalKey(c,x,y);
            return;
    }
    glApp::postRedisplay();
}

/// callback for shift-click, with unprojected mouse position
void processMouseClick(int mx, int my, const Vector3 & a, int)
{
    View& view = glApp::currentView();
    view.loadView();
    view.setPickProjection(mx, my, 16, 16);
    pam.selected = mesh.pick();
    glApp::flashText("selected %i", pam.selected);
}

/// callback for shift-drag, with unprojected mouse positions
void processMouseDrag(int, int, Vector3 & a, const Vector3 & b, int)
{
    glApp::flashText("drag %.1f %.1f %.1f", b.XX, b.YY, b.ZZ);
    glApp::postRedisplay();
}

/// timer callback
void timer(int value)
{
    if ( animate )
        load(file_index+1);
    glutTimerFunc(pam.delay, timer, 1);
}


int display(View& view)
{
    view.openDisplay();
    if ( pam.point_style && mesh.nbPoints() > 0 )
    {
        mesh.drawPoints(pam.point_size, pam.point_color);
    }
    if ( pam.face_style && mesh.nbFaces() > 0 )
    {
        Vector3 V = view.depthAxis();
        float axis[4] = { float(V.XX), float(V.YY), float(V.ZZ), 0 };
        mesh.drawFaces(axis, pam.face_color, pam.selected);
    }
    view.closeDisplay();
    return 0;
}

//------------------------------------------------------------------------------

void help()
{
    printf("Simple mesh viewer by Francois Nedelec, Copyright EMBL 2015\n"
           "Command-line options:\n"
           "   parameters  display list of parameters\n"
           "   help        display this help\n"
           "   keys        display list of keyboard controls\n"
           "   size=###    set window size\n"
           "   P=###       set value of parameter `P'\n");
}

void help_keys()
{
    printf("Keyboard commands:\n"
           "   SPACE   reset view\n");
}

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
    Glossary arg;
    if ( arg.read_strings(argc-1, argv+1) )
        return EXIT_FAILURE;
    
    if ( arg.use_key("parameters") )
    {
        pam.write(std::cout);
        return EXIT_SUCCESS;
    }
    
    if ( arg.use_key("help") )
    {
        help();
        help_keys();
        return EXIT_SUCCESS;
    }
    
    if ( pam.config.size() )
        arg.read_file(pam.config);
    
    pam.read(arg);

    if ( arg.has_key("directory") )
        pam.path = arg.value("directory");
    
    file_list = FilePath::list_dir(pam.path.c_str(), "rec");

#if 0
    printf("%s:\n", pam.path.c_str());
    for ( std::string& f : file_list )
        printf("    %s\n", f.c_str());
#endif
    
    glApp::initialize();
    glApp::views[0].read(arg);
    arg.print_warnings(stderr, 1, "\n");

    glutInit(&argc, argv);
    glApp::setDimensionality(3);
    glApp::normalKeyFunc(processNormalKey);
    glApp::actionFunc(processMouseClick);
    glApp::actionFunc(processMouseDrag);
    glApp::newWindow(display);
    glApp::attachMenu();
    glApp::setScale(3);
    gle::initialize();

    timer(0);
    load(file_index);
    glutMainLoop();
}
