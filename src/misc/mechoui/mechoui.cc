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
#include "mesh.h"


MechouiParam pam;

Mesh mesh;

bool animate = false;
unsigned file_index = 0;
std::vector<std::string> file_list;

//------------------------------------------------------------------------------

void load(unsigned index)
{
    if ( index < file_list.size() )
    {
        file_index = index;
        if ( pam.directory != "." )
            pam.file = pam.directory + "/" + file_list[index];
        else
            pam.file = file_list[index];
        mesh.read(pam.file.c_str());
        glApp::flashText("%s: %i points", pam.file.c_str(), mesh.nbPoints());
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
        case ' ':
            animate = ! animate;
            break;
        case 'c':
            pam.selected = ( pam.selected + 1 ) % 17;
            break;
        case 'p':
        {
            if ( pam.delay < 6 )
                pam.delay = 3;
            else
                pam.delay /= 2;
            glApp::flashText("delay %i ms", pam.delay);
        } break;
        case 'o':
            pam.delay = ( pam.delay < 1024 ) ? 2*pam.delay : 2048;
            glApp::flashText("delay %i ms", pam.delay);
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
        case 'P':
            pam.point_style = ! pam.point_style;
        default:
            glApp::processNormalKey(c,x,y);
            return;
    }
    glApp::postRedisplay();
}

/// callback for shift-click, with unprojected mouse position
void  processMouseClick(int mx, int my, const Vector3 & a, int)
{
    glApp::flashText("click %i %i", mx, my);
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);

    View& view = glApp::currentView();
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPickMatrix(mx, my, 2, 2, viewport);
    view.setProjection();
    pam.selected = mesh.pick();
    glMatrixMode(GL_MODELVIEW);
}

/// callback for shift-drag, with unprojected mouse positions
void  processMouseDrag(int, int, Vector3 & a, const Vector3 & b, int)
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
    pam.face_color = view.front_color;
    mesh.draw(pam);
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
    int ax = 1;
    
    if ( argc > 1 && FilePath::is_dir(argv[1]) )
    {
        FilePath::change_dir(argv[1]);
        ++ax;
    }
    
    if ( arg.read_strings(argc-ax, argv+ax) )
        return EXIT_FAILURE;

    if ( arg.use_key("help") )
    {
        help();
        help_keys();
        return EXIT_SUCCESS;
    }
    
    if ( arg.use_key("parameters") )
    {
        pam.write(std::cout);
        return EXIT_SUCCESS;
    }
    
    glutInit(&argc, argv);
    glApp::setDimensionality(3);
    glApp::normalKeyFunc(processNormalKey);
    glApp::actionFunc(processMouseClick);
    glApp::actionFunc(processMouseDrag);

    if ( pam.config.size() )
        arg.read_file(pam.config);
    
    pam.read(arg);
    glApp::currentView().read(arg);
    arg.print_warnings(std::cerr, 1, "\n");
    glApp::newWindow(display);
    glApp::attachMenu();
    glApp::setScale(2);
    gle::initialize();
    
    file_list = FilePath::list_dir(pam.directory.c_str(), "rec");
    load(file_index);
    timer(0);
    
    glutMainLoop();
}
