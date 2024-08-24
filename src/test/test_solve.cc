// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
// This was created by F. Nedelec on April 2003

#include <cstdlib>
#include <cstdio>
#include "glut.h"
#include "gym_matrix.h"

#include "vector2.h"
#include "matrix22.h"
#include "random.h"


Vector2 pg(5,0,0), px(1,0,0), pf, pn, pc, pp;

//rate of calculation:
int delay = 100;
int normalization = 1;
int correct = 1;
int mobile = 1;

real radius = 1;
real time_step = 1./16;
real stiff = 10;
real noise = 1;

void timerFunction(int value)
{
    real s, h = stiff * time_step;
    Matrix22 I, P, D, C;
    I = Matrix22::one();
    real xn = px.norm();
    Vector2 pxn = px.normalized();
    
    P(0,0) = 1 - pxn[0] * pxn[0];
    P(1,0) =   - pxn[1] * pxn[0];
    P(0,1) =   - pxn[0] * pxn[1];
    P(1,1) = 1 - pxn[1] * pxn[1];
    
    Vector2 rhs, f, fx;
    Vector2 random(RNG.gauss()*noise/stiff, RNG.gauss()*noise/stiff, 0);
    
    
    printf("\nstiff*time_step %2f noise %f correct %i norm %i :\n",
           stiff*time_step, noise, correct, normalization);
    printf("px    : "); px.println();
    //P.println();
    
    //free point obtained without the constraints:
    pf = ( px + h * ( pg + random ) ) / ( 1.0 + h );
    printf("free  : "); pf.println();
    
    //projection of the free point on the constraints
    pp = radius * pf.normalized();
    printf("proj  : "); pp.println();
    
    //adding the projector in the dynamic matrix
    rhs = px + h * ( P * ( pg + random ));
    pc = (I + h * P).inverted() * rhs;
    printf("const : "); pc.println();
    
    f = ( pg - px );
    s = dot(pxn, f) / xn;
    fx = 2 * s * pxn - f / xn;
    
    //projector with its corrections due to the derivative of the constraints:
    C(0,0) = pxn[0] * fx[0] - s;
    C(1,0) = pxn[1] * fx[0];
    C(0,1) = pxn[0] * fx[1];
    C(1,1) = pxn[1] * fx[1] - s;
    
    //C.println();
    rhs = px + h * ( P * ( pg + random - C * px ));
    D = I + h * ( P * ( I - C ));
    //D.println();
    //D.inverse();
    
    pn = D.inverted() * rhs;
    printf("new   : ");   pn.println();
    
    if ( mobile )        px = correct ? pn : pc;
    if ( normalization ) px = px.normalized(radius);
    
    glutPostRedisplay();
    glutTimerFunc(delay, timerFunction, 1);
}


void display()
{
    glClear(GL_COLOR_BUFFER_BIT);
    
    glColor3f(0,0,1.0);
    glBegin(GL_LINE_LOOP);
    for ( real ii = 0 ; ii < 6.28; ii+=0.0314 )
        glVertex2d( radius*std::cos(ii), radius*std::sin(ii) );
    glEnd();
    
    glPointSize(7.0);
    glBegin(GL_POINTS);
    
    glColor3f(1.f, 1.f, 1.f);
    glVertex2d(px[0], px[1]);
    
    glColor3f(0.f, 0.f, 1.f);
    glVertex2d(pg[0], pg[1]);
    
    glColor3f(1.f, 0.f, 0.f);
    glVertex2d(pf[0], pf[1]);
    
    glColor3f(0.f, 1.f, 0.f);
    glVertex2d(pc[0], pc[1]);
    
    glColor3f(0.f, 0.f, 1.f);
    glVertex2d(pn[0], pn[1]);
    
    glColor3f(1.f, 0.f, 1.f);
    glVertex2d(pp[0], pp[1]);
    
    glEnd();
    
    glColor3f(1.0,1.0,1.0);
    glBegin(GL_LINE_LOOP);
    glVertex2d(px[0], px[1]);
    glVertex2d(pg[0], pg[1]);
    glEnd();
    
    glFlush();
}


//------------------------------------------------------------------------------

//size of viewing box:
GLdouble vsize[] = { 10, 10 };

//the zoom factor:
real zoomSaved, zoom = 1;
int mouseAction, mouseX, mouseY;
GLint viewport[4];


//----------------matrices to compute the inverse projection of mouse locations

#define MOUSE_ZOOM    GLUT_RIGHT_BUTTON
#define MOUSE_SET     GLUT_LEFT_BUTTON
#define MENU_BUTTON   GLUT_MIDDLE_BUTTON

void setModelView()
{
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glScalef(zoom, zoom, zoom);
    glutPostRedisplay();
}


void processNormalKey(unsigned char c, int x, int y)
{
    switch (c) {
        case 27:
        case 'q':
            exit(EXIT_SUCCESS);
            
        case 'o':
            delay *= 2; break;
        case 'p':
            if ( delay > 1 ) delay /= 2;
            break;
        case 'K':
            stiff *= 2; break;
        case 'J':
            stiff /= 2; break;
        case 'k':
            stiff += 1; break;
        case 'j':
            stiff -= 1; break;
            
        case 'i':
            noise *= 2; break;
        case 'u':
            noise /= 2; break;
            
        case 'n':
            normalization = !normalization;
            break;
            
        case 'c':
            correct = !correct;
            break;
            
        case 'm':
            mobile = !mobile;
            break;
            
        case 'z':
            px.set( 1, 0, 0);
            pg.set( 5, 0, 0 );
            break;
            
        default:
            printf("normal key %c %i %i\n", c, x, y);
    }
}

enum MENUS_ID { MENU_QUIT };

void processMenu(int item)
{
    if ( item == MENU_QUIT )
        exit(EXIT_SUCCESS);
}

void windowReshaped(int w, int h)
{
    glViewport(0, 0, w, h);
    glGetIntegerv(GL_VIEWPORT, viewport);
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    double ratio = w * vsize[1] / double( vsize[0] * h );
    
    if ( ratio > 1 )
        glOrtho(-vsize[0], vsize[0], -vsize[1]/ratio, vsize[1]/ratio, 0, 1);
    else
        glOrtho(-vsize[0]*ratio, vsize[0]*ratio, -vsize[1], vsize[1], 0, 1);
}

//------------------------------------------------------------------------------

void unproject(float x, float y)
{
    GLfloat mat_model[16];
    GLfloat mat_proj[16];
    float pt[4] = { x, viewport[3]-y, 0, 1 };
    glGetFloatv(GL_MODELVIEW_MATRIX, mat_model);
    glGetFloatv(GL_PROJECTION_MATRIX, mat_proj);
    gym::unproject(pt, mat_model, mat_proj, viewport);
    pg.set(pt[0], pt[1], pt[2]);
    glutPostRedisplay();
}


void processMouse(int button, int state, int x, int y)
{
    //  printf("button %i %i %i %i\n", button, state, x, y);
    if ( state != GLUT_DOWN ) return;
    mouseAction = button;
    
    mouseX = x;
    mouseY = y;
    
    switch( mouseAction )
    {
        case MOUSE_ZOOM:
            zoomSaved = zoom;
            break;
            
        case MOUSE_SET:
            unproject(x, y);
            break;
    }
}


void processMotion(int x, int y)
{
    real d;
    switch( mouseAction )
    {
        case MOUSE_ZOOM:
            
            d = 1 + 4 * real( x - mouseX ) / viewport[2];
            if ( d > 0 ) zoom = zoomSaved * d;
                setModelView();
            break;
            
        case MOUSE_SET:
            unproject(x, y);
            break;
    }
}


void initGLUT()
{
    glClearColor(0.f, 0.f, 0.f, 0.f);
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_LINE_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glHint(GL_POINT_SMOOTH_HINT, GL_FASTEST);
    glHint(GL_LINE_SMOOTH_HINT, GL_FASTEST);
    
    glutCreateMenu(processMenu);
    glutAddMenuEntry("Quit", MENU_QUIT);
    glutAttachMenu(MENU_BUTTON);
    
    setModelView();
    glutTimerFunc(100, timerFunction, 1);
}

int main(int argc, char* argv[])
{
    RNG.seed();
    glutInit(&argc, argv);
    
    glutInitDisplayMode( GLUT_SINGLE | GLUT_RGBA );
    glutInitWindowSize(400, 400);
    glutInitWindowPosition(50, 50);
    glutCreateWindow(argv[0]);
    
    initGLUT();
    
    glutDisplayFunc(display);
    glutReshapeFunc(windowReshaped);
    glutMouseFunc(processMouse);
    glutMotionFunc(processMotion);
    glutKeyboardFunc(processNormalKey);
    
    glutMainLoop();
}
