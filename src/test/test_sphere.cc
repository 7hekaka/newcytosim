// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
// Test for SphericalCode

#include <pthread.h>

#include "spherical_code.h"
#include "glapp.h"
#include "glut.h"
#include "gle.h"
#include "gym_text.h"
#include "gle_flute.h"

using namespace gle;

SphericalCode S, T;
SphericalCode * front = &S;

pthread_t       thread;
pthread_mutex_t lock;

int n_points = 12;

//------------------------------------------------------------------------------

void batch(unsigned long nbp, int repeat)
{
    size_t iterations = 0;
    real energy = INFINITY;
    real distance = 0;
    
    printf("%4lu pts :", nbp);
    printf(" %6.4f :", S.expectedDistance(nbp));
    
    for ( int m=0; m < repeat; ++m )
    {
        iterations += S.distributePoints(nbp, 1e-4, 1<<14);
        
        if ( S.finalEnergy() < energy )
            energy = S.finalEnergy();
        
        if ( S.minimumDistance() > distance )
            distance = S.minimumDistance();
    }
    
    printf("   distance %9.6f",    distance);
    printf("   energy %14.5f",     energy);
    printf("   iterations %12lu\n", iterations);
}

//------------------------------------------------------------------------------

void* calculateSphere(void * arg)
{
    glApp::setMessage("Calculating...");
    glApp::postRedisplay();

    if ( front == &S )
    {
        T.distributePoints(n_points, 1e-4, 1<<14);
        front = &T;
    }
    else
    {
        S.distributePoints(n_points, 1e-4, 1<<14);
        front = &S;
    }
    
    glApp::setMessage("");
    glApp::postRedisplay();
        
    pthread_mutex_unlock(&lock);
    pthread_exit(0);
}

//------------------------------------------------------------------------------
void processNormalKey(unsigned char c, int x, int y)
{
    switch (c)
    {
        case 'r': n_points-=256;  break;
        case 't': n_points-=32;   break;
        case 'y': n_points+=1;    break;
        case 'u': n_points+=4;    break;
        case 'i': n_points+=16;   break;
        case 'o': n_points+=128;  break;
        case 'p': n_points+=1024; break;
        case 'q' : exit(1);
            
        default:
            glApp::processNormalKey(c,x,y);
            return;
    }
    if ( n_points < 1 )
        n_points = 1;
    
    if ( 0 == pthread_mutex_trylock(&lock) )
    {
        pthread_create(&thread, 0, &calculateSphere, (void *)1);
    }
    else
    {
        glApp::flashText("already calculating...");
    }
}

//------------------------------------------------------------------------------
void drawVertices()
{
    size_t cnt = front->nbPoints();
    flute3* flu = gle::mapBufferV3(4*front->nbPoints());
    for ( size_t i = 0; i < cnt; ++i )
        flu[i] = { Vector3(front->addr(i)) };
    glPointSize(5);
    gle::unmapBufferV3();
    glDrawArrays(GL_POINTS, 0, cnt);
}

void nameVertices()
{
    char tmp[128];
    for ( size_t i=0; i < front->nbPoints(); ++i )
    {
        GLfloat x, y, z;
        front->putPoint(&x, &y, &z, i);
        snprintf(tmp, sizeof(tmp), "%lu", i);
        gym::drawText(Vector3(x, y, z), tmp, GLUT_BITMAP_8_BY_13, 0.5);
    }
}


void display(View& view, int)
{
    view.openDisplay();
    glPointSize(7);
    
    if ( front == &T )
        glColor3f(0.f, 0.f, 1.0);
    else
        glColor3f(0.f, 0.7, 0.f);

    drawVertices();
    nameVertices();
    
#if ( 1 )
    const real E = 0.1;
    flute8* flu = gle::mapBufferC4V4(4*front->nbPoints());
    gle_color col(1,1,1), lor(0,0,0);
    size_t n = 0;
    for ( size_t i = 0; i < front->nbPoints(); ++i )
    {
        Vector3 b, c, a(front->addr(i));
        a.orthonormal(b, c);
        flu[n++] = { col, a };
        flu[n++] = { lor, a+E*b };
        flu[n++] = { col, a };
        flu[n++] = { lor, a+E*c };
    }
    glLineWidth(3);
    gle::unmapBufferC4V4();
    glDrawArrays(GL_LINES, 0, n);
#endif
    
    if ( 0 )
    {
        glColor4f(0.3, 0.3, 0.3, 0.5);
        glDepthMask(GL_FALSE);
        glutSolidSphere(0.98,30,30);
        glDepthMask(GL_TRUE);
    }
    view.closeDisplay();
}

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
    RNG.seed();
    
    if ( argc == 3 ) 
    {
        size_t min = strtoul(argv[1], 0, 10);
        size_t max = strtoul(argv[2], 0, 10);
        for ( size_t nbp = min; nbp < max; nbp += 7)
            batch(nbp, 16);
        return EXIT_SUCCESS;
    }
    
    if ( argc == 2 ) 
        n_points = strtoul(argv[1], 0, 10);
    
    pthread_mutex_init(&lock, 0);
    front->distributePoints(n_points, 1e-4, 1<<14);
    
    glutInit(&argc, argv);
    glApp::setDimensionality(3);
    glApp::attachMenu(GLUT_RIGHT_BUTTON);
    glApp::normalKeyFunc(processNormalKey);
    glApp::newWindow(display);
    glApp::setScale(3);
    gle::initialize();
    glutMainLoop();
}
