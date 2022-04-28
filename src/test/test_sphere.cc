// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
// Test for SphericalCode

#include <pthread.h>

#include "spherical_code.h"
#include "glapp.h"
#include "glut.h"
#include "gle.h"
#include "gym_flute.h"
#include "gym_draw.h"
#include "gym_text.h"


SphericalCode S, T;
SphericalCode * front = &S;

pthread_t thread;
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
void drawVertices(const float col[4])
{
    size_t cnt = front->nbPoints();
    flute8* flu = gym::mapBufferC4V4(4*front->nbPoints());
    for ( size_t i = 0; i < cnt; ++i )
        flu[i] = { col, Vector3(front->addr(i)) };
    gym::unmapBufferC4V4();
    gym::drawPoints(5, 0, cnt);
}

void nameVertices(const float col[4])
{
    char tmp[128];
    for ( size_t i=0; i < front->nbPoints(); ++i )
    {
        real x, y, z;
        front->putPoint(&x, &y, &z, i);
        snprintf(tmp, sizeof(tmp), "%lu", i);
        gym::drawText(Vector3(x, y, z), BITMAP_8_BY_13, col, tmp, 0.5);
    }
}

void drawTangents(const float col[4], const float lor[4])
{
    const real E = 0.1;
    flute8* flu = gym::mapBufferC4V4(4*front->nbPoints());
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
    gym::unmapBufferC4V4();
    gym::drawLines(3, 0, n);
}


int display(View& view)
{
    gle_color col(1,1,1), lor(0,0,0);
    view.openDisplay();
    glPointSize(7);
    
    if ( front == &T )
        col[1] = 0.5f;
    drawVertices(col);
    nameVertices(col);
    drawTangents(col, lor);
    view.closeDisplay();
    return 0;
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
    glApp::normalKeyFunc(processNormalKey);
    glApp::newWindow(display);
    glApp::attachMenu();
    glApp::setScale(3);
    gle::initialize();
    glutMainLoop();
}
