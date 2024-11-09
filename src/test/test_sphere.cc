// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
// Test for SphericalCode

#include <pthread.h>

#include "spherical_code.h"
#include "glapp.h"
#include "glut.h"
#include "gle.h"
#include "gym_flute.h"
#include "gym_draw.h"
#include "gym_view.h"
#include "gym_cap.h"
#include "gym_check.h"
#include "fg_stroke.h"


SphericalCode S, T;
SphericalCode * front = &S;

pthread_t thread;
pthread_mutex_t lock;

size_t n_points = 26;

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
    pthread_exit(nullptr);
}

//------------------------------------------------------------------------------
void processNormalKey(unsigned char c, int x, int y)
{
    switch (c)
    {
        case 'q' : exit(1);
        case 'r': n_points-=256; break;
        case 't': n_points-=32;  break;
        case 'y': n_points+=1;   break;
        case 'u': n_points+=4;   break;
        case 'i': n_points+=16;  break;
        case 'o': n_points+=64;  break;
        case 'p': n_points+=256; break;
        case 's': front = (front==&S?&T:&S); return;
        default: glApp::processNormalKey(c,x,y); return;
    }
    if ( n_points < 1 )
        n_points = 1;
    
    if ( 0 == pthread_mutex_trylock(&lock) )
    {
        pthread_create(&thread, nullptr, &calculateSphere, (void *)1);
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
    flute3* flu = gym::mapBufferV3(4*front->nbPoints());
    for ( size_t i = 0; i < cnt; ++i )
    {
        real const* ptr = front->addr(i);
        flu[i] = { (float)ptr[0], (float)ptr[1], (float)ptr[2] };
    }
    gym::unmapBufferV3();
    gym::drawPoints(5, 0, cnt);
}

void nameVertices(View& view)
{
    gym::disableLighting();
    gym::disableAlphaTest();
    gym::cancelRotation();
    char tmp[32] = { 0 };
    for ( size_t i=0; i < front->nbPoints(); ++i )
    {
        snprintf(tmp, sizeof(tmp), "%lu", i);
        real const* ptr = front->addr(i);
        view.strokeString(ptr[0], ptr[1], ptr[2], tmp);
    }
    gym::restoreAlphaTest();
    gym::restoreLighting();
    gym::load_ref();
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
    gym::cleanupCV();
}


int display(View& view)
{
    gym_color col(1,1,1), lor(0,0,0);
    view.openDisplay();
    glPointSize(7);
    
    if ( front == &T )
        col[1] = 0.5f;
    gym::color(col);
    drawVertices();
    gym::color(1,1,1);
    nameVertices(view);
    //drawTangents(col, lor);
    view.closeDisplay();
    return 0;
}

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
    RNG.seed();
    
    if ( argc == 3 ) 
    {
        size_t min = strtoul(argv[1], nullptr, 10);
        size_t max = strtoul(argv[2], nullptr, 10);
        for ( size_t nbp = min; nbp < max; nbp += 7)
            batch(nbp, 16);
        return EXIT_SUCCESS;
    }
    
    if ( argc == 2 ) 
        n_points = strtoul(argv[1], nullptr, 10);
    
    pthread_mutex_init(&lock, nullptr);
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
