// Cytosim was created by Francois Nedelec. Copyright Cambridge University 2022
// This is a visual test for the rasterizer used for attachment of Hands in Cytosim
// Created by FJN in October 2002, 10.12.2024

#include <ctime>
#include "gle.h"
#include "glut.h"
#include "glapp.h"
#include "real.h"
#include "random.h"
#include "rasterizer.h"
#include "vector2.h"
#include "vector3.h"

#include "gym_flute.h"
#include "gym_draw.h"
#include "gym_view.h"
#include "gym_cap.h"

// select between 2D or 3D mode:
#define FLAT_RASTERIZER 0


#if FLAT_RASTERIZER
typedef Vector2 Vector;
#else
typedef Vector3 Vector;
#endif


const int SIZE = 10;  // dimension of grid
const unsigned MAX = 64;  // max nb of points

unsigned n_pts = 2;
real radius = 2;

Vector shift(0, 0, 0);
Vector delta(1, 1, 1);
Vector pts[MAX];

#if FLAT_RASTERIZER
int hit[2*SIZE+1][2*SIZE+1];
#else
int hit[2*SIZE+1][2*SIZE+1][2*SIZE+1];
#endif

extern bool rasterizer_draws;

//------------------------------------------------------------------------------

void drawDot(const Vector P, float size, gym_color col)
{
    flute8 * flu = gym::mapBufferC4V4(1);
    flu[0] = { col, P };
    gym::unmapBufferC4V4();
    gym::ref_view();
    gym::drawSquarePoints(size, 0, 1);
    gym::cleanupCV();
}

void drawDots(unsigned N, const Vector P[], float size, gym_color col)
{
    flute8 * flu = gym::mapBufferC4V4(N);
    unsigned i = 1;
    flu[0] = { col.mix(gym_color(1,1,1,1)), P[0] };
    for ( ; i < N; ++i )
        flu[i] = { col, P[i] };
    gym::unmapBufferC4V4();
    gym::ref_view();
    gym::drawPoints(size, 0, i);
    gym::cleanupCV();
}

void drawLines(unsigned N, const Rasterizer::Vertex2 P[], float size)
{
    gym_color col(1, 0, 1);
    flute8 * flu = gym::mapBufferC4V4(N);
    unsigned i = 0;
    for ( ; i < N; ++i )
        flu[i] = { col, P[i].XX, P[i].YY, 0.f };
    flu[0] = { gym_color(0, 1, 0), P[0].XX, P[0].YY, 0.f };
    gym::unmapBufferC4V4();
    gym::drawLineStrip(size, 0, i);
    gym::cleanupCV();
}

void newPoints()
{
    for ( unsigned i = 0; i < MAX ; ++i )
        pts[i] = (SIZE-1) * Vector::randS();
}


void clearHits()
{
    for ( int i = 0; i <= 2*SIZE; ++i )
    for ( int j = 0; j <= 2*SIZE; ++j )
#if FLAT_RASTERIZER
        hit[i][j] = 0;
#else
    for ( int k = 0; k <= 2*SIZE; ++k )
        hit[i][j][k] = 0;
#endif
}


void punch(int x_inf, int x_sup, int y, int z, void*)
{
    //printf(" { %i %i %i } ", x_inf, x_sup, y);
    if ( std::abs(y) > SIZE ) return;
    if ( std::abs(z) > SIZE ) return;
    
    for ( int x = x_inf; x <= x_sup; ++x )
    {
        //printf(" { %i %i }", x, y);

        if ( -SIZE <= x && x <= SIZE )
        {
#if FLAT_RASTERIZER
            ++hit[x+SIZE][y+SIZE];
#else
            ++hit[x+SIZE][y+SIZE][z+SIZE];
#endif
        }
    }
}

void paint(int x_inf, int x_sup, int y, int z, void*)
{
    gym_color col(0, 1, 0, 0.5);
    flute8 * flu = gym::mapBufferC4V4(2);
    flu[0] = { col, (float)x_inf, (float)y, (float)z };
    flu[1] = { col, (float)x_sup, (float)y, (float)z };
    gym::unmapBufferC4V4();
    gym::drawLines(1, 0, 2);
    gym::drawPoints(4, 0, 1);
    gym::cleanupCV();
}

void rasterize(Vector P, Vector Q, void (func)(int, int, int, int, void*))
{
    real iPQ = 1 / ( P - Q ).norm();
#if FLAT_RASTERIZER
    Rasterizer::paintRectangle(func, nullptr, P, Q, iPQ, radius, shift, delta);
    //Rasterizer::paintRectangle(func, nullptr, P, Q, iPQ, radius);
    //Rasterizer::paintBox2D(func, nullptr, P, Q, radius, shift, delta);
#else
    //Rasterizer::paintHexagonalPrism(func, nullptr, P, Q, iPQ, radius, shift, delta);
    Rasterizer::paintCuboid(func, nullptr, P, Q, iPQ, radius, shift, delta);
    //Rasterizer::paintBox3D(func, nullptr, P, Q, radius, shift, delta);
#endif
}

//------------------------------------------------------------------------------

/// check if 'x' is within distance 'radius' from segment [pq]
int inCylinder(Vector const& P, Vector const& Q, Vector X)
{
    X -= P;
    const Vector PQ = Q - P;
    const real len = PQ.normSqr();
    real abs = dot(PQ, X) / len;
    abs = max_real(0, abs);
    abs = min_real(1, abs);
    return ( X - PQ*abs ).normSqr() <= radius * radius;
}


bool checkPoint(Vector const& P, Vector const& Q, int i, int j, int k)
{
    Vector V(i,j,k);
    int in = inCylinder(P,Q,V);
    
#if FLAT_RASTERIZER
    int ht = hit[i+SIZE][j+SIZE];
#else
    int ht = hit[i+SIZE][j+SIZE][k+SIZE];
#endif
    
    if ( ht != in )
    {
        if ( ht )
            drawDot(V, 7, gym_color(1,0,0,0.5));
        else
            drawDot(V, 7, gym_color(0.5,0.5,0.5));
    }
    
    return ( ht != in );
}


bool testFatLine(Vector P, Vector Q)
{
    bool res = false;
    clearHits();
    rasterize(P, Q, punch);
    for ( int i = -SIZE; i <= SIZE; ++i )
    for ( int j = -SIZE; j <= SIZE; ++j )
#if FLAT_RASTERIZER
        res |= checkPoint(P, Q, i, j, 0);
#else
    for ( int k = -SIZE; k <= SIZE; ++k )
        res |= checkPoint(P, Q, i, j, k);
#endif
    return res;
}


void manyTest()
{
    size_t nb_trials = 1<<10;
    do {
        pts[0] = (SIZE-1) * Vector::randS();
        pts[1] = (SIZE-1) * Vector::randS();
        if ( testFatLine(pts[0], pts[1]) )
            break;
    } while ( nb_trials-- > 0 );
}

//------------------------------------------------------------------------------


void processNormalKey(unsigned char c, int x=0, int y=0)
{
    switch (c) {
        case 27:
        case 'q':
            exit(EXIT_SUCCESS);
        case ' ':
            newPoints();
            break;
        case '0':
            glApp::resetView();
            break;
        case 'p': n_pts = std::min(n_pts+1, MAX); break;
        case 'o': n_pts = std::max(n_pts-1, 1U); break;
        case 'P': n_pts = std::min(n_pts+8, MAX); break;
        case 'O': n_pts = std::max(n_pts-8, 1U); break;
        case '2': n_pts = 2; break;
        case '3': n_pts = 3; break;
        case '4': n_pts = 4; break;
        case '5': n_pts = 5; break;
        case 'r':
            manyTest();
            break;
        case 'd':
            rasterizer_draws = !rasterizer_draws;
            break;
        case 'h':
            printf("keyboard commands:\n"
                   " space : draw a new random distribution\n"
                   " p     : increase number of points\n"
                   " o     : decrease number of points\n"
                   " r     : perform many tests as fast as possible\n");
        default:
            glApp::processNormalKey(c, 0, 0);
    }
    glApp::postRedisplay();
}

//------------------------------------------------------------------------------

void drawGridPoints()
{
    size_t S = 2 * SIZE + 1;
    gym_color col(.5f, .5f, .5f);
    flute8 * flu = gym::mapBufferC4V4(S*S);
    unsigned n = 0;
    for ( int i = -SIZE; i <= SIZE; i += 1)
    for ( int j = -SIZE; j <= SIZE; j += 1)
        flu[n++] = { col, (float)i, (float)j, 0.f };
    gym::unmapBufferC4V4();
    gym::drawPoints(1, 0, n);
    gym::cleanupCV();
}

void drawGrid()
{
    float S = (float)SIZE;
    gym_color col(.5f, .5f, .5f);
    flute8 * flu = gym::mapBufferC4V4(2*SIZE+1);
    unsigned n = 0;
    for ( int i = -SIZE; i <= SIZE; i += 5 )
    {
        float I = (float)i;
        flu[n++] = { col, I, -S, 0.f };
        flu[n++] = { col, I,  S, 0.f };
        flu[n++] = { col, -S, I, 0.f };
        flu[n++] = { col,  S, I, 0.f };
    }
    gym::unmapBufferC4V4();
    gym::drawLines(0.5, 0, n);
    gym::cleanupCV();
}


int display(View& view)
{
    view.openDisplay();

    gym::disableLighting();
#if FLAT_RASTERIZER
    drawGrid();
    drawGridPoints();
#endif
    
    // draw original points:
    drawDots(n_pts, pts, 16, gym_color(1,0,1));
    
#if FLAT_RASTERIZER
    if ( n_pts > 2 )
    {
        Rasterizer::Vertex2 vex[MAX];
        for ( size_t i = 0; i < n_pts; ++i )
            vex[i] = pts[i];
        
        // compute and draw convex-hull
        int nb = Rasterizer::convexHull2D(n_pts, vex);
        drawLines(nb, vex, 2);

        Rasterizer::paintPolygon2D(paint, nullptr, nb, vex, 0);
    }
    else
#endif
    {
        Vector P = pts[0];
        Vector Q = pts[1];
        testFatLine(P,Q);
        rasterize(P,Q,paint);
        view.closeDisplay();
    }
    return 0;
}

/* 
 This only work if rasterizer does not issue openGL commands */
void speedTest(size_t cnt)
{
    clearHits();
    
    for ( size_t n = 1; n < n_pts; ++n )
    {
        Vector P = pts[n-1];
        Vector Q = pts[n];
        real iPQ = 1 / ( P - Q ).norm();
        for ( size_t c = 0; c < cnt; ++c )
        {
#if FLAT_RASTERIZER
            Rasterizer::paintRectangle(punch, nullptr, P, Q, iPQ, radius, shift, delta);
#else
            Rasterizer::paintCuboid(punch, nullptr, P, Q, iPQ, radius, shift, delta);
#endif
        }
    }
}

//------------------------------------------------------------------------------

int main(int argc, char* argv[])
{
    RNG.seed();
    newPoints();
    if ( argc > 1 )
    {
        rasterizer_draws = false;
        speedTest(strtoul(argv[1], nullptr, 10));
    }
    else
    {
        rasterizer_draws = true;
        glutInit(&argc, argv);
        glApp::setDimensionality(3);
        glApp::normalKeyFunc(processNormalKey);
        glApp::setScale(2*(SIZE+radius+1));
        glApp::newWindow(display);
        glApp::attachMenu();
        gle::initialize();

        glutMainLoop();
    }
}

