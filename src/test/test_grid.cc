// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University
// Started FJN on January 2008

#ifndef DIM
#  define DIM 3
#endif

#include "vector.h"
#include "random.h"
#include "glapp.h"
#include "glut.h"
#include "real.h"
#include "timer.h"
#include "gle.h"
#include "gym_flute.h"
#include "gym_flute_dim.h"
#include "gym_draw.h"
#include "gym_cap.h"

#include "grid.h"
#include "grid_display.h"

// grid definitions:
const int range = 1;
real inf[] = {-range,-range, 0 };
real sup[] = { range, range, 1 };
index_t num[] = { 10, 7, 1 };


Grid<real, DIM> grid;

index_t cell_indx;
int coord[4] = { 0 };
Vector3 mouse(0,0,0);
Vector3 node(0,0,0);


bool roundRegions = false;
real regionRadius = 0.125;

//------------------------------------------------------------------------------
void throwMarbles(int cnt)
{
    grid.setValues(0);
    real w[3] = { 0, 0, 0 };
    for( int n = 0; n < cnt; ++n )
    {
        w[0] = range * RNG.sreal();
        w[1] = range * RNG.sreal();
        w[2] = RNG.preal();
        grid(w) += RNG.pint32(4);
    }
}


void processNormalKey(unsigned char c, int x=0, int y=0)
{
    switch (c)
    {
        case 'p':
            for ( int d = 0; d < DIM; ++d )
                grid.setPeriodic(d, !grid.isPeriodic(d));
            break;
            
        case 'i':
            //decrease region-radius
            if ( regionRadius > 1 )
                regionRadius /= M_SQRT2;
            if ( roundRegions )
                grid.createRoundRegions(regionRadius);
            else
                grid.createSquareRegions(regionRadius);
            glApp::flashText("radius = %f", regionRadius);
            break;

        case 'o':
            // increase region-radius
            regionRadius *= M_SQRT2;
            if ( roundRegions )
                grid.createRoundRegions(regionRadius);
            else
                grid.createSquareRegions(regionRadius);
            glApp::flashText("radius = %f", regionRadius);
            break;

        case 'r':
            if ( roundRegions )
                grid.createRoundRegions(regionRadius);
            else
                grid.createSquareRegions(regionRadius);
            glApp::flashText("radius = %f", regionRadius);
            break;

        case 's':
            grid.createSideRegions(regionRadius);
            break;

        case 'S':
            grid.createSideRegions(1);
            break;

        case 'h':
            printf("Shift-click to position test-point\n");
            return;

        case ' ': 
            throwMarbles(6);
            glApp::currentView().reset();
            break;

        default:
            glApp::processNormalKey(c,x,y);
            return;
    }
    glApp::postRedisplay();
}


//------------------------------------------------------------------------------

///set callback for shift-click, with unprojected click position
void processMouseClick(int, int, const Vector3 & a, int)
{
    mouse = a;
    cell_indx = grid.index(mouse);
    grid.setPositionFromIndex(node, cell_indx, 0);
    grid.setCoordinatesFromIndex(coord, cell_indx);
    
    char str[32];

    if ( grid.hasRegions() )
    {
        real val = grid.sumValuesInRegion(cell_indx);
        snprintf(str, sizeof(str), "cell %u : coord %i %i : %.0f marbles",
                 cell_indx, coord[0], coord[1], val);
    } 
    else
    {
        snprintf(str, sizeof(str), "cell %u : coord %i %i", cell_indx, coord[0], coord[1]);
    }
    
    glApp::setMessage(str);
    glApp::postRedisplay();
}

///set callback for shift-drag, with unprojected mouse and click positions
void processMouseDrag(int mx, int my, Vector3 & a, const Vector3 & b, int m)
{
    processMouseClick(mx, my, b, m);
}

//------------------------------------------------------------------------------
#if ( DIM >= 3 )
static gym_color field_color(int, const real& val, Vector3 const&)
{
    return gym_color(val/5.0, 0, 0);
}
#else
static gym_color field_color(int, const real& val, Vector2 const&)
{
    return gym_color(val/5.0, 0, 0);
}
#endif


int display(View& view)
{
    view.openDisplay();
    gym::disableLighting();

#if ( DIM >= 3 )
    Vector3 dir = view.depthAxis();
    drawValues(grid, field_color, 0, dir, 0);
#elif ( DIM > 1 )
    drawValues(grid, field_color, 0);
#endif
    
    float gray[4] = {0.5f,0.5f,0.5f,1.f};
    float white[4] = {1,1,1,1};
    float yellow[4] = {1,1,0,1};
    float blue[4] = {0,0,1,1};

    //--------------draw a grid in gray:
    
    gym::color(gray);
    drawBoundaries(grid, 1);

    //--------------draw content of cells
    const real gold = 2.0 / ( 1.0 + sqrt(5) );
    flute4D * flu = gym::mapBufferC4VD(16*grid.nbCells()+2);
    flute4D * ptr = flu;

    for ( index_t c = 0 ; c < grid.nbCells(); ++c )
    {
        int cnt = grid.icell(c);
        // use Fibonacci's spiral:
        if ( cnt > 0 )
        {
            Vector x, y;
            grid.setPositionFromIndex(x, c, 0.2);
            grid.setPositionFromIndex(y, c, 0.8);
            for ( int u = 0; u < std::min(16, cnt); ++u )
            {
                Vector off(fmod(u*gold,1.0), float(u)/cnt, 0.5);
                *ptr++ = { yellow, x+(y-x).e_mul(off) };
            }
        }
    }

    //-------------draw selected-cell
    *ptr++ = { white, mouse };
    *ptr++ = { blue, node };
    gym::unmapBufferC4VD();
    gym::drawPoints(10, 0, ptr-flu);
    gym::cleanupCV();

    //-------------draw region
    if ( grid.hasRegions() )
    {
        char str[32];
        int const* offset = nullptr;
        size_t nR = grid.getRegion(offset, cell_indx);
        for ( size_t j = 0; j < nR; ++j )
        {
            Vector pos;
            grid.setPositionFromIndex(pos, cell_indx+offset[j], 0.4);
            snprintf(str, sizeof(str), "%d", offset[j]);
            view.drawText(pos, white, str, 0);
        }
    }
    if ( 1 )
    {
        char str[16];
        real val = grid.interpolate(mouse);
        snprintf(str, sizeof(str), "cell %u %f", cell_indx, val);
        glApp::setMessage(str);
    }
    view.closeDisplay();
    return 0;
}


void speedTest()
{
    printf("Speed test...");

    real L[] = { 0, 0, 0};
    real R[] = { 1, 1, 1};
    index_t S[] = { 10, 10, 10};

    Grid<float, 3> tmp;
    tmp.setDimensions(L, R, S);
    tmp.createCells();
    tmp.setValues(0);

    real w[3], A = 0.5;
    for ( int cc=0; cc<10000; ++cc )
    {
        w[0] = RNG.preal();
        w[1] = RNG.preal();
        w[2] = RNG.preal();
        for ( int x = 0; x < 1024; ++x )
        {
            ++tmp(w);
            ++tmp(w+Vector(A,0,0));
            ++tmp(w+Vector(0,A,0));
            ++tmp(w+Vector(A,A,0));
            ++tmp(w+Vector(0,0,A));
            ++tmp(w+Vector(A,0,A));
            ++tmp(w+Vector(0,A,A));
            ++tmp(w+Vector(A,A,A));
        }
    }

    FILE* test = fopen("testgrid.out","w");
    tmp.printValues(test, 0);
    fclose(test);
    printf("wrote file testgrid.out\n");
}


void testInterpolate(unsigned CNT)
{
    real L[] = { 0.0, 0.0, 0.0 };
    real R[] = { 1.0, 1.0, 1.0 };
    index_t S[] = { 100, 100, 100 };
    
    const unsigned MAX = 1 << 14;
    real  rand[MAX+3] = { 0 };
    for ( size_t i = 0; i < MAX+3; ++i )
        rand[i] = RNG.preal();

    Grid<double, 3> map;
    map.setDimensions(L, R, S);
    map.createCells();
    map.setValues(0);
    
    for ( size_t i = 0; i < CNT; ++i )
    {
        real W[] = { RNG.preal(), RNG.preal(), RNG.preal() };
        ++map(W);
    }

    real ** vec = new real*[CNT];
    for ( size_t i = 0; i < CNT; ++i )
        vec[i] = rand + RNG.pint32(MAX);

    tick();
    real sum = 0;
    for ( size_t r = 0; r < 100; ++r )
    for ( size_t i = 0; i < CNT; ++i )
        sum += map.interpolate3D(vec[i]) + map.interpolate3D(vec[i]);
    printf("  interpolate3D  CPU %7.3f  sum = %f\n", tock(CNT), sum);
    
    tick();
    real som = 0;
    for ( size_t r = 0; r < 100; ++r )
    for ( size_t i = 0; i < CNT; ++i )
        som += map.interpolate(vec[i]) + map.interpolate(vec[i]);
    printf("  interpolate    CPU %7.3f  sum = %f\n", tock(CNT), som);

    delete[] vec;
}


int main(int argc, char* argv[])
{
    RNG.seed();

    if ( argc > 1 )
    {
        testInterpolate(1<<20);
        return 0;
    }

    //initialize the grid:
    grid.setDimensions(inf, sup, num);
    grid.createCells();
    //grid.setPeriodic(0, true);
    //grid.setPeriodic(1, true);
    throwMarbles(1<<DIM);

    glutInit(&argc, argv);
    glApp::setDimensionality(DIM);
    glApp::actionFunc(processMouseClick);
    glApp::actionFunc(processMouseDrag);
    glApp::normalKeyFunc(processNormalKey);
    glApp::newWindow(display);
    glApp::attachMenu();
    glApp::setScale(2*range+1);
    gle::initialize();
    glutMainLoop();
}
