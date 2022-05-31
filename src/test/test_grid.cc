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

#include "grid.h"
#include "grid_display.h"

//area of the grid
const int range = 5;
real     left[] = {  -range,  -range,  0 };
real    right[] = {   range,   range,  0.5 };
size_t n_cell[] = { 2*range,   range,  2 };


typedef Grid<real, DIM> grid_type;
grid_type myGrid;

size_t cell_indx;
int coord[4] = { 0 };
Vector3 pos(0,0,0);
Vector3 nod(0,0,0);


#define TEST_REGIONS
real regionRadius = 1.5;

//------------------------------------------------------------------------------
void throwMarbles(int cnt)
{
    myGrid.setValues(0);
    real w[3] = { 0, 0, 0 };
    for( int n = 0; n < cnt; ++n )
    {
        w[0] = range * RNG.sreal();
        w[1] = range * RNG.sreal();
        w[2] = range * RNG.sreal();
        ++myGrid(w);
    }
}


void processNormalKey(unsigned char c, int x=0, int y=0)
{
    switch (c)
    {
        case 'p':
            for ( size_t d = 0; d < DIM; ++d )
                myGrid.setPeriodic(d, !myGrid.isPeriodic(d));
            break;
            
        case 'i':
            //decrease region-radius
            if ( regionRadius > 1 )
                regionRadius -= 0.25;
            myGrid.createRoundRegions(regionRadius);
            //myGrid.createSquareRegions(regionRadius);
            glApp::flashText("radius = %f", regionRadius);
            break;

        case 'o':
            // increase region-radius
            regionRadius += 0.25;
            myGrid.createRoundRegions(regionRadius);
            //myGrid.createSquareRegions(regionRadius);
            glApp::flashText("radius = %f", regionRadius);
            break;

        case 'r':
            myGrid.createRoundRegions(regionRadius);
            glApp::flashText("radius = %f", regionRadius);
            break;

        case 's':
            myGrid.createSideRegions(regionRadius);
            break;
            
        case 'h':
            printf("Shift-click to position test-point\n");
            return;

        case ' ': 
            throwMarbles(32);
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
    pos = a;
    cell_indx = myGrid.index(pos);
    myGrid.setPositionFromIndex(nod, cell_indx, 0);
    myGrid.setCoordinatesFromIndex(coord, cell_indx);
    
    char str[32];

    if ( myGrid.hasRegions() )
    {
        real num = myGrid.sumValuesInRegion(cell_indx);
        snprintf(str, sizeof(str), "cell %lu : coord %i %i : %.0f marbles",
                 cell_indx, coord[0], coord[1], num);
    } 
    else
    {
        snprintf(str, sizeof(str), "cell %lu : coord %i %i", cell_indx, coord[0], coord[1]);
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
static gle_color field_color(int, const real& val, Vector3 const&)
{
    return gle_color(val/5.0, 0, 0);
}
#else
static gle_color field_color(int, const real& val, Vector2 const&)
{
    return gle_color(val/5.0, 0, 0);
}
#endif


int display(View& view)
{
    view.openDisplay();

#if ( DIM >= 3 )
    Vector3 dir = view.depthAxis();
    drawValues(myGrid, field_color, 0, dir, 0);
#elif ( DIM > 1 )
    drawValues(myGrid, field_color, 0);
#endif
    
    float gray[4] = {1,1,1,.5f};
    float white[4] = {1,1,1,1};
    float yellow[4] = {0,1,1,1};

    //--------------draw a grid in gray:
    
    gym::color(gray);
    drawBoundaries(myGrid, 1);

    //--------------draw content of cells
    const real gold = 2.0 / ( 1.0 + sqrt(5) );
    flute4D * flu = gym::mapBufferC4VD(16*myGrid.nbCells()+2);
    flute4D * ptr = flu;

    for ( size_t c = 0 ; c < myGrid.nbCells(); ++c )
    {
        int cnt = myGrid.icell(c);
        // use Fibonacci's spiral on the sphere:
        if ( cnt > 0 )
        {
            Vector x, y;
            myGrid.setPositionFromIndex(x, c, 0.3);
            myGrid.setPositionFromIndex(y, c, 0.7);
            for ( int u = 0; u < std::min(16, cnt); ++u )
            {
                Vector off(fmod(u*gold,1.0), float(u)/cnt, 0.5);
                *ptr++ = { yellow, x+(y-x).e_mul(off) };
            }
        }
    }

    //-------------draw selected-cell
    gle_color lor(1,1,0);
    *ptr++ = { lor, Vector2(pos) };
    *ptr++ = { lor, Vector2(nod) };
    gym::unmapBufferC4VD();
    gym::drawPoints(12, 0, ptr-flu);
    gym::cleanup();

    //-------------draw region
    if ( myGrid.hasRegions() )
    {
        char str[16];
        int * offset = nullptr;
        size_t nb = myGrid.getRegion(offset, cell_indx);
        for ( size_t j = 0; j < nb; ++j )
        {
            Vector x;
            myGrid.setPositionFromIndex(x, cell_indx+offset[j], 0.4);
            snprintf(str, sizeof(str), "%lu", j);
            view.drawText(x, white, str, 0);
        }
    }
    else
    {
        char str[16];
        real val = myGrid.interpolate(pos);
        snprintf(str, sizeof(str), "cell %lu %f", cell_indx, val);
        glApp::setMessage(str);
    }
    view.closeDisplay();
    return 0;
}


void speedTest()
{
    printf("Real test...");

    real   L[] = { 0, 0, 0};
    real   R[] = { 1, 1, 1};
    size_t S[] = { 10, 10, 10};

    Grid<float, 3> map;
    map.setDimensions(L, R, S);
    map.createCells();
    map.setValues(0);

    real w[3];
    for ( int cc=0; cc<10000; ++cc )
    {
        w[0] = RNG.preal();
        w[1] = RNG.preal();
        w[2] = RNG.preal();
        for ( int x = 0; x < 1000; ++x )
        {
            ++map( w );
            ++map( w );
            ++map( w );
            ++map( w );
            ++map( w );
            ++map( w );
            ++map( w );
            ++map( w );
            ++map( w );
            ++map( w );
        }
    }

    FILE* test = fopen("testgrid.out","w");
    map.printValues(test, 0);
    fclose(test);
    printf("wrote file testgrid.out\n");
}


void testInterpolate(unsigned CNT)
{
    real   L[] = { 0.0, 0.0, 0.0 };
    real   R[] = { 1.0, 1.0, 1.0 };
    size_t S[] = { 100, 100, 100 };
    
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
    myGrid.setDimensions(left, right, n_cell);
    myGrid.createCells();
    //myGrid.setPeriodic(0, true);
    //myGrid.setPeriodic(1, true);
    throwMarbles(8*(1<<DIM));

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
