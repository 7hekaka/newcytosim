// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
/**
 This is a program to export simulation objects, in a format readable by Blender
 
 It stores the points of the filaments as they are generated.
 Attachment positions of Couple are relocated to the closest filament-point.

 F. Nedelec, 20.09.2017 and 29.01.2018 to 3.02.2018, 1.02.2019
*/

#include <fstream>
#include "frame_reader.h"
#include "iowrapper.h"
#include "glossary.h"
#include "simul.h"

#include <vector>
#include <map>

bool binary = 0;
FILE * file = nullptr;

typedef std::vector<Vector3> PointList;
typedef std::map<ObjectID, PointList> PointMap;

PointMap points;

// list corresponding to current fiber
PointList * record = 0;


void help(std::ostream& os)
{
    os << "Syntax:\n";
    os << "    cymart CLASS [style=filament|actin|microtubule] [INPUTFILE] [binary=1]\n";
    os << "Cymart exports Cytosim's objects to a Blender-friendly format\n";
    os << "\n";
    os << "`CLASS` names the property of the fiber (use `all` to export all fibers)\n";
    os << "Select `style=filament` or `style=actin` or `style=microtubule`.\n";
    os << "Each frame of the trajectory file is sent to a separate file.\n";
    os << "By default these files are plain ASCII text files\n";
    os << "This is " << DIM << "D\n";
}


FILE * openFile(const char base[], unsigned inx)
{
    char name[256] = { 0 };
    FILE * f = nullptr;
    
    if ( binary )
    {
        snprintf(name, sizeof(name), "%s%04u.sweet", base, inx);
        f = fopen(name, "wb");
    }
    else
    {
        snprintf(name, sizeof(name), "%s%04u.txt", base, inx);
        f = fopen(name, "w");
    }
    if ( !f )
        std::cerr << "Error: could not create `" << name << "'\n";
    if ( f && ferror(f) )
    {
        fclose(f);
        std::cerr << "Error making `" << name << "'\n";
        f = nullptr;
    }
    return f;
}


//------------------------------------------------------------------------------
#pragma mark - write data chunk

void writeBinary(FILE* file, uint16_t a, uint32_t b, uint16_t c, Vector3 const& pos)
{
    fwrite(&a, 2, 1, file);
    fwrite(&b, 4, 1, file);
    fwrite(&c, 2, 1, file);
    float vec[3] = { 0 };
    pos.store(vec);
    fwrite(vec, 3, sizeof(float), file);
}


void writeText(FILE* file, uint16_t a, uint32_t b, uint16_t c, Vector3 const& pos)
{
    fprintf(file, "%i %i %i ", a, b, c);
    fprintf(file, "%.6f %.6f %.6f\n", pos.XX, pos.YY, pos.ZZ);
}

void writePoint(uint16_t a, uint32_t b, uint16_t c, Vector3 const& pos)
{
    if ( binary )
        writeBinary(file, a, b, c, pos);
    else
        writeText(file, a, b, c, pos);
}

//------------------------------------------------------------------------------
#pragma mark - write data chunk

Vector3 closestPoint(Fiber const* fib, const Vector3 pos)
{
    PointList const& pts = points[fib->identity()];
    
    if ( pts.empty() )
        return Vector3(0,0,0);
    
    Vector3 res = pts[0];
    real dis = distanceSqr(pos, res);
    
    for ( Vector3 const& vec : pts )
    {
        real d = distanceSqr(pos, vec);
        if ( d < dis )
        {
            res = vec;
            dis = d;
        }
    }
    return res;
}


Vector3 closestPoint(Fiber const* fib, const Vector2 pos)
{
    return closestPoint(fib, Vector3(pos));
}


void clearPoints()
{
    for ( PointMap::iterator i = points.begin(); i != points.end(); ++i )
        i->second.clear();
    points.clear();
}

//------------------------------------------------------------------------------
#pragma mark -

/// save Microtubules protofilament
void drawFilament(Fiber const& fib)
{
    uint16_t a = 1;
    uint16_t b = 0;
    uint32_t c = fib.identity();

    // axial translation between two sucessive monomers:
    const real dab = 0.004;
    
    real ab = 0;
    int cnt = 0;
    // increment until we reach the minus end
    while ( ab <= fib.abscissaM() )
    {
        ++cnt;
        ab += dab;
    }
    // draw the monomers until the plus end:
    while ( ab < fib.abscissaP() )
    {
        // alternate colors:
        b = 1 + ( cnt & 1 );
        writePoint(a, b, c, Vector3(fib.posP(ab)));
        ab += dab;
    }
}


/// save double helical filament
void drawActin(Fiber const& fib)
{
    uint16_t a = 1;
    uint16_t b = 0;
    uint32_t c = fib.identity();

    // axial translation between two sucessive monomers:
    const real dab = 0.00275;
    // distance from central axis to center of monomers
    real off = 0.0045 - dab;
    
    // rotation angle between consecutive monomers
    const real dan = -166 * M_PI / 180;
    const real cs = std::cos(dan);
    const real sn = std::sin(dan);
    
    real ab = 0;
    Vector3 p(0,0,0); //position of monomer;
    Vector3 d(fib.dirEndM());   // unit tangent to centerline
    Vector3 n = fib.adjustedNormal(d);
    //std::clog << fib.reference() << " " << n << "    " << n.normSqr() << " " << n*d << '\n';
    
    int cnt = 0;
    // rotate normal until we reach the minus end
    while ( ab < fib.abscissaM() )
    {
        ++cnt;
        n = d.rotateOrtho(n, cs, sn);
        ab += dab;
    }
    
    // draw the monomers until the plus end
    while ( ab < fib.abscissaP() )
    {
        d = Vector3(fib.dir(ab));
        n = d.rotateOrtho(n, cs, sn);
        p = Vector3(fib.pos(ab)) + off * n;

        // use different tones to individualize the two strands:
        b = 1 + ( cnt & 1 );
        
        // change color for the last 2 monomers:
        if ( ab + dab > fib.abscissaP() )
            b += 2;

        writePoint(a, b, c, p);
        ab += dab;
        ++cnt;
    }
}


/// save tubular structure made of 13-protofilaments
void drawMicrotubule(Fiber const& fib)
{
    real da[] = {0,0.000923,0.001846,0.002769,0.003692,0.004615,0.005538,0.006461,0.007384,0.008308,0.009231,0.010154,0.011077};
    real dx[] = {0.8855,0.5681,0.1205,-0.3546,-0.7485,-0.9709,-0.9709,-0.7485,-0.3546,0.1205,0.5681,0.8855,1.0000};
    real dy[] = {-0.4647,-0.8230,-0.9927,-0.9350,-0.6631,-0.2393,0.2393,0.6631,0.9350,0.9927,0.8230,0.4647,0};

    uint16_t a = 1;
    uint16_t b = 0;
    uint32_t c = fib.identity();

    // axial translation between two sucessive monomers:
    const real dab = 0.004;
    // enlarged radius of monomers makes them overlap slighlty
    const real rad = 0.7 * dab;
    // distance from central axis to center of monomers
    real off = 0.025 / 2 - rad;
    
    real ab = dab * std::ceil( fib.abscissaM() / dab );
    Vector3 d(fib.dir(ab));   // unit tangent vector
    Vector3 n = fib.adjustedNormal(d);
    
    const real abmax = fib.abscissaP();

    int cnt = 0;
    while ( ab < abmax )
    {
        d = Vector3(fib.dir(ab));
        Vector3 p(fib.pos(ab));
        
        // adjust 'n' to keep it orthogonal to 'd':
        n = d.orthogonal(n, 1.0);
        
        // set two vectors orthogonal to 'd' of norm 'off':
        Vector3 e = n * off;
        Vector3 f = cross(d, e);
        
        // color of alpha subunit:
        b = 1;
        
        real up = abmax - ( cnt & 1 ? 0 : dab );
  
        for ( int i = 0; i < 13; ++i )
        {
            if ( ab + da[i] < up )
            {
                if ( cnt & 1 )
                {
                    // change color for beta-tubulin near the end:
                    if ( ab + da[i] + 4 * dab > abmax )
                        b = 3;
                    else
                        b = 2;
                }
                writePoint(a, b, c, p+dx[i]*e+dy[i]*f+da[i]*d);
            }
        }
        
        ab += dab;
        ++cnt;
    }
}


void drawLink(Couple const* cop)
{
    int   dat[3] = { 7, (int)cop->identity(), (int)cop->prop->number() };
    float vec[6] = { 0 };
    
    closestPoint(cop->fiber1(), cop->posHand2()).store(vec);
    closestPoint(cop->fiber2(), cop->posHand1()).store(vec+3);

    if ( binary )
    {
        fwrite(dat, 3, sizeof(int), file);
        fwrite(vec, 6, sizeof(float), file);
    }
    else
    {
        fprintf(file, "%i %i %i ", dat[0], dat[1], dat[2]);
        fprintf(file, "%.6f %.6f %.6f %.6f %.6f %.6f\n", vec[0], vec[1], vec[2], vec[3], vec[4], vec[5]);
    }
}

//------------------------------------------------------------------------------
#pragma mark -


int main(int argc, char* argv[])
{
    if ( argc < 2 || strstr(argv[1], "help") )
    {
        help(std::cout);
        return EXIT_SUCCESS;
    }
    
    std::string fiber_type = argv[1];
    Property * selected = 0;
    
    Simul simul;
    Glossary arg;
    FrameReader reader;

    if ( arg.read_strings(argc-2, argv+2) )
        return EXIT_FAILURE;

    int style = 1;
    arg.set(binary, "binary");
    arg.set(style, "style", {{"filament", 1}, {"actin", 2}, {"microtubule", 3}});

    std::string input = Simul::TRAJECTORY;
    arg.set(input, ".cmo") || arg.set(input, "input");

    unsigned frame = 0;
    RNG.seed();

    try
    {
        FiberSet const& fibers = simul.fibers;
        CoupleSet const& couples = simul.couples;
        simul.loadProperties();
        if ( fiber_type != "all" )
            selected = simul.properties.find_or_die("fiber", fiber_type);
        
        reader.openFile(input);
        
        // process all frames in the file:
        while ( 0 == reader.loadNextFrame(simul) )
        {
            switch ( style )
            {
                case 1: file = openFile("filament", frame);    break;
                case 2: file = openFile("actin", frame);       break;
                case 3: file = openFile("microtubule", frame); break;
            }
            
            clearPoints();
            
            // save fibers in the natural order:
            for ( Fiber const* fib = fibers.firstID(); fib; fib = fibers.nextID(fib) )
            {
                if ( selected == 0 || fib->prop == selected )
                {
                    //if ( ! binary ) fprintf(file, "%% fiber %li\n", fib->identity());
                    record = &points[fib->identity()];
                    
                    switch ( style )
                    {
                        case 1: drawFilament(*fib);    break;
                        case 2: drawActin(*fib);       break;
                        case 3: drawMicrotubule(*fib); break;
                    }
                    
                    //std::clog << "   fiber " << fib->reference() << " has " << record->size() << " points\n";
                }
            }
            fclose(file);

            std::clog << "frame " << frame << " has " << points.size() << " fibers\n";
            
            // save links from Couple in natural order:
            file = openFile("link", frame);
            for( Couple const* cop = couples.firstID(); cop; cop = couples.nextID(cop) )
            {
                Fiber const* f1 = cop->fiber1();
                Fiber const* f2 = cop->fiber2();
                if ( f1 && f2 && ( selected==0 || ( f1->prop==selected && f2->prop==selected )))
                    drawLink(cop);
            }
            fclose(file);
            ++frame;
        }
    }
    catch( Exception & e )
    {
        std::cerr << e.brief() << '\n';
        return EXIT_FAILURE;
    }
}
