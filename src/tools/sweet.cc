// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
/**
 This is a program to export simulation objects, in a compact and efficient format
 
 It generated a binary or a text file following the sweet16 format:
 
 [2 bytes unsigned integer: class]
 [4 bytes unsigned integer: identification]
 [2 bytes unsigned integer: info/color]
 [4 bytes float: X coordinate]
 [4 bytes float: Y coordinate]
 [4 bytes float: Z coordinate]

 F. Nedelec, 20.11.2020
*/

#include <fstream>
#include "frame_reader.h"
#include "iowrapper.h"
#include "glossary.h"
#include "simul.h"

#include <vector>
#include <map>

FILE * file = nullptr;

void help(std::ostream& os)
{
    os << "\n";
    os << "Sweet exports Cytosim's object coordinates to files\n";
    os << "\n";
    os << "Syntax:\n";
    os << "       sweet class [INPUTFILE] [binary=1]\n";
    os << "\n";
    os << "`class` names the objects to be exported (use `all` to export all fibers)\n";
    os << "Each frame of the trajectory file is sent to a separate file.\n";
    os << "By default these files are plain ASCII text files\n";
}


FILE * openFile(const char base[], unsigned inx, bool binary)
{
    char name[256] = { 0 };
    FILE * f = nullptr;
    
    if ( binary )
    {
        snprintf(name, sizeof(name), "%s%04u.bin", base, inx);
        f = fopen(name, "wb");
    }
    else
    {
        snprintf(name, sizeof(name), "%s%04u.txt", base, inx);
        f = fopen(name, "w");
    }
    if ( !f )
    {
        std::cerr << "Aborted: could not create file " << name << '\n';
        exit(EXIT_FAILURE);
    }
    if ( ferror(f) )
    {
        fclose(f);
        std::cerr << "Aborted: could not create file " << name << '\n';
        exit(EXIT_FAILURE);
    }
    return f;
}

//------------------------------------------------------------------------------
#pragma mark -


void writeBinary(const unsigned int dat[3], Vector3 const& pos)
{
    uint16_t a = dat[0];
    uint32_t b = dat[1];
    uint16_t c = dat[2];
    if ( a != dat[0] || b != dat[1] || c != dat[2] )
    {
        fprintf(stderr, "sweet20 data format overflow\n");
        exit(EXIT_FAILURE);
    }
    fwrite(&a, 1, 2, file);
    fwrite(&b, 1, 4, file);
    fwrite(&c, 1, 2, file);
    float vec[3] = { 0 };
    pos.store(vec);
    fwrite(vec, 3, sizeof(float), file);
}


void writeText(const unsigned int dat[3], Vector3 const& pos)
{
    fprintf(file, "%i %i %i ", dat[0], dat[1], dat[2]);
    fprintf(file, "%.6f %.6f %.6f\n", pos.XX, pos.YY, pos.ZZ);
}


//------------------------------------------------------------------------------
#pragma mark -


int main(int argc, char* argv[])
{
    if ( argc > 1 && strstr(argv[1], "help") )
    {
        help(std::cout);
        return EXIT_SUCCESS;
    }
    
    Glossary arg;
    if ( arg.read_strings(argc-1, argv+1) )
        return EXIT_FAILURE;

    bool binary = 0;
    arg.set(binary, "binary");
    std::string input = TRAJECTORY;
    arg.set(input, ".cmo") || arg.set(input, "input");

    Simul simul;
    simul.loadProperties();
    FrameReader reader;
    reader.openFile(input);
    
    if ( !reader.good() )
    {
        fprintf(stderr, "could not open input file\n");
        return EXIT_FAILURE;
    }
    
    try
    {
        unsigned frame = 0;
        // process all frames in the file:
        while ( 0 == reader.loadNextFrame(simul) )
        {
            file = openFile("fil", frame++, binary);
            // process fibers in the natural order:
            for ( Fiber const* fib = simul.fibers.firstID(); fib; fib = simul.fibers.nextID(fib) )
            {
                unsigned int data[3] = { 0 };
                data[0] = fib->prop->number();
                data[1] = fib->identity();
                if ( binary )
                {
                    for ( size_t i = 0; i < fib->nbPoints(); ++i )
                        writeBinary(data, fib->posPoint(i));
                }
                else
                {
                    for ( size_t i = 0; i < fib->nbPoints(); ++i )
                        writeText(data, fib->posPoint(i));
                }
            }
            fclose(file);
        }
    }
    catch( Exception & e )
    {
        std::cerr << e.brief() << '\n';
        return EXIT_FAILURE;
    }
    
    return EXIT_SUCCESS;
}
