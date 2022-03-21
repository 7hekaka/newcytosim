// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

/**
 This is a program to analyse simulation results:
 it reads a trajectory-file, and print some data from it.
*/

#include "splash.h"
#include "frame_reader.h"
#include "iowrapper.h"
#include "glossary.h"
#include "messages.h"
#include "parser.h"
#include "simul.h"
#include "simul_prop.h"


void help(std::ostream& os)
{
    os << "Synopsis: generate reports/statistics about cytosim's objects\n";
    os << "\n";
    os << "Syntax:\n";
    os << "       analyse WHAT [frame=INTEGER]\n";
    os << "\n";
    os << "Analyse will generate the same reports as Simul::report()\n";
    os << "See the HTML documentation of Simul::report() for a list of possible values for WHAT\n";
    os << "\n";
    os << "The report is send to the standard output";
    os << "\n";
}

//------------------------------------------------------------------------------


int main(int argc, char* argv[])
{
    if ( argc < 2 || strstr(argv[1], "help") )
    {
        help(std::cout);
        return EXIT_SUCCESS;
    }
    
    if ( strstr(argv[1], "info") || strstr(argv[1], "version")  )
    {
        splash(std::cout);
        print_version(std::cout);
        return EXIT_SUCCESS;
    }

    std::string what = argv[1];
    
    Glossary arg;
    Simul simul;
    FrameReader reader;

    if ( arg.read_strings(argc-2, argv+2) )
        return EXIT_FAILURE;

    try
    {
        simul.loadProperties();
        reader.openFile(simul.prop.system_file);
    }
    catch( Exception & e )
    {
        std::clog << e.brief() << '\n';
        return EXIT_FAILURE;
    }

    Cytosim::silent();
    size_t frm = 0;
    
    // multiple frame indices can be specified:
    if ( arg.set(frm, "frame") )
    {
        size_t s = 0;
        do {
            // try to load the specified frame:
            if ( 0 == reader.loadFrame(simul, frm) )
            {
                simul.poly_report(std::cout, what, arg, frm);
            }
            else
            {
                std::cerr << "Error: missing frame " << frm << '\n';
                return EXIT_FAILURE;
            }
            ++s;
        } while ( arg.set(frm, "frame", s) );
    }
    else
    {
        // process every frame in the file:
        while ( 0 == reader.loadNextFrame(simul) )
        {
            simul.poly_report(std::cout, what, arg, frm);
        }
    }
}
