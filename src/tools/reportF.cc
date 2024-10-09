// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
/**
 This is a program to analyse simulation results:
 it reads a trajectory-file, and print some data from it.
*/

#include <fstream>
#include <sstream>

#include "stream_func.h"
#include "frame_reader.h"
#include "iowrapper.h"
#include "glossary.h"
#include "messages.h"
#include "splash.h"
#include "parser.h"
#include "simul.h"

char root[256] = "report";

void help(std::ostream& os)
{
    os << "Cytosim-reportF\n";
    os << "       generates reports/statistics from a trajectory file\n";
    os << "       this tool is simular to 'report', but generates one file per frame\n";
    os << " Syntax:\n";
    os << "       reportF WHAT [OPTIONS]\n";
    os << " Options:\n";
    os << "       precision=INTEGER\n";
    os << "       column=INTEGER\n";
    os << "       verbose=0\n";
    os << "       frame=INTEGER\n";
    os << "       period=INTEGER\n";
    os << "       input=FILE_NAME\n";
    os << "       root=STRING\n";
    os << " This tool must be invoked in a directory containing the simulation output,\n";
    os << " and it will generate reports by calling Simul::report(). The only required\n";
    os << " argument `WHAT` determines what sort of data will be generated. Many options are\n";
    os << " available, please check the HTML documentation for a list of possibilities.\n\n";
    os << " By default, all frames in the file are processed in order, but a initial frame\n";
    os << " and a periodicity can be specified.\n";
    os << " The input trajectory file is `objects.cmo` unless otherwise specified.\n";
    os << " The result is sent to standard output unless a file is specified as `output=NAME`\n";
    os << " Attention: there should be no whitespace in any of the option.\n";
    os << "\n";
    os << "of the trajectory to a different file. These files are named:\n";
    os << "    ROOT####.txt\n";
    os << "where #### is the frame number and ROOT can be specified.\n";
    os << "Made with format version " << Simul::currentFormatID << " and DIM=" << DIM << "\n";
}

//------------------------------------------------------------------------------

void report(Simul const& sim, std::string const& what, Glossary& opt, size_t frm)
{
    char filename[512];
    snprintf(filename, sizeof(filename), "%s%04lu.txt", root, frm);
    std::ofstream os(filename);
    sim.poly_report(os, what, opt, frm);
}

//------------------------------------------------------------------------------


int main(int argc, char* argv[])
{
    if ( argc < 2 || strstr(argv[1], "help") )
    {
        help(std::cout);
        return EXIT_SUCCESS;
    }
    
    if ( strstr(argv[1], "info") || strstr(argv[1], "--version")  )
    {
        splash(std::cout);
        print_version(std::cout);
        return EXIT_SUCCESS;
    }
    
    Simul simul;
    Glossary arg;
    FrameReader reader;

    std::string input = Simul::TRAJECTORY;
    std::string str, what = argv[1];

    if ( arg.read_strings(argc-2, argv+2) )
        return EXIT_FAILURE;

    size_t frame = 0;
    size_t period = 1;

    arg.set(input, ".cmo") || arg.set(input, "input");
    if ( arg.use_key("-") ) arg.define("verbose", "0");

    RNG.seed();

    if ( arg.set(str, "root") )
        strncpy(root, str.c_str(), sizeof(root));
    
    try
    {
        simul.loadProperties();
        reader.openFile(input);
        
        arg.set(frame, "frame");
        if ( arg.set(period, "period") )
            period = std::max(1UL, period);
    }
    catch( Exception & e )
    {
        std::clog << e.brief() << '\n';
        return EXIT_FAILURE;
    }
    
    Cytosim::silent();
    size_t cnt = 0;
    
    try
    {
        if ( reader.loadFrame(simul, frame) )
        {
            std::cerr << "Error: missing frame " << frame << '\n';
            return EXIT_FAILURE;
        }
        if ( DIM != reader.vectorSize() )
        {
            std::cerr << "Error: dimensionality missmatch between `report` and file\n";
            return EXIT_FAILURE;
        }
        
        do {
            report(simul, what, arg, frame);
            frame += period;
            ++cnt;
        } while ( 0 == reader.loadFrame(simul, frame) );
    }
    catch( Exception & e )
    {
        std::cerr << e.brief() << "\n";
        return EXIT_FAILURE;
    }
    
    arg.print_warnings(stderr, cnt, "\n");
}
