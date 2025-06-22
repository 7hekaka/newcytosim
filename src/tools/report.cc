// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University
/**
 This is a program to analyse simulation results:
 it reads a trajectory-file, and print some data from it.
*/

#include <fstream>
#include <sstream>

#include "stream_func.h"
#include "frame_reader.h"
#include "iowrapper.h"
#include "filepath.h"
#include "glossary.h"
#include "messages.h"
#include "splash.h"
#include "parser.h"
#include "simul.h"
#include "slice.h"

Simul simul;
FrameReader reader;
int prefix = 0;


void help(std::ostream& os)
{
    os << "Cytosim-report\n";
    os << "       generates reports/statistics from a trajectory file\n";
    os << " Syntax:\n";
    os << "       report [time] WHAT [OPTIONS]\n";
    os << " Options:\n";
    os << "       precision=INTEGER\n";
    os << "       column=INTEGER\n";
    os << "       verbose=0\n";
    os << "       frame=INTEGER[,INTEGER[,INTEGER[,INTEGER]]]\n";
    os << "       frame=INTEGER:INTEGER:INTEGER\n";
    os << "       input=FILE_NAME\n";
    os << "       output=FILE_NAME\n\n";
    os << "  This tool must be invoked in a directory containing the simulation output,\n";
    os << "  and it will generate reports by calling Simul::report(). The only required\n";
    os << "  argument `WHAT` determines what sort of data will be generated. Many options are\n";
    os << "  available, please check the HTML documentation for a list of possibilities.\n\n";
    os << "  By default, all frames in the file are processed in order, but a frame index,\n";
    os << "  or multiple indices can be specified (the first frame has index 0).\n";
    os << "  A periodicity can also be specified (ignored if multiple frames are specified).\n\n";
    os << "  The input trajectory file is `objects.cmo` unless otherwise specified.\n";
    os << "  The result is sent to standard output unless a file is specified as `output=NAME`\n";
    os << "  Attention: there should be no whitespace in any of the option.\n\n";
    os << "Examples:\n";
    os << "       report fiber:points\n";
    os << "       report fiber:points frame=10 > fibers.txt\n";
    os << "       report fiber:points frame=10,20 > fibers.txt\n";
    os << "       report fiber:points period=8 > fibers.txt\n";
    os << "Made with format version " << Simul::currentFormatID << " and DIM=" << DIM << "\n";
}

//------------------------------------------------------------------------------

void report_prefix(std::ostream& os, Simul const& sim, std::string const& what, Glossary& opt, size_t frm)
{
    char str[256] = { 0 };
    size_t str_len = 0;
    
    if ( prefix & 1 )
        str_len += snprintf(str, sizeof(str), "%9.3f ", sim.time());
    
    if ( prefix & 2 )
        str_len += snprintf(str+str_len, sizeof(str)-str_len, "%9lu ", frm);
    
    std::stringstream ss;
    sim.poly_report(ss, what, opt, -1);
    StreamFunc::prefix_lines(os, ss, str, 0, '%');
}


int report(std::ostream& os, std::string const& what, Glossary& opt, size_t frm)
{
    try
    {
        if ( 0 == reader.loadFrame(simul, frm) )
        {
            if ( prefix )
                report_prefix(os, simul, what, opt, frm);
            else
                simul.poly_report(os, what, opt, frm);
            return 0;
        }
        return 1;
    }
    catch( Exception & e )
    {
        std::cerr << e.brief() << "\n";
        return 2;
    }
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
    
    Glossary arg;

    std::string input = Simul::TRAJECTORY;
    std::string str, what;
    std::ofstream ofs;
    std::ostream out(std::cout.rdbuf());

    // check for prefix:
    int ax = 1;
    while ( argc > ax+1 )
    {
        if ( 0 == strncmp(argv[ax], "time", 4) )
            prefix |= 1;
        else if ( 0 == strncmp(argv[ax], "frame", 5) )
            prefix |= 2;
        else
            break;
        ++ax;
    }
    
    what = argv[ax++];
    if ( arg.read_strings(argc-ax, argv+ax) )
        return EXIT_FAILURE;

    // change working directory if specified:
    if ( arg.has_key("directory") )
    {
        FilePath::change_dir(arg.value("directory"));
        std::clog << "Cytosim working directory is " << FilePath::get_cwd() << '\n';
    }

#if BACKWARD_COMPATIBILITY < 50
    if ( arg.set(str, "prefix") && str=="time" )
        prefix = 1;
#endif
    
    arg.set(input, ".cmo") || arg.set(input, "input");
    if ( arg.use_key("-") ) arg.define("verbose", "0");

    RNG.seed();
    
    try
    {
        simul.loadProperties();
        reader.openFile(input);

        if ( arg.set(str, "output") )
        {
            ofs.open(str.c_str());
            if ( ofs.is_open() )
                out.rdbuf(ofs.rdbuf());
            else
            {
                std::clog << "Cannot open output file\n";
                return EXIT_FAILURE;
            }
        }
    }
    catch( Exception & e )
    {
        std::clog << e.brief() << '\n';
        return EXIT_FAILURE;
    }
    
    Cytosim::silent();
    Slice slice;

    // process first record, and check things:
    size_t frm = slice.first();
    if ( reader.loadFrame(simul, frm) )
    {
        std::cerr << "Error: missing frame " << frm << '\n';
        return EXIT_FAILURE;
    }
    if ( DIM != reader.vectorSize() )
    {
        std::cerr << "Error: dimensionality missmatch between `report` and file\n";
        return EXIT_FAILURE;
    }

    int cnt = 0;
    // multiple ranges of indices can be specified:
    if ( arg.num_values("frame") > 0 )
    {
        int inp = 0;
        while ( arg.set(str, "frame", inp++) )
        {
            slice.parse(str.c_str());
            for ( frm = slice.first(); frm <= slice.last(); frm += slice.increment() )
            {
                if ( 0 == report(out, what, arg, frm) )
                    ++cnt;
                else
                    break;
            }
        }
    }
    else
    {
        // if 'frame' is not specified, a 'period' can still be specified
        int p = 1;
        if ( arg.set(p, "period") )
            slice.increment(p);
        for ( frm = slice.first(); frm <= slice.last(); frm += slice.increment() )
        {
            if ( 0 == report(out, what, arg, frm) )
                ++cnt;
            else
                break;
        }
    }
    out << '\n';
    arg.print_warnings(stderr, cnt, "\n");
}
