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
    os << "       period=INTEGER\n";
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

void report_prefix(Simul const& sim, std::ostream& os, std::string const& what, size_t frm, Glossary& opt)
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


void report(Simul const& sim, std::ostream& os, std::string const& what, size_t frm, Glossary& opt)
{
    try
    {
        if ( prefix )
            report_prefix(sim, os, what, frm, opt);
        else
            sim.poly_report(os, what, opt, frm);
    }
    catch( Exception & e )
    {
        std::cerr << e.brief() << '\n';
        exit(EXIT_FAILURE);
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
    
    Simul simul;
    Glossary arg;
    FrameReader reader;

    std::string input = Simul::TRAJECTORY;
    std::string str, what;
    std::ofstream ofs;
    std::ostream out(std::cout.rdbuf());

    // check for prefix:
    int ax = 1;
    while ( argc > ax+1 )
    {
        if ( strstr(argv[ax], "time") )
            prefix |= 1;
        else if ( strstr(argv[ax], "frame") )
            prefix |= 2;
        else
            break;
        ++ax;
    }
    
    what = argv[ax++];
    if ( arg.read_strings(argc-ax, argv+ax) )
        return EXIT_FAILURE;

#if BACKWARD_COMPATIBILITY < 50
    if ( arg.set(str, "prefix") && str=="time" )
        prefix = 1;
#endif
    
    size_t frame = 0;
    size_t period = 1;

    arg.set(input, ".cmo") || arg.set(input, "input");
    if ( arg.use_key("-") ) arg.define("verbose", 0);

    RNG.seed();

    try
    {
        simul.loadProperties();
        reader.openFile(input);
        
        // get arguments:
        if ( arg.set(frame, "frame") )
            period = 0;
        arg.set(period, "period");

        if ( arg.set(str, "output") )
        {
            try {
                ofs.open(str.c_str());
            }
            catch( ... )
            {
                std::clog << "Cannot open output file\n";
                return EXIT_FAILURE;
            }
            out.rdbuf(ofs.rdbuf());
        }
    }
    catch( Exception & e )
    {
        std::clog << e.brief() << '\n';
        return EXIT_FAILURE;
    }
    
    Cytosim::silent();
    
    // process first record, at index 'frame':
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

    report(simul, out, what, frame, arg);
    size_t cnt = 1;

    if ( arg.num_values("frame") > 1 )
    {
        // multiple record indices were specified:
        size_t s = 1;
        while ( arg.set(frame, "frame", s) )
        {
            // try to load the specified frame:
            if ( 0 == reader.loadFrame(simul, frame) )
            {
                report(simul, out, what, frame, arg);
                ++cnt;
            }
            else
            {
                std::cerr << "Error: missing frame " << frame << '\n';
                return EXIT_FAILURE;
            }
            ++s;
        }
    }
    else if ( period > 0 )
    {
        // process every 'period' record:
        size_t f = frame;
        while ( 0 == reader.loadNextFrame(simul)  )
        {
            ++f;
            if ( f % period == frame % period )
            {
                report(simul, out, what, f, arg);
                ++cnt;
            }
        }
    }
    
    out << '\n';
    arg.print_warnings(std::cerr, cnt, "\n");
}
