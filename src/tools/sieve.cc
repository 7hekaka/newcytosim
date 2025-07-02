// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "simul.h"
#include "parser.h"
#include "glossary.h"
#include "iowrapper.h"
#include "exceptions.h"
#include "simul_prop.h"


void help(std::ostream& os)
{
    os << "Cytosim-sieve:\n";
    os << "   Sieve reads a cytosim trajectory file, loading frames into memory,\n";
    os << "   and writes them to a new or existing file using the latest cytosim format.\n";
    os << "   The output can be generated in either binary or text format.\n";
    os << "   A category of objects can be removed by specifying `skip=WHAT`.\n";
    os << "   If the specified output file already exists, data is appended to it.\n";
    os << "Usage:\n";
    os << "    sieve input_file output_file [options]\n\n";
    os << "Options:\n";
    os << "    dim=INT                  process files with specified dimensionality\n";
    os << "    binary=1                 use binary format (default)\n";
    os << "    binary=0                 use text format (default if output ends with .txt\n";
    os << "    skip=WHAT                remove all objects of class WHAT\n";
    os << "    simul.skip_free_single=1 do not read unbound singles\n";
    os << "    simul.skip_free_couple=1 do not read unbound couples\n";
    os << "    frame=INDEX              process only specified frame\n\n";
    os << "Examples:\n";
    os << "    sieve objects.cmo objects.txt\n";
    os << "    sieve objects.cmo objects.txt skip=couple\n";
    os << "Made with format version " << Simul::currentFormatID << " and DIM=" << DIM << "\n";
}


int main(int argc, char* argv[])
{
    int binary = 1;
    unsigned dim = DIM;
    std::string skip;
    Simul simul;
    Glossary arg;
    ObjectSet * skip_set = nullptr;

    if ( argc < 3 )
    {
        help(std::cout);
        return EXIT_SUCCESS;
    }
    
    // check extension of output file:
    char const* ext = strrchr(argv[2], '.');
    if ( ext && 0 == memcmp(ext, ".txt", 4) )
        binary = 0;
    
    if ( arg.read_strings(argc-3, argv+3) )
        return EXIT_FAILURE;
    
    if ( arg.set(skip, "skip") )
       skip_set = simul.findSet(skip);
    
    arg.set(dim, "dim");
    arg.set(binary, "binary");

    Inputter in(DIM);
    try {
        simul.loadProperties();
        in.open(argv[1], "rb");
    }
    catch( Exception & e ) {
        std::cerr << e.brief() << '\n';
        return EXIT_FAILURE;
    }
    
    // read parameters from command line
    for ( int n = 3; n < argc; ++n )
    {
        if ( simul.readParameter(argv[n]) )
            std::clog << ">>>>>> set " << argv[n] << "\n";
    }
    
    std::clog << ">>>>>> Sieve `" << argv[1] << "' ---> `" << argv[2];
    std::clog << "'  binary: " << binary << "\n";
    
    size_t frm = 0, frame = ~0;
    // the index for a single frame can be specified:
    bool has_frame = arg.set(frame, "frame");
    
    while ( in.good() && frm <= frame )
    {
        try {
            if ( simul.reloadObjects(in) )
                return EXIT_SUCCESS;
        }
        catch( Exception & e ) {
            std::cerr << "Frame " << frm << ":" << e.brief() << '\n';
        }
        
        if ( has_frame && frm < frame )
            continue;
        
        if ( in.vectorSize() != dim )
        {
            std::cerr << "Aborting due to dimensionality mismatch (file is " << in.vectorSize() << "D)\n";
            return 1;
        }

        if ( skip_set )
            skip_set->erase();
        
        //simul.reportInventory(std::cout);
        //std::clog << "\r" << std::setw(5) << frm << "   ";
        
        try {
            simul.writeObjects(argv[2], true, binary);
        }
        catch( Exception & e ) {
            std::cerr << e.brief() << '\n';
            return EXIT_FAILURE;
        }
        ++frm;
    }
}
