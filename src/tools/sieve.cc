// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "simul.h"
#include "parser.h"
#include "glossary.h"
#include "iowrapper.h"
#include "exceptions.h"
#include "simul_prop.h"


void help()
{
    printf("Cytosim-sieve:\n\n");
    printf("   Sieve reads a cytosim trajectory file, loading frames into memory,\n");
    printf("   and writes them to a new file using the latest cytosim format.\n");
    printf("   The output can be generated in either binary or text format.\n");
    printf("   A category of objects can be removed by specifying `skip=WHAT`.\n");
    printf("   If the specified output file already exists, data is appended to it.\n");
    printf("   This writes %iD files with format %i (real = %lu bytes)\n", DIM, Simul::currentFormatID, sizeof(real));
    printf("   (built on %s)\n\n", __DATE__);
    printf("Usage:\n");
    printf("    sieve input_file output_file [options]\n\n");
    printf("Options:\n");
    printf("    dim=INT            process files with specified dimensionality\n");
    printf("    binary=1           use binary format (default)\n");
    printf("    binary=0           use text format (default if output ends with .txt\n");
    printf("    skip=WHAT          remove all objects of class WHAT\n");
    printf("    skip_free_single=1 remove unbound singles\n");
    printf("    skip_free_couple=1 remove unbound couples\n");
    printf("    frame=INDEX        process only specified frame\n\n");
    printf("Examples:\n");
    printf("    sieve objects.cmo objects.txt\n");
    printf("    sieve objects.cmo objects.txt skip=couple\n");
}


int main(int argc, char* argv[])
{
    unsigned dim = DIM;
    Simul simul;
    Glossary arg;
    bool binary = true;
    ObjectSet * skip_set = nullptr;
    std::string skip;

    if ( argc < 3 )
    {
        help();
        return EXIT_SUCCESS;
    }
    
    // check extension of output file:
    char const* ext = strrchr(argv[2], '.');
    if ( ext )
    {
        if ( 0 == strncmp(ext, ".txt", 4) )
            binary = 0;
    }
    
    if ( arg.read_strings(argc-3, argv+3) )
        return EXIT_FAILURE;
    
    if ( arg.set(skip, "skip") )
       skip_set = simul.findSet(skip);
    
    arg.set(dim, "dim");
    arg.set(binary, "binary");
    arg.set(simul.prop->skip_free_single, "skip_free_single");
    arg.set(simul.prop->skip_free_couple, "skip_free_couple");

    Inputter in(DIM);
    try {
        simul.loadProperties();
        in.open(argv[1], "rb");
    }
    catch( Exception & e ) {
        std::cerr << e.brief() << '\n';
        return EXIT_FAILURE;
    }
    
    std::clog << ">>>>>> Sieve `" << argv[1] << "' -> `" << argv[2] << "'\n";
    
    size_t frm = 0, frame = 0;
    
    // a frame index can be specified:
    bool has_frame = arg.set(frame, "frame");
    
    while ( in.good() )
    {
        try {
            if ( simul.reloadObjects(in) )
                return EXIT_SUCCESS;
        }
        catch( Exception & e ) {
            std::cerr << "Frame " << frm << ":" << e.brief() << '\n';
        }
        
        if ( in.vectorSize() != dim )
        {
            std::cerr << "Abort: dimensionality mismatch (file is " << in.vectorSize() << "D)\n";
            return 1;
        }

        if ( skip_set )
            skip_set->erase();
        
        /*
        simul.reportInventory(std::cout);
        std::clog << "\r" << std::setw(5) << cnt << "   ";
        */
        
        try {
            if ( !has_frame || ( frm == frame ) )
                simul.writeObjects(argv[2], true, binary);
        }
        catch( Exception & e ) {
            std::cerr << e.brief() << '\n';
            return EXIT_FAILURE;
        }
        if ( has_frame && frm == frame )
            break;
        ++frm;
    }
}
