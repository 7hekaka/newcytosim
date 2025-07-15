// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University

#include "exceptions.h"
#include "glossary.h"


int main(int argc, char* argv[])
{
    Glossary arg;
    std::string s;

    if ( arg.read_strings(argc-1, argv+1) )
        return EXIT_FAILURE;

    try {
        
        // read file (recognized by its extension) if provided on command line:
        if ( arg.set(s, ".cym") )
            arg.read_file(s);
        
        // print content of Glossary:
        printf("%lu keys:\n", arg.num_keys());
        arg.print(std::cout, "    > ");
        
        // extract values from Glossary:
        int i = 0;
        float f = 0;

        if ( arg.set(f, "f", "float") ) printf("float : %f\n", f);
        if ( arg.set(s, "s", "string") ) printf("string : %s\n", s.c_str());
        if ( arg.set(i, "i", "integer") ) printf("integer : %i\n", i);
    }
    catch ( Exception& e )
    {
        std::cout << e.brief() << '\n';
        return EXIT_FAILURE;
    }
}
