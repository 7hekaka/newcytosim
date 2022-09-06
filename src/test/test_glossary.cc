// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University

#include "exceptions.h"
#include "glossary.h"


int main(int argc, char* argv[])
{
    Glossary arg;
    
    if ( arg.read_strings(argc-1, argv+1) )
        return EXIT_FAILURE;

    try {
        
        int i = 0;
        float f = 0;
        std::string s;
        
        // read file if provided on command line:
        // the file name is recognized by its extension
        if ( arg.set(s, ".cym") )
            arg.read_file(s);
        
        // print content of Glossary:
        printf("%lu keys:\n", arg.num_keys());
        arg.print(std::cout, "    > ");
        
        // extract values from Glossary:
        if ( arg.set(i, "integer") )  printf("integer : %i\n", i);
        if ( arg.set(f, "float") )    printf("float : %f\n", f);
        if ( arg.set(s, "string") )   printf("string : %s\n", s.c_str());
    }
    catch ( Exception& e )
    {
        std::cout << e.brief() << '\n';
        return EXIT_FAILURE;
    }
}
