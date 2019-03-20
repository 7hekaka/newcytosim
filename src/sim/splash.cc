// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "splash.h"
#include "assert_macro.h"
#include "real.h"


void splash(std::ostream& os)
{
    os << " ------------------------------------------------------------- \n";
    os << "|  CytoSIM  -  www.cytosim.org  -  version PI  -  Jan  2019   |\n";
    os << " ------------------------------------------------------------- \n";
}


void print_version(std::ostream& os)
{
    os << "    Precision: " << sizeof(real) << " bytes, " << REAL_EPSILON << "\n";
    os << "    Built on " <<__DATE__<< " at " <<__TIME__<< "\n";
    
#ifdef COMPILER_VERSION
    os << "    with " << COMPILER_VERSION << "\n";
#else
    os << "    with unknown compiler\n";
#endif
    os << "    C++ version " << __cplusplus << "\n";

#ifdef CODE_VERSION
    os << "    Code version " << CODE_VERSION << "\n";
#else
    os << "    Code version unknown\n";
#endif
    
#ifdef NDEBUG
    os << "    (no assertions)\n";
#else
    os << "    with assertions\n";
#endif
}
