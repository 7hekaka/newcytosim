// Cytosim was created by Francois Nedelec. Copyright 2023 Cambridge University

#include <iostream>
#include "assert_macro.h"
#include "hand_prop.h"
#include "modulo.h"
#include "fiber.h"
#include "real.h"
#include "dim.h"


/// print header line identifying the project
inline void splash(std::ostream& os)
{
    os << " ------------------------------------------------------------- \n";
    os << "|  CytoSIM " <<DIM<<"D  -  www.cytosim.org  -  version PI  -  02.2022  |\n";
    os << " ------------------------------------------------------------- \n";
}


/// print general info about the program
inline void print_version(std::ostream& os)
{
    os << "   Dimension " << DIM;
    os << "  Periodic " << ENABLE_PERIODIC_BOUNDARIES;
    os << "  AttachPool " << POOL_UNATTACHED;
    os << "  Precision " << sizeof(real) << "\n";

    os << "   Fiber has lattice " << FIBER_HAS_LATTICE;
    os << " density " << FIBER_HAS_DENSITY;
    os << " family " << FIBER_HAS_FAMILY;
    os << " glue " << FIBER_HAS_GLUE;
    os << " bind_closest " << BIND_CLOSEST_FIBER << "\n";

    os << "   Built " <<__DATE__<< " " <<__TIME__<< " with " <<__VERSION__<< "\n";
    
#ifdef CODE_VERSION
    os << "   Code version " << CODE_VERSION;
#else
    os << "   Code version unknown";
#endif
    
#ifdef NDEBUG
    os << " (no assertions)\n";
#else
    os << " with assertions\n";
#endif
}

