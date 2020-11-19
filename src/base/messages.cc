// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "messages.h"

namespace Cytosim
{
    /// alias to standard output
    Output out(std::cout);
    
    /// alias to standard log
    Output log(std::clog);
    
    /// alias to standard error
    Output warn(std::cerr, 32U, "WARNING: ");
    
    /// supress all output
    void silent()
    {
        Cytosim::out.silent();
        Cytosim::log.silent();
        Cytosim::warn.silent();
    }
}
