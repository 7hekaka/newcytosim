// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University

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
