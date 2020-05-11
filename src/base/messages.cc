// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "messages.h"

namespace Cytosim
{
    /// alias to standard output
    Output out(std::cout);
    
    /// for logs
    Output log(std::clog);
    
    /// for warnings
    Output warn(std::cerr, 32U, "WARNING: ");
    
    /// turn all output off
    void silence_all()
    {
        Cytosim::out.silent();
        Cytosim::log.silent();
        Cytosim::warn.silent();
    }
}
