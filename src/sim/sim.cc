// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include <csignal>

#include "simul.h"
#include "parser.h"
#include "messages.h"
#include "glossary.h"
#include "exceptions.h"
#include "print_color.h"
#include "backtrace.h"
#include "filepath.h"
#include "splash.h"
#include "tictoc.h"
#include "unistd.h"


void help(std::ostream& os)
{
    os << "sim [OPTIONS] [FILE]\n";
    os << "  FILE    run specified config file, if ending with `.cym'\n";
    os << "  +       redirect output to terminal instead of `messages.cmo'\n";
    os << "  -       suppress output\n";
    os << "  info    print build options\n";
    os << "  help    print this message\n";
}


void handle_signal(int sig)
{
    /*
     A signal handler is restricted to call only async-signal-safe-functions
     practically speaking, most syscalls(2) only
     */
    char str[128] = { 0 };
    strncpy(str, "\nCytosim received signal   \n", 128);
    str[26] = (char)('0' + ( sig     % 10));
    str[25] = (char)('0' + ((sig/10) % 10));
    (void) write(STDERR_FILENO, str, 29);
    print_backtrace();
    abort();
    //_exit(sig);
}

void handle_abort(int sig)
{
    // this will prevent crash reports to be generated
    write(STDERR_FILENO, "\nCytosim catched abort\n", 23);
    _exit(sig);
}

void handle_interrupt(int sig)
{
    Cytosim::out << "killed " << sig << "\n" << std::endl;
    _exit(sig);
}

//------------------------------------------------------------------------------
//=================================  MAIN  =====================================
//------------------------------------------------------------------------------

/// for normal output:
using std::cout;
/// where errors are printed:
using std::cerr;

int main(int argc, char* argv[])
{
    // register callback to catch interrupting signals:
    if ( signal(SIGINT, handle_interrupt) )
        cerr << "Could not register SIGINT handler\n";
    if ( signal(SIGTERM, handle_interrupt) )
        cerr << "Could not register SIGTERM handler\n";
#if 0
    if ( signal(SIGSEGV, handle_signal) )
        cerr << "Could not register SIGSEGV handler\n";
    if ( signal(SIGILL,  handle_signal) )
        cerr << "Could not register SIGILL handler\n";
    if ( signal(SIGABRT, handle_abort) )
        cerr << "Could not register SIGABRT handler\n";
#endif
    // catch division by zero and Nan
    //feenableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);

    Glossary arg;

    //parse the command line:
    if ( arg.read_strings(argc-1, argv+1) )
        return 1;

    if ( arg.use_key("help") || arg.use_key("--help") )
    {
        splash(cout);
        help(cout);
        return 0;
    }

    if ( arg.use_key("info") || arg.use_key("--version")  )
    {
        splash(cout);
        print_version(cout);
        return 0;
    }
    
    if ( arg.use_key("-") )
    {
        Cytosim::out.silent();
        Cytosim::log.silent();
        Cytosim::warn.silent();
    }
    else if ( ! arg.use_key("+") )
    {
        Cytosim::out.open("messages.cmo");
        Cytosim::log.redirect(Cytosim::out);
        Cytosim::warn.redirect(Cytosim::out);
    }

    // change working directory if specified:
    if ( arg.has_key("directory") )
    {
        FilePath::change_dir(arg.value("directory"));
        std::clog << "Cytosim working directory is " << FilePath::get_cwd() << '\n';
    }

    Cytosim::out << "% " << TicToc::date() << "\n";
    print_version(Cytosim::out);
    
    Simul simul;
    try {
        simul.initialize(arg);
    }
    catch( Exception & e ) {
        print_magenta(cerr, e.brief());
        cerr << e.info() << '\n';
        return 2;
    }
    catch(...) {
        print_red(cerr, "Error: an exception occurred during initialization\n");
        return 3;
    }
    
    arg.print_warnings(cerr, 1, " on command line\n");
    time_t sec = TicToc::seconds_since_1970();
    
    try {
        Parser(simul, 1, 1, 1, 1, 1).readConfig();
    }
    catch( Exception & e ) {
        print_magenta(cerr, e.brief());
        cerr << e.info() << '\n';
        return 4;
    }
    catch(...) {
        print_red(cerr, "Error: an exception occurred\n");
        return 5;
    }
    
    Cytosim::out << "% " << TicToc::date() << "\n";
    sec = TicToc::seconds_since_1970() - sec;
    Cytosim::out << "end  " << sec << " s ( " << (float)(sec)/3600 << " h )\n";
    Cytosim::out.close();
}
