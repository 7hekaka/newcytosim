// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include <csignal>
#include <unistd.h>

#include "backtrace.h"
#include "signal_handlers.h"
#include "print_color.h"
#include "exceptions.h"
#include "time_date.h"
#include "filepath.h"
#include "glossary.h"
#include "messages.h"
#include "simul.h"
#include "parser.h"
#include "splash.h"

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
    ssize_t __attribute__((unused)) u;
    u = write(STDERR_FILENO, str, 29);
    print_backtrace();
    abort();
    //_exit(sig);
}

void handle_abort(int sig)
{
    // this will prevent crash reports to be generated
    ssize_t __attribute__((unused)) u;
    u = write(STDERR_FILENO, "\nCytosim catched abort\n", 23);
    _exit(sig);
}

void handle_interrupt(int sig)
{
    Cytosim::out("killed ", sig, "\n");
    _exit(sig);
}

//------------------------------------------------------------------------------

int main(int argc, char* argv[])
{
    register_signal_handlers();
    // register callback to catch interrupting signals:
    if ( signal(SIGINT, handle_interrupt) )
        std::cerr << "Could not register SIGINT handler\n";
    if ( signal(SIGTERM, handle_interrupt) )
        std::cerr << "Could not register SIGTERM handler\n";
    // catch division by zero and Nan
    //feenableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);

    //parse the command line:
    Glossary arg;
    if ( arg.read_strings(argc-1, argv+1) )
        return 1;

    if ( arg.use_key("help") || arg.use_key("--help") )
    {
        splash(std::cout);
        help(std::cout);
        return 0;
    }

    if ( arg.use_key("info") || arg.use_key("--version")  )
    {
        splash(std::cout);
        print_version(std::cout);
        return 0;
    }
    
    if ( arg.use_key("-") )
    {
        Cytosim::out.silent();
        Cytosim::log.silent();
        Cytosim::warn.silent();
    }
    else if ( arg.use_key("+") )
    {
        Cytosim::log.redirect(std::cerr);
        Cytosim::warn.redirect(std::cerr);
    }
    else
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

    Cytosim::out("% ", TimeDate::date_string(), " PID ", getpid(), '\n');
    print_version(Cytosim::out);
    
    Simul simul;
    try {
        simul.prop.read(arg);
        arg.print_warnings(stderr, 1, " on command line\n");
        simul.initCytosim();
        time_t sec = TimeDate::seconds_since_1970();
        // read and execute default config file:
        Parser(&simul, 1, 1, 1, 1, 1, 1).readConfig();
        Cytosim::out("% ", TimeDate::date_string(), '\n');
        sec = TimeDate::seconds_since_1970() - sec;
        Cytosim::out << "end  " << sec << " s ( " << (float)(sec)/3600 << " h )\n";
    }
    catch( Exception & e ) {
        print_magenta(stderr, e.brief());
        fputs(e.info().c_str(), stderr);
        return 4;
    }
    catch(...) {
        print_red(stderr, "Error: an exception occurred\n");
        return 5;
    }
}
