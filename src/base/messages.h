// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University

#ifndef  MESSAGES_H
#define  MESSAGES_H

#include <iostream>
#include <sstream>
#include <fstream>


/// This facility provides some control over output
/** F. Nedelec, 16.03.2018 */
namespace Cytosim
{
    /// a class representing an output stream
    class Output
    {
        /// prefix to all messages
        std::string pref_;

        /// pointer to the current destination of output
        std::ostream out_;
        
        /// file stream, if open
        std::ofstream ofs_;

        /// alias to /dev/null
        std::ofstream nul_;

        /// number of output events still permitted
        size_t cnt_;

    public:
        
        /// create stream directed to given stream with `max_output` allowed
        Output(std::ostream& os, size_t sup = 0x1p31, std::string const& p = "")
        : pref_(p), out_(std::cout.rdbuf()), cnt_(sup)
        {
            out_.rdbuf(os.rdbuf());
            nul_.open("/dev/null");
        }
        
        /// redirect output to given file
        void open(std::string const& filename)
        {
            ofs_.open(filename.c_str());
            out_.rdbuf(ofs_.rdbuf());
        }
        
        /// close file
        void close()
        {
            if ( ofs_.is_open() )
                ofs_.close();
            out_.rdbuf(std::cout.rdbuf());
        }
        
        /// return current output
        operator std::ostream&()
        {
            return out_;
        }
        
        /// true if output is /dev/null
        bool is_silent()
        {
            return out_.rdbuf() == nul_.rdbuf();
        }

        /// direct output to /dev/null
        void silent()
        {
            out_.rdbuf(nul_.rdbuf());
        }
        
        /// direct output to given stream
        void redirect(Output const& x)
        {
            out_.rdbuf(x.out_.rdbuf());
        }
        
        /// std::ostream style output operator
        template < typename T >
        std::ostream& operator << (T const& x)
        {
            if ( out_.good() && cnt_ )
            {
                --cnt_;
                if ( pref_.size() > 0 )
                    out_.write(pref_.data(), pref_.size());
                out_ << x;
                return out_;
            }
            return nul_;
        }
        
        /// output string with flush
        void operator()(const char* str)
        {
            operator<<(str);
            out_.flush();
        }

        /// C-style `printf()` syntax followed by flush
        template < typename Arg1, typename... Args >
        void operator()(const char* fmt, Arg1 arg1, Args&&... args)
        {
            char str[1024] = { 0 };
            snprintf(str, sizeof(str), fmt, arg1, args...);
            operator<<(str);
            out_.flush();
        }

    };
    
    /// for normal output
    extern Output out;

    /// for logging
    extern Output log;

    /// for warnings
    extern Output warn;

    /// suppress all output
    void silent();
}


/// a macro to print some text only once
#define LOG_ONCE(a) { static bool V=1; if (V) { V=0; Cytosim::log << a; } }

#endif
