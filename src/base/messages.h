// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University

#ifndef  MESSAGES_H
#define  MESSAGES_H

#include <iostream>
#include <sstream>
#include <fstream>


/// Different Output instances provide some control over output
namespace Cytosim
{
    /// a class representing an output stream
    class Output
    {
        /// prefix to all messages
        std::string pref_;
        
        /// the last message that was printed
        std::string last_;

        /// destination of output
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
        
        /// destructor closes the file
        ~Output()
        {
            close();
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
        
        /// output string if different from last line, and flush
        std::ostream& operator << (std::string const& arg)
        {
            //std::clog << "[" << last_ << "][" << arg << "]\n";
            if ( last_ != arg )
            {
                last_ = arg;
                if ( out_.good() && cnt_ )
                {
                    --cnt_;
                    out_ << pref_ << arg;
                    out_.flush();
                    return out_;
                }
                return nul_;
            }
            return out_;
        }
        
        template <typename Arg1>
        void stringify(std::ostringstream& oss, Arg1 arg1)
        {
            oss << arg1;
        }

        template <typename Arg1, typename... Args>
        void stringify(std::ostringstream& oss, Arg1 arg1, Args... args)
        {
            oss << arg1;
            stringify(oss, args...);
        }

        
        /// write `s`
        void operator()(const std::string& arg)
        {
            operator<<(arg);
        }

        template <typename Arg1, typename... Args>
        void operator()(Arg1 arg1, Args... args)
        {
            std::ostringstream oss;
            stringify(oss, arg1, args...);
            operator<<(oss.str());
        }

        /// C-style `printf()` syntax followed by flush
        template < typename Arg1, typename... Args >
        void print(const char* fmt, Arg1 arg1, Args&&... args)
        {
            char tmp[1024] = { 0 };
            snprintf(tmp, sizeof(tmp), fmt, arg1, args...);
            operator<<(tmp);
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
