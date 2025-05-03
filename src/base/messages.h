// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University

#ifndef  MESSAGES_H
#define  MESSAGES_H

#include <iostream>
#include <sstream>
#include <fstream>


/// create std::string using the C-style `printf()` syntax
template < typename Arg1, typename... Args >
std::string format(const char* fmt, Arg1 arg1, Args&&... args)
{
    char tmp[1024] = { 0 };
    static_assert(std::is_trivial<Arg1>::value, "non-trivial type");
    snprintf(tmp, sizeof(tmp), fmt, arg1, args...);
    return (tmp);
}


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
        
        /// write `s`
        void operator()(const std::string& arg)
        {
            operator<<(arg);
        }

        /// write `s` followed by `a`
        template <typename A>
        void operator()(const std::string& s, const A& a)
        {
            std::ostringstream oss;
            oss << s << a;
            operator<<(oss.str());
        }
        
        /// write `s` followed by `a` and `b`
        template <typename A, typename B>
        void operator()(const std::string& s, const A& a, const B& b)
        {
            std::ostringstream oss;
            oss << s << a << b;
            operator<<(oss.str());
        }
        
        /// write `s` followed by `a`, `b` and `c`
        template <typename A, typename B, typename C>
        void operator()(const std::string& s, const A& a, const B& b, const C& c)
        {
            std::ostringstream oss;
            oss << s << a << b << c;
            operator<<(oss.str());
        }
        
        /// write `s` followed by `a`, `b`, `c` and `d`
        template <typename A, typename B, typename C, typename D>
        void operator()(const std::string& s, const A& a, const B& b, const C& c, const D& d)
        {
            std::ostringstream oss;
            oss << s << a << b << c << d;
            operator<<(oss.str());
        }
        
        /// write `s` followed by `a`, `b`, `c`, `d` and `e`
        template <typename A, typename B, typename C, typename D, typename E>
        void operator()(const std::string& s, const A& a, const B& b, const C& c, const D& d, const E& e)
        {
            std::ostringstream oss;
            oss << s << a << b << c << d << e;
            operator<<(oss.str());
        }
        
        /// write `s` followed by `a`, `b`, `c`, `d`, `e` and `f`
        template <typename A, typename B, typename C, typename D, typename E, typename F>
        void operator()(const std::string& s, const A& a, const B& b, const C& c, const D& d, const E& e, const F& f)
        {
            std::ostringstream oss;
            oss << s << a << b << c << d << e << f;
            operator<<(oss.str());
        }
        
        /// write `s` followed by `a`, `b`, `c`, `d`, `e`, `f` and `g`
        template <typename A, typename B, typename C, typename D, typename E, typename F, typename G>
        void operator()(const std::string& s, const A& a, const B& b, const C& c, const D& d, const E& e, const F& f, const G& g)
        {
            std::ostringstream oss;
            oss << s << a << b << c << d << e << f << g;
            operator<<(oss.str());
        }
        
        /// write `s` followed by `a`, `b`, `c`, `d`, `e`, `f`, `g` and `h`
        template <typename A, typename B, typename C, typename D, typename E, typename F, typename G, typename H>
        void operator()(const std::string& s, const A& a, const B& b, const C& c, const D& d, const E& e, const F& f, const G& g, const H& h)
        {
            std::ostringstream oss;
            oss << s << a << b << c << d << e << f << g << h;
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
