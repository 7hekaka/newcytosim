// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
//Some error conditions are handled by throwing exceptions.
//here we define a very primite Exception class for cytosim

#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H

#include "assert_macro.h"
#include <string>
#include <sstream>


/// A mechanism to handle errors (see C++ manual)
/**
 The exception carry a 'message' and associated 'info', which are both strings.
 The message is set by the constructor, and the info is set by the << operator.
 
 Usage: Throw an Exception (not a pointer), and catch a reference to an exception.
 This ensures proper memory managment (coordinated calls of constructor / destructor)
*/
class Exception 
{
protected:
    
    /// brief description of the issue
    std::string msg_;

    /// background information
    std::string info_;
    
    /// file in which the exception occurred
    std::string file_;

public:
    
    /// Creator with empty message
    Exception()
    {
    }
    
    /// constructor with given message
    Exception(std::string const& m)
    {
        msg_ = m;
        //printf("Exception(%s)\n", msg.c_str());
    }
    
    /// return the message
    std::string brief() const
    {
        if ( file_.empty() )
            return "Error, " + msg_;
        else
            return "Error in `" + file_ + "' : " + msg_;
    }
    
    /// return size of supplementary message
    size_t has_info() const
    {
        return info_.size();
    }
    
    /// return supplementary message
    std::string info() const
    {
        if ( info_.size() > 0 )
             return ": " + info_;
        return "\n";
    }

    /// return C-string pointer to message
    char const* what() const
    {
        return msg_.c_str();
    }
    
    /// change the message
    std::string message() const
    {
        return msg_;
    }

    /// change the message
    void message(const std::string& m)
    {
        msg_ = m;
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

    template <typename... Args>
    Exception(const std::string& s, Args... args)
    {
        std::ostringstream oss(s);
        stringify(oss, args...);
        msg_ = oss.str();
    }
    
    /// define file
    void set_file(const std::string& f)
    {
        file_ = f;
    }

    /// append string to info
    Exception& operator << (const std::string& arg)
    {
        info_.append(arg);
        return *this;
    }
    
    /// append C-string to info
    Exception& operator << (const char arg[])
    {
        info_.append(arg);
        return *this;
    }

    /// append string-representation of `x` to info
    template<typename T>
    Exception& operator << (const T& x)
    {
        std::ostringstream oss;
        oss << x;
        *this << oss.str();
        return *this;
    }
};


//------------------------------------------------------------------------------
/// This class is thrown if a parameter value is invalid
class InvalidParameter : public Exception 
{
    
public:
    
    /// constructor
    InvalidParameter() : Exception()
    {
        //printf("new InvalidParameter [%s]\n", m.c_str());
    }
    
    /// constructor
    InvalidParameter(const std::string m) : Exception(m)
    {
        //printf("new InvalidParameter [%s]\n", m.c_str());
    }
    
    template <typename... Args>
    InvalidParameter(const std::string& s, Args... args) : Exception(s, args...) {}
};


//------------------------------------------------------------------------------
/// InvalidSyntax is thrown while parsing config file
class InvalidSyntax : public Exception 
{
    
public :
    
    /// constructor
    InvalidSyntax(std::string const& m) : Exception(m)
    {
        //printf("new InvalidSyntax [%s]\n", m.c_str());
    }
};

//------------------------------------------------------------------------------
/// InvalidIO is thrown during file Input/Output
class InvalidIO : public Exception 
{
    
public :
    
    /// constructor
    InvalidIO(const std::string m) : Exception(m)
    {
        //printf("new InvalidIO [%s]\n", m.c_str());
    }
};


#endif
