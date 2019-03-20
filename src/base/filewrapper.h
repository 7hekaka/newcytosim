// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
///F. Nedelec, EMBL, October 2006. nedelec@embl.de

#ifndef  FILEWRAPPER_H
#define  FILEWRAPPER_H

/**
 The keyword _FILE_OFFSET_BITS affects the code in <cstdio> etc.
 Defining this keywords allow us to use 64bits fpos_t, 
 and thus to support files above 2GB
 */
#define _FILE_OFFSET_BITS  64


#include "assert_macro.h"
#include <cstdio>
#include <sys/types.h>
#include <string>


/// A simple wrapper around a C-FILE
/**
 The FileWrapper has a cast-operator to FILE*,
 and it can thus be used directly in the functions of the C-library.
 */
class FileWrapper
{    
protected:
    
    /// the C-file descriptor
    FILE*       mFile;
    
    /// the name of the file or some other information:
    std::string mPath;
    
public:
    
    /// constructor - no file
    explicit FileWrapper();
    
    /// constructor which opens a file
    FileWrapper(FILE* , const char * path = nullptr);
    
    /// constructor which opens a file
    FileWrapper(const char* name, const char* mode);

    /// destructor 
    virtual ~FileWrapper();
    
    /// constructor from an already opened file
    void operator =(FILE *);
    
    /// automatic conversion to a FILE *
    operator FILE*()                     { return mFile; }
    
    /// open a file
    int     open(const char* name, const char* mode);
    
    /// rewind file
    void    rewind()                     { if ( mFile ) { std::clearerr(mFile); std::rewind(mFile); } }

    /// rewind file
    void    clearerr()                   { if ( mFile ) std::clearerr(mFile); }

    /// close file
    void    close();
    
    /// return the file pointer
    FILE*   file()                       { return mFile; }
    
    /// the path of the file, or of the last attempt to open a file
    const char * path()          const   { return mPath.c_str(); }
    
    /// true if output goes to stdout
    bool    std()               const    { return mFile==stdout; }
    
    /// true if at end of file
    bool    eof()               const    { return mFile && std::feof(mFile); }
    
    /// return the value of ferror()
    int     error()             const    { return std::ferror(mFile); }

    /// true if file is good for writing / reading
    bool    good()              const    { return mFile && !std::ferror(mFile); }

    /// set `p` to current reading position of file, relative to beggining
    long    pos()               const    { if ( mFile ) return ftell(mFile); else return 0; }

    /// set `p` to current reading position of file, relative to beggining
    int     get_pos(fpos_t& p)  const    { return fgetpos(mFile, &p); }

    /// change current reading position to `p` specified from the beginning of the file
    void    set_pos(const fpos_t& p)     { fsetpos(mFile, &p); }
    
    
    /// put a C-string
    void    put_line(const char * line, bool end = 0);

    /// put a C++ string
    void    put_line(const std::string&, bool end = 0);

    /// read until character `end` is found and set `line`, including terminating character
    void    get_line(std::string& line, char end='\n');
    
    /// read stream until given string is found
    void    skip_until(const char * str);
    
    /// lock file for current thread
    void    lock()                       { flockfile(mFile); }
    
    /// unlock file for current thread
    void    unlock()                     { funlockfile(mFile); }
    
    /// report next character to be read
    int     peek()                       { int c=getc_unlocked(mFile); if ( c != EOF ) ungetc(c, mFile); return c; }
    
    /// read a character
    int     get_char()                   { return getc_unlocked(mFile); }

    /// read a byte
    uint8_t get_byte()                   { return (uint8_t)getc_unlocked(mFile); }

    /// unget character from input
    void    unget(int c)                 { ungetc(c, mFile); }

    /// write a character
    int     put_char(int c)              { return putc_unlocked(c, mFile); }

    /// flush
    void    flush()                      { std::fflush(mFile); }

};

#endif

