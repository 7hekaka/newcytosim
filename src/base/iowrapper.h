// Cytosim was created by Francois Nedelec. Copyright Cambridge University 2020

#ifndef  IOWRAPPER_H
#define  IOWRAPPER_H

#include <cstdio>
#include <stdint.h>
#include "filewrapper.h"

/// Input in text or binary mode with automatic byte-swapping for endianess compatibility
class Inputter : public FileWrapper
{
private:
        
    /// The format ID of the input: this allow backward compatibility with older formats
    unsigned format_;
    
    /// The dimensionality of vectors stored in the file
    unsigned vecsize_;
    
    /** if the state is stored in a binary format, binary_
        is set to 1 or 2. with 2, byte order is swapped automatically
        this occurs for example when reading a simulation calculated 
        on PC from mac, or vice et versa.
        */
    int      binary_;

public:
    
    /// set defaults (not-binary)
    void     reset();
    
    /// Constructor
    Inputter(unsigned d) : FileWrapper(nullptr) { reset(); vecsize_=d; }
    
    /// Constructor
    Inputter(unsigned d, FILE* f, const char* path=nullptr) : FileWrapper(f, path) { reset(); vecsize_=d; }
    
    /// constructor which opens a file
    Inputter(unsigned d, const char* name, bool bin) : FileWrapper(name, bin?"rb":"r") { reset(); vecsize_=d; }

    /// return dimensionnally of vectors
    unsigned vectorSize()     const { return vecsize_; }
    
    /// Set dimentionnality of vectors
    void     vectorSize(unsigned d) { vecsize_ = d; }
    
    /// returns the type of input
    unsigned formatID()       const { return format_; }

    /// returns the type of input
    void     formatID(unsigned f)   { format_ = f; }

    /// Returns 1 for native binary format, 2 for non-native binary format, and 0 if not binary
    int      binary()         const { return binary_; }
    
    /// initialize the automatic swapping of bytes in the binary format
    void     setEndianess(const char[2]);
    
    /// Read integer on 2 bytes
    int16_t  readInt16();
    /// Read integer on 4 bytes
    int32_t  readInt32();

    /// Read unsigned integer on 1 byte
    uint8_t  readUInt8();
    /// Read unsigned integer on 2 bytes
    uint16_t readUInt16();
    /// Read unsigned integer on 4 bytes
    uint32_t readUInt32();
    /// Read unsigned integer on 8 bytes
    uint64_t readUInt64();
    
    /// Reads one float on 4 bytes
    float    readFloat();
    /// Reads one double on 8 bytes
    double   readDouble();
    
    /// Reads one vector, returning D coordinates in the array of size D
    void     readFloats(float[], unsigned D);
    /// Reads one vector, returning D coordinates in the array of size D
    void     readFloats(double[], unsigned D);
    
    /// Reads `n` vector, returning D coordinates for each, in the array of size n*D
    void     readFloats(float[], size_t n, unsigned D);
    /// Reads `n` vector, returning D coordinates for each, in the array of size n*D
    void     readFloats(double[], size_t n, unsigned D);

};


#pragma mark -


/// Output in text or binary mode in the native endianess
class Outputter : public FileWrapper
{
    
private:
        
    /// Flag for binary output
    bool binary_;

public:

    /// constructor
    Outputter();
    
    /// constructor which opens a file
    Outputter(FILE* f, bool b) : FileWrapper(f, nullptr), binary_(b) {};

    /// constructor which opens a file where `a` specifies append and `b` binary mode.
    Outputter(const char* name, bool a, bool b=false);
    
    /// Open a file where `a` specifies append and `b` binary mode.
    int  open(const char* name, bool a, bool b=false);
    
    /// Sets to write in binary format
    void binary(bool b) { binary_ = b; }
    
    /// Return the current binary format
    bool binary() const { return binary_; }

    /// Puts given string, and '01' or '10', to specify the byte order 
    void writeEndianess();
        
    /// Inserts a new line symbol, but only in text output mode
    void writeSoftNewline();
    
    /// Inserts `N` space(s), but only in text output mode
    void writeSoftSpace(size_t N = 1);
    
    /// Write integer on 1 byte 
    void writeInt8(int, char before=' ');
    /// Write integer on 2 bytes
    void writeInt16(int, char before=' ');
    /// Write integer on 4 bytes
    void writeInt32(int, char before=' ');
    
    /// Write unsigned integer on 1 byte  
    void writeUInt8(unsigned, char before=' ');
    /// Write unsigned integer on 2 bytes
    void writeUInt16(unsigned, char before=' ');
    /// Write unsigned integer on 4 bytes
    void writeUInt32(unsigned, char before=' ');
    /// Write unsigned integer on 4 bytes
    void writeUInt64(unsigned long, char before=' ');

    /// Write value on 4 bytes
    void writeFloat(float);
    /// Write value on 4 bytes
    void writeFloat(double x) { writeFloat((float)x); }
    /// Write multiple values using 4 bytes per value, and possibly a character before
    void writeFloats(const float*, size_t, char before=0);
    /// Write multiple values using 4 bytes per value, and possibly a character before
    void writeFloats(const double*, size_t, char before=0);

    /// Write value on 8 bytes
    void writeDouble(double);
    /// Write multiple values using 8 bytes per value
    void writeDoubles(const double*, size_t, char before=0);

    int  writeChar(int c, int b)
    {
        if ( binary_ )
            return putc_unlocked(c|b, mFile);
        else
            return putc_unlocked(c, mFile);
    }

};

#endif
