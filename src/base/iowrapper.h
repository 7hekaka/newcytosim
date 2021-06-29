// Cytosim was created by Francois Nedelec. Copyright Cambridge University 2020

#ifndef  IOWRAPPER_H
#define  IOWRAPPER_H

#include <cstdio>
#include <stdint.h>
#include "filewrapper.h"


/// the highest bit is not used by ASCII codes
constexpr uint8_t HIGH_BIT = 128;
constexpr uint8_t LOW_BITS = 127;


/// Input in text or binary mode with automatic byte-swapping for endianess compatibility
class Inputter : public FileWrapper
{
private:
        
    /// The format ID of the input: this allow backward compatibility with older formats
    size_t format_;
    
    /// The dimensionality of vectors stored in the file
    size_t vecsize_;
    
    /** if the state is stored in a binary format, binary_
        is set to 1 or 2. with 2, byte order is swapped automatically
        this occurs for example when reading a simulation calculated 
        on PC from mac, or vice et versa.
        */
    int binary_;
    
public:
    
    /// set defaults (not-binary)
    void     reset();
    
    /// Constructor
    Inputter(size_t d) : FileWrapper(nullptr) { reset(); vecsize_=d; }
    
    /// Constructor
    Inputter(size_t d, FILE* f, const char* path=nullptr) : FileWrapper(f, path) { reset(); vecsize_=d; }
    
    /// constructor which opens a file
    Inputter(size_t d, const char* name, bool bin) : FileWrapper(name, bin?"rb":"r") { reset(); vecsize_=d; }

    /// return dimensionnally of vectors
    size_t   vectorSize()     const { return vecsize_; }
    
    /// Set dimentionnality of vectors
    void     vectorSize(size_t d) { vecsize_ = d; }
    
    /// returns the type of input
    size_t   formatID()       const { return format_; }

    /// returns the type of input
    void     formatID(size_t f)   { format_ = f; }

    /// Returns 1 for native binary format, 2 for non-native binary format, and 0 if not binary
    int      binary()         const { return binary_; }
    
    /// initialize the automatic swapping of bytes in the binary format
    void     setEndianess(const char[2]);
    
    /// Read ASCII integer
    int      readInt();
    /// Read integer on 2 bytes
    int16_t  readInt16();
    /// Read integer on 4 bytes
    int32_t  readInt32();

    /// Read ASCII integer
    unsigned readUInt();
    /// Read unsigned integer on 1 byte
    uint8_t  readUInt8();
    /// Read unsigned integer on 2 bytes
    uint16_t readUInt16();
    /// Read unsigned integer on 4 bytes
    uint32_t readUInt32();
    /// Read unsigned integer on 8 bytes
    uint64_t readUInt64();
    
    /// Read unsigned integer on 2 bytes
    uint16_t readUInt16bin();
    /// Read unsigned integer on 4 bytes
    uint32_t readUInt32bin();

    /// Reads float in [0, 1] stored on 2 bytes
    float    readFixed();
    /// Read angle on 2 bytes
    float    readAngle();
    /// Read angle on 2 bytes
    float    readPositiveAngle();
    /// Reads one float on 4 bytes
    float    readFloat();
    /// Reads one double on 8 bytes
    double   readDouble();
    
    /// Reads one vector, setting `cnt` coordinates in the array
    void     readFloats(float[], size_t dim);
    /// Reads one vector, setting `cnt` coordinates in the array
    void     readFloats(double[], size_t dim);

    /// Reads one vector, setting `cnt` coordinates in the array
    void     readFloats(size_t cnt, float[], size_t dim);
    /// Reads one vector, setting `cnt` coordinates in the array
    void     readFloats(size_t cnt, double[], size_t dim);

    /// Reads one vector, setting `cnt` coordinates in the array
    void     readDoubles(double[], size_t D);

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
    
    /// Write integer in ASCII
    void writeInt(int);
    /// Write integer on 1 byte
    void writeInt8(int);
    /// Write integer on 2 bytes
    void writeInt16(int);
    /// Write integer on 4 bytes
    void writeInt32(int);
    
    /// Write unsigned integer in ASCII
    void writeUInt(unsigned);
    /// Write unsigned integer on 1 byte
    void writeUInt8(unsigned);
    /// Write unsigned integer on 2 bytes
    void writeUInt16(unsigned);
    /// Write unsigned integer on 4 bytes
    void writeUInt32(unsigned);
    /// Write unsigned integer on 4 bytes
    void writeUInt64(unsigned long);

    /// Write unsigned integer in ASCII
    void writeUInt(unsigned, char before);
    /// Write unsigned integer on 1 byte
    void writeUInt16(unsigned, char before);
    /// Write unsigned integer on 2 bytes
    void writeUInt32(unsigned, char before);

    /// check if x would overflow the fixed format
    static bool overflowFixed(float x) { int16_t i=int16_t(x*2048.f); return i != int16_t(i); }
    /// store float in [0, 1] using 2 bytes
    void writeFixed(float);
    
    /// store an angle in [-PI, PI] using 2 bytes
    void writeAngle(float);
    /// store an angle in [0, PI] using 2 bytes
    void writePositiveAngle(float);

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

       
    /// Add new line symbol, but only in text output mode
    void writeSoftNewline() { if ( !binary_ ) put_char('\n'); }
    
    /// Add a space, but only in text output mode
    void writeSoftSpace() { if ( !binary_ ) put_char(' '); }
    
    /// put a C++ string
    void writeLine(const std::string& arg) { put_line(arg); }

    /// put character
    void writeChar(const int c) { putc_unlocked(c, mFile); }

};

#endif
