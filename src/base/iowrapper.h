// Cytosim was created by Francois Nedelec. Copyright Cambridge University 2020

#ifndef  IOWRAPPER_H
#define  IOWRAPPER_H

#include <cstdio>
#include <cstdint>
#include "filewrapper.h"



/// Input in text or binary mode with automatic byte-swapping for endianess compatibility
class Inputter : public FileWrapper
{
private:
        
    /// The format ID of the input: this allow backward compatibility with older formats
    int format_;
    
    /// The dimensionality of vectors stored in the file
    unsigned vecsize_;
    
    /** if the state is stored in a binary format, binary_ is 1 or 2.
        with 'binary_==2', byte order is swapped automatically
        to allow reading/writing accross big- and little-endian machines.
        */
    int binary_;
    
public:
    
    /// set defaults (not-binary)
    void reset();
    
    /// Constructor
    Inputter(unsigned d) : FileWrapper(nullptr) { reset(); vecsize_=d; }
    
    /// Constructor
    Inputter(unsigned d, FILE* f, const char* path=nullptr) : FileWrapper(f, path) { reset(); vecsize_=d; }
    
    /// constructor which opens a file
    Inputter(unsigned d, const char* name, bool bin) : FileWrapper(name, bin?"rb":"r") { reset(); vecsize_=d; }

    /// return dimensionnally of vectors
    unsigned vectorSize() const { return vecsize_; }
    
    /// Set dimentionnality of vectors
    void vectorSize(unsigned d) { vecsize_ = d; }
    
    /// return file format version identification number
    int formatID() const { return format_; }

    /// set file format version identification number
    void setFormatID(int f) { format_ = f; }

    /// Returns 1 for native binary format, 2 for non-native binary format, and 0 if not binary
    int binary() const { return binary_; }
    
    /// initialize the automatic swapping of bytes in the binary format
    void setEndianess(const char[2]);
    
    /// Read ASCII integer
    int readInt();
    /// Read integer on 2 bytes
    int readInt16();
    /// Read integer on 4 bytes
    int readInt32();

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
    float readFixed();
    /// Read angle on 2 bytes
    float readAngle();
    /// Read two angles on 4 bytes
    void readEulerAngles(float&, float&);
    /// Reads one float on 4 bytes
    float readFloatBinary();
    /// Reads one float on 4 bytes
    float readFloat();
    /// Reads one double on 8 bytes
    double readDouble();
    
    /// Reads one vector, setting `cnt` coordinates in the array
    void readFloats(float[], unsigned dim);
    /// Reads one vector, setting `cnt` coordinates in the array
    void readFloats(double[], unsigned dim);

    /// Reads one vector, setting `cnt` coordinates in the array
    void readFloats(unsigned cnt, float[], unsigned dim);
    /// Reads one vector, setting `cnt` coordinates in the array
    void readFloats(unsigned cnt, double[], unsigned dim);

    /// Reads one vector, setting `cnt` coordinates in the array
    void readDoubles(double[], unsigned D);

};


#pragma mark -


/// Output in text or binary mode in the native endianess
class Outputter : public FileWrapper
{
    
private:
        
    /// Flag for binary output
    int binary_;

public:

    /// constructor
    Outputter();
    
    /// constructor which opens a file, in binary mode if 'b==true'
    Outputter(FILE* f, int b) : FileWrapper(f, nullptr), binary_(b) {};

    /// constructor which opens a file where `a` specifies append and `b` binary mode.
    Outputter(const char* name, bool a, int b=0);
    
    /// Open a file where `a` specifies append and `b` binary mode.
    int open(const char* name, bool a, int b=0);
    
    /// Sets to write in binary format
    void binary(int b) { binary_ = b; }
    
    /// Return the current binary format
    int binary() const { return binary_; }

    /// Puts given string, and '01' or '10', to specify the byte order 
    void writeEndianess();
    
    /// Write integer in ASCII
    void writeInt(int, char before);
    /// Write unsigned integer in ASCII
    void writeUInt(unsigned);
    /// Write unsigned integer in ASCII
    void writeUInt(unsigned, char before);
    /// Write unsigned integer in ASCII
    void writeUInt(unsigned long, char before);

    /// Write integer on 1 byte
    void writeInt8(int);
    /// Write integer on 2 bytes
    void writeInt16(int);
    /// Write integer on 4 bytes
    void writeInt32(int);
    
    /// Write unsigned integer on 1 byte
    void writeUInt8(unsigned);
    /// Write unsigned integer on 2 bytes
    void writeUInt16(unsigned);
    /// Write unsigned integer on 4 bytes
    void writeUInt32(unsigned);
    /// Write unsigned integer on 4 bytes
    void writeUInt64(unsigned long);

    /// Write unsigned integer on 1 byte
    void writeUInt16(unsigned, char before);
    /// Write unsigned integer on 2 bytes
    void writeUInt32(unsigned, char before);
    /// Write unsigned integer on 2 bytes
    void writeUInt16Binary(unsigned);

    /// store float in [0, 1] using 2 bytes
    void writeSignedFixed(float);
    /// store float in [-1, 1] using 2 bytes
    void writePositiveFixed(float);

    /// store an angle in [-PI, PI] using 2 bytes
    void writeAngle(float);
    /// store `a` in [-PI, PI] and `b` in [0, PI], using 2 bytes for each
    void writeEulerAngles(float a, float b);

    /// Write value on 4 bytes
    void writeFloat(float);
    /// Write value on 4 bytes
    void writeFloat(double x) { writeFloat((float)x); }
    /// disable any implicit conversion
    template <typename T> void writeFloat(T x) = delete;
    
    /// Write value on 4 bytes
    void writeFloatBinary(float);
    /// Write multiple values using 4 bytes per value, and possibly a character before
    void writeFloatsBinary(const float*, size_t);
    /// Write multiple values using 4 bytes per value, and possibly a character before
    void writeFloatsBinary(const double*, size_t);

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
    void write(const std::string&);

    /// put character
    void writeChar(const int c) { putc_unlocked(c, mFile); }

};

#endif
