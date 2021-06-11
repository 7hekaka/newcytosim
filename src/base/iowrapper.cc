// Cytosim was created by Francois Nedelec. Copyright Cambridge University 2020

#include "iowrapper.h"
#include "exceptions.h"
#include "byteswap.h"


/// check the size of some types that are baked in the code
static void sanityCheck()
{
    bool okay = true;
    okay &= ( 2 == sizeof(uint16_t) );
    okay &= ( 4 == sizeof(uint32_t) );
    okay &= ( 8 == sizeof(uint64_t) );
    okay &= ( 4 == sizeof(float) );
    okay &= ( 8 == sizeof(double) );
    if ( ! okay )
    {
        fprintf(stderr, "Error: non-standard types in Inputter\n");
        exit(EXIT_FAILURE);
    }
}


//==============================================================================
#pragma mark - INPUT

void Inputter::reset()
{
    format_  = 0;
    vecsize_ = 3;
    binary_  = 0;
    sanityCheck();
}


/**
 Reads a short and compares with the native storage, to set
 binary_=1, for same-endian or binary_ = 2, for opposite endian
*/
void Inputter::setEndianess(const char data[2])
{
    char native[3] = { 0 };
    *((uint16_t*)native) = 12592U;
    //binary_ = 1 for same-endianess, 2 for opposite-endianess:
    binary_ = 1 + ( data[0] != native[0] );
}


int Inputter::readInt()
{
    int i;
    if ( 1 != fscanf(mFile, " %i", &i) )
        throw InvalidIO("readInt failed");
    return i;
}


int16_t Inputter::readInt16()
{
    if ( ! binary_ )
        return readInt();
    
    int16_t v;
    if ( 1 != fread(&v, 2, 1, mFile) )
        throw InvalidIO("readInt16() failed");
    if ( binary_ == 2 )
        v = byteswap(v);
    return v;
}


int32_t Inputter::readInt32()
{
    if ( ! binary_ )
        return readInt();
    
    int32_t v;
    if ( 1 != fread(&v, 4, 1, mFile) )
        throw InvalidIO("readInt32() failed");
    if ( binary_ == 2 )
        v = byteswap(v);
    return v;
}


unsigned Inputter::readUInt()
{
    unsigned u;
    if ( 1 != fscanf(mFile, " %u", &u) )
        throw InvalidIO("readUInt failed");
    return u;
}


uint8_t Inputter::readUInt8()
{
    if ( ! binary_ )
        return readInt();
    
    return get_byte();
}


uint16_t Inputter::readUInt16bin()
{
    uint16_t v;
    if ( 1 != fread(&v, 2, 1, mFile) )
        throw InvalidIO("readUInt16() failed");
    if ( binary_ == 2 )
        v = byteswap(v);
    return v;
}


uint16_t Inputter::readUInt16()
{
    if ( ! binary_ )
        return readInt();
    
    uint16_t v;
    if ( 1 != fread(&v, 2, 1, mFile) )
        throw InvalidIO("readUInt16() failed");
    if ( binary_ == 2 )
        v = byteswap(v);
    return v;
}


uint32_t Inputter::readUInt32bin()
{
    uint32_t v;
    if ( 1 != fread(&v, 4, 1, mFile) )
        throw InvalidIO("readUInt32() failed");
    if ( binary_ == 2 )
        v = byteswap(v);
    return v;
}


uint32_t Inputter::readUInt32()
{
    if ( ! binary_ )
        return readInt();
    
    uint32_t v;
    if ( 1 != fread(&v, 4, 1, mFile) )
        throw InvalidIO("readUInt32() failed");
    if ( binary_ == 2 )
        v = byteswap(v);
    return v;
}


uint64_t Inputter::readUInt64()
{
    if ( ! binary_ )
        return readInt();
    
    uint64_t v;
    if ( 1 != fread(&v, 8, 1, mFile) )
        throw InvalidIO("readUInt64() failed");
    if ( binary_ == 2 )
        v = byteswap(v);
    return v;
}


float Inputter::readFixed()
{
    uint16_t i;
    if ( 1 != fread(&i, 2, 1, mFile) )
        throw InvalidIO("readFixed() failed");
    if ( binary_ == 2 )
        i = byteswap(i);
    return float(i) * 0x1p-11;
}


float Inputter::readFloat()
{
    float v;
    if ( binary_ )
    {
        if ( 1 != fread(&v, 4, 1, mFile) )
            throw InvalidIO("readFloat() failed");
        if ( binary_ == 2 )
            v = byteswap(v);
    }
    else
    {
        if ( 1 != fscanf(mFile, " %f", &v) )
            throw InvalidIO("readFloat failed");
    }
    return v;
}


double Inputter::readDouble()
{
    double v;
    if ( binary_ )
    {
        if ( 1 != fread(&v, 8, 1, mFile) )
            throw InvalidIO("readDouble() failed");
        if ( binary_ == 2 )
            v = byteswap(v);
    }
    else
    {
        if ( 1 != fscanf(mFile, " %lf", &v) )
            throw InvalidIO("readDouble failed");
    }
    return v;
}


/**
 This will read `vecsize_` floats, and set `dim` values in ptr[], filling in with zeros.
 The default vector size can be changed by calling `vectorSize(INT)`
 */
void Inputter::readFloats(float ptr[], const size_t dim)
{
    size_t stop = std::min(vecsize_, dim);
    size_t d = 0;
    while ( d < stop )
        ptr[d++] = readFloat();
    while ( d < dim )
        ptr[d++] = 0.0f;
    for ( d = stop; d < vecsize_; ++d )
        readFloat();
}


/**
 This will read `vecsize_` floats, and set `dim` values in ptr[], filling in with zeros.
 */
void Inputter::readFloats(double ptr[], const size_t dim)
{
    size_t stop = std::min(vecsize_, dim);
    size_t d = 0;
    while ( d < stop )
        ptr[d++] = readFloat();
    while ( d < dim )
        ptr[d++] = 0.0;
    for ( d = stop; d < vecsize_; ++d )
        readFloat();
}


/**
This will read `n * vecsize_` floats, and store `n * dim` values in ptr[].
*/
void Inputter::readFloats(const size_t cnt, float ptr[], const size_t dim)
{
    if ( dim < vecsize_ || ! binary_ )
    {
        for ( size_t i = 0; i < cnt ; ++i )
            readFloats(ptr+dim*i, dim);
        return;
    }

    size_t n = cnt * vecsize_;
    if ( n != fread(ptr, 4, n, mFile) )
        throw InvalidIO("readFloats(D) failed");

    if ( binary_ == 2 )
    {
        for ( size_t i = 0; i < n; ++i )
            ptr[i] = byteswap(ptr[i]);
    }
    if ( vecsize_ < dim )
    {
        size_t u = cnt;
        while ( u-- > 0 )
        {
            size_t i = dim;
            while ( i-- > vecsize_ )
                ptr[u*dim+i] = 0.f;
            while ( i-- > 0 )
                ptr[u*dim+i] = ptr[u*vecsize_+i];
        }
    }
}



/**
 This will read `vecsize_` doubles, and set `cnt` values in ptr[], filling in with zeros.
 */
void Inputter::readDoubles(double ptr[], const size_t cnt)
{
    size_t stop = std::min(vecsize_, cnt);
    size_t d = 0;
    while ( d < stop )
        ptr[d++] = readDouble();
    while ( d < cnt )
        ptr[d++] = 0.0;
    for ( d = stop; d < vecsize_; ++d )
        readDouble();
}

//==============================================================================
#pragma mark - OUTPUT

Outputter::Outputter()
: FileWrapper(stdout) 
{
    binary_ = false;
    sanityCheck();
}


Outputter::Outputter(const char* name, const bool a, const bool b)
{
    open(name, a, b);
    sanityCheck();
}


int Outputter::open(const char* name, const bool a, const bool b)
{
    binary_ = b;
    
    //create a 'mode' string appropriate for Windows OS
    char m[3] = { 0 };
    
    if ( a )
        m[0] = 'a';
    else
        m[0] = 'w';
    
    if ( b )
        m[1] = 'b';
        
    return FileWrapper::open(name, m);
}


void Outputter::writeEndianess()
{
    //the value corresponds to the ASCII code of "01"
    uint16_t x = 12592U;
    if ( 2 != fwrite(&x, 1, 2, mFile) )
        throw InvalidIO("writeEndianess failed");
}


void Outputter::writeInt(const int n, char before)
{
    if ( 2 > fprintf(mFile, "%c%i", before, n) )
        throw InvalidIO("writeInt failed");
}


void Outputter::writeInt8(const int n, char before)
{
    if ( !binary_ )
        return writeInt(n, before);
    
    int8_t v = (int8_t)n;
    
    if ( n != v )
        throw InvalidIO("value out of range for writeInt8()");
    if ( 1 != fwrite(&v, 1, 1, mFile) )
        throw InvalidIO("writeInt8() failed");
}


void Outputter::writeInt16(const int n, char before)
{
    if ( !binary_ )
        return writeInt(n, before);

    int16_t v = (int16_t)n;
    
    if ( n != v )
        throw InvalidIO("value out of range for writeInt16()");

    if ( 2 != fwrite(&v, 1, 2, mFile) )
        throw InvalidIO("writeInt16() failed");
}


void Outputter::writeInt32(const int n, char before)
{
    if ( !binary_ )
        return writeInt(n, before);

    int32_t v = (int32_t)n;
    
    if ( n != v )
        throw InvalidIO("value out of range for writeInt32()");
    
    if ( 4 != fwrite(&v, 1, 4, mFile) )
        throw InvalidIO("writeInt32() failed");
}


void Outputter::writeUInt(const unsigned n, char before)
{
    if ( before )
    {
        if ( 2 > fprintf(mFile, "%c%u", before, n) )
            throw InvalidIO("writeUInt8 failed");
    }
    else {
        if ( 1 > fprintf(mFile, "%u", n) )
            throw InvalidIO("writeUInt8 failed");
    }
}

void Outputter::writeUInt8(const unsigned n, char before)
{
    if ( !binary_ )
        return writeUInt(n, before);

    uint8_t v = (uint8_t)n;
    
    if ( n != v )
        throw InvalidIO("value out of range for writeUInt8()");
    
    if ( 1 != fwrite(&v, 1, 1, mFile) )
        throw InvalidIO("writeUInt8() failed");
}


void Outputter::writeUInt16(const unsigned n, char before)
{
    if ( !binary_ )
        return writeUInt(n, before);

    uint16_t v = (uint16_t)n;
    
    if ( n != v )
        throw InvalidIO("value out of range for writeUInt16()");

    if ( 2 != fwrite(&v, 1, 2, mFile) )
        throw InvalidIO("writeUInt16() failed");
}


void Outputter::writeUInt32(const unsigned n, char before)
{
    if ( !binary_ )
        return writeUInt(n, before);

    uint32_t v = (uint32_t)n;
    
    if ( n != v )
        throw InvalidIO("value out of range for writeUInt32()");
    
    if ( 4 != fwrite(&v, 1, 4, mFile) )
        throw InvalidIO("writeUInt32() failed");
}


void Outputter::writeUInt64(const unsigned long n, char before)
{
    if ( !binary_ )
        return writeUInt(n, before);

    uint64_t v = (uint64_t)n;
    
    if ( n != v )
        throw InvalidIO("value out of range for writeUInt64()");
    
    if ( 8 != fwrite(&v, 1, 8, mFile) )
        throw InvalidIO("writeUInt64() failed");
}


void Outputter::writeFixed(const float x)
{
    uint16_t i = uint16_t(x * 2048.f);
    if ( 2 != fwrite(&i, 1, 2, mFile) )
        throw InvalidIO("writeFixed() failed");
}


void Outputter::writeFloat(const float x)
{
    if ( binary_ )
    {
        if ( 4 != fwrite(&x, 1, 4, mFile) )
            throw InvalidIO("writeFloat() failed");
    }
    else
    {
        if ( 6 > fprintf(mFile, " %.6f", x) )
            throw InvalidIO("writeFloat failed");
    }
}


void Outputter::writeFloats(const float* a, const size_t n, char before)
{
    if ( before && !binary_ )
        putc(before, mFile);
    
    for ( size_t d = 0; d < n; ++d )
        writeFloat(a[d]);
}


void Outputter::writeFloats(const double* a, const size_t n, char before)
{
    if ( before && !binary_ )
        putc(before, mFile);
    
    for ( size_t d = 0; d < n; ++d )
        writeFloat(a[d]);
}


void Outputter::writeDouble(const double x)
{
    if ( binary_ )
    {
        if ( 8 != fwrite(&x, 1, 8, mFile) )
            throw InvalidIO("writeDouble() failed");
    }
    else
    {
        if ( 10 > fprintf(mFile, " %.8lf", x) )
            throw InvalidIO("writeDouble failed");
    }
}


void Outputter::writeDoubles(const double* a, const size_t n, char before)
{
    if ( before && !binary_ )
        putc(before, mFile);
    
    for ( size_t d = 0; d < n; ++d )
        writeDouble(a[d]);
}


void Outputter::writeSoftNewline()
{
    if ( !binary_ )
        putc('\n', mFile);
    fflush(mFile);
}


void Outputter::writeSoftSpace(size_t N)
{
    if ( !binary_ )
    {
        while ( N > 0 ) {
            fprintf(mFile, " ");
            N--;
        }
    }
}

