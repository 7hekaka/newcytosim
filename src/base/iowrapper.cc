// Cytosim was created by Francois Nedelec. Copyright Cambridge University 2020

#include <cmath>
#include "iowrapper.h"
#include "exceptions.h"


template<typename T> static inline T byteswap16(T& x)
{
    static_assert(sizeof(T)==2, "2 != sizeof()");
    return static_cast<T>(__builtin_bswap16(x));
}

template<typename T> static inline T byteswap32(T& x)
{
    static_assert(sizeof(T)==4, "4 != sizeof()");
    return static_cast<T>(__builtin_bswap32(x));
}

template<typename T> static inline T byteswap64(T& x)
{
    static_assert(sizeof(T)==8, "8 != sizeof()");
    return static_cast<T>(__builtin_bswap64(x));
}

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
    if ( 1 != fscanf(mFile, " %d", &i) )
        throw InvalidIO("readInt failed");
    return i;
}


int Inputter::readInt16()
{
    if ( ! binary_ )
        return readInt();
    
    int16_t v;
    if ( 1 != fread(&v, 2, 1, mFile) )
        throw InvalidIO("readInt16() failed");
    if ( binary_ == 2 )
        v = byteswap16(v);
    return v;
}


int Inputter::readInt32()
{
    if ( ! binary_ )
        return readInt();
    
    int32_t v;
    if ( 1 != fread(&v, 4, 1, mFile) )
        throw InvalidIO("readInt32() failed");
    if ( binary_ == 2 )
        v = byteswap32(v);
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
        return readUInt();
    
    return (uint8_t)get_char();
}


uint16_t Inputter::readUInt16bin()
{
#if 1
    union { uint16_t u; uint8_t c[2]; } u16;
    if ( binary_ == 2 )
    {
        u16.c[1] = get_char();
        u16.c[0] = get_char();
    }
    else
    {
        u16.c[0] = get_char();
        u16.c[1] = get_char();
    }
    return u16.u;
#else
    uint16_t v;
    if ( 1 != fread(&v, 2, 1, mFile) )
        throw InvalidIO("readUInt16bin() failed");
    if ( binary_ == 2 )
        v = byteswap16(v);
    return v;
#endif
}


uint16_t Inputter::readUInt16()
{
    if ( ! binary_ )
        return readUInt();
#if 1
    union { uint16_t u; uint8_t c[2]; } u16;
    if ( binary_ == 2 )
    {
        u16.c[1] = get_char();
        u16.c[0] = get_char();
    }
    else
    {
        u16.c[0] = get_char();
        u16.c[1] = get_char();
    }
    return u16.u;
#else
    uint16_t v;
    if ( 1 != fread(&v, 2, 1, mFile) )
        throw InvalidIO("readUInt16() failed");
    if ( binary_ == 2 )
        v = byteswap16(v);
    return v;
#endif
}


uint32_t Inputter::readUInt32bin()
{
    uint32_t v;
    if ( 1 != fread(&v, 4, 1, mFile) )
        throw InvalidIO("readUInt32bin() failed");
    if ( binary_ == 2 )
        v = byteswap32(v);
    return v;
}


uint32_t Inputter::readUInt32()
{
    if ( ! binary_ )
        return readUInt();
    
    uint32_t v;
    if ( 1 != fread(&v, 4, 1, mFile) )
        throw InvalidIO("readUInt32() failed");
    if ( binary_ == 2 )
        v = byteswap32(v);
    return v;
}


uint64_t Inputter::readUInt64()
{
    if ( ! binary_ )
        return readUInt();
    
    uint64_t v;
    if ( 1 != fread(&v, 8, 1, mFile) )
        throw InvalidIO("readUInt64() failed");
    if ( binary_ == 2 )
        v = byteswap64(v);
    return v;
}


float Inputter::readFixed()
{
    if ( binary_ )
    {
        uint16_t i;
        if ( 1 != fread(&i, 2, 1, mFile) )
            throw InvalidIO("readFixed() failed");
        if ( binary_ == 2 )
            i = byteswap16(i);
        constexpr float F = 1.f / 65535.f;
        return float(i) * F;
    }
    else
    {
        float v;
        if ( 1 != fscanf(mFile, " %f", &v) )
            throw InvalidIO("readAngle failed");
        return v;
    }
}


float Inputter::readAngle()
{
    assert_true( binary_ );
    int16_t i;
    if ( 1 != fread(&i, 2, 1, mFile) )
        throw InvalidIO("readAngle() failed");
    if ( binary_ == 2 )
        i = byteswap16(i);
    return float(i) * 1.e-4f;
}


void Inputter::readEulerAngles(float& a, float& b)
{
    assert_true( binary_ );
    uint16_t iu[2];
    if ( 2 != fread(&iu, 2, 2, mFile) )
        throw InvalidIO("readEulerAngles() failed");
    if ( binary_ == 2 )
    {
        iu[0] = byteswap16(iu[0]);
        iu[1] = byteswap16(iu[1]);
    }
    a = float(*((int16_t*)iu)) * 1.e-4f;
    b = float(iu[1]) * 5.e-5f;
}


float Inputter::readFloatBinary()
{
    assert_true( binary_ );
    union { float f; uint8_t i[4]; } u;
#if 1
    u.i[0] = get_char();
    u.i[1] = get_char();
    u.i[2] = get_char();
    u.i[3] = get_char();
#else
    if ( 1 != fread(&u, sizeof(float), 1, mFile) )
        throw InvalidIO("readFloat failed");
#endif
    if ( binary_ == 2 )
        u.f = byteswap32(u.f);
    return u.f;
}


float Inputter::readFloat()
{
    union { float f; uint8_t i[4]; } u;
    if ( binary_ )
    {
#if 1
        u.i[0] = get_char();
        u.i[1] = get_char();
        u.i[2] = get_char();
        u.i[3] = get_char();
#else
        if ( 1 != fread(&u, 4, 1, mFile) )
            throw InvalidIO("readFloat() failed");
#endif
        if ( binary_ == 2 )
            u.f = byteswap32(u.f);
    }
    else
    {
        if ( 1 != fscanf(mFile, " %f", &u.f) )
            throw InvalidIO("readFloat failed");
    }
    return u.f;
}


double Inputter::readDouble()
{
    double v;
    if ( binary_ )
    {
        if ( 1 != fread(&v, 8, 1, mFile) )
            throw InvalidIO("readDouble() failed");
        if ( binary_ == 2 )
            v = byteswap64(v);
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
void Inputter::readFloats(float flt[], const unsigned dim)
{
    unsigned stop = std::min(vecsize_, dim);
    unsigned d = 0;
    if ( binary_ )
    {
        if ( stop != fread(flt, sizeof(float), stop, mFile) )
            throw InvalidIO("readFloats(float) failed");
        if ( binary_ == 2 ) {
            for ( int i = 0; i < stop; ++i )
                flt[i] = byteswap32(flt[i]);
        }
        d = stop;
    } else {
        while ( d < stop )
            flt[d++] = readFloat();
    }
    while ( d < dim )
        flt[d++] = 0.0f;
    for ( d = stop; d < vecsize_; ++d )
        readFloat();
}


/**
 This will read `vecsize_` floats, and set `dim` values in ptr[], filling in with zeros.
 */
void Inputter::readFloats(double ptr[], const unsigned dim)
{
    unsigned stop = std::min(vecsize_, dim);
    unsigned d = 0;
    if ( binary_ && dim < 4 )
    {
        float flt[4] = { 0 };
        size_t cnt = fread(flt, sizeof(float), stop, mFile);
        if ( cnt < stop )
        {
            for ( size_t d = cnt; d < stop; ++d )
                ptr[d] = 0;
            throw InvalidIO("readFloats(double) failed");
        }
        if ( binary_ == 2 ) {
            for ( ; d < stop; ++d )
                ptr[d] = byteswap32(flt[d]);
        } else {
            for ( ; d < stop; ++d )
                ptr[d] = flt[d];
        }
    }
    else
    {
        while ( d < stop )
            ptr[d++] = readFloat();
    }
    while ( d < dim )
        ptr[d++] = 0.0;
    for ( d = stop; d < vecsize_; ++d )
        readFloat();
}


/// unpack and zero-out data to match dimensionality, eg { X, Y, 0 } if vsz==2
static void unpack_floats(int cnt, float flt[], const int vsz, const int dim)
{
    assert_true(vsz < dim);
    //printf("\n/"); for ( int i = 0; i < cnt*vsz; ++i ) printf(" %6.2f", flt[i]);
    int u = cnt;
    while ( u-- > 0 )
    {
        int i = dim-1;
        do
            flt[u*dim+i] = 0.f;
        while ( --i > vsz );
        do
            flt[u*dim+i] = flt[u*vsz+i];
        while ( --i >= 0 );
    }
    //printf("\nL"); for ( int i = 0; i < cnt*dim; ++i ) printf(" %6.2f", flt[i]);
}


/// unpack and zero-out data to match dimensionality, eg { X, Y, 0 } if vsz==2
static void pack_floats(int cnt, float flt[], const int vsz, const int dim)
{
    assert_true(dim < vsz);
    //printf("\n/"); for ( int i = 0; i < cnt*vsz; ++i ) printf(" %6.2f", flt[i]);
    for ( int i = 0; i < cnt ; ++i )
    {
        for ( int d = 0; d < dim; ++d )
            flt[i*dim+d] = flt[i*vsz+d];
    }
    //printf("\nL"); for ( int i = 0; i < cnt*dim; ++i ) printf(" %6.2f", flt[i]);
}


/**
This will read `n * vecsize_` floats, and store `n * dim` values in ptr[].
*/
void Inputter::readFloats(const unsigned cnt, float ptr[], const unsigned dim)
{
    if ( ! binary_ )
    {
        for ( unsigned i = 0; i < cnt ; ++i )
        {
            for ( unsigned d = 0; d < vecsize_; ++d )
            {
                float f = 0.f;
                if ( 1 != fscanf(mFile, " %f", &f) )
                    throw InvalidIO("readFloat failed");
                ptr[dim*i+d] = ( d < dim ? (double)f : 0.f );
            }
        }
        return;
    }

    // read all values in one call to fread()
    unsigned n = cnt * vecsize_;
    if ( n != fread(ptr, sizeof(float), n, mFile) )
        throw InvalidIO("readFloats(F) failed");

    if ( binary_ == 2 )
    {
        for ( unsigned i = 0; i < n; ++i )
            ptr[i] = byteswap32(ptr[i]);
    }
    if ( vecsize_ < dim )
        unpack_floats(cnt, ptr, vecsize_, dim);
    else if ( dim < vecsize_ )
        pack_floats(cnt, ptr, vecsize_, dim);
}



/**
This will read `n * vecsize_` floats, and store `n * dim` values in ptr[].
*/
void Inputter::readFloats(const unsigned cnt, double ptr[], const unsigned dim)
{
    if ( vecsize_ < dim )
    {
        for ( unsigned i = 0; i < cnt*dim; ++i )
            ptr[i] = 0.0;
    }

    if ( ! binary_ )
    {
        for ( unsigned i = 0; i < cnt ; ++i )
        {
            for ( unsigned d = 0; d < vecsize_; ++d )
            {
                float f = 0.f;
                if ( 1 != fscanf(mFile, " %f", &f) )
                    throw InvalidIO("readFloat failed");
                ptr[dim*i+d] = ( d < dim ? (double)f : 0.0 );
            }
        }
        return;
    }

    size_t n = cnt * vecsize_;
    float * flt = new float[n];
    
    // read all values in one call to fread()
    if ( n != fread(flt, sizeof(float), n, mFile) )
        throw InvalidIO("readFloats(D) failed");

    if ( binary_ == 2 )
    {
        for ( unsigned i = 0; i < n; ++i )
            flt[i] = byteswap32(flt[i]);
    }
    unsigned top = std::min(dim, vecsize_);
    for ( unsigned i = 0; i < cnt; ++i )
    {
        for ( unsigned d = 0; d < top; ++d )
            ptr[i*dim+d] = (double)flt[i*vecsize_+d];
    }
    delete[] flt;
}



/**
 This will read `vecsize_` doubles, and set `cnt` values in ptr[], filling in with zeros.
 */
void Inputter::readDoubles(double ptr[], const unsigned cnt)
{
    unsigned stop = std::min(vecsize_, cnt);
    unsigned d = 0;
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


Outputter::Outputter(const char* name, const bool a, const int b)
{
    open(name, a, b);
    sanityCheck();
}


int Outputter::open(const char* name, const bool a, const int b)
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
    if ( 1 != fwrite(&x, 2, 1, mFile) )
        throw InvalidIO("writeEndianess failed");
}


void Outputter::write(const std::string& arg)
{
    fwrite(arg.c_str(), 1, arg.size(), mFile);
}


void Outputter::writeInt(const int n, char before)
{
    if ( 2 > fprintf(mFile, "%c%i", before, n) )
        throw InvalidIO("writeInt failed");
}


void Outputter::writeUInt(const unsigned n)
{
    if ( 1 > fprintf(mFile, "%u", n) )
        throw InvalidIO("writeUInt failed");
}


void Outputter::writeUInt(const unsigned n, char before)
{
    if ( 2 > fprintf(mFile, "%c%u", before, n) )
        throw InvalidIO("writeUInt failed");
}


void Outputter::writeUInt(const unsigned long n, char before)
{
    if ( 2 > fprintf(mFile, "%c%lu", before, n) )
        throw InvalidIO("writeUInt failed");
}


void Outputter::writeInt8(const int n)
{
    if ( !binary_ )
        return writeInt(n, ' ');
    
    int8_t v = (int8_t)n;
    
    if ( n != v )
        throw InvalidIO("value out of range for writeInt8()");
    if ( 1 != fwrite(&v, 1, 1, mFile) )
        throw InvalidIO("writeInt8() failed");
}


void Outputter::writeInt16(const int n)
{
    if ( !binary_ )
        return writeInt(n, ' ');

    int16_t v = (int16_t)n;
    
    if ( n != v )
        throw InvalidIO("value out of range for writeInt16()");

    if ( 1 != fwrite(&v, 2, 1, mFile) )
        throw InvalidIO("writeInt16() failed");
}


void Outputter::writeInt32(const int n)
{
    if ( !binary_ )
        return writeInt(n, ' ');

    int32_t v = (int32_t)n;
    
    if ( n != v )
        throw InvalidIO("value out of range for writeInt32()");
    
    if ( 1 != fwrite(&v, 4, 1, mFile) )
        throw InvalidIO("writeInt32() failed");
}


void Outputter::writeUInt8(const unsigned n)
{
    if ( !binary_ )
        return writeUInt(n, ' ');

    uint8_t v = (uint8_t)n;
    
    if ( n != v )
        throw InvalidIO("value out of range for writeUInt8()");
    
    if ( 1 != fwrite(&v, 1, 1, mFile) )
        throw InvalidIO("writeUInt8() failed");
}


void Outputter::writeUInt16Binary(const unsigned n)
{
    assert_true( binary_ );
    uint16_t v = (uint16_t)n;
    
    if ( n != v )
        throw InvalidIO("value out of range for writeUInt16()");

    if ( 1 != fwrite(&v, 2, 1, mFile) )
        throw InvalidIO("writeUInt16() failed");
}

void Outputter::writeUInt16(const unsigned n)
{
    if ( !binary_ )
        return writeUInt(n, ' ');

    uint16_t v = (uint16_t)n;
    
    if ( n != v )
        throw InvalidIO("value out of range for writeUInt16()");

    if ( 1 != fwrite(&v, 2, 1, mFile) )
        throw InvalidIO("writeUInt16() failed");
}


void Outputter::writeUInt32(const unsigned n)
{
    if ( !binary_ )
        return writeUInt(n, ' ');

    uint32_t v = (uint32_t)n;
    
    if ( n != v )
        throw InvalidIO("value out of range for writeUInt32()");
    
    if ( 1 != fwrite(&v, 4, 1, mFile) )
        throw InvalidIO("writeUInt32() failed");
}


void Outputter::writeUInt64(const unsigned long n)
{
    if ( !binary_ )
        return writeUInt(n, ' ');

    uint64_t v = (uint64_t)n;
    
    if ( n != v )
        throw InvalidIO("value out of range for writeUInt64()");
    
    if ( 1 != fwrite(&v, 8, 1, mFile) )
        throw InvalidIO("writeUInt64() failed");
}


void Outputter::writeUInt16(const unsigned n, char before)
{
    if ( !binary_ )
        return writeUInt(n, before);

    uint16_t v = (uint16_t)n;
    
    if ( n != v )
        throw InvalidIO("value out of range for writeUInt16()");

    if ( 1 != fwrite(&v, 2, 1, mFile) )
        throw InvalidIO("writeUInt16() failed");
}


void Outputter::writeUInt32(const unsigned n, char before)
{
    if ( !binary_ )
        return writeUInt(n, before);

    uint32_t v = (uint32_t)n;
    
    if ( n != v )
        throw InvalidIO("value out of range for writeUInt32()");
    
    if ( 1 != fwrite(&v, 4, 1, mFile) )
        throw InvalidIO("writeUInt32() failed");
}


void Outputter::writePositiveFixed(const float x)
{
    if ( binary_ )
    {
        int32_t i = static_cast<int32_t>(std::round(x * 65535.f));
        uint16_t u = (uint16_t)std::max(0, std::min(i, 65535));
        if ( u != i )
            fprintf(stderr, "writePositiveFixed(%f) out-of-range\n", x);
        if ( 1 != fwrite(&u, 2, 1, mFile) )
            throw InvalidIO("writePositiveFixed() failed");
    }
    else
    {
        if ( 6 > fprintf(mFile, " %.6f", x) )
            throw InvalidIO("writePositiveFixed failed");
    }
}


void Outputter::writeSignedFixed(const float x)
{
    if ( binary_ )
    {
        int32_t i = static_cast<int32_t>(std::round(x * 32767.f));
        int16_t u = (int16_t)std::max(-32768, std::min(i, 32767));
        if ( u != i )
            fprintf(stderr, "writeSignedFixed(%f) out-of-range\n", x);
        if ( 1 != fwrite(&u, 2, 1, mFile) )
            throw InvalidIO("writeSignedFixed() failed");
    }
    else
    {
        if ( 6 > fprintf(mFile, " %.6f", x) )
            throw InvalidIO("writeSignedFixed failed");
    }
}


/*
 Since the angle is within [-PI, PI], using int16_t (max 32767), by scaling by 10000,
 expanding to the range [-31416, 31416]. The delta is about ~ 10-4 radian
 */
void Outputter::writeAngle(const float x)
{
    assert_true( binary_ );
    int16_t i = int16_t(x * 10000.f);
    if ( 1 != fwrite(&i, 2, 1, mFile) )
        throw InvalidIO("writeAngle() failed");
}

/**
 Store `a` in [-PI, PI], using int16_t (max 32767), by scaling by 10000,
 expanding to the range [-31416, 31416]. The delta is about ~ 10-4 radian

 Store `b` in [0, PI] using uint16_t (max 65535), by scaling by 20000,
 expanding to the range [0, 62832]. The delta is about ~ 5.10-5 radian
 */
void Outputter::writeEulerAngles(const float a, const float b)
{
    assert_true( binary_ );
    const float sup = 3.1999f;
    bool valid = ( std::fabs(a) < sup ) & ( 0 <= b ) & ( b < sup );
    uint16_t iu[2];
    *(int16_t*)iu = int16_t(a * 10000.f);
    iu[1] = uint16_t(b * 20000.f);
    if ( !valid || 2 != fwrite((int16_t*)iu, 2, 2, mFile) )
        throw InvalidIO("writeEulerAngles() failed");
}


void Outputter::writeFloatBinary(const float x)
{
    assert_true( binary_ );
    if ( 1 != fwrite(&x, sizeof(float), 1, mFile) )
        throw InvalidIO("writeFloat() failed");
}


void Outputter::writeFloat(const float x)
{
    if ( binary_ )
    {
        if ( 1 != fwrite(&x, sizeof(float), 1, mFile) )
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

void Outputter::writeFloatsBinary(const float* a, const size_t n)
{
    assert_true( binary_ );
    if ( n != fwrite(a, sizeof(float), n, mFile) )
        throw InvalidIO("writeFloat() failed");
    for ( size_t d = 0; d < n; ++d )
        writeFloat(a[d]);
}

void Outputter::writeFloatsBinary(const double* a, const size_t n)
{
    assert_true( binary_ );
    float * f = new float[n];
    for ( size_t i = 0; i < n; ++i )
        f[i] = static_cast<float>(a[i]);
    size_t res = fwrite(f, sizeof(float), n, mFile);
    delete[] f;
    if ( res != n )
        throw InvalidIO("writeFloatsBinary(double*) failed");
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
        if ( 8 != fwrite(&x, 8, 1, mFile) )
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

