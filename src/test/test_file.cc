// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University

#include <stdint.h>
#include <cstdio>

/// used to address the highest bit which is not used by ASCII codes
const uint8_t HIGH_BIT = 128;
/// bit mask for those bits which are used by ASCII codes
const uint8_t LOW_BITS = 127;


FILE * openFile(const char * filename, const char* code)
{
    if ( filename[0] == '\0' )
        return 0;
    FILE * f = fopen(filename, code);
    if ( f )
    {
        if ( ferror(f) )
            fclose(f);
        else
            return f;
    }
    return 0;
}


int main()
{
    uint8_t g = 'x';
    int id = 0x70101;
    char filename[] = "test_file";
    FILE * F = openFile(filename, "wb");
    if ( F ) {
        //uint32_t u = id | (uint32_t(g|HIGH_BIT)<<24);
        uint8_t u[4] = { g, uint8_t((id>>16)&0xFF), uint8_t((id>>8)&0xFF), uint8_t(id&0xFF) };
        fwrite(&u, 4, 1, F);
        fclose(F);
        printf(" %i ", id);
    }
    
    FILE * G = openFile(filename, "rb");
    {
        uint8_t u[4] = { 0, 0, 0, 0 };
        fread(u, 4, 1, G);
        for ( int i = 0; i < 4; ++i )
            printf("|%x|", u[i]);
        id = ( int(u[1]) << 16 ) + ( int(u[2]) << 8 ) + int(u[3]);
        printf(" %i", id);
    }
    printf("|done\n");
}
