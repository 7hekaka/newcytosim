// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

/**
 'Frametool' is a simple utility that can read and extract frames in
 Cytosim's trajectory files, usually called "objects.cmo".
 
 It only uses the START and END tags of frames, and does not care
 about the organization of the data contained between these tags.
 
 'Frametool' can extract frames from the file, which is useful
 for example to reduce the size of 'objects.cmo' by dropping some frames.

 You can reduce the file size by half by dropping every odd frame:
 > frametool objects.cmo 0:2: > o.cmo
 > mv o.cmo objects.cmo
 
 Another tool 'sieve' can be used to read/write object-files,
 allowing finer manipulation of the simulation frames.
 
 Last updated 11.12.2019
*/

#include <errno.h>
#include <cstdio>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>


enum { COUNT, COPY, LAST, SIZE, EPID, SPLIT };
enum { UNKNOWN, FRAME_START, FRAME_SECTION, FRAME_END };

FILE * output = stdout;
const size_t buf_size = 128;
char buf[buf_size];

unsigned long frame_pid = 0;
double frame_time = 0;


FILE * openfile(char name[], char const* mode)
{
    FILE * file = fopen(name, mode);
    
    if ( file==nullptr )
    {
        fprintf(stderr, "Could not open file `%s'\n", name);
        return nullptr;
    }
    
    if ( ferror(file) )
    {
        fclose(file);
        fprintf(stderr, "Error opening file `%s'\n", name);
        return nullptr;
    }
    
    return file;
}


/**
 read a line, and returns a code indicating if this is the start
 or the end of a cytosim frame
 */
int whatline(FILE* in, FILE* out)
{
    char *const end = buf + buf_size - 1;
    char * ptr = buf;

    int c = 0;
    do {
        c = getc_unlocked(in);
        
        if ( c == EOF )
            return EOF;

        if ( ptr < end )
            *ptr++ = (char)c;
        
        if ( out )
            putc_unlocked(c, out);
        
    } while ( c != '\n' );
    
    // fill-in with zeros:
    if ( ptr < end )
        *ptr = 0;
    else
        *end = 0;
    
    if ( *buf == '#' )
    {
        if ( 0 == strncmp(buf, "#frm ", 5) )     return FRAME_START;
        if ( 0 == strncmp(buf, "#frame ", 7) )   return FRAME_START;
        if ( 0 == strncmp(buf, "#Cytosim ", 9) )
        {
            frame_pid = strtoul(buf+10, nullptr, 10);
            return FRAME_START;
        }
        if ( 0 == strncmp(buf, "#end ", 5) )     return FRAME_END;
        if ( 0 == strncmp(buf, " #end ", 6) )    return FRAME_END;
        if ( 0 == strncmp(buf, "#section ", 9) ) return FRAME_SECTION;
        if ( 0 == strncmp(buf, "#time ", 6) )
            frame_time = strtod(buf+6, nullptr);
    }
    return UNKNOWN;
}


//=============================================================================


/// Slice represents a regular subset of indices
class Slice
{
    size_t s; ///< start
    size_t i; ///< increment
    size_t e; ///< end
    
public:
    
    Slice()
    {
        s =  0;
        i =  1;
        e = ~0UL;
    }
    
    Slice(const char arg[])
    {
        s =  0;
        i =  1;
        e = ~0UL;

        int c = 0;
        char * str;
        errno = 0;
        s = strtoul(arg, &str, 10);
        if ( errno ) goto finish;
        if ( *str == ':' )
        {
            c = 1;
            ++str;
            i = strtoul(str, &str, 10);
            if ( errno ) goto finish;
            if ( *str == ':' )
            {
                ++str;
                if ( *str == 0 ) return;
                c = 2;
                e = strtoul(str, &str, 10);
                if ( errno ) goto finish;
            }
        }
        if ( *str ) goto finish;

        if ( c == 0 )
            e = s;
        if ( c == 1 )
        { e = i; i = 1; }

        //fprintf(stderr, "frametool slice %lu:%lu:%lu\n", s, i, e);
        return;
    finish:
        fprintf(stderr, "syntax error in `%s', expected START:INCREMENT:END\n", arg);
        exit(EXIT_FAILURE);
    }
    
    bool match(size_t n)
    {
        if ( n < s )
            return false;
        if ( e < n )
            return false;
        return 0 == ( n - s ) % i;
    }
    
    size_t last()
    {
        return e;
    }
};

//=============================================================================

void countFrame(const char str[], FILE* in)
{
    size_t frm = 0;
    int code = 0;
    do {
        code = whatline(in, nullptr);
        frm += ( code == FRAME_END );
    } while ( code != EOF );
    
    printf("%40s: %lu frames\n", str, frm);
}


void sizeFrame(FILE* in)
{
    long pos = 0;
    long cnt = 0, old = 0;
    size_t frm = 0;

    while ( !ferror(in) )
    {
        switch(whatline(in, nullptr))
        {
            case FRAME_END: {
                size_t kb = ( ftell(in) - pos ) >> 10;
                printf("pid %lu   frame %6lu   time: %10.5f %6lu kB %7li lines (%+li)\n",
                       frame_pid, frm, frame_time, kb, cnt, cnt-old);
                old = cnt;
                ++frm;
            }
            case FRAME_START:
                pos = ftell(in);
                cnt = 0;
                break;
            case EOF:
                return;
            default:
                ++cnt;
        }
    }
}


void extract(FILE* in, FILE* file, Slice sli)
{
    size_t frm = 0;
    FILE * out = sli.match(0) ? file : nullptr;

    while ( !ferror(in) )
    {
        switch(whatline(in, out))
        {
            case FRAME_START:
                if ( frm > sli.last() )
                    return;
                out = sli.match(frm) ? file : nullptr;
                break;
            case FRAME_END:
                if ( ++frm > sli.last() )
                    return;
                out = sli.match(frm) ? file : nullptr;
                break;
            case EOF:
                return;
        }
    }
}


void extract_pid(FILE* in, unsigned long pid)
{
    FILE* out = nullptr;
    
    while ( !ferror(in) )
    {
        switch(whatline(in, nullptr))
        {
            case FRAME_START:
            case FRAME_END:
                if ( pid == frame_pid )
                    out = output;
                else
                    out = nullptr;
                break;
            case EOF:
                return;
        }
    }
}


void extractLast(FILE* in)
{
    fpos_t pos, start;
    fgetpos(in, &start);

    int code = 0;
    while ( code != EOF )
    {
        code = whatline(in, nullptr);
        if ( code == FRAME_END )
        {
            start = pos;
            fgetpos(in, &pos);
        }
    }
    
    clearerr(in);
    fsetpos(in, &start);
    
    int c = 0;
    while ( 1 )
    {
        c = getc_unlocked(in);
        if ( c == EOF )
            break;
        putchar(c);
    }
    
    putchar('\n');
}


void split(FILE * in)
{
    size_t frm = 0;
    char name[128] = { 0 };
    snprintf(name, sizeof(name), "objects%04lu.cmo", frm);
    FILE * out = fopen(name, "w");
    
    while ( !ferror(in) )
    {
        switch(whatline(in, out))
        {
            case FRAME_START:
                snprintf(name, sizeof(name), "objects%04lu.cmo", frm);
                out = openfile(name, "w");
                if ( out ) {
                    flockfile(out);
                    fprintf(out, "#Cytosim\n");
                }
                break;
            case FRAME_END:
                if ( out ) {
                    funlockfile(out);
                    fclose(out);
                    out = nullptr;
                }
                ++frm;
                break;
            case EOF:
                return;
        }
    }
}

//=============================================================================

void help()
{
    printf("Synopsis:\n");
    printf("    `frametool` can list the frames present in a trajectory file,\n");
    printf("     and extract selected ones\n");
    printf("Syntax:\n");
    printf("    frametool FILENAME \n");
    printf("    frametool FILENAME size\n");
    printf("    frametool FILENAME pid=PID\n");
    printf("    frametool FILENAME INDICES\n");
    printf("    frametool FILENAME split\n");
    printf(" where INDICES specifies an integer or a range of integers as:\n");
    printf("        INDEX\n");
    printf("        START:END\n");
    printf("        START:\n");
    printf("        START:INCREMENT:END\n");
    printf("        START:INCREMENT:\n");
    printf("        last\n");
    printf(" The option 'split' will create one file for each frame in the input\n");
    printf("Examples:\n");
    printf("    frametool objects.cmo 0:2:\n");
    printf("    frametool objects.cmo 0:10\n");
    printf("    frametool objects.cmo last\n");
    printf("    frametool objects.cmo split\n");
}


bool is_file(const char path[])
{
    struct stat s;
    if ( 0 == stat(path, &s) )
        return S_ISREG(s.st_mode);
    return false;
}

bool is_dir(const char path[])
{
    struct stat s;
    if ( 0 == stat(path, &s) )
        return S_ISDIR(s.st_mode);
    return false;
}


int main(int argc, char* argv[])
{
    int has_file = 0;
    int mode = COUNT;
    char cmd[256] = "";
    char filename[256] = "objects.cmo";
    char outputname[256] = { 0 };
    unsigned long pid = 0;

    for ( int i = 1; i < argc ; ++i )
    {
        char * arg = argv[i];
        if ( 0 == strncmp(arg, "help", 4) )
        {
            help();
            return EXIT_SUCCESS;
        }
        char *dot = strrchr(arg, '.');
        if ( is_file(arg) )
        {
            if ( has_file++ )
            {
                fprintf(stderr, "error: only one input file can be specified\n");
                return EXIT_SUCCESS;
            }
            strncpy(filename, arg, sizeof(filename));
        }
        else if ( dot && 0==strcmp(dot, ".cmo" ) )
        {
            if ( outputname[0] )
            {
                fprintf(stderr, "error: only one output file can be specified\n");
                return EXIT_SUCCESS;
            }
            strncpy(outputname, arg, sizeof(outputname));
        }
        else if ( is_dir(arg) )
            snprintf(filename, sizeof(filename), "%s/objects.cmo", arg);
        else
        {
            strncpy(cmd, argv[i], sizeof(cmd));
            
            if ( isdigit(*cmd) )
                mode = COPY;
            else if ( 0 == strncmp(cmd, "last", 4) )
                mode = LAST;
            else if ( 0 == strncmp(cmd, "split", 5) )
                mode = SPLIT;
            else if ( 0 == strncmp(cmd, "size", 4) || *cmd == '+' )
                mode = SIZE;
            else if ( 0 == strncmp(cmd, "count", 5) )
                mode = COUNT;
            else if ( 0 == strncmp(cmd, "pid=", 4) )
            {
                mode = EPID;
                errno = 0;
                pid = strtoul(cmd+4, nullptr, 10);
                if ( errno )
                {
                    fprintf(stderr, "syntax error");
                    return EXIT_FAILURE;
                }
            }
            else
            {
                fprintf(stderr, "unexpected command (for help, invoke `frametool help`)\n");
                return EXIT_FAILURE;
            }
        }
    }
    
    //----------------------------------------------
    
    FILE * file = openfile(filename, "r");
    if ( !file || ferror(file) )
    {
        fprintf(stderr, "failed to open input file\n");
        return EXIT_FAILURE;
    }
    flockfile(file);

    if ( mode == COUNT )
        countFrame(filename, file);
    else if ( mode == SIZE )
        sizeFrame(file);
    else
    {
        if ( *outputname )
        {
            output = fopen(outputname, "wb");
            if ( ! output || ferror(output) )
            {
                fprintf(stderr, "failed to open output file\n");
                return EXIT_FAILURE;
            }
        }
        if ( mode == COPY )
            extract(file, output, Slice(cmd));
        else if ( mode == LAST )
            extractLast(file);
        else if ( mode == EPID )
            extract_pid(file, pid);
        else if ( mode == SPLIT )
            split(file);
        if ( output != stdout )
        {
            fprintf(stderr, "> %s\n", outputname);
            fclose(output);
        }
    }
    
    funlockfile(file);
    fclose(file);
    return EXIT_SUCCESS;
}
