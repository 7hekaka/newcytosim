// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

/**
 `Frametool` is a simple utility that can read and extract frames in
 Cytosim's trajectory files, usually called "objects.cmo".
 
 Frametool only recognizes the START and END tags of frames, and will copy
 all the data contained between these tags verbatim. It does not load the
 objects stored in the data, and is not able to modify them.
 
 Frametool can be used to extract one frame from the file, or multiple frames,
 specified using numerical indices: START:PERIOD:END
 These indices simply reflect the order in which the frame appear in the file,
 starting with index 0 for the first frame.

 A typical use case is to reduce the file by dropping some intermediate frames.
 Example, this will drop odd frames:
    > frametool objects.cmo 0:2: o.cmo
    > mv o.cmo objects.cmo
 
 Another tool `sieve` can be used to read/write object-files,
 allowing for a finer manipulation of the simulation frames.
 
 FJN, last updated 12.12.2023, 18.12.2024
*/

#include <unistd.h>
#include <errno.h>
#include <cstdio>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "slice.h"

enum { COUNT, COPY, LAST, SIZE, EPID, SPLIT };
enum { UNKNOWN, FRAME_START, FRAME_SECTION, FRAME_END, TIME_LINE };

const size_t buf_size = 128;
char buf[buf_size];

FILE * output = stdout;
unsigned long frame_pid = 0;
double frame_time = 0;
double time_added = 0;

/** Open a C-style File */
FILE * openFile(char name[], char const* mode)
{
    FILE * file = fopen(name, mode);
    
    if ( file==nullptr )
    {
        fprintf(stderr, "Could not open `%s'\n", name);
        return nullptr;
    }
    
    if ( ferror(file) )
    {
        fclose(file);
        fprintf(stderr, "Error opening `%s'\n", name);
        return nullptr;
    }
    
    return file;
}


/**
 read a line, and send every characters to 'out', and returns a code
 indicating if this line was the start or the end of a cytosim frame
 */
int whatLine(FILE* in, FILE* out)
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
    
    // terminate buffer with zeros:
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
        if ( 0 == strncmp(buf, "#time ", 6) ) {
            frame_time = strtod(buf+6, nullptr);
            // add a second line to change the time:
            if ( out && time_added > 0 )
                fprintf(out, "#time %.6f sec\n", frame_time+time_added);
            return TIME_LINE;
        }
    }
    return UNKNOWN;
}

//=============================================================================

void countFrame(const char str[], FILE* in)
{
    double time = -1;
    size_t frm = 0;
    int code = 0;
    do {
        code = whatLine(in, nullptr);
        if ( code == FRAME_END )
        {
            time = frame_time;
            ++frm;
        }
    } while ( code != EOF );
    
    printf("%40s: %lu frames, time: %10.3f", str, frm, time);
    if ( frame_time != time )
        printf(" + incomplete frame at %10.3f\n", frame_time);
    else
        printf("\n");
}


void sizeFrame(FILE* in, int details)
{
    long pos = 0, sec = 0;
    char str[buf_size] = { 0 };
    size_t frm = 0;

    while ( !ferror(in) )
    {
        switch(whatLine(in, nullptr))
        {
            case FRAME_SECTION:
                if ( details ) {
                    size_t bytes = ftell(in) - sec;
                    if ( str[0] )
                    {
                        if ( bytes > 1024 )
                            printf(" %32s : %6lu kB\n", str, bytes>>10);
                        else
                            printf(" %32s : %7lu B\n", str, bytes);
                    }
                    strncpy(str, buf, sizeof(str));
                    if ( isspace(str[strlen(buf)-1]) )
                        str[strlen(buf)-1] = 0;
                }
                sec = ftell(in);
                break;
            case FRAME_END: {
                size_t kb = ( ftell(in) - pos ) >> 10;
                printf("pid %lu   frame %6lu   time: %10.5f %6lu kB\n",
                       frame_pid, frm, frame_time, kb);
                sec = 0;
                ++frm;
            } // intentional fallthrough
            case FRAME_START:
                pos = ftell(in);
                *str = 0;
                break;
            case EOF:
                if ( sec )
                {
                    size_t kb = ( ftell(in) - pos ) >> 10;
                    printf("+INCOMPLETE frame %6lu   time: %10.5f %6lu kB\n",
                           frm, frame_time, kb);
                }
                return;
        }
    }
}


void extract(FILE* in, FILE* file, Slice const& sli)
{
    size_t frm = 0;
    FILE * out = sli.match(0) ? file : nullptr;

    while ( !ferror(in) )
    {
        switch(whatLine(in, out))
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


void extractPID(FILE* in, unsigned long pid)
{
    FILE* out = nullptr;
    
    while ( !ferror(in) )
    {
        switch(whatLine(in, nullptr))
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
        code = whatLine(in, nullptr);
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


void separateFrames(FILE * in)
{
    size_t frm = 0;
    char name[128] = { 0 };
    snprintf(name, sizeof(name), "objects%04lu.cmo", frm);
    FILE * out = fopen(name, "w");
    
    while ( !ferror(in) )
    {
        switch(whatLine(in, out))
        {
            case FRAME_START:
                snprintf(name, sizeof(name), "objects%04lu.cmo", frm);
                out = openFile(name, "w");
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
    printf("    frametool FILENAME size+\n");
    printf("    frametool FILENAME pid=PID\n");
    printf("    frametool FILENAME time+=FLOAT\n");
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
    int details = 0;
    char slice[256] = "";
    char filename[256] = "objects.cmo";
    char outputname[256] = { 0 };
    unsigned long pid = 0;

    for ( int i = 1; i < argc ; ++i )
    {
        char * arg = argv[i];
        char * dot = strrchr(argv[i], '.');
        if ( 0 == strncmp(arg, "help", 4) )
        {
            help();
            return EXIT_SUCCESS;
        }
        if ( is_file(arg) )
        {
            if ( 0 == has_file++ )
                strncpy(filename, arg, sizeof(filename));
            else
            {
                if ( outputname[0] )
                {
                    fprintf(stderr, "error: only one input file can be specified\n");
                    return EXIT_SUCCESS;
                }
                strncpy(outputname, arg, sizeof(outputname));
            }
        }
        else if ( is_dir(arg) )
        {
            snprintf(filename, sizeof(filename), "%s/objects.cmo", arg);
            if ( ! is_file(filename) )
            {
                fprintf(stderr, "error: missing file %s/objects.cmo\n", arg);
                return EXIT_SUCCESS;
            }
        }
        else if ( dot && 0==memcmp(dot, ".cmo", 5) )
        {
            if ( outputname[0] )
            {
                fprintf(stderr, "error: only one output file can be specified\n");
                return EXIT_SUCCESS;
            }
            strncpy(outputname, arg, sizeof(outputname));
        }
        else
        {
            if ( isdigit(*arg) || *arg == ':' )
            {
                mode = COPY;
                strncpy(slice, argv[i], sizeof(slice));
            }
            else if ( 0 == strncmp(arg, "last", 4) )
                mode = LAST;
            else if ( 0 == strncmp(arg, "split", 5) )
                mode = SPLIT;
            else if ( 0 == strncmp(arg, "size", 4) )
            {
                mode = SIZE;
                if ( arg[4] == '+' ) details = 1;
            }
            else if ( 0 == strncmp(arg, "+", 1) )
            {
                mode = SIZE;
                details = 1;
            }
            else if ( 0 == strncmp(arg, "count", 5) )
                mode = COUNT;
            else if ( 0 == strncmp(arg, "pid=", 4) )
            {
                mode = EPID;
                errno = 0;
                pid = strtoul(arg+4, nullptr, 10);
                if ( errno )
                {
                    fprintf(stderr, "syntax error in `%s`\n", arg);
                    return EXIT_FAILURE;
                }
            }
            else if ( 0 == strncmp(arg, "time+=", 6) )
            {
                errno = 0;
                time_added = strtod(arg+6, nullptr);
                if ( errno )
                {
                    fprintf(stderr, "syntax error in `%s`\n", arg);
                    return EXIT_FAILURE;
                }
            }
            else
            {
                fprintf(stderr, "unexpected command `%s` (for help, invoke `frametool help`)\n", arg);
                return EXIT_FAILURE;
            }
        }
    }
    
    //----------------------------------------------
    
    FILE * file = openFile(filename, "r");
    if ( !file )
        return EXIT_FAILURE;
    flockfile(file);

    if ( mode == COUNT )
        countFrame(filename, file);
    else if ( mode == SIZE )
        sizeFrame(file, details);
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
        
        if ( output == stdout && isatty(1) )
            fprintf(stderr, "Error: cannot send output to terminal!\n");
        else
        {
            if ( mode == COPY )
                extract(file, output, Slice(slice));
            else if ( mode == LAST )
                extractLast(file);
            else if ( mode == EPID )
                extractPID(file, pid);
            else if ( mode == SPLIT )
                separateFrames(file);
        }
        
        if ( output != stdout )
        {
            fprintf(stderr, "> %s\n", outputname);
            fclose(output);
        }
    }
    
    funlockfile(file);
    fclose(file);
}
