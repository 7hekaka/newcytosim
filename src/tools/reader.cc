// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University
/*
 This is mostly a test for the class defined in frame_reader.h
 but it can be used to navigate frames inside a trajectory file 'objects.cmo'
*/

#include <cstring>
#include <cctype>
#include <cstdlib>

#include "filepath.h"
#include "glossary.h"
#include "messages.h"
#include "iowrapper.h"
#include "stream_func.h"
#include "frame_reader.h"
#include "simul.h"
#include "parser.h"

//------------------------------------------------------------------------------

void help(std::ostream& os)
{
    os << "The `reader` reads Cytosim's trajectory file\n";
    os << "Syntax: reader [OPTIONS] [DIRECTORY] INPUT_FILE_NAME output=FILE_NAME\n";
    os << "\n";
    os << "OPTIONS:\n";
    os << "     help       display this message\n";
    os << "     binary=0   write text coordinates in `file_out'\n";
    os << "     binary=1   write binary coordinates in `file_out'\n";
    os << "     verbose=?  set the verbose level\n";
    os << "Made with format version " << Simul::currentFormatID << " and DIM=" << DIM << "\n";
}

void instructions(std::ostream& os = std::cout)
{
    os << "Commands understood at prompt:\n";
    os << "  'q'      quit\n";
    os << "  'n'      read next frame\n";
    os << "  'w'      write frame\n";
    os << "  'c'      clear buffer without changing positions\n";
    os << "  'r'      rewind\n";
    os << "  'e'      erase state\n";
    os << " INTEGER   read specified frame if possible\n";
}

void inventory(Simul& simul, size_t frame)
{
    printf("loaded frame %li @ %.3f:", frame, simul.time());
    std::stringstream ss;
    simul.reportInventory(ss);
    StreamFunc::prefix_lines(std::cout, ss, "    ", 0, 0);
    printf("\n");
}

//------------------------------------------------------------------------------

char * get_input(char * str, size_t len)
{
    printf("? ");
    return fgets(str, len, stdin);
}

int main(int argc, char* argv[])
{
    Simul simul;
    Glossary arg;
    FrameReader reader;
    size_t frm = 0;

    if ( arg.read_strings(argc-1, argv+1) )
        return EXIT_FAILURE;

    if ( arg.use_key("help") )
    {
        help(std::cout);
        instructions();
        return EXIT_SUCCESS;
    }
    
    if ( arg.has_key("directory") )
        FilePath::change_dir(arg.value("directory"));

    std::string input = Simul::TRAJECTORY;
    std::string output = "system.cmo";
    
    arg.set(output, "output");
    arg.set(input, "input") || arg.set(input, ".cmo");

    int binary = 1;
    arg.set(binary, "binary");
    
    try {
        simul.loadProperties();
        reader.openFile(input);
    }
    catch( Exception & e )
    {
        std::clog << e.brief() << '\n';
        return EXIT_FAILURE;
    }
    
    if ( reader.badFile() )
    {
        printf("File could not be oppened\n");
        return EXIT_FAILURE;
    }
    
    if ( arg.set(frm, "frame") )
    {
        if ( reader.loadFrame(simul, frm) )
            printf("Frame %lu could not be loaded\n", frm);
        else
            printf("loaded frame %li, time %.3f:\n", reader.currentFrame(), simul.time());
    }
    
    instructions();
    char cmd[1024] = "\0";
    while ( get_input(cmd, sizeof(cmd)) )
    {
        if ( isdigit(cmd[0]))
        {
            char * end = nullptr;
            size_t f = strtoul(cmd, &end, 10);
            if ( errno )
                printf("Reader: syntax error in `%s'\n", cmd);
            else if ( end > cmd )
            {
                try {
                    if ( 0 != reader.loadFrame(simul, f) )
                    {
                        printf("Reader: frame %lu not found.", f);
                        reader.clear();
                    }
                    else
                        frm = reader.currentFrame();
                }
                catch( Exception & e ) {
                    printf("reader.loadFrame(%lu) exception: %s\n", f, e.what());
                }
            }
        }
        else switch( cmd[0] )
        {
            case '\n':
            case 'n':
                try {
                    int err = reader.loadNextFrame(simul);
                    if ( err ) printf("reader.loadNextFrame error: %i\n", err);
                }
                catch( Exception & e ) {
                    printf("reader.loadNextFrame exception: %s\n", e.what());
                }
                break;
            case 'w':
                simul.writeObjects(output, true, binary);
                printf("%s <--- frame %lu\n", output.c_str(), reader.currentFrame());
                continue;
            case 'e':
                simul.eraseObjects();
                simul.eraseProperties();
                break;
            case 'b':
                binary = !binary;
                printf("Reader: binary = %i\n", binary);
                continue;
            case 'c':
                reader.clearPositions();
                continue;
            case 'r':
                reader.rewind();
                break;
            case 'q': case 'Q': case 27:
                return 0;
            default:
                printf("I do not understand `%c`\n", cmd[0]);
                break;
        }
        inventory(simul, reader.currentFrame());
    }
}
