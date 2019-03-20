/**
 Cytosim control using a USB-connected MIDI board, with several scenarios
 For "Nuit Blanche" in Paris, 6.10.2018
 Gaelle LETORT and Francois NEDELEC
*/

#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <errno.h>
#include <math.h>
#include <signal.h>
#include <unistd.h>
#include <vector>

// Using Real-Time MIDI by Gary Scavone:
#include "RtMidi.h"

/// which scenario is running [0 ... 4]
unsigned mode = 0;

/// child process id:
pid_t child = 0;

/// file descriptor for a pipe: [0]=READ, [1]=WRITE
int fds[2] = { 0 };


void start(const char* path, char *const command[])
{
    // create a unidirectional pipe:
    if ( pipe(fds) < 0 )
    {
        perror("pipe");
        exit(1);
    }
    
    // create a child process
    child = fork();
    
    if ( child == -1 )
    {
        perror("fork");
        exit(1);
    }
    
    if ( child == 0 )
    {
        // this code executed by the child process
        // child closes pipe exit:
        close(fds[1]);
        // map the standard-input to the pipe exit:
        while ((dup2(fds[0], STDIN_FILENO) == -1) && (errno == EINTR)) {}
        // run executable (this should not return, except if error occurred)
        execv(path, command);
        // the command failed, and error is indicated by 'errno':
        perror("execl");
        write(STDERR_FILENO, "while executing command: ", 25);
        write(STDERR_FILENO, path, strlen(path));
        for ( int i = 1; command[i]; ++i )
        {
            write(STDERR_FILENO, " ", 1);
            write(STDERR_FILENO, command[i], strlen(command[i]));
        }
        write(STDERR_FILENO, "\n", 1);
        _exit(1);
    }
    
    // this code executed by the parent process
    // close pipe entry:
    close(fds[0]);
}


void stop(int)
{
    close(fds[1]);
    kill(child, SIGTERM);
    fds[0] = 0;
    fds[1] = 0;
    child = 0;
}


// start cytosim using the pipe
void start(unsigned new_mode)
{
    const size_t LEN = 1024;
    char mem[6*LEN] = { 0 };
    char* args[6] = { mem, mem+LEN, mem+2*LEN, mem+3*LEN, mem+4*LEN, mem+5*LEN };

    // keep mode in [0, 4]
    mode = new_mode % 5;
    
    snprintf(args[0], LEN, "play%i", mode);
    snprintf(args[1], LEN, "config%i.cym", mode);
    snprintf(args[2], LEN, "live");
    snprintf(args[3], LEN, "fullscreen=1");
    args[4] = 0;
    args[5] = 0;
    
#if ( 0 )
    //This works when executables are located in the current working directory:
    start(args[0], args);
#else
    // use a full path to spawn child
    char path[LEN] = { 0 };
    char cwd[LEN] = { 0 };

    // get current working directory
    getcwd(cwd, LEN);
    
    // build the full path to executable:
    snprintf(path, LEN, "%s/%s", cwd, args[0]);
    
    if ( 0 )
    {
        printf("    > in %s\n", cwd);
        printf("    > %s %s %s %s\n", path, args[1], args[2], args[3]);
    }
    
    start(path, args);
#endif
}


void restart(unsigned arg)
{
    stop(0);
    start(arg);
}

//------------------------------------------------------------------------------

/// max value of the sliders
const float maxValue = 127.0;

/// Converts MIDI value (from 0 to 127) to [min, max] range, linear dependency
float linear(int val, float min, float max)
{
    return ( min + (max-min) * ( (float)val / maxValue) );
}

/// Converts MIDI value (from 0 to 127) to [min, max] range, x^2 dependency
float quadratic(int val, float min, float max)
{
    return ( min + (max-min) * float(val*val)/(maxValue*maxValue) );
}


int makeCommand(char * str, size_t len, int slider, int value)
{
    switch ( slider )
    {
        case 1:    case 95:
            return snprintf(str, len, "change system { viscosity=%f }\n", linear(value, 0.1, 1.0));
        case 2:    case 96:
            return snprintf(str, len, "change filament { rigidity=%f }\n", quadratic(value, 0.1, 20));
        case 3:    case 97:
            return snprintf(str, len, "change filament { growing_speed=%f }\n", linear(value, 0, 0.2));
        case 4:    case 100:
            return snprintf(str, len, "change motor { unloaded_speed=%f }\n", linear(value, -0.8, 0.8));
        case 5:    case 101:
            return snprintf(str, len, "change dynein { unloaded_speed=%f }\n", linear(value, -0.8, 0.8));
        case 6:    case 102:
            return snprintf(str, len, "change centrosome { stiffness=1000, %f }\n", linear(value, 1, 1000));
        case 7:    case 103:
            return snprintf(str, len, "change complex { stiffness=%f }\n", linear(value, 0, 256));
        case 8:    case 106:
        {
            float w = linear(value, 0.54, 1.0);
            float h = ( 0.85 * 0.85 * 0.21 ) / ( w * w );
            return snprintf(str, len, "change all space { dimensions= %f %f %f }\n", w, w, h);
        }
        default:
            printf("Slider %i (value %i) not implemented\n", slider, value);
            break;
    }
    return 0;
}


void printMessage( std::vector<unsigned char>& mes )
{
    for ( size_t i = 0; i < mes.size(); ++i )
        printf("\t %i",(int) mes[i]);
    printf("\n");
}


void goLive(RtMidiIn& midi, unsigned mode)
{
    std::vector<unsigned char> msg;
    char cmd[4096];

    start(mode);

    while ( child )
    {
        midi.getMessage(&msg);
        // has a recieved a message from MIDI input, of length 3
        if ( msg.size() == 3 )
        {
            //printMessage(msg);
            // first number ?

            // which slider is changed
            int slider = (int) msg[1];
            // slider's value
            int value = (int) msg[2];

            /**
             Check buttons dedicated to restart the simulation, directly with
             different scenario. These buttons send two MIDI signals, as
             the key is pressed down and when it is released.
             */
            if ( slider == 48 ) { if ( !value ) restart(0); continue; }
            if ( slider == 50 ) { if ( !value ) restart(1); continue; }
            if ( slider == 52 ) { if ( !value ) restart(2); continue; }
            if ( slider == 53 ) { if ( !value ) restart(3); continue; }
            if ( slider == 55 ) { if ( !value ) restart(4); continue; }
            
            // make the command suitable to cytosim
            int len = makeCommand(cmd, sizeof(cmd), slider, value);
            if ( len > 0 )
            {
                //write(STDOUT_FILENO, "> ", 2); write(STDOUT_FILENO, cmd, sizeof(cmd));
                
                // send string through pipe entry:
                ssize_t s = write(fds[1], cmd, len);
                if ( s != len )
                {
                    printf("Error: pipe is broken\n");
                    break;
                }
            }
        }
        else
        {
            printf("Unexpected MIDI message:\n");
            printMessage(msg);
            break;
        }

        // Sleep for 5 milliseconds.
        usleep(5000);
    }
}

//------------------------------------------------------------------------------

void usage()
{
    printf("\nUsage: cytomaster PORT MODE\n");
    printf("     PORT = the MIDI device to use\n");
    printf("     MODE = [ 0, 1, 2, 3, 4 ]\n\n");
    printf(" use `cytomaster 999' to list all known MIDI ports\n\n");
}


void namePorts(RtMidiIn& midi)
{
    unsigned int np = midi.getPortCount();
    
    if ( np < 1 )
        printf("no MIDI port detected!\n");
    else
    {
        for ( unsigned int p = 0; p < np; ++p )
        {
            midi.openPort(p);
            printf("   MIDI port %u is [ %s ]\n", p, midi.getPortName().c_str());
            midi.closePort();
        }
    }
}


int main( int argc, char *argv[] )
{
    if ( argc < 2 || 3 < argc )
    {
        usage();
        return 1;
    }
    
    RtMidiIn midi(RtMidi::MACOSX_CORE);
    
    // Check available ports vs. specified.
    unsigned int port = (unsigned)atoi(argv[1]);
    unsigned int nPorts = midi.getPortCount();
    
    if ( nPorts < 1 )
    {
        printf("no MIDI port detected!\n");
        return 1;
    }
    if ( port >= nPorts )
    {
        printf("Invalid port specified (only %u ports detected)\n", nPorts);
        namePorts(midi);
        return 1;
    }
    
    try
    {
        midi.openPort(port);
    }
    catch ( RtMidiError &error )
    {
        error.printMessage();
        return 0;
    }
    
    // Don't ignore sysex, timing, or active sensing messages.
    midi.ignoreTypes(false, false, false);
    
    // Install an interrupt handler function.
    signal(SIGINT, stop);
    
    // Periodically check input queue.
    printf("Listening to [ %s ] ... quit with Ctrl-C.\n", midi.getPortName().c_str());
    
    if ( 2 < argc )
        mode = (unsigned)atoi(argv[2]);
    
    goLive(midi, mode);
    
    return 0;
}
