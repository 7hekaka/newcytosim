// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

/**
 This code is outdated, please refer to report.cc
 */

#include "simul.h"
#include "frame_reader.h"
#include "iowrapper.h"
#include "messages.h"
#include "exceptions.h"
#include "parser.h"

Simul sim;

enum Stage { INIT, COUNT, PRINT };

//=====================================================================
void myAnalysis(const Stage action, const int frame=0)
{
    switch( action )
    {
        case INIT:
        {
            //init your variables
            //frame is not defined
        } break;
            
        case COUNT:
        {
            //analyse, and sum up
            //frame is the index of the current frame in the file
        } break;
            
        case PRINT:
        {
            //print a report
            //nbFrames indicates the total number of frames given
        } break;
    }
}


//===================================================================

void countMTtips(const Stage action, const int frame=0)
{
    static Grid <real, 1> density;
    static Grid <real, 1> denR, denL;

    switch( action )
    {
        case INIT:
        {
            real xdim = sim.spaces.master()->max_extension();
            density.create1D(-xdim, xdim, 100);
            denR.create1D(-xdim, xdim, 100);
            denL.create1D(-xdim, xdim, 100);
        } break;
            
        case COUNT:
        {
            denR.setValues(0);
            denL.setValues(0);
            for ( Fiber* M=sim.fibers.first(); M ; M=M->next() )
            {
                real X = M->posEndM().XX;
                ++density(&X);
                if ( M->dirEndM().XX > 0 )
                    ++denR(&X);
                else
                    ++denL(&X);
            }
        } break;
            
        case PRINT: {
            density.printValuesWithRange(stdout);
        } break;
    }
}

//===================================================================
void countMTends(const Stage action, const int frame=0)
{
    static Grid<real, 1> density;
    
    switch( action )
    {
        case INIT:
        {
            density.create1D(-6, 6, 24);
        } break;
            
        case COUNT:
        {
            for ( Fiber* mt=sim.fibers.first(); mt ; mt=mt->next() ) {
                Vector v = mt->posEndM();
                ++density(&v.XX);
            }
        } break;
            
        case PRINT:
        {
            density.printValuesWithRange(stdout);
        } break;
    }
}

//===================================================================
void countBentFibers(const Stage action, const int frame=0)
{
    static int res = 0;
    
    switch( action )
    {
        case INIT:
        {
            res = 0;
        } break;
            
        case COUNT:
        {
            for ( Fiber* mt=sim.fibers.first(); mt ; mt=mt->next() ) {
                if ( dot(mt->dirEndP(), mt->dirEndM()) < 0 )
                    ++res;
            }
        } break;
            
        case PRINT:
        {
            printf("%% curling: %i events in %i frames\n", res, frame);
            printf("%.4f\n", res/float(frame));
        } break;
    }
}


//=====================================================================
void countMTlength(const Stage action, const int frame=0)
{
    static Grid<real, 1> mtlength;

    switch( action )
    {
        case INIT:
        {
            mtlength.create1D(0, 40, 100);
            mtlength.setValues(0);
        } break;
            
        case COUNT:
        {
            for ( Fiber * mt=sim.fibers.first(); mt ; mt=mt->next() )
            {
                real x = mt->length();
                ++mtlength( &x );
            }
        } break;

        case PRINT:
        {
            FILE* outfile = fopen("mtl.out","w");
            mtlength.printValuesWithRange(outfile);
            fclose(outfile);
        } break;
    }
}

//=====================================================================

void countCouples(const Stage action, const int frame=0)
{
    const int mx = 4;
    static int cnt[3][mx];
    
    switch( action ) {
        case INIT: {
            for ( int ii = 0; ii < 3; ++ii )
                for ( int jj = 0; jj < mx; ++jj )
                    cnt[ii][jj] = 0;
        } break;
            
        case COUNT: {
            //add by state and type:
            Couple * cxi;
            
            for ( cxi=sim.couples.firstFF(); cxi ; cxi = cxi->next() ) {
                unsigned pi = cxi->prop->number();
                if ( pi < mx ) ++cnt[0][pi];
            }
            
            for ( cxi=sim.couples.firstAF(); cxi ; cxi = cxi->next() ) {
                unsigned pi = cxi->prop->number();
                if ( pi < mx ) ++cnt[1][pi];
            }

            for ( cxi=sim.couples.firstFA(); cxi ; cxi = cxi->next() ) {
                unsigned pi = cxi->prop->number();
                if ( pi < mx ) ++cnt[1][pi];
            }
            
            for ( cxi=sim.couples.firstAA(); cxi ; cxi = cxi->next() ) {
                unsigned pi = cxi->prop->number();
                if ( pi < mx ) ++cnt[2][pi];
            }
        } break;
            
        case PRINT: {
            for ( int t=0; t<mx; ++t ) {
                printf("%i  %10i %10i %10i\n", t, cnt[0][t], cnt[1][t], cnt[2][t]);
            }
        } break;
    }
}

//=====================================================================

void countSingle(const Stage action, const int frame=0)
{
    switch( action ) {
        case INIT: {
        } break;
            
        case COUNT: {
        } break;
            
        case PRINT: {
        } break;
    }
}

//=====================================================================

void help()
{
    printf("Cytosim / analyse:\n");
    printf("  Prints cumulative statistics on the objects\n");
    printf("USAGE: analyse [ -f(index) ] [-f(index1-index2)] [ command_name1 command_name2 ... ]\n");
    printf("\n");
    printf("Commands supported:\n");
    printf(" fibers           report Fiber lengths and states\n");
    printf(" singles          report Single numbers in each state\n");
    printf(" couples          report Couples numbers in each state\n");
    printf(" bent             report curling fibers\n");
    printf("\n");
    printf("analyse2  DIM=%i, compiled %s, %s\n", DIM, __TIME__, __DATE__ );
}

//=====================================================================

std::string commands;
unsigned int first_frame = 0;
unsigned int last_frame  = 1e6;


int parseCommandLine(char arg[])
{
    if ( '-' == arg[0] ) {
        switch( arg[1] ) {
            case 'h': help(); exit(EXIT_SUCCESS);
            case 'f':   //frame index
                if ( arg[2] ) {
                    if ( 1 == sscanf(arg+2, "%u-%u", &first_frame, &last_frame )) {
                        if ( arg[strlen(arg)-1] != '-' )
                            last_frame = first_frame;
                    }
                }
                return 0;
        }
    }
    return 1;
}

//=====================================================================
int main(int argc, char* argv[])
{
    if ( argc < 2 ) {
        help();
        return EXIT_SUCCESS;
    }
    
    std::string cmd;
    for ( int ii = 1; ii < argc; ++ii ) {
        if ( 0 != parseCommandLine(argv[ii]) ) {
            cmd.append(" ");
            cmd.append(argv[ii]);
            cmd.append(" ");
        }
    }
    
    Cytosim::all_silent();

    FrameReader reader;
    try {
        sim.loadProperties();
        reader.openFile(TRAJECTORY);
    }
    catch( Exception & e ) {
        std::cerr << "Aborted: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    
    const char * commands = cmd.c_str();
    //initialization:
    if ( strstr(commands, " stat " ))         myAnalysis(INIT);
    if ( strstr(commands, " fibers  " ))      countMTlength(INIT);
    if ( strstr(commands, " couples " ))      countCouples(INIT);
    if ( strstr(commands, " singles " ))      countSingle(INIT);
    if ( strstr(commands, " minus " ))        countMTtips(INIT);
    if ( strstr(commands, " ends " ))         countMTends(INIT);
    if ( strstr(commands, " bent " ))         countBentFibers(INIT);
    
    //counting:
    int nbFrames = 0;
    for ( unsigned int frame = first_frame; frame <= last_frame ; ++frame )
    {
        try
        {
            if ( 0 == reader.loadFrame(sim, frame) )
            {
                ++nbFrames;
                //printf("frame %3i\n", frame);
                if ( strstr(commands, " stat " ))      myAnalysis(COUNT, frame);
                if ( strstr(commands, " fibers " ))    countMTlength(COUNT, frame);
                if ( strstr(commands, " couples " ))   countCouples(COUNT, frame);
                if ( strstr(commands, " singles " ))   countSingle(COUNT, frame);
                if ( strstr(commands, " minus " ))     countMTtips(COUNT, frame);
                if ( strstr(commands, " ends " ))      countMTends(COUNT, frame);
                if ( strstr(commands, " bent " ))      countBentFibers(COUNT, frame);
            }
            else
            {
                if ( ! reader.good() ) break;
                std::cerr << "Error: missing frame " << frame << std::endl;
            }

        }
        catch( Exception & e )
        {
            std::cerr << "Error: " << e.message() << std::endl;
        }
    }
    
    //finishing:
    printf("%% %i frames were analyzed:\n", nbFrames);
    if ( strstr(commands, " stat " ))       myAnalysis(PRINT, nbFrames);
    if ( strstr(commands, " fibers " ))     countMTlength(PRINT, nbFrames);
    if ( strstr(commands, " couples " ))    countCouples(PRINT, nbFrames);
    if ( strstr(commands, " singles " ))    countSingle(PRINT, nbFrames);
    if ( strstr(commands, " minus " ))      countMTtips(PRINT, nbFrames);
    if ( strstr(commands, " ends " ))       countMTends(PRINT, nbFrames);
    if ( strstr(commands, " bent " ))       countBentFibers(PRINT, nbFrames);

    
    return EXIT_SUCCESS;
}
