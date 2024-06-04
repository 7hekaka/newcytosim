// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University
/*
 This is a test for controlling cytosim's play by sending commands via the standard input
 This program initiates a child process to run 'play' and send commands via a Pipe
 Author: F. Nedelec, 6.03.2018 -- 19.05.2018
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <signal.h>
#include <string.h>
#include <errno.h>
#include <cmath>
#include <sys/param.h>

#include "gym_color.h"

/// child process id:
pid_t child = 0;

/// file descriptor for a pipe: [0]=READ, [1]=WRITE
int fds[2];


int start(const char* path, char *const command[])
{
    // create communication pipe:
    if ( pipe(fds) < 0 )
    {
        perror("pipe");
        errno = 0;
        return 1;
    }
    
    // create a child process
    child = fork();
    
    if ( child == -1 )
    {
        perror("fork");
        errno = 0;
        return 1;
    }
    
    if ( child == 0 )
    {
        // this code executed by the child process
        // child closes pipe exit:
        close(fds[1]);
        // map the standard-input to the pipe exit:
        while ((dup2(fds[0], STDIN_FILENO) == -1) && (errno == EINTR)) {}
        // execv() should not return, except in case of error
        execv(path, command);
        // the command failed, and error is indicated by 'errno':
        perror("execv");
        ssize_t __attribute__((unused)) u;
        u = write(STDERR_FILENO, "while executing command:", 24);
        for ( int i = 0; command[i]; ++i )
        {
            u = write(STDERR_FILENO, " ", 1);
            u = write(STDERR_FILENO, command[i], strlen(command[i]));
        }
        u = write(STDERR_FILENO, "\n", 1);
        errno = 0;
        _exit(1);
    }
    
    // this code executed by the parent process
    // close pipe entry:
    close(fds[0]);
    return 0;
}


void stop()
{
    close(fds[1]);
    kill(child, SIGTERM);
    child = 0;
}


// this one might give you a headache:
int command_zoom(char cmd[], size_t len, int num)
{
    float z = 1.0f + 0.25f * cosf(num * M_PI / 180.0);
    return snprintf(cmd, len, "change all simul display { zoom=%.3f; }\n", z);
}

// build command suitable to cytosim
int command(char cmd[], size_t len, int num)
{
    float angle = num * M_PI / 180.0;
    if ( num & 1 )
    {
        // change zoom, to test smoothness of response
        float z = 1.0f + 0.25f * cosf(angle*0.25f);
        return snprintf(cmd, len, "change all simul display { zoom=%.3f; }\n", z);
    }
    else
    switch (( num >> 1 ) % 6 )
    {
        case 0:
        {
            gym_color col = gym_color::hue_color(angle);
            std::string str = col.to_string();
            return snprintf(cmd, len, "change all simul display {back_color=%s}\n", str.c_str());
        }
        case 1:
            return snprintf(cmd, len, "change kinesin { unloaded_speed=%.2f }\n", std::cos(angle));
        case 2:
        {
            int s = round(2+14*std::fabs(cosf(angle)));
            return snprintf(cmd, len, "change all hand display { size=%i; }\n", s);
        }
        case 3:
        {
            gym_color col = gym_color::hue_color(M_PI+angle);
            std::string str = col.to_string();
            return snprintf(cmd, len, "change all fiber display { color=%s }\n", str.c_str());
        }
        case 4:
        {
            float s = 0.5*round(16+15*cosf(angle));
            return snprintf(cmd, len, "change all fiber display { line=%.1f; }\n", s);
        }
        case 5:
        {
            //return snprintf(cmd, len, "change all simul display {save_images=1}\n");
            float s = 0.5*round(16+15*cosf(angle+M_PI));
            return snprintf(cmd, len, "change all fiber display { points=%.1f, 1; };\n", s);
        }
    }
    return 0;
}


int main(int argc, char* argv[])
{
    char path[PATH_MAX] = { 0 };
    char cmd[1024] = { 0 };
    ssize_t len = 0;

    if ( argc < 2 )
    {
        printf("Invoke 'test_pipe' with executable and arguments, eg:");
        printf("   test_pipe bin/play config.cym\n");
        return 1;
    }
    
    // resolve full path of executable:
    if ( 0 == realpath(argv[1], path) )
    {
        printf("error: could not resolve executable `%s'\n", argv[1]);
        return 1;
    }
    
    int err = start(path, argv+1);
    if ( err ) return err;
    sleep(1);

    // start controlling:
    for ( int a = 0; a < 100*360; ++a )
    {
        len = command(cmd, sizeof(cmd), a);
    
        if ( a == 360*2 )
        {
            // restart child once
            stop();
            start(path, argv+1);
        }
        
        if ( len > 0 )
        {
            //write(STDOUT_FILENO, "SENDING ", 8);
            //write(STDOUT_FILENO, cmd, sizeof(cmd));
            
            // send string through pipe entry:
            ssize_t s = write(fds[1], cmd, len);
            if ( s != len )
            {
                printf("Error: pipe is broken\n");
                break;
            }
        }
        usleep(20000);
    }
    
    //stop();
    printf("controller terminated\n");
}

