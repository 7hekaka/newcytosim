// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University
/*
 This controls cytosim's play by sending commands via a direct pipe
 Author: FJ Nedelec, 6.03.2018 -- 01.04.2022
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <signal.h>
#include <string.h>
#include <errno.h>
#include <cmath>
#include <sys/param.h>


/// child process id:
pid_t child = 0;

/// file descriptor for a pipe: [0]=READ, [1]=WRITE
int fds[2];


void signal_handler(int sig)
{
    ssize_t __attribute__((unused)) u;
    u = write(STDERR_FILENO, "\n* * * * *\n", 11);
    psignal(sig, "CYMWIZ");
    u = write(STDERR_FILENO, "* * * * *\n", 10);
}


int start(const char* path, char *const command[])
{
    // create a unidirectional pipe:
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
        close(fds[1]);  // close pipe exit
        // map the standard-input to the pipe exit:
        while ((dup2(fds[0], STDIN_FILENO) == -1) && (errno == EINTR)) {}
        // execv() should not return, except if error occurred
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
        _exit(1);
    }
    
    // this code executed by the parent process
    close(fds[0]);    // close pipe entry
    return 0;
}


void stop()
{
    close(fds[1]);
    kill(child, SIGTERM);
    child = 0;
}


// build command suitable to cytosim
int command(char cmd[], size_t len, int num)
{
    float A = num * M_PI / 180.0;
    float R = 3.0;
    float X = R * cos(A);
    float Y = R * sin(A);
    return snprintf(cmd, len, "set handle clamp {%.2f %.2f}\n", X, Y);
}


int main(int argc, char* argv[])
{
    char path[PATH_MAX] = { 0 };
    char cmd[1024] = { 0 };
    ssize_t len = 0;

    if ( argc < 2 )
    {
        printf("Invoke 'cymwiz' with executable and arguments, eg:");
        printf("   cymwiz bin/play config.cym\n");
        return 1;
    }
    
    // resolve full path of executable:
    if ( 0 == realpath(argv[1], path) )
    {
        printf("error: could not resolve executable `%s'\n", argv[1]);
        return 1;
    }
    
    if ( signal(SIGPIPE, signal_handler) )
        write(STDERR_FILENO, "no SIGPIPE handler\n", 19);

    if ( 0 == start(path, argv+1) )
    {
        sleep(1);
        
        // start controlling:
        for ( int a = 0; a < 100*360; ++a )
        {
            len = command(cmd, sizeof(cmd), a);
            
            if ( len > 0 )
            {
                //printf(">>> %s", cmd);
                // send string through pipe entry:
                ssize_t s = write(fds[1], cmd, len);
                if ( s != len )
                {
                    printf("Error: pipe is broken\n");
                    break;
                }
            }
            usleep(50000);
        }
    }
    printf("controller ended\n");
}

