// Cytosim was created by Francois Nedelec. Copyright Cambridge University 2024
// Created by Francois Nedelec on 03/09/2014; +ncurses 08/10/2024.


#include "print_color.h"
#include <cstdio>
#include <cstdlib>
#include <unistd.h>


#ifndef PRINT_IN_COLOR
#   define PRINT_IN_COLOR 1
#endif

#if PRINT_IN_COLOR
#   include <curses.h>
#endif

// these are ANSI escape sequence for UNIX-based systems
#define KNRM  "\x1B[0m"
#define KBLD  "\x1B[1m"
#define KUND  "\x1B[4m"
#define KREV  "\x1B[7m"

#define KBLK  "\x1B[30m"
#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define KYEL  "\x1B[33m"
#define KBLU  "\x1B[34m"
#define KMAG  "\x1B[35m"
#define KCYN  "\x1B[36m"
#define KWHT  "\x1B[37m"

#define KBLDRED  "\x1B[1m\x1B[31m"
#define KBLDGRN  "\x1B[1m\x1B[32m"
#define KBLDYEL  "\x1B[1m\x1B[33m"
#define KBLDBLU  "\x1B[1m\x1B[34m"
#define KBLDMAG  "\x1B[1m\x1B[35m"
#define KBLDCYN  "\x1B[1m\x1B[36m"
#define KBLDWHT  "\x1B[1m\x1B[37m"

// using a macro to concatenate strings:
#define FORMAT_STRING(COLOR) COLOR "%s" KNRM


/* Returns the number of colors supported by the terminal */
int num_colors(int fd)
{
    int res = 0;
#if PRINT_IN_COLOR
    if ( !isatty(fd) )
    {
        //printf("[FILE %i is not a TTY]", fd);
        return 0;
    }
    stdscr = initscr();
    if ( stdscr != nullptr )
    {
        start_color();
        if ( has_colors() )
            res = COLOR_PAIRS;
    }
    refresh();
    endwin();
    //printf("[FILE %i has %i colors]", fd, res);
#endif
    return res;
 }


/* Check if terminal supports color printing */
inline static bool terminal_has_colors(FILE * file)
{
    return num_colors(fileno(file)) > 7;
}


void print_red(FILE * file, std::string const& str)
{
#if PRINT_IN_COLOR
    if ( terminal_has_colors(file) )
        fprintf(file, FORMAT_STRING(KBLDRED), str.c_str());
    else
#endif
        fputs(str.c_str(), file);
}

void print_green(FILE * file, std::string const& str)
{
#if PRINT_IN_COLOR
    if ( terminal_has_colors(file) )
        fprintf(file, FORMAT_STRING(KBLDGRN), str.c_str());
    else
#endif
        fputs(str.c_str(), file);
}

void print_yellow(FILE * file, std::string const& str)
{
#if PRINT_IN_COLOR
    if ( terminal_has_colors(file) )
        fprintf(file, FORMAT_STRING(KBLDYEL), str.c_str());
    else
#endif
        fputs(str.c_str(), file);
}

void print_blue(FILE * file, std::string const& str)
{
#if PRINT_IN_COLOR
    if ( terminal_has_colors(file) )
        fprintf(file, FORMAT_STRING(KBLDBLU), str.c_str());
    else
#endif
        fputs(str.c_str(), file);
}

void print_magenta(FILE * file, std::string const& str)
{
#if PRINT_IN_COLOR
    if ( terminal_has_colors(file) )
        fprintf(file, FORMAT_STRING(KBLDMAG), str.c_str());
    else
#endif
        fputs(str.c_str(), file);
}

void print_cyan(FILE * file, std::string const& str)
{
#if PRINT_IN_COLOR
    if ( terminal_has_colors(file) )
        fprintf(file, FORMAT_STRING(KBLDCYN), str.c_str());
    else
#endif
        fputs(str.c_str(), file);
}

void print_bold(FILE * file, std::string const& str)
{
#if PRINT_IN_COLOR
    if ( terminal_has_colors(file) )
        fprintf(file, FORMAT_STRING(KBLD), str.c_str());
    else
#endif
        fputs(str.c_str(), file);
}


