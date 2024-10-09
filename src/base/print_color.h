// Cytosim was created by Francois Nedelec. Copyright Cambridge University 2024
// FJN 03.09.2014.

/**
 This defines functions to print text in color on terminals that support it
 It is based on the original ANSI specification that defined 8 standard colors:
 https://en.wikipedia.org/wiki/ANSI_escape_code
 
 This can be disabled by setting `PRINT_IN_COLOR` to 0 in `print_color.cc`
 
 To detect the number of colors supported by the terminal, this uses `ncurses`:
 https://en.wikipedia.org/wiki/Ncurses
 */

#include <cstdio>
#include <string>


/// return number of colors supported by terminal associated with file descriptor.
int num_colors(int fildes);

/// print text in red
void print_red(FILE*, std::string const&);

/// print text in green
void print_green(FILE*, std::string const&);

/// print text in yellow
void print_yellow(FILE*, std::string const&);

/// print text in blue
void print_blue(FILE*, std::string const&);

/// print text in magenta
void print_magenta(FILE*, std::string const&);

/// print text in cyan
void print_cyan(FILE*, std::string const&);

/// print text in bold
void print_bold(FILE*, std::string const&);
