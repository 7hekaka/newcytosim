// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef TOKENIZER_H
#define TOKENIZER_H

#include "assert_macro.h"
#include <iostream>
#include <string>
#include <vector>

/// elementary tokenizer
/**
 The Tokenizer is used to cut a character stream into words,
 and to interpret Cytosim's configuration file
 */
namespace Tokenizer
{
    /// return closing delimiter corresponding to `arg`, or 0 if this is not a known delimiter
    char block_delimiter(char arg);
    
    /// same as block_delimiter()
    //inline char block_delimiter(int c) { return block_delimiter((char)c); }
    
    /// discard characters until the next new-line character
    int skip_line(std::istream&);

    /// skip space and new-line if `eat_line==true`, returning the next character
    int skip_space(std::istream&, bool eat_line);

    /// skip upcomming space-like characters, including new-line if `eat_line`==true
    int get_character(std::istream&, bool eat_space=true, bool eat_line=false);
    
    /// read one signed integer, and throw exception if spacing character does not follow
    bool get_integer(std::istream&, long&);
    
    /// read one unsigned integer, and throw exception if spacing character does not follow
    bool get_integer(std::istream&, unsigned long&);
    
    /// split string using the given separator.
    std::vector<std::string> split(std::string&, char sep, bool get_empty_fields);
    
    /// read symbol from string if possible, truncating string if successful
    std::string split_symbol(std::string&);

    /// split a polysymbol 'symbol:num' or 'work:num' such as 'filament1'
    bool split_polysymbol(std::string&, long& num);

    /// read integer from string if possible, truncating string if successful
    bool split_integer(long&, std::string&);
    
    /// read unsigned integer from string if possible, truncating string if successful
    bool split_integer(unsigned long&, std::string&);
    
    /// read unsigned integer from string if possible, truncating string if successful
    bool split_integer(unsigned&, std::string&);
    
    /// read multiple forms of integer numbers
    std::string get_integer(std::istream&);

    /// return next token if it looks like a number
    std::string get_real(std::istream&);

    /// return next token if it looks like a hexadecimal number
    std::string get_hexadecimal(std::istream&);

    /// return next token if it looks like a variable name
    std::string get_symbol(std::istream&, bool eat_line=false);

    /// return next token if it looks like a variable name
    std::string get_symbol(std::string const&, size_t&);

    /// return 'true' if next symbol matches 'arg'
    bool has_symbol(std::string const&, size_t&, std::string const& arg, bool eat_line=false);

    /// return next token composed of multiple symbols concatenated with ':'
    std::string get_polysymbol(std::istream&, bool eat_line=false);

    /// return next token if it looks like a file name
    std::string get_path(std::istream&, bool eat_line=false);
    
    /// return next token on the same line (this does not eat new-line characters)
    std::string get_token(std::istream&);
    
    /// accumulate characters until the next new-line character
    std::string get_line(std::istream&);

    /// read text until `c_out` is encountered, assuming `c_in` was already read
    std::string get_block_text(std::istream&, char c_in, char c_out);
    
    /// read text until `c_out` is encountered, assuming `c_in` was already read
    std::string get_blocked_text(std::string const&, size_t&, char c_in, char c_out);

    /// skip spaces and read a block delimited by `c_in`, or return empty string if `c_in` is not found
    std::string get_block(std::istream&, char c_in, bool or_die=false);
    
    /// read a delimited set of characters, return block with delimiters included
    std::string get_block(std::istream&);
    
    /// remove `c_in` and matching delimiters from string, return false if not present
    bool strip_block(std::string&, char c_in);

    /// return text read until `what` is found, stoping immediately before
    std::string get_until(std::istream&, std::string what);

    /// remove characters present in `ws` found at both edges of `str`
    std::string trim(std::string const& str, const std::string& ws = " \t\n");
    
}

#endif


