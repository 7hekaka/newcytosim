// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
// Created by Francois Nedelec on 1/7/2009.

#include "assert_macro.h"
#include "tokenizer.h"
#include "exceptions.h"
#include <errno.h>

// Use second definition to trace execution
#define VLOG(ARG) ((void) 0)
//#define VLOG(ARG) std::clog << ARG;


char Tokenizer::block_delimiter(char c)
{
    switch (c)
    {
        case '(': return ')';
        case '{': return '}';
        case '[': return ']';
        case '"': return '"';
    }
    return 0;
}

//------------------------------------------------------------------------------

std::string Tokenizer::get_line(std::istream& is)
{
    std::string res;
    res.reserve(1024);
    std::getline(is, res);
    return res;
}


int Tokenizer::skip_line(std::istream& is)
{
    int c;
    do {
        c = is.get();
    } while ( c != '\n' && c != EOF );
    return is.peek();
}


int Tokenizer::skip_space(std::istream& is, bool eat_line)
{
    int c = is.peek();
    while ( isspace(c) )
    {
        if (( c == '\n' ) & !eat_line )
            break;
        is.get();
        c = is.peek();
    }
    return c;
}


int Tokenizer::get_character(std::istream& is, bool eat_space, bool eat_line)
{
    int c = 0;
    do {
        c = is.get();
        if ( c == EOF )
            return EOF;
        if (( c == '\n' ) & ! eat_line )
            break;
    } while ( eat_space && isspace(c) );
    return c;
}

//------------------------------------------------------------------------------
#pragma mark -


int valid_symbol(int c)
{
    return isalnum(c) | ( c == '_' );
}

int valid_path_start(int c)
{
    return isalpha(c) | (c=='/') | (c=='\\') | (c=='.');
}

int valid_path(int c)
{
    return isalnum(c) | (c=='_') | (c=='-') | (c=='/') | (c=='\\') | (c=='.') | (c==':');
}

int valid_hexadecimal(int c)
{
    return isxdigit(c) | ( c=='x' );
}

std::string get_stuff(std::istream& is, int (*valid)(int))
{
    std::string res;
    int c = is.peek();
    while ( valid(c) )
    {
        res.push_back(is.get());
        c = is.peek();
    }
    return res;
}

std::string get_stuff(std::string const& arg, size_t& sci, int (*valid)(int))
{
    size_t pos = sci;
    while ( valid(arg[sci]) )
        ++sci;
    return arg.substr(pos, sci-pos);
}

//------------------------------------------------------------------------------
#pragma mark -

/**
 get_symbol() reads consecutive characters:
 - starting with a alpha-character,
 - followed by alpha-numerical characters or '_'
 .
 */
std::string Tokenizer::get_symbol(std::istream& is, bool eat_line)
{
    std::string res;
    int c = skip_space(is, eat_line);
    
    if ( isalpha(c) )
        res = get_stuff(is, valid_symbol);
    
    VLOG("Tokenizer:SYMBOL |" << res << "|\n");
    return res;
}

/**
 get_symbol() reads consecutive characters:
 - starting with a alpha-character,
 - followed by alpha-numerical characters or '_'
 .
 */
std::string Tokenizer::get_symbol(std::string const& str, size_t& sci)
{
    std::string res;
    int c = str[sci];
    while ( isspace(c) )
    {
        if ( c == '\n' )
            break;
        c = str[++sci];
    }

    if ( isalpha(c) )
        res = get_stuff(str, sci, valid_symbol);
    
    VLOG("Tokenizer:SYMBOL |" << res << "|\n");
    return res;
}


/**
 This will split 'arg' into 'symbol|stuff'.
 It will return 'symbol' and truncate arg to 'stuff'.
 */
std::string Tokenizer::split_symbol(std::string& arg)
{
    size_t i = 0;
    while ( isspace(arg[i]) )
        ++i;
    size_t s = i;
    
    while ( valid_symbol(arg[i]) )
        ++i;
    size_t e = i - s;
    
    std::string res = arg.substr(s, e-s);
    
    while ( isspace(arg[i]) )
        ++i;

    // remove extracted part:
    arg.erase(0, i);
    return res;
}


bool Tokenizer::has_symbol(std::string const& str, size_t& sci, std::string const& arg, bool eat_line)
{
    if ( arg[sci] )
    {
        size_t pos = sci;
        if ( get_symbol(str, sci) == arg )
            return true;
        sci = pos;
    }
    return false;
}


/**
 get_polysymbol() reads multiple symbols concatenated with ':'
 */
std::string Tokenizer::get_polysymbol(std::istream& is, bool eat_line)
{
    std::string res = get_symbol(is, eat_line);
    while ( is.peek() == ':' )
    {
        res += (char)is.get();
        res += get_stuff(is, valid_symbol);
    }
    VLOG("Tokenizer:POLYSYMBOL |" << res << "|\n");
    return res;
}


/**
 Separates a word and a number, eg. 'cell1' into 'cell' and '1'
 It also splits 'cell:1' into 'cell' and '1'
 Negative numbers are coded as 'cell~1', splited as 'cell' and '-1'
 */
bool Tokenizer::split_polysymbol(std::string& arg, long& num)
{
    long val = 0;
    std::istringstream is(arg);
    std::string res = get_symbol(is, false);
    size_t n = res.size();
    int c = is.good() ? is.peek() : 0;
    if ( c == ':' )
    {
        is.get();
        // spliting as symbol:number
        is >> val;
        if ( val <= 0 )
            return false;
    }
    else if ( c == '~' )
    {
        is.get();
        // spliting as symbol:negative_number
        is >> val;
        if ( val <= 0 )
            return false;
        val = -val;
    }
    else
    {
        // spliting as word|number, concatenated
        is.clear();
        is.seekg(0);
        n = get_stuff(is, isalpha).size();
        is >> val;
    }
    // check if splitting was successful:
    if ( is.fail() )
        return false;
    num = val;
    arg.resize(n);
    return true;
}


/**
 Return next token if it looks likes a path to a file name
 */
std::string Tokenizer::get_path(std::istream& is, bool eat_line)
{
    int c = skip_space(is, eat_line);
    
    if ( c == '*' )
    {
        is.get();
        return "*";
    }
    
    if ( !valid_path_start(c) )
        return "";
    
    std::string res = get_stuff(is, valid_path);
    
    VLOG("Tokenizer: FILENAME |" << res << "|\n");
    
    return res;
}


//------------------------------------------------------------------------------
#pragma mark - Numbers

/**
 Read multiple forms of integer:
 - integer-constant: digit-sequence [exponent-part]
 - digit-sequence:   [digit-sequence] digit
 - digit         :   one of [0123456789]
 - exponent-part: ('e' or 'E') [sign] digit-sequence
 */
std::string Tokenizer::get_integer(std::istream& is)
{
    bool accept_expon = true;
    std::string res;
    
    if ( ! is.good() )
        return "";
    
    int c = is.get();
    
    if (( c == '+' ) | ( c == '-' ))
    {
        int d = is.peek();
        if (( d == EOF )| !isdigit(d) )
        {
            is.unget();
            return "";
        }
        res.push_back((char)c);
        c = is.get();
    }
    
    while ( c != EOF )
    {
        if ( isdigit(c) )
            res.push_back((char)c);
        else if ((( c == 'e' ) | ( c == 'E' )) & accept_expon )
        {
            accept_expon = false;
            int d = is.peek();
            // accept an optional sign character
            if ( isdigit(d) )
                res.push_back((char)c);
            else if (( d == '+' ) | ( d == '-' ))
            {
                res.push_back((char)c);
                res.push_back((char)is.get());
            }
            else
                break;
        }
        else
            break;
        c = is.get();
    }
    
    if ( c != EOF )
        is.unget();
    
    VLOG("Tokenizer: INTEGER |" << res << "|\n");
    return res;
}

/**
 Read an integer, or return false if this was not possible.
 This does not skip newline characters.
 Note that the value of `var` will not change, if the input fails.
 */
bool Tokenizer::get_integer(std::istream& is, unsigned long& var)
{
    int c = skip_space(is, false);
    if ( c == '-' )
        throw InvalidParameter("a non-negative integer is expected");
    std::streampos isp = is.tellg();
    unsigned long num = 0;
    is >> std::noskipws >> num;
    if ( is.fail() )
    {
        is.clear();
        is.seekg(isp);
        return false;
    }
    if ( is.peek() == '.' )
    {
        // declare error, if a 'dot' sign follows
        is.seekg(isp);
        throw InvalidParameter("an integer is expected");
    }
    var = num;
    return true;
}

/**
 Read an integer, or return false if that was not possible.
 This does not skip newline characters.
 Note that the value of `var` will not change, if the input fails.
 */
bool Tokenizer::get_integer(std::istream& is, long& var)
{
    skip_space(is, false);
    std::streampos isp = is.tellg();
    long num = 0;
    is >> std::noskipws >> num;
    if ( is.fail() )
    {
        is.clear();
        is.seekg(isp);
        return false;
    }
    if ( is.peek() == '.' )
    {
        // declare error, if a 'dot' sign follows
        is.seekg(isp);
        throw InvalidParameter("an integer is expected");
    }
    var = num;
    return true;
}


std::vector<std::string> Tokenizer::split(std::string& str, char sep, bool get_empty_fields)
{
    std::vector<std::string> res;
    std::istringstream iss(str);
    while (!iss.eof())
    {
        std::string s;
        getline(iss, s, sep);
        if ( get_empty_fields || !s.empty() )
            res.push_back(s);
    }
    return res;
}

/**
 Split string `arg` into an integer, a space, and the remaining string.
 Any space after the integer is discarded. `arg` is truncated in the process.
 */
bool Tokenizer::split_integer(long& val, std::string& arg)
{
    long num;
    size_t i = 0;
    try {
        num = std::stol(arg, &i, 10);
    }
    catch (...)
    {
        return false;
    }
    if ( i > 0 && isspace(arg[i]) )
    {
        val = num;
        // skip any additional space-like characters:
        while ( isspace(arg[i]) )
            ++i;
        // remove consumed characters:
        arg.erase(0, i);
        return true;
    }
    return false;
}

/**
 Split string `arg` into an integer, a space, and the remaining string.
 Any space after the integer is discarded. `arg` is truncated in the process.
 */
bool Tokenizer::split_integer(unsigned long& val, std::string& arg)
{
    unsigned long num;
    size_t i = 0;
    try {
        num = std::stoul(arg, &i, 10);
    }
    catch (...)
    {
        return false;
    }
    if ( i > 0 && isspace(arg[i]) )
    {
        val = num;
        // skip any additional space-like characters:
        while ( isspace(arg[i]) )
            ++i;
        // remove consumed characters:
        arg.erase(0, i);
        return true;
    }
    return false;
}


/// like the other split_integer, but check for range
bool Tokenizer::split_integer(unsigned& val, std::string& arg)
{
    unsigned long ul = val;
    if ( split_integer(ul, arg) )
    {
        val = static_cast<unsigned>(ul);
        if ( static_cast<unsigned long>(val) != ul )
            throw InvalidParameter("out-of-range split_integer(unsigned)");
        return true;
    }
    return false;
}

/**
 Read a number specified in the standard (US) format with an optional '.':
 floating-point-constant: 
             [sign] fractional-constant [exponent-part]
             [sign] digit-sequence [exponent-part]
 fractional-constant: digit-sequence.digit-sequence
 digit-sequence:     [digit-sequence] digit
 digit         : one of [0123456789]
 exponent-part : ('e' or 'E') [sign] digit-sequence
 sign          : '+' or '–'
*/
std::string Tokenizer::get_real(std::istream& is)
{
    bool accept_point = true;
    bool accept_expon = true;
    std::string res;
    
    if ( ! is.good() )
        return "";

    int c = is.get();
    
    if (( c == '+' ) | ( c == '-' ))
    {
        int d = is.peek();
        if ( (d == EOF) | !isdigit(d) )
        {
            is.unget();
            return "";
        }
        res.push_back((char)c);
        c = is.get();
    }

    while ( c != EOF )
    {
        if ( isdigit(c) )
            res.push_back((char)c);
        else if ( (c == '.') & accept_point )
        {
            res.push_back((char)c);
            accept_point = false;
        }
        else if ((( c == 'e' ) | ( c == 'E' )) & accept_expon )
        {
            accept_expon = false;
            // only accept integer within exponent
            accept_point = false;
            int d = is.peek();
            // accept an optional sign character
            if ( isdigit(d) )
                res.push_back((char)c);
            else if (( d == '+' ) | ( d == '-' ))
            {
                res.push_back((char)c);
                res.push_back((char)is.get());
            }
            else
                break;
        }
        else
            break;
        c = is.get();
    }
    
    if ( c != EOF )
        is.unget();
    
    VLOG("Tokenizer: REAL |" << res << "|\n");
    return res;
}


std::string Tokenizer::get_hexadecimal(std::istream& is)
{
    return get_stuff(is, valid_hexadecimal);
}

//------------------------------------------------------------------------------
#pragma mark -


/**
 get_token() reads the next 'token' in stream.
 It can be:
 - a symbol starting with an alpha character
 - a block enclosed by delimiters: {}, () or ""
 - a number 
 - the new line character ('\n') if 'eat_line=true'
 - a single character, except '\n' if 'eat_line=false'
 */

std::string Tokenizer::get_token(std::istream& is)
{
    int c = skip_space(is, false);

    if ( isalpha(c) )
        return get_symbol(is);

    if ( block_delimiter(c) )
        return get_block_text(is, (char)is.get(), block_delimiter(c));
    
    c = is.get();
    int d = is.peek();

    if ( d == EOF )
        return std::string(1,(char)c);
    
    if (( c == '0' ) & ( d == 'x' ))
    {
        is.unget();
        return get_hexadecimal(is);
    }
    
    if ( isdigit(c) | ((( c == '-' ) | ( c == '+' )) & isdigit(d)) )
    {
        is.unget();
        return get_real(is);
    }
    
    // anything else is one character long:
    VLOG("Tokenizer: ASCII |" << c << "|\n");
    return std::string(1,(char)c);
}

//------------------------------------------------------------------------------

/**
 This will read a block, assuming that the opening delimiter has been read already.
 It will read characters until the given closing delimiter `c_out` is found.
 if `c_in` is not zero, the block is returned with `c_in` and `c_out` at the
 first and last positions.
 if `c_in` is null, the character `c_out` is read but not appended to the result.
 This throws an exception if `c_out` is not found.
 */
std::string Tokenizer::get_block_text(std::istream& is, char c_in, const char c_out)
{
    assert_true(c_out);
    std::string res;
    res.reserve(2032);
    char c = 0;
    
    if ( c_in )
        res.push_back(c_in);
    is.get(c);
    
    while ( is.good() )
    {
        if ( c == c_out )
        {
            if ( c_in )
                res.push_back(c_out);
            return res;
        }
        else if ( block_delimiter(c) )
            res.append( get_block_text(is, c, block_delimiter(c)) );
        else if (( c == ')' ) | ( c == '}' ))
            throw InvalidSyntax("unclosed delimiter '"+std::string(1,c_in)+"'");
        else
            res.push_back(c);
        is.get(c);
    }
    
    throw InvalidSyntax("missing '"+std::string(1,c_out)+"'");
    return "";
}


/**
 This will read a block, assuming that the opening delimiter has been read already.
*/
std::string Tokenizer::get_blocked_text(std::string const& arg, size_t& sci, char c_in, const char c_out)
{
    assert_true(c_out);
    size_t start = sci;
    char c;
    
    do {
        c = arg[sci++];
        if ( c == c_out )
            return arg.substr(start, sci-start);
        else if ( block_delimiter(c) )
            get_blocked_text(arg, sci, c, block_delimiter(c));
        else if ( c == ')' || c == '}' || c == ']' )
            throw InvalidSyntax("unclosed delimiter '"+std::string(1,c_in)+"'");
    } while ( c );
    
    throw InvalidSyntax("missing '"+std::string(1,c_out)+"'");
    return "";
}

/**
 This will skip spaces and new-lines until a character is found.
 If this character is equal to `c_in`, then the block is read and returned.
 Otherwise returns empty string "".
 
 @returns content of the block without delimiters
 */
std::string Tokenizer::get_block(std::istream& is, char c_in, bool or_die)
{
    assert_true(c_in);
    
    int c = skip_space(is, true);
    
    if ( c == '%' )
        c = skip_line(is);
    
    if ( c == c_in )
    {
        is.get();
        std::string res = get_block_text(is, 0, block_delimiter(c_in));
        VLOG("Tokenizer:BLOCK |" << res << "|\n");
        return res;
    }

    if ( or_die )
        throw InvalidSyntax("expected {}-delimited block `{...}'");

    return "";
}


std::string Tokenizer::get_block(std::istream& is)
{
    int c = skip_space(is, true);
    
    if ( block_delimiter(c) )
        return get_block_text(is, (char)c, block_delimiter(c));
    
    return "";
}


bool Tokenizer::strip_block(std::string& arg, char c_in)
{
    size_t s = arg.size();
    
    char c_out = block_delimiter(c_in);

    if ( s < 2 || arg[0] != c_in || arg[s-1] != c_out )
        return false;
    
    arg = arg.substr(1, s-2);
    return true;
}


//------------------------------------------------------------------------------
#pragma mark -


std::string Tokenizer::get_until(std::istream& is, std::string what)
{
    std::string res;
    res.reserve(16384);
    size_t d = 0;
    char c = 0;
    is.get(c);
    
    while ( is.good() )
    {
        if ( c == what[d] )
        {
            ++d;
            if ( what[d] == '\0' )
                break;
        }
        else
        {
            if ( d == 0 )
            {
                res.push_back(c);
            }
            else
            {
                res.push_back(what[0]);
                if ( d > 1 ) {
                    is.seekg(-(int)d, std::ios_base::cur);
                    d = 0;
                } else {
                    if ( c == what[0] )
                        d = 1;
                    else {
                        res.push_back(c);
                        d = 0;
                    }
                }
            }
        }
        is.get(c);
    }
    //std::clog << "get_until|" << res << '\n';
    return res;
}


std::string Tokenizer::trim(std::string const& str, const std::string& ws)
{
    std::string::size_type s = str.find_first_not_of(ws);
    if ( s == std::string::npos )
        return std::string();
    std::string::size_type e = str.find_last_not_of(ws);
    return str.substr(s, 1+e-s);
}

