// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include "glossary.h"
#include "tokenizer.h"
#include "stream_func.h"
#include "print_color.h"
#include "filepath.h"
#include <fstream>
#include <cctype>
#include <iomanip>

/// This is used to align text in error messages
const char PREF[] = "          : ";


// Use the second definition to get some reports:
#if 0
#  define VLOG0(ARG) std::clog << ARG;
#  define VLOG1(ARG) std::clog << ARG;
#  define VLOG2(ARG) std::clog << ARG;
#else
#  define VLOG0(ARG) ((void) 0)
#  define VLOG1(ARG) ((void) 0)
#  define VLOG2(ARG) ((void) 0)
#endif

//------------------------------------------------------------------------------

Glossary::Glossary()
{
}

Glossary::Glossary(std::istream& in)
{
    read(in);
}

Glossary::Glossary(const std::string& arg)
{
    std::istringstream iss(arg);
    read(iss);
}

Glossary::Glossary(const std::string& key, const std::string& val)
{
    define(key, val);
}


/** If successful, `var` is set to the content of the block without delimiters */
int Glossary::set_block(std::string & var, char c_in, key_type const& key, size_t inx) const
{
    rec_type const* rec = values(key);
    
    if ( rec && inx < rec->size() )
    {
        val_type const& val = rec->at(inx);
        
        if ( val.defined_ )
        {
            std::string const& arg = val.value_;
            size_t s = arg.size();
            if ( s < 2 || arg[0] != c_in )
                return 0;
            if ( arg[s-1] != Tokenizer::block_delimiter(c_in) )
                return 0;
            var = arg.substr(1, s-2);
            ++val.read_;
            return 1;
        }
    }
    return 0;
}

//------------------------------------------------------------------------------
#pragma mark - Formating

std::string format_value(std::string const& str)
{
    if ( std::string::npos != str.find_first_of(" ;") )
        return '(' + str + ')';
    else
        return str;
}


std::string format_count(unsigned c)
{
    if ( c == 0 )
        return " (unread)";
    if ( c == 1 )
        return " (read)";
    return " (read " + std::to_string(c) + "x)";
}
    

std::string format(Glossary::pair_type const& pair)
{
    std::string res = pair.first;
    if ( pair.second.size() > 0 )
    {
        res += " = " + format_value(pair.second[0].value_);
        for ( unsigned i = 1; i < pair.second.size(); ++i )
            res += ", " + format_value(pair.second[i].value_);
        res += ";";
    }
    return res;
}

/**
 Write the usage-counter for each value.
 The width of each record will match what is printed by Glossary::print()
 */
std::string format_counts(Glossary::pair_type const& pair)
{
    std::string res = pair.first;
    if ( pair.second.size() > 0 )
    {
        res += " = " + pair.second[0].value_ + format_count(pair.second[0].read_);
        for ( unsigned i = 1; i < pair.second.size(); ++i )
            res += ", " + pair.second[i].value_ + format_count(pair.second[i].read_);
        res += ";";
    }
    return res;
}


//------------------------------------------------------------------------------
#pragma mark - Access

bool Glossary::has_key(key_type const& k) const
{
    return ( mTerms.end() != mTerms.find(k) );
}


/**
 @return 'false' if the key is not defined and `true` if the key is defined without a value,
 and finally the value if the key was defined with a value that can be read as a boolean.
 */
bool Glossary::use_key(key_type const& k)
{
    map_type::iterator w = mTerms.find(k);
    
    if ( w != mTerms.end() )
    {
        rec_type& rec = w->second;
        bool var = true;

        // if a value is defined, we return this value:
        if ( rec.size() > 0 )
        {
            val_type const& val = rec.at(0);
            if ( val.defined_ )
            {
                std::istringstream iss(val.value_);
                iss >> var;
                if ( iss.fail() )
                {
                    print_magenta(stderr, "warning: unexpected value `"+val.value_+"' specified for `"+w->first+"'\n");
                    var = false;
                }
            }
        }
        mTerms.erase(w);
        return var;
    }
    return false;
}


void Glossary::clear(key_type const& key)
{
    map_type::iterator w = mTerms.find(key);
    
    if ( w != mTerms.end() )
        mTerms.erase(w);
}


void Glossary::clear_except(key_type const& key)
{
    map_type::iterator w = mTerms.find(key);
    
    if ( w == mTerms.end() )
        clear();
    else
    {
        rec_type rec = w->second;
        clear();
        mTerms[key] = rec;
    }
}


void Glossary::clear_reads()
{
    for ( auto& i : mTerms )
    {
        for ( val_type const& v : i.second )
            v.read_ = 0;
    }
}


unsigned Glossary::num_reads(std::string const& key) const
{
    unsigned res = 0;
    map_type::const_iterator w = mTerms.find(key);
    if ( w != mTerms.end() )
    {
        for ( val_type const& v : w->second )
            res += v.read_;
    }
    return res;
}


/** transfer read counts from 'opt' to 'this' */
void Glossary::add_reads(Glossary const& opt)
{
    for ( auto& i : mTerms )
    {
        map_type::const_iterator w = opt.mTerms.find(i.first);
        if ( w != opt.mTerms.end() )
        {
            size_t sup = std::min(i.second.size(), w->second.size());
            //std::clog << " Glossary::add_reads `" << w->first << "' " << sup << "\n";
            for ( size_t n = 0; n < sup; ++n )
                i.second[n].read_ = std::max(i.second[n].read_, w->second[n].read_);
        }
    }
}


Glossary Glossary::get_term(key_type const& key) const
{
    Glossary res;
    map_type::const_iterator w = mTerms.find(key);
    
    if ( w != mTerms.end() )
        res.mTerms[key] = w->second;
    
    return res;
}


Glossary Glossary::unused_terms() const
{
    Glossary res;
    for ( auto const& i : mTerms )
    {
        bool used = false;
        for ( val_type const& v : i.second )
            if ( v.read_ == 0 )
                used = true;
        
        if ( !used )
            res.mTerms[i.first] = i.second;
    }
    return res;
}

//------------------------------------------------------------------------------
#pragma mark - Values access

size_t Glossary::num_values(key_type const& k) const
{
    map_type::const_iterator w = mTerms.find(k);
    if ( w != mTerms.end() )
        return w->second.size();
    else
        return 0;
}


bool Glossary::has_value(key_type const& key, size_t inx) const
{
    map_type::const_iterator w = mTerms.find(key);
    if ( w != mTerms.end() )
        return (inx < w->second.size()) && w->second[inx].value_.size();
    return false;
}


Glossary::rec_type * Glossary::values(key_type const& key)
{
    map_type::iterator w = mTerms.find(key);
    return ( w == mTerms.end() ) ? nullptr : &( w->second );
}


Glossary::rec_type const* Glossary::values(key_type const& key) const
{
    map_type::const_iterator w = mTerms.find(key);
    return ( w == mTerms.end() ) ? nullptr : &( w->second );
}


std::string Glossary::value(key_type const& key, size_t inx) const
{
    map_type::const_iterator w = mTerms.find(key);
    if ( w != mTerms.end() )
    {
        if ( inx < w->second.size() )
        {
            w->second[inx].read_++;
            return w->second[inx].value_;
        }
    }
    return "";
}


/**
 This is equivalement to value(key, inx) == val, except
 that the counter is incremented only if there is an exact match
 */
bool Glossary::value_is(key_type const& key, size_t inx, std::string const& val) const
{
    map_type::const_iterator w = mTerms.find(key);
    if ( w != mTerms.end() )
    {
        if ( inx < w->second.size() )
        {
            if ( w->second[inx].value_ == val )
            {
                w->second[inx].read_++;
                return true;
            }
        }
    }
    return false;
}

//------------------------------------------------------------------------------
#pragma mark - Value definitions

/**
 This reads a KEY followed by the assignement operator
 @returns
 - 0 if no valid key is found
 - 1 if the key is immediately followed by the '=' sign
 - 2 otherwise
*/
int Glossary::read_key(Glossary::pair_type& res, std::istream& is)
{
    std::string k = Tokenizer::get_symbol(is, false);

    if ( k.empty() )
        return 1;

    res.first = k;
    
    int op = Tokenizer::get_character(is);

    if ( op == '=' )
    {
        VLOG2("Glossary:KEY     |" << res.first << "| = \n");
        return 0;
    }
    
    VLOG2("Glossary:KEY     |" << res.first << "|\n");
    return 2;
}


/**
 push value at the back of `res.second`
 */
void Glossary::add_rhs_value(Glossary::key_type const& key, Glossary::rec_type& rec, std::string& str, bool def)
{
    //remove any space at the end of the string:
    std::string val = Tokenizer::trim(str);
    
    VLOG2("Glossary:SET" << std::setw(20) << key << "[" << rec.second.size() << "] = |" << val << "|\n");

    rec.emplace_back(val, def);
}


/**
 return true if `c` can constitute a value.
 space are allowed since vectors are read as sets of space-separated values
*/
bool valid_value(const int c)
{
   return isalnum(c) || c==' ' || c=='/' || c == '#' || c==':' || c=='\t'
    || c=='_' || c=='.' || c=='+' || c=='-' || c=='~';
}


/**
 read one right-hand side value of an assignment
 return 1 if parsing should continue with the same key
 */
int Glossary::read_rhs(Glossary::key_type const& key, Glossary::rec_type& rec, std::istream& is)
{
    // skip spaces, but do not eat lines
    char c = Tokenizer::get_character(is);
    char d = Tokenizer::block_delimiter(c);
    bool delimited = 0;
    
    std::string k;
    if ( d )
    {
        delimited = true;
        if ( c == '(' )
            k = Tokenizer::get_block_text(is, 0, d); // parentheses are stripped
        else
            k = Tokenizer::get_block_text(is, c, d);
        c = Tokenizer::get_character(is);
    }
    else
    {
        while ( valid_value(c) || c == '(' )
        {
            if ( c == '(' )
                k.append(Tokenizer::get_block_text(is, c, ')'));
            else
                k.push_back(char(c));
            c = is.get();
        }
    }
    //std::clog << (char)c << "|";
        
    if ( c == EOF || c == '\n' || c == '\r' )
    {
        if ( k.size() || delimited )
            add_rhs_value(key, rec, k, true);
        return 0;
    }
    
    if ( c == ';' )
    {
        add_rhs_value(key, rec, k, k.size() || delimited);
        return 0;
    }

    if ( c == ',' )
    {
        add_rhs_value(key, rec, k, k.size() || delimited);
        return 1;
    }

    if ( c == '%' )
    {
        if ( k.size() || delimited )
            add_rhs_value(key, rec, k, true);
        Tokenizer::skip_line(is);
        return 0;
    }
    
    if ( c == '\\' )
    {
        if ( k.size() || delimited )
            add_rhs_value(key, rec, k, true);
        // go to next line:
        Tokenizer::skip_space(is, true);
        return 1;
    }
    
    is.unget();
    throw InvalidSyntax("syntax error: unexpected `"+std::string(1,(char)c)+"'");
}


/**
 The value of `mode' will determine how to handle duplicate values:
 - `0`, previous values are erased without warning,
 - `1`, a prexisting symbol in not altered, but no exception is thrown
 - `2`, an exception is thrown for any duplicate symbol
 */
void Glossary::add_entry(Glossary::pair_type const& pair, int no_overwrite)
{
    VLOG0("Glossary:ENTRY" << pair.second.size() << "   " << pair << '\n');
    
    map_type::iterator w = mTerms.find(pair.first);
    
    if ( w == mTerms.end() )
    {
        mTerms.insert(pair);
    }
    else
    {
        // the key already exists:
        rec_type & rec = w->second;
        // check every value of the argument:
        for ( unsigned i = 0; i < pair.second.size(); ++i )
        {
            if ( rec.size() <= i )
                rec.push_back(pair.second[i]);
            else
            {
                // the record already exists:
                if ( !rec[i].defined_  ||  !no_overwrite )
                    rec[i] = pair.second[i];
                else if ( pair.second[i].value_ != rec[i].value_  &&  no_overwrite > 1 )
                    throw InvalidSyntax("conflicting definition: "+format(*w)+" and "+format(pair));
            }
        }
    }
}


/// define one value for the key at specified index: `key[inx]=rhs`.
void Glossary::define(key_type const& key, const std::string& rhs, size_t inx)
{
    std::string val = Tokenizer::trim(rhs);
    map_type::iterator w = mTerms.find(key);

    if ( w == mTerms.end() )
    {
        if ( inx > 0 )
            throw InvalidSyntax("index out of range in Glossary::define");
        VLOG1("Glossary:DEFINE    "+key+" = |"+rhs+"|\n");
        mTerms[key].emplace_back(val, true);
    }
    else
    {
        VLOG1("Glossary:DEFINE    "+key+"[" << inx << "] = |" << val << "|\n");
        rec_type & rec = w->second;
        if ( rec.size() > inx )
            rec[inx] = val_type(val, true);
        else if ( rec.size() == inx )
            rec.emplace_back(val, true);
        else
            throw InvalidSyntax("index out of range in Glossary::define");
    }
}


/// define possibly muliple values for the key: `key = rhs`.
void Glossary::define_rhs(key_type const& key, const std::string& rhs)
{
    std::istringstream iss(rhs);
    map_type::iterator w = mTerms.find(key);
    VLOG1("Glossary:DEFINE_RHS "+key+" = |"+rhs+"|\n");

    if ( w == mTerms.end() )
    {
        pair_type pair;
        pair.first = Tokenizer::trim(key);
        while ( read_rhs(pair.first, pair.second, iss) );
        mTerms.insert(pair);
    }
    else
    {
        // replace all values:
        w->second.clear();
        while ( read_rhs(key, w->second, iss) );
    }
}


/**
 This should be equivalent to read('k = rhs')
 */
void Glossary::add_value(key_type const& key, const std::string& arg)
{
    std::string val = Tokenizer::trim(arg);
    map_type::iterator w = mTerms.find(key);
    // define a new key, or add value to existing key:
    if ( w == mTerms.end() )
        mTerms[key].emplace_back(val, true);
    else
        w->second.emplace_back(val, true);
}


//------------------------------------------------------------------------------
#pragma mark - Output

void Glossary::print_counts(std::ostream& os, std::string const& prefix) const
{
    for ( auto const& i : mTerms )
    {
        os << prefix << format_counts(i);
        std::endl(os);
    }
}


void Glossary::print(std::ostream& os, std::string const& prefix) const
{
    for ( auto const& i : mTerms )
    {
        os << prefix << format(i);
        std::endl(os);
    }
}


std::string Glossary::to_string() const
{
    std::stringstream ss;
    ss << "{ ";
    for ( auto const& i : mTerms )
    {
        ss << i.first;
        if ( i.second.size() > 0 )
        {
            ss << '=' << i.second[0].value_;
            for ( unsigned x = 1; x < i.second.size(); ++x )
                ss << ", " + i.second[x].value_;
        }
        ss << "; ";
    }
    ss << '}';
    return ss.str();
}


std::ostream& operator << (std::ostream& os, Glossary::pair_type const& arg)
{
    os << "  " << format(arg);
    return os;
}


std::ostream& operator << (std::ostream& os, Glossary const& arg)
{
    arg.print(os, "  ");
    return os;
}

//------------------------------------------------------------------------------
#pragma mark - Input

/**
 @copydetails Glossary::add_entry
 */
void Glossary::read_entry(std::istream& is, int no_overwrite)
{
    int c = Tokenizer::skip_space(is, true);

    if ( c == EOF )
        return;
        
    // skip comments:
    if ( c == '%' )
    {
        Tokenizer::skip_line(is);
        return;
    }
    
    std::streampos isp = is.tellg();
    pair_type pair;
    try
    {
        int code = read_key(pair, is);
        if ( code == 1 )
            throw InvalidParameter("unexpected token");
        if ( code )
            throw InvalidParameter("syntax error");
        
        while ( read_rhs(pair.first, pair.second, is) );
        
        if ( pair.second.empty() )
            throw InvalidSyntax("expected value after `"+pair.first+"=`");
    }
    catch( Exception& e )
    {
        e << StreamFunc::marked_line(is, isp, PREF);
        throw;
    }
    
    add_entry(pair, no_overwrite);
}


/**
 @copydetails Glossary::add_entry
 */
void Glossary::read(std::istream& is, int no_overwrite)
{
    std::istream::sentry s(is);
    if ( s )
        while ( is.good() )
            read_entry(is, no_overwrite);
}


void Glossary::read(std::string const& str, int no_overwrite)
{
    VLOG2("Glossary:READ STR |" << str << "|\n");
    std::istringstream iss(str);
    read(iss, no_overwrite);
}


void Glossary::read_file(const char path[], int no_overwrite)
{
    std::ifstream is(path, std::ifstream::in);
    if ( is.good() )
        read(is, no_overwrite);
    else
        throw InvalidIO("could not open Glossary file");
}


/**
 This is useful to parse the command-line strings given to main().
 
 The following syntax will be accepted:
 FILE.EXT
 and recorded as:
 EXT = FILE.EXT
 
 Strings corresponding to existing directories
 */
void Glossary::read_string(const char arg[], int no_overwrite)
{
    pair_type pair;
    VLOG0("Glossary:ARG      |" << arg << "|\n");
    if ( strchr(arg, '=') )
    {
        /*
         Here is a key specified with one of more value:
         */
        std::istringstream iss(arg);
        if ( 0 == read_key(pair, iss) )
        {
            while ( read_rhs(pair.first, pair.second, iss) );
            if ( pair.second.empty() )
                throw InvalidSyntax("expected value after `"+pair.first+"=`");
            add_entry(pair, no_overwrite);
        }
    }
    else
    {
        // Handle a key specified without any value
        if ( FilePath::is_dir(arg) )
            add_value("directory", arg);
        else
        {
            /*
             Find last occurence of '.' to identify a potential file name,
             while anything else is treated as a orphan string */
            char const* c = strrchr(arg, '.');
            if ( c )
                add_value(c, arg);
            else
            {
                pair.first = arg;
                add_entry(pair, no_overwrite);
            }
        }
    }
}


/**
 This is useful to parse the command-line strings given to main().
 
 The following syntax will be accepted:
 FILE.EXT
 and recorded as:
 EXT = FILE.EXT
 
 Strings corresponding to existing directories
 */
int Glossary::read_strings(int argc, char* argv[], int no_overwrite)
{
    int res = 0;
    for ( int i = 0; i < argc; ++i )
    {
        try
        {
            read_string(argv[i], no_overwrite);
        }
        catch( Exception & e )
        {
            print_magenta(stderr, e.brief());
            fputs(e.info().c_str(), stderr);
            res = 1;
        }
    }
    return res;
}


std::istream& operator >> (std::istream& is, Glossary& glos)
{
    glos.read(is);
    return is;
}

//------------------------------------------------------------------------------
#pragma mark -

/**
 builds the warning message in `msg'.
 
 @returns:
 - 4 if the parameter was not read
 - 2 if one of the value was not read
 - 1 if some of the value were used multiple times
 - 0 otherwise
 .
 */
int Glossary::warning(Glossary::pair_type const& pair, std::string& msg, size_t threshold)
{
    int code = 4;
    const rec_type& rec = pair.second;
        
    for ( size_t i = 0; i < rec.size(); ++i )
    {
        val_type const& val = rec[i];
        if ( val.read_ > 0 )
            code &= ~4; // one value used
        if ( !val.read_ && val.defined_ )
            code |= 2;  // one value not used
        else if ( val.read_ > threshold )
            code |= 1;  // one value overused
    }
    
    if ( code & 4 )
        msg = "Warning, the parameter `" + format(pair) + "' was ignored";
    else if ( code & 2 )
        msg = "Warning, a value was ignored: " + format_counts(pair);
    else if ( code & 1 )
        msg = "Warning, some value might have been overused: " + format_counts(pair);
    
    return code;
}


/**
 @returns the type of warning associated with the entire set of terms
 If the return value is not zero, something was added to 'msg', and
 this might be printed by the calling function.
 */
int Glossary::has_warning(std::string& msg, size_t threshold) const
{
    int res = 0;
    std::string war;
    for ( auto const& i : mTerms )
    {
        if ( i.first.empty() )
            continue;
        int val = warning(i, war, threshold);
        if ( val )
        {
            if ( res ) msg.push_back('\n');
            msg.append(war);
            res |= val;
        }
    }
    return res;
}


void Glossary::print_warnings(FILE * file, size_t threshold, std::string const& msg) const
{
    std::string war;
    if ( has_warning(war, threshold) )
    {
        print_yellow(file, war);
        fputs(msg.c_str(), file);
    }
}

//------------------------------------------------------------------------------

/**
 This copies the string, removing spaces
*/
template <>
void Glossary::set_value(std::string& var, key_type const& key, std::string const& val)
{
    //var = Tokenizer::trim(val);
    // strip double quotes, ASCII(") == 34:
    if ( val.size() > 3 && val.front() == 34 && val.back() == 34 )
        var = val.substr(1, val.size()-2);
    else
        var = val;

    VLOG2("Glossary:SET STR   " << key << " = |" << var << "|\n");
}


/**
 This reads a floating point value,
 also accepting 'inf', '+inf' and '-inf' for INFINITY values
 */
template <>
void Glossary::set_value(float& var, key_type const& key, std::string const& val)
{
/*
    // Infinite values are normally handled by std::strtof()
    if ( val == "inf" || val == "+inf" )
    {
        var = INFINITY;
        return;
    }
    
    if ( val == "-inf" )
    {
        var = -INFINITY;
        return;
    }
*/
    char const* ptr = val.c_str();
    char * end = nullptr;
    float num = strtof(ptr, &end);
    if ( end == ptr || not_space(end) )
        throw InvalidSyntax("could not set scalar value `"+key+"' from `"+val+"'");
    var = num;
}

/**
 This reads a floating point value,
 also accepting 'inf', '+inf' and '-inf' for INFINITY values
 */
template <>
void Glossary::set_value(double& var, key_type const& key, std::string const& val)
{
/*
    // Infinite values are normally handled by std::strtod()
    if ( val == "inf" || val == "+inf" )
    {
        var = INFINITY;
        return;
    }
    
    if ( val == "-inf" )
    {
        var = -INFINITY;
        return;
    }
*/
    char const* ptr = val.c_str();
    char * end = nullptr;
    double num = strtod(ptr, &end);
    if ( end == ptr || not_space(end) )
        throw InvalidSyntax("could not set scalar value `"+key+"' from `"+val+"'");
    var = num;
}
