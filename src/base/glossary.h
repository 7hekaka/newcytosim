// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#ifndef GLOSSARY_H
#define GLOSSARY_H

#include "exceptions.h"
#include "assert_macro.h"
#include <iostream>
#include <sstream>
#include <errno.h>
#include <string>
#include <vector>
#include <map>


/// This is used to align text in the error messages
extern const char PREF[];


/// Glossary holds a list of (key, values) where both key and values are strings
/** 
 This class is used for reading configuration files:
 - Reads a std::istream to builds a std::map of <key, record>
 - Simple syntax based on = ( ) { } “ , ; % focused on flexible value setting.
 - each `record` is a list of values
 - Provides values upon requests with function set(key, index) .
 - A counter records the usage of the values.
 .
 
 Notes:
-# There can be an arbitrary number of Keys, and an abitrary number of values for each key.
-# Values are kept as strings, and are converted at request by templated functions:

       template <typename T> int set(T & ptr, std::string key)

   - `key` is the name under which the value appeared,
   - The template argument `T` defines the type of the parameter,
     and the value-string is interpreted accordingly,
   - The interpreted value is stored in `ptr`,
   - Returns 1 if `ptr` was set successfully, 0 otherwise
   .
-# The method warning() can report values that have not been used, 
 or that have been used more than once.
.
 
 Class reviewed by Andre Clapson on 10.03.2011.
 
 @todo parse and instantiate values like 'random.uniform()' or 'PI*30' or '0.1/60'?
*/

class Glossary
{
public:
    
    /// a key is defined by its name, a string
    typedef std::string key_type;
    
    /// an string ASCII-encoded value with metadata
    struct val_type 
    {
        /// the value specified as a string
        std::string value_;

        /// true if this value has been explicitly defined by the user
        bool defined_;
        
        /// the number of times this value has been read
        mutable unsigned read_;

        /// constructor
        val_type() { defined_=false; read_=0; }
        
        /// constructor with initialization
        val_type(std::string const& s, bool d) { value_=s; defined_=d; read_=0; }
    };
    
    /// a record is a set of values associated with a key
    typedef std::vector<val_type>         rec_type;
    
    /// type for the list of (key, record)
    typedef std::map<key_type, rec_type>  map_type;
    
    /// type of a pair (key, record)
    typedef std::pair<key_type, rec_type> pair_type;

    /// type for a dictionary of terms given to set(T&, ...)
    template < typename T >
    using dict_type = std::initializer_list<std::pair<Glossary::key_type, T> >;

private:
    
    /// ordered list of key-values pairs
    map_type mTerms;
    
    //-------------------------------------------------------------------------------
    #pragma mark -

    /// write a value, adding enclosing parenthesis if it contains space characters
    static std::string format_value(const std::string&);
    
    /// write the number of time parameter was used
    static std::string format_count(size_t);
    
    /// report unused values and values used more than `threshold` times
    static int warning(pair_type const&, std::string& msg, size_t threshold);

    /// read key and assignement operator
    static int read_key(pair_type&, std::istream&);
    
    /// read one right-hand-side entry of an assignement
    static int read_rhs(key_type const&, rec_type&, std::istream&);
    
    /// add right-hand-side entry to pair.second
    static void add_rhs_value(key_type const&, rec_type&, std::string&, bool);

    /// register a new pair into the dictionnary
    void add_entry(pair_type const&, int no_overwrite);
    
    /// disabled
    void define(key_type const& key, void*);
    
    //-------------------------------------------------------------------------------
    
    /// returns first non-space character in null-terminated C-string
    static char const* not_space(const char s[])
    {
        while(*s)
        {
            if ( isspace(*s) )
                ++s;
            else
                return s;
        }
        return nullptr;
    }
    
    /// true if `str` is composed of alpha characters and '_'
    static int is_alpha(std::string const& str)
    {
        if ( str.empty() )
            return 0;
        for ( char c : str )
        {
            if ( ! isalpha(c) &&  c != '_' )
                return 0;
        }
        return 1;
    }
    
    /// check if string could be a number
    static int is_number(std::string const& str)
    {
        char const* ptr = str.c_str();
        char * end = nullptr;
        
        errno = 0;
        long i = strtol(ptr, &end, 10);
        
        if ( !errno && end > ptr && !not_space(end) )
            return 2 + ( i >= 0 );
        
        errno = 0;
        double d = strtod(ptr, &end);
        if ( !errno && end > ptr && !not_space(end) )
            return 4 + ( d >= 0 );
        
        return 0;
    }
    
    /// set `var` from string `val`
    template <typename T>
    static void set_value(T & var, key_type const& key, std::string const& val)
    {
        if ( val.empty() )
            throw InvalidSyntax("could not set `"+key+"' from empty string");
        
        std::istringstream iss(val);
        
        iss >> var;
        
        if ( iss.fail() )
            throw InvalidSyntax("could not set `"+key+"' from `"+val+"'");
        
        // check if all characters were used:
        if ( ! iss.eof() )
        {
            char const* chr = val.c_str() + iss.tellg();
            if ( not_space(chr) )
            {
                std::cerr << "Warning: ignored trailing `" << chr;
                std::cerr << "' in `" << key << " = " << val << "'\n";
            }
        }
    }
    
    /// set enum of type T using a dictionary of correspondances
    template <typename T>
    static void set_value(T & var, key_type const& key, std::string const& val,
                          dict_type<T> const& dict, bool accept_all)
    {
        for ( auto const& kv : dict )
        {
            // accept both the named-value and the representation of the value
            if ( val == kv.first || val == std::to_string(kv.second) )
            {
                var = kv.second;
                //std::clog << "KeyList::set  " << key << " -> " << val << '\n';
                return;
            }
        }
        
        if ( accept_all && std::is_enum<T>::value )
        {
            // accept a numeric value
            size_t idx;
            int i = std::stoi(val, &idx, 10);
            T t = static_cast<T>(i);
            if ( idx == val.size() && (int)t == i )
            {
                std::clog << "Warning: accepting unknown value " << t << " for `" << key << "'\n";
                var = t;
                return;
            }
        }
        
        InvalidParameter e("could not set `"+key+"' from `"+val+"'");
        e << "\nKnown values are:\n";
        for ( auto const& kv : dict )
            e << PREF << kv.first << " = " << kv.second << '\n';
        throw e;
    }
    
    //-------------------------------------------------------------------------------
    #pragma mark -

public:
    
    /// initialize
    explicit Glossary();

    /// this constructor calls read(in)
    explicit Glossary(std::istream& in);

    /// this constructor calls read(string)
    explicit Glossary(const std::string&);

    /// this constructor calls define(key, val)
    explicit Glossary(const std::string& key, const std::string& val);

    //-------------------------------------------------------------------------------

    /// true if no key were set
    bool empty() const { return mTerms.empty(); }

    /// total number of keys
    size_t num_keys() const { return mTerms.size(); }
    
    /// return `true` if key is present, even if no value is associated with it
    bool has_key(key_type const&) const;
    
    /// return `true` if key is present, and delete key
    bool use_key(key_type const&);
    
    /// remove given key
    void clear(key_type const&);
    
    /// clear all entries
    void clear() { mTerms.clear(); }
    
    /// remove all keys except the given one
    void clear_except(key_type const&);

    /// clear usage counts for all entries
    void clear_reads();
    
    /// return number of times some value was read
    unsigned num_reads(key_type const&) const;

    /// add usage counts from another Glossary
    void add_reads(Glossary const&);

    /// ordered list of key-values pairs; needed for python access
    map_type terms() const { return mTerms; };

    /// create a new Glossary with only the given key
    Glossary get_term(key_type const&) const;
    
    /// create a new Glossary with terms that were not used
    Glossary unused_terms() const;

    /// return number of values associated with a key
    size_t num_values(key_type const&) const;
    
    /// return true if key is present and a value was set for given index
    bool has_value(key_type const&, size_t inx = 0) const;
    
    /// gives a pointer to the values corresponding to a key, or null if the key is not present
    rec_type * values(key_type const&);

    /// gives a const pointer to the values corresponding to a key, or null if the key is not present
    rec_type const* values(key_type const&) const;
    
    /// return copy of value corresponding to `key[inx]`, or empty string if this value is not present
    std::string value(key_type const&, size_t inx = 0) const;
    
    /// returns true if `key[inx]==val`, or false otherwise. Counter is incremented in case of match
    bool value_is(key_type const& key, size_t inx, std::string const& val) const;
    
    /// set message if values were unused or used multiple times; return warning code
    int has_warning(std::string&, size_t threshold = 1) const;
    
    /// print message about unused values and values used multiple times
    void print_warnings(FILE*, size_t threshold, std::string const&) const;

    //-------------------------------------------------------------------------------
    #pragma mark -
    
    /// define one value for the key at specified index: `key[inx] = val`.
    void define(key_type const& key, std::string const& val, size_t inx);
    
    /// define one value from class T, for the key: `key[inx] = to_string(val)`.
    template <typename T>
    void define(key_type const& key, T const& val, size_t inx = 0)
    {
        std::ostringstream oss;
        oss << val;
        define(key, oss.str(), inx);
    }
    
    /// define possibly multiple values for the key at specified index: `key[inx] = val`.
    void define_rhs(key_type const& key, std::string const& val);

    /// add value to key
    void add_value(key_type const&, std::string const& val);

    /// update the glossary to include one assignment read from stream
    void read_entry(std::istream&, int no_overwrite = 2);

    /// update the glossary to include assignments stored in a stream
    void read(std::istream&, int no_overwrite = 2);

    /// update the glossary to include assignments stored in a string
    void read(const std::string&, int no_overwrite = 2);
    
    /// read file specified in path
    void read_file(const char path[], int no_overwrite = 2);
    
    /// read a file specified by name
    void read_file(std::string const& str, int no = 2) { read_file(str.c_str(), no); }
    
    /// read a C-style argument
    void read_string(const char arg[], int no_overwrite = 2);

    /// read C-style command-line arguments, return 0 if success
    int read_strings(int argc, char* argv[], int no_overwrite = 2);

    /// write all [key, values]
    void print(std::ostream&, std::string const& prefix = "") const;
    
    /// write all [key, values] with number of access times
    void print_counts(std::ostream&, std::string const& prefix = "") const;
    
    /// write all [key, values] on a single line
    std::string to_string() const;

    //-------------------------------------------------------------------------------
    #pragma mark -

    /// try to set `var` from `key[inx]`, without recording that the parameter was read.
    template <typename T>
    int peek(T & var, key_type const& key, size_t inx = 0) const
    {
        rec_type const* rec = values(key);
        
        if ( rec && inx < rec->size() )
        {
            val_type const& val = rec->at(inx);
            
            if ( val.defined_ )
            {
                set_value(var, key, val.value_);
                return 1;
            }
        }
        return 0;
    }

    /// try to set `var` from `key[inx]`. @return 1 if `var` was set, 0 otherwise
    /** An internal counter is incremented to record that the value was read */
    template <typename T>
    int set(T & var, key_type const& key, size_t inx = 0) const
    {
        rec_type const* rec = values(key);
        
        if ( rec && inx < rec->size() )
        {
            val_type const& val = rec->at(inx);
            
            if ( val.defined_ )
            {
                set_value(var, key, val.value_);
                ++val.read_;
                return 1;
            }
        }
        return 0;
    }
    
    /// try to set `var` from `key[inx]` or `alt[jax]`. @return 1 if `var` was set, 0 otherwise
    /** An internal counter is incremented to record that the value was read */
    template <typename T>
    int set(T & var, key_type const& key, key_type const& alt) const
    {
        int res = set(var, key, 0);
        if ( !res )
            res = set(var, alt, 0);
        return res;
    }
    
    /// try to set `var` from `key1[0]`, `key2[0]` or `key3[0]`. @return 1 if `var` was set, 0 otherwise
    /** An internal counter is incremented to record that the value was read */
    template <typename T>
    int set(T & var, key_type const& key1, key_type const& key2, key_type const& key3) const
    {
        int res = set(var, key1, 0);
        if ( !res )
            res = set(var, key2, 0);
        if ( !res && !key3.empty() )
            res = set(var, key3, 0);
        return res;
    }

    /// try to set `var` from `key[inx]` or `alt[jax]`. @return 1 if `var` was set, 0 otherwise
    /** An internal counter is incremented to record that the value was read */
    template <typename T>
    int set(T & var, key_type const& key, size_t inx, key_type const& alt, size_t jax) const
    {
        int res = set(var, key, inx);
        if ( !res )
            res = set(var, alt, jax);
        return res;
    }
    
    
    /// try to set `var` from `key[inx]`, if this is a block starting by 'c_in'. @return 1 if `var` was set
    /** If successful, `var` is set to the content of the block without delimiters */
    int set_block(std::string & var, char c_in, key_type const& key, size_t inx = 0) const;

    
    /// try to set `cnt` values in `ptr[]`, starting at `key[0]`. @return number of values set
    /** The internal counters are incremented to record that the values were read */
    template <typename T>
    int set(T * ptr, size_t cnt, key_type const& key) const
    {
        rec_type const* rec = values(key);
        
        if ( !rec )
            return 0;
        
        int set = 0;
        for ( size_t inx = 0; inx < rec->size() && inx < cnt; ++inx )
        {
            val_type const& val = rec->at(inx);
            
            if ( val.defined_ )
            {
                set_value(ptr[inx], key, val.value_);
                ++rec->at(inx).read_;
                ++set;
            }
        }
        return set;
    }
   

    /// try to set `var` from `key[inx]`, using `dict`. @return 1 if `var` was set, 0 otherwise
    /** An internal counter is incremented to record that the value was read */
    template <typename T>
    int set(T & var, key_type const& key, unsigned inx, dict_type<T> const& dict, bool accept_all = false) const
    {
        rec_type const* rec = values(key);
        
        if ( rec && inx < rec->size() )
        {
            val_type const& val = rec->at(inx);
            
            if ( val.defined_ )
            {
                set_value(var, key, val.value_, dict, accept_all);
                ++val.read_;
                return 1;
            }
        }
        return 0;
    }
    
    /// try to set `var` from `key[inx]` or `alt[jax]`, using the dictionary `dict`
    /** An internal counter is incremented to record that the value was read */
    template <typename T>
    int set(T & var, key_type const& key, unsigned inx, key_type const& alt, unsigned jax,
            dict_type<T> const& dict, bool accept_all = false) const
    {
        int res = set(var, key, inx, dict, accept_all);
        if ( !res )
            res = set(var, alt, jax, dict, accept_all);
        return res;
    }

    /// try to set `var` from `key[0]`, using the dictionary `dict`
    /** An internal counter is incremented to record that the value was read */
    template <typename T>
    int set(T & var, key_type const& key, dict_type<T> const& dict, bool accept_all = false) const
    {
        return set(var, key, 0, dict, accept_all);
    }
    
    //-------------------------------------------------------------------------------
    #pragma mark -

    /// try to set `var` from `key[inx]`. @return 1 if `var` was set, 0 otherwise
    /** An internal counter is incremented to record that the value was read */
    static unsigned least_used_index(rec_type const* rec)
    {
        unsigned i = 0U;
        unsigned r = ~0U;
        for ( unsigned v = 0; v < rec->size(); ++v )
        {
            if ( rec->at(v).defined_ && rec->at(v).read_ < r )
            {
                r = rec->at(v).read_;
                i = v;
            }
        }
        return i;
    }
    
    std::string least_used_value(key_type const& key)
    {
        rec_type const* rec = values(key);
        if ( rec  &&  0 < rec->size() )
        {
            unsigned i = least_used_index(rec);
            val_type const& val = rec->at(i);
            ++val.read_;
            return val.value_;
        }
        return "";
    }
    
    /// try to set `var` from `key[inx]`. @return 1 if `var` was set, 0 otherwise
    /** An internal counter is incremented to record that the value was read */
    template <typename T>
    int set_from_least_used_value(T & var, key_type const& key) const
    {
        rec_type const* rec = values(key);
        
        if ( rec  &&  0 < rec->size() )
        {
            unsigned i = least_used_index(rec);
            //std::clog << key << " least used " << i << "\n";
            val_type const& val = rec->at(i);
            if ( val.defined_ )
            {
                set_value(var, key, val.value_);
                ++val.read_;
                return 1;
            }
        }
        return 0;
    }

    /// true if value of `key[inx]` is composed of alpha characters and '_'
    int is_alpha(key_type const& key, unsigned inx) const
    {
        rec_type const* rec = values(key);
        if ( !rec || inx >= rec->size() )
            return 0;
        return is_alpha(rec->at(inx).value_);
    }
    
    
    /// check if value of `key[inx]` could be a number
    /**
     @returns a bitwise field with individual bits set as follows:
     - 1 if not negative (that is >= 0 )
     - 2 if integer
     - 4 if float
     .
     eg. result is 0 if not a number; result is 3 for a positive integer
     */
    int is_number(key_type const& key, unsigned inx) const
    {
        rec_type const* rec = values(key);
        if ( !rec || inx >= rec->size() )
            return 0;
        return is_number(rec->at(inx).value_);
    }
    
    int is_positive_integer(key_type const& key, unsigned inx) const
    {
        return is_number(key, inx) == 3;
    }
    
    template <typename T>
    int set_positive_integer(T & var, key_type const& key, unsigned inx = 0) const
    {
        rec_type const* rec = values(key);
        if ( rec && inx < rec->size() )
        {
            val_type const& val = rec->at(inx);
            
            if ( val.defined_ && is_number(rec->at(inx).value_)==3 )
            {
                set_value(var, key, val.value_);
                ++val.read_;
                return 1;
            }
        }
        return 0;
    }

};


#pragma mark -


/// special function for std::string arguments.
template <>
void Glossary::set_value(std::string& var, key_type const&, std::string const&);

/// special function for float
template <>
void Glossary::set_value(float& var, key_type const&, std::string const&);

/// special function for double
template <>
void Glossary::set_value(double& var, key_type const&, std::string const&);


/// input from stream
std::istream& operator >> (std::istream&, Glossary&);

/// output of one value
std::ostream& operator << (std::ostream&, const Glossary::pair_type&);

/// output of all values
std::ostream& operator << (std::ostream&, const Glossary&);

#endif


