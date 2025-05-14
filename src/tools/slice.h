// Cytosim was created by Francois Nedelec. Copyright 2025 Cambridge University

#include <cstdio>
#include <cstring>

/// Slice represents a regular subset of positive indices
/**
 A slice can be constructed from a string,
 where A, B and S are positive integers:
 
 A       just A
 A:B     from A to B, B included
 :B      from 0 to B, B included
 A:      from A to MAX_SIZE_T (a large number)
 A:S:B   from A to B, with interval S
 A:S:    from A to MAX_SIZE_T, with interval S
 */
class Slice
{
    size_t s; ///< start
    size_t i; ///< increment
    size_t e; ///< end
    
public:
    
    Slice()
    {
        s = 0;
        i = 1;
        e = ~0UL;
    }
    
    /// read string specification with given separator
    void parse(const char arg[], const char sep = ':')
    {
        s = 0;
        i = 1;
        e = ~0UL;

        errno = 0;
        char * str = nullptr;
        if ( *arg == sep ) {
            s = 0;
            e = strtoul(arg+1, &str, 10);
            if ( errno ) goto finish;
        } else if ( *arg ) {
            s = strtoul(arg, &str, 10);
            if ( errno ) goto finish;
            if ( *str == sep ) {
                e = strtoul(str+1, &str, 10);
                if ( errno ) goto finish;
            }
            else
                e = s;
        } else return;
        if ( *str == sep )
        {
            i = e;
            if ( str[1] == 0 )
                e = ~0UL;
            else
            {
                e = strtoul(str+1, &str, 10);
                if ( errno ) goto finish;
                if ( *str ) goto finish;
            }
        }
        return;
    finish:
        fprintf(stderr, "syntax error in `%s', expected START:INCREMENT:END\n", arg);
    }
    
    /// constructor from string
    Slice(const char arg[], const char sep = ':')
    {
        parse(arg, sep);
    }
    
    size_t first() const
    {
        return s;
    }

    size_t last() const
    {
        return e;
    }

    size_t increment() const
    {
        return i;
    }
    
    void first(size_t x)
    {
        s = x;
    }

    void last(size_t x)
    {
        e = x;
    }

    void increment(size_t x)
    {
        i = x;
    }

    size_t count() const
    {
        size_t res = 0;
        for ( size_t x = s; x <= e; x += i )
            ++res;
        return res;
    }
        
    bool match(size_t n) const
    {
        if ( n < s )
            return false;
        if ( e < n )
            return false;
        return 0 == ( n - s ) % i;
    }

    void print(FILE * file = stderr) const
    {
        if ( e < s )
            fprintf(file, "empty ");
        fprintf(file, "Slice(%lu:%lu:%lu)\n", s, i, e);
    }
};

