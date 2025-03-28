// Cytosim was created by Francois Nedelec. Copyright 2019 Cambridge University.
//
//  evaluator.h
//
//  Created by Francois Nedelec on 08/02/2019.
//  Copyright 2019 Cambridge University. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string>
#include <vector>

#include "exceptions.h"
#include "random.h"

/// a minimal math expression evaluator
/**
 This can evaluate boolean expressions like `X^2 + (Y-3)^3 < 4'
 We could use instead [tinyexpr](https://codeplea.com/tinyexpr)
*/
class Evaluator
{
    /// variable names must be a single letter
    typedef std::pair<std::string, real> variable;

    /// list of variables
    typedef std::vector<variable> variable_list;
    
private:

    /// pointer
    mutable char const* ptr;
    
    /// list of variables
    variable_list variables_;

    static void print_variables(std::ostream& os, variable_list const& list)
    {
        os << "Known variables:\n";
        for ( variable const& v : list )
            os << "   " << v.first << " = " << v.second << "\n";
    }
    
    static int function_(std::string arg)
    {
        if ( arg == "cos" ) return 1;
        if ( arg == "sin" ) return 2;
        if ( arg == "tan" ) return 3;
        if ( arg == "cosd" ) return 4;
        if ( arg == "sind" ) return 5;
        if ( arg == "tand" ) return 6;
        if ( arg == "cosh" ) return 8;
        if ( arg == "sinh" ) return 9;
        if ( arg == "tanh" ) return 10;
        if ( arg == "exp" ) return 16;
        if ( arg == "expm1" ) return 17;
        if ( arg == "log" ) return 18;
        if ( arg == "log2" ) return 19;
        if ( arg == "sqrt" ) return 20;
        if ( arg == "abs" ) return 21;
        if ( arg == "floor" ) return 32;
        if ( arg == "round" ) return 33;
        if ( arg == "trunc" ) return 34;
        if ( arg == "frac" ) return 35;
        if ( arg == "rand" ) return 40;
        if ( arg == "prand" ) return 41;
        if ( arg == "srand" ) return 42;
        return 0;
    }
    
    static real eval_function_(int fun, real x)
    {
        switch ( fun )
        {
            case 1: return cos(x);
            case 2: return sin(x);
            case 3: return tan(x);
            case 4: return cos(x*(M_PI/180));
            case 5: return sin(x*(M_PI/180));
            case 6: return tan(x*(M_PI/180));
            case 8: return cosh(x);
            case 9: return sinh(x);
            case 10: return tanh(x);
            case 16: return exp(x);
            case 17: return expm1(x);
            case 18: return log(x);
            case 19: return log2(x);
            case 20: return sqrt(x);
            case 21: return abs(x);
            case 32: return floor(x);
            case 33: return round(x);
            case 34: return trunc(x);
            case 35: return x-trunc(x);
            case 40: return RNG.preal();
            case 41: return RNG.preal();
            case 42: return RNG.sreal();
        }
        return x;
    }
    
    void skip_space() const
    {
        while ( isspace(*ptr) )
            ++ptr;
    }

    std::string token_() const
    {
        char const* s = ptr++;
        while ( isalnum(*ptr) )
            ++ptr;
        std::string tok;
        tok.append(s, ptr-s);
        //std::clog << "token: " << tok << "\n";
        return tok;
    }
    
    real value_(std::string var) const
    {
        for ( variable const& v : variables_ )
        {
            if ( var == v.first )
                return v.second;
        }
        print_variables(std::clog, variables_);
        throw InvalidSyntax("unknown variable '"+var+"'");
        return 0;
    }
    
    real number_() const
    {
        errno = 0;
        char * end = nullptr;
        real d = strtod(ptr, &end);
        if ( errno )
            throw InvalidSyntax("expected a number");
        ptr = end;
        //std::clog << "number: " << d << "  " << ptr << "\n";
        return d;
    }
    
    real factor_() const
    {
        //std::clog << "factor: " << ptr << "\n";
        skip_space();
        char c = *ptr;
        if ( '0' <= c && c <= '9' )
            return number_();
        if ( c == '(' )
        {
            ++ptr; // '('
            //std::clog << " (   " << ptr << "\n";
            real res = expression_();
            skip_space();
            if ( *ptr != ')' )
                throw InvalidSyntax("missing closing parenthesis");
            ++ptr; // ')'
            return res;
        }
        if ( c == '+' )
        {
            ++ptr;
            return factor_();
        }
        if ( c == '-' )
        {
            ++ptr;
            return -factor_();
        }
        if ( isalpha(*ptr) )
        {
            std::string tok = token_();
            int fun = function_(tok);
            if ( fun && *ptr == '(' )
            {
                ++ptr; // '('
                real val = NAN;
                if ( *ptr != ')' )
                    val = expression_();
                if ( *ptr != ')' )
                    throw InvalidSyntax("missing closing parenthesis");
                ++ptr; // ')'
                return eval_function_(fun, val);
            }
            return value_(tok);
        }
        throw InvalidSyntax("unexpected syntax");
    }
    
    real term_() const
    {
        //std::clog << "term: " << ptr << "\n";
        real result = factor_();
        while ( 1 )
        {
            skip_space();
            char c = *ptr;
            if ( c == '*' )
            {
                ++ptr;
                result *= factor_();
            }
            else if ( c == '/' )
            {
                ++ptr;
                result /= factor_();
            }
            else if ( c == '^' )
            {
                ++ptr;
                result = std::pow(result, factor_());
            }
            else
                return result;
        }
    }
    
    real expression_() const
    {
        //std::clog << "expression: " << ptr << "\n";
        real result = term_();
        while ( 1 )
        {
            skip_space();
            char c = *ptr;
            if ( c == '+' )
            {
                ++ptr;
                result += term_();
            }
            else if ( c == '-' )
            {
                ++ptr;
                result -= term_();
            }
            else if ( c == '<' )
            {
                ++ptr;
                return ( result < term_() );
            }
            else if ( c == '>' )
            {
                ++ptr;
                return ( result > term_() );
            }
            else if ( c == '&' )
            {
                ++ptr;
                real rhs = term_();
                return (( result != 0 ) && ( rhs != 0 ));
            }
            else if ( c == '|' )
            {
                ++ptr;
                real rhs = term_();
                return (( result != 0 ) || ( rhs != 0 ));
            }
            else
                return result;
        }
    }
    
public:
    
    Evaluator(std::initializer_list<variable> const& v) : variables_(v)
    {
        variables_.push_back({"pi", M_PI});
        //print_variables(std::clog, variables_);
    }
    
    /// evaluate expression in string argument
    real eval(char const* str) const
    {
        real res = NAN;
        try
        {
            ptr = str;
            res = expression_();
            //std::clog << "eval(" << str << ") = " << res << "\n";
        }
        catch( Exception& e )
        {
            e << " while evaluating `" << str << "'";
            throw;
        }
        return res;
    }
    
    /// evaluate expression in string argument
    real eval(std::string const& str) const
    {
        return eval(str.c_str());
    }
    
    /// return string-representation of evaluated expression
    std::string eval_(std::string const& str) const
    {
        real x = eval(str.c_str());
        long i = (long) x;
        if ( i == x )
            return std::to_string(i);
        return std::to_string(x);
    }
};
