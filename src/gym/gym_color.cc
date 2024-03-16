// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#include "gym_color.h"
#include "gym_color_list.h"
#include "stream_func.h"
#include "exceptions.h"
#include "random.h"
#include <iomanip>
#include <cctype>
#include <cmath>


//-----------------------------------------------------------------------
#pragma mark Color Input/Output

static const char hexadecimal_digit[] = "0123456789ABCDEF";

/**
 Write '0xRRGGBBAA' to 'str[]'
 */
void gym_color::hexadecimal(char * str) const
{
    *str++ = '0';
    *str++ = 'x';
    for ( int i = 0; i < 4; ++i )
    {
        uint8_t d = uint8_t(255*col_[i]);
        if ( i == 3 && d == 255 ) break;
        *str++ = hexadecimal_digit[d&15];
        *str++ = hexadecimal_digit[d>>4];
    }
    *str = 0;
}


std::string gym_color::to_string() const
{
    char str[16] = { 0 };
    hexadecimal(str);
    return std::string(str);
}


uint8_t hex2byte(int c)
{
    if ( '0' <= c && c <= '9' )
        return (uint8_t)(c - '0');
    else if ( 'A' <= c && c <= 'F' )
        return (uint8_t)(c - ( 'A' - 10 ));
    else if ( 'a' <= c && c <= 'f' )
        return (uint8_t)(c - ( 'a' - 10 ));
    throw InvalidSyntax("invalid hexadecimal digit");
}

uint8_t hex2byte(int a, int b)
{
    return (uint8_t)( hex2byte(a) << 4 ) | hex2byte(b);
}

/**
 A color is composed of 4 components (Red, Green, Blue, Alpha),
 and can be specified in different ways:
 -# with a hexadecimal integer: 0xFF0000FF  or  0xff0000ff
 -# with 3 or 4 floats: (1 0 0) or (1 0 0 1)
 -# with a name:  red
 -# with a number: #1
 .
*/
std::istream& operator >> (std::istream& is, gym_color& col)
{
    col.set_white();

    std::istream::sentry s(is);
    if ( s )
    {
        int c = is.get();
        int d = is.peek();
        
        if ( isalpha(c) )
        {
            is.unget();
            std::string name;
            is >> name;
            if ( name == "random" )
            {
                c = RNG.pint32(gym::nb_alt_color());
                col = gym::alt_color(c);
                return is;
            }
            try {
                col = gym::std_color(name);
            }
            catch ( InvalidSyntax & e )
            {
                is.setstate(std::istream::failbit);
                std::cerr << e.brief() << " using white\n";
                gym::print_std_colors(std::cerr);
                col = 0xFFFFFFFF;
            }
        }
        else if ( '#'==c  &&  isdigit(d) )
        {
            unsigned int nb;
            is >> nb;
            col = gym::alt_color(nb);
        }
        else if ( '0'==c  &&  'x'==d )
        {
            is.get();
            uint8_t u[4] = { 0xFF, 0xFF, 0xFF, 0xFF };
            int i = 0;
            while ( i < 4 )
            {
                c = is.get();
                d = is.get();
                if ( c == EOF || d == EOF )
                {
                    if ( i < 3 )
                        throw InvalidSyntax("incomplete hexadecimal color specification");
                    is.clear();
                    break;
                }
                u[i++] = hex2byte(c, d);
            }
            col.set_bytes(u[0], u[1], u[2], u[3]);
        }
        else if ( isdigit(c) )
        {
            is.unget();
            gym_color::COLOF r, g, b, a=1;
            if ( is >> r >> g >> b )
            {
                is >> a;
                is.clear();
                col.set(r,g,b,a);
            }
        }
    }
    return is;
}


std::ostream& operator << (std::ostream& os, gym_color const& arg)
{
    char str[12] = { 0 };
    arg.hexadecimal(str);
    os << str;
    return os;
}


//-----------------------------------------------------------------------
#pragma mark - Color making primitives


/**
 r,g,b values are from 0 to 1
 h = [0, 360], s = [0,1], v = [0,1]
 if s == 0, then h = -1 (undefined)
*/

void gym_color::RGB2HSV(const COLOF r, const COLOF g, const COLOF b, COLOF* h, COLOF* s, COLOF* v)
{
    COLOF mn, mx, delta;
    mn = std::min(r, std::min(g, b));
    mx = std::max(r, std::max(g, b));
    *v = mx;
    delta = mx - mn;
    if ( mx != 0 )
        *s = delta / mx;
    else {
        *s = 0;
        *h = -1;
        return;
    }
    if ( r == mx )
        *h = ( g - b ) / delta;       // between yellow & magenta
    else if ( g == mx )
        *h = 2 + ( b - r ) / delta;   // between cyan & yellow
    else
        *h = 4 + ( r - g ) / delta;   // between magenta & cyan
    *h *= 60;                         // degrees
    if ( *h < 0 )
        *h += 360;
}


void gym_color::HSV2RGB(const COLOF h, const COLOF s, const COLOF v, COLOF* r, COLOF* g, COLOF* b)
{
    int i;
    COLOF f, p, q, t;
    if ( s == 0 ) {
        // achromatic (gray)
        *r = *g = *b = v;
        return;
    }
    COLOF hc = h/60;               // sector 0 to 5
    i = (int)std::floor(hc);
    f = hc - (COLOF)i;             // fractional part of h
    p = v * ( 1 - s );
    q = v * ( 1 - s * f );
    t = v * ( 1 - s * ( 1 - f ) );
    switch( i )
    {
        case 0:  *r = v; *g = t; *b = p; break;
        case 1:  *r = q; *g = v; *b = p; break;
        case 2:  *r = p; *g = v; *b = t; break;
        case 3:  *r = p; *g = q; *b = v; break;
        case 4:  *r = t; *g = p; *b = v; break;
        case 5:  *r = v; *g = p; *b = q; break;
        default: *r = 1; *g = 1; *b = 1; break;
    }
}


/**
 set a RGB color as a function of a Hue value `a` in [-PI, PI].
 The colors follow in this order: red, green, blue, red ...
*/
void gym_color::set_hue_components(COLOF& r, COLOF& g, COLOF& b, const COLOF h)
{
    COLOF x = 3 * COLOF( h * M_1_PI + 1 );
    int i = (int)std::floor(x);
    COLOF f = x-(COLOF)i;
    switch( i % 6 )
    {
        case 0: r = 1;   g = f;   b = 0;   break;
        case 1: r = 1-f; g = 1;   b = 0;   break;
        case 2: r = 0;   g = 1;   b = f;   break;
        case 3: r = 0;   g = 1-f; b = 1;   break;
        case 4: r = f;   g = 0;   b = 1;   break;
        case 5: r = 1;   g = 0;   b = 1-f; break;
        default: r = 1;  g = 0;   b = 0;   break; // never executed;
    }
}


/**
 set a RGB color as a function of a 3D vector with components in [-1, 1],
 with alpha-component equal to `a`.
 Two opposite vectors gives approximately complementary colors.
 */
gym_color gym_color::radial_color(const COLOF x, const COLOF y, const COLOF z, const COLOF a)
{
    COLOF pX = std::max(0.0f,-x), nX = std::min(0.0f,-x);
    COLOF pY = std::max(0.0f, y), nY = std::min(0.0f, y);
    COLOF pZ = std::max(0.0f, z), nZ = std::min(0.0f, z);
    return gym_color(pX-0.5f*(nY+nZ), pY-0.5f*(nX+nZ), pZ-0.5f*(nX+nY), a);
}

/**
 set a RGB color as a function of a 2D vector, using the angle in the XY plane
*/
gym_color gym_color::radial_colorXY(const COLOF x, const COLOF y, const COLOF a)
{
    COLOF r, g, b;
    set_hue_components(r, g, b, atan2f(y, x));
    return gym_color(r, g, b, a);
}

/**
 set a RGB color as a function of a Hue value `h` in [-PI, PI],
 with alpha-component equal to `a`.
 The colors follow in this order: red, green, blue, red ...
 */
gym_color gym_color::hue_color(const COLOF h, const COLOF a)
{
    COLOF r, g, b;
    set_hue_components(r, g, b, h);
    return gym_color(r, g, b, a);
}

/**
set a RGB color as a function of a value `h` in [0, 4] with transparency `a`.
The result varies from dark-blue, blue, cyan, yellow, orange to red:
- 0 : green
- 1 : blue
- 2 : red
*/
gym_color gym_color::bgr_color(const COLOF h, const COLOF a)
{
    COLOF r = clamp(h-1);
    COLOF g = clamp(1-h);
    COLOF b = clamp(h);
    return gym_color(r, g, b, a);
}

/**
set a RGB color as a function of a value `h` in [0, 4] with transparency `a`.
The result varies from dark-blue, blue, cyan, yellow, orange to red:
- 0 : dark-blue
- 1 : blue
- 2 : cyan
- 3 : red
- 4 : full red
*/
gym_color gym_color::jet_color(const COLOF h, const COLOF a)
{
    COLOF r, g, b;
    if ( h <= 0.4 )
    {
        r = 0;
        g = 0;
        b = 0.4f;
    }
    else if ( h < 4.0 )
    {
        int i = (int)std::floor(h);
        COLOF f = h-(COLOF)i;
        switch( i )
        {
            case 0:  r = 0;   g = 0;   b = f;   break;
            case 1:  r = 0;   g = f;   b = 1;   break;
            case 2:  r = f;   g = 1;   b = 1-f; break;
            case 3:  r = 1;   g = 1-f; b = 0;   break;
            default: r = 1;   g = 0;   b = 0;   break;
        }
    }
    else
    {
        r = 1;
        g = 0;
        b = 0;
    }
    return gym_color(r, g, b, a);
}

/**
set a RGB color as a function of a value h in [0, 4] with transparency `a`.
The result varies from black to white:
- 0 : black
- 1 : blue
- 2 : green
- 3 : red
- 4 : yellow
- 5 : white
*/
gym_color gym_color::dark_jet_color(const COLOF h, const COLOF a)
{
    COLOF r, g, b;
    if ( h <= 0.1 )
    {
        r = 0;
        g = 0;
        b = 0.1f;
    }
    else if ( h < 5.0 )
    {
        int i = (int)std::floor(h);
        COLOF f = h-(COLOF)i;
        switch( i )
        {
            case 0:  r = 0;   g = 0;   b = f;   break;
            case 1:  r = 0;   g = f;   b = 1-f; break;
            case 2:  r = f;   g = 1-f; b = 0;   break;
            case 3:  r = 1;   g = f;   b = 0;   break;
            case 4:  r = 1;   g = 1;   b = f;   break;
            default: r = 1;   g = 1;   b = 1;   break;
        }
    }
    else
    {
        r = 1;
        g = 1;
        b = 1;
    }
    return gym_color(r, g, b, a);
}


/**
set a RGB color as a function of a value h in [0, 4].
The result vary from black, blue, cyan, yellow, orange to red:
- 0 : black
- 1 : blue
- 2 : green
- 3 : red
- 4 : yellow
- 5 : white
*/
gym_color gym_color::jet_color_alpha(const COLOF h)
{
    COLOF r, g, b, a;
    if ( h <= 0 )
    {
        r = 0;
        g = 0;
        b = 0;
        a = 0;
    }
    else if ( h < 5.0 )
    {
        int i = (int)std::floor(h);
        COLOF f = h-(COLOF)i;
        switch( i )
        {
            case 0:  r = 0;   g = 0;   b = 1;   a = f; break;
            case 1:  r = 0;   g = f;   b = 1-f; a = 1; break;
            case 2:  r = f;   g = 1-f; b = 0;   a = 1; break;
            case 3:  r = 1;   g = f;   b = 0;   a = 1; break;
            case 4:  r = 1;   g = 1;   b = f;   a = 1; break;
            default: r = 1;   g = 1;   b = 1;   a = 1; break;
        }
    }
    else
    {
        r = 1;
        g = 1;
        b = 1;
        a = 1;
    }
    return gym_color(r, g, b, a);
}

