// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#ifndef GYM_COLOR_H
#define GYM_COLOR_H

#include <string>
#include <iostream>
#include <cstdint>

/**
 gym_color implements colors with 4-components:
 - Red
 - Green
 - Blue
 - Alpha = transparency
 .
 
 This class implements the `RGBA` format using an array of 4 floats.
 
 F. Nedelec -- Merged two older color classes on 23.08.2015
 */
/// Color with 4 components: red, green, blue, alpha (RGBA)
class gym_color
{
#pragma mark - Static methods
public:

    /// type used to quantify color components (should be float, normally)
    typedef float COLOF;

private:
    
    /// color components in [0, 1]: Red, Green, Blue, Alpha
    COLOF col_[4];

    /// concatenate 4 bytes into a 32-bit integer
    static uint32_t combine(uint32_t R, uint32_t G, uint32_t B, uint32_t A)
    {
        constexpr uint32_t K = 0xFF;
        return (R&K) << 24 | (G&K) << 16 | (B&K) << 8 | (A&K);
    }
    
    /// concatenate 4 bytes into an integer
    static uint32_t combine(COLOF R, COLOF G, COLOF B, COLOF A)
    {
        return combine(uint32_t(255*R), uint32_t(255*G), uint32_t(255*B), uint32_t(255*A));
    }
    
    /// return value clamped to [0, 1]
    static COLOF clamp(COLOF s) { return std::max(COLOF(0), std::min(s, COLOF(1))); }

#pragma mark - Private methods

    /// return integer representation of color ('col_')
    uint32_t rgba() const
    {
        return combine(col_[0], col_[1], col_[2], col_[3]);
    }
    
    /// update 'col_' to match RGBA values encoded in integer argument
    void update_float(uint32_t arg)
    {
        col_[0] = COLOF( 0xFF & ( arg >> 24 ) ) / 255;
        col_[1] = COLOF( 0xFF & ( arg >> 16 ) ) / 255;
        col_[2] = COLOF( 0xFF & ( arg >>  8 ) ) / 255;
        col_[3] = COLOF( 0xFF & arg ) / 255;
    }
    
#pragma mark - Public methods

public:
    
    /// set to white
    void set_white()
    {
        col_[0] = 1;
        col_[1] = 1;
        col_[2] = 1;
        col_[3] = 1;
    }
    
    /// set to black
    void set_black()
    {
        col_[0] = 0;
        col_[1] = 0;
        col_[2] = 0;
        col_[3] = 1;
    }

    /// specify floating point components
    void set(COLOF r, COLOF g, COLOF b)
    {
        col_[0] = clamp(r);
        col_[1] = clamp(g);
        col_[2] = clamp(b);
        col_[3] = 1.0;
    }

    /// specify floating point components
    void set(COLOF r, COLOF g, COLOF b, COLOF a)
    {
        col_[0] = clamp(r);
        col_[1] = clamp(g);
        col_[2] = clamp(b);
        col_[3] = clamp(a);
    }
    
    /// export floating point components
    void store(COLOF& r, COLOF& g, COLOF& b, COLOF& a) const
    {
        r = col_[0];
        g = col_[1];
        b = col_[2];
        a = col_[3];
    }
    
    /// export floating point components to array
    void store(COLOF c[4]) const
    {
        c[0] = col_[0];
        c[1] = col_[1];
        c[2] = col_[2];
        c[3] = col_[3];
    }

    /// specify components with bytes
    void set_bytes(uint8_t r, uint8_t g, uint8_t b, uint8_t a)
    {
        col_[0] = COLOF(r) / 255;
        col_[1] = COLOF(g) / 255;
        col_[2] = COLOF(b) / 255;
        col_[3] = COLOF(a) / 255;
    }

    /// export components as bytes
    void put_bytes(uint8_t& r, uint8_t& g, uint8_t& b, uint8_t& a) const
    {
        uint32_t col = rgba();
        r = 0xFF & (uint8_t)( col >> 24 );
        g = 0xFF & (uint8_t)( col >> 16 );
        b = 0xFF & (uint8_t)( col >> 8 );
        a = 0xFF & (uint8_t)( col );
    }
    
#pragma mark - Constructors

    /// default constructor
    gym_color()
    {
        col_[0] = 0;
        col_[1] = 0;
        col_[2] = 0;
        col_[3] = 0;
    }
    
    /// constructor
    gym_color(const uint32_t& u)
    {
        update_float(u);
    }
    
    /// constructor from RGB values, with Alpha component = 1.0
    gym_color(const COLOF& r, const COLOF& g, const COLOF& b)
    {
        set(r,g,b,1.0f);
    }

    /// constructor from RGBA components
    gym_color(const COLOF& r, const COLOF& g, const COLOF& b, const COLOF& a)
    {
        set(r,g,b,a);
    }
    
#pragma mark - Public methods

    void operator = (const uint32_t& arg)
    {
        update_float(arg);
    }
    
    /// true if colors are roughly equal
    bool operator ==(const gym_color col) const { return rgba() == col.rgba(); }
    
    /// true if colors are significantly different
    bool operator !=(const gym_color col) const { return rgba() != col.rgba(); }
    
    COLOF const* colors() const { return col_; }

    // conversion operator to float[4]
    operator COLOF const*() const { return col_; }

    /// access to float components
    COLOF& operator [] (int i) { return col_[i]; }

    COLOF r() const { return col_[0]; }
    COLOF g() const { return col_[1]; }
    COLOF b() const { return col_[2]; }
    COLOF a() const { return col_[3]; }
    
    uint8_t alpha_() const { return uint8_t(255*col_[3]); }

    COLOF red()   const { return col_[0]; }
    COLOF green() const { return col_[1]; }
    COLOF blue()  const { return col_[2]; }
    COLOF alpha() const { return col_[3]; }

    void set_red  (COLOF s) { col_[0] = clamp(s); }
    void set_green(COLOF s) { col_[1] = clamp(s); }
    void set_blue (COLOF s) { col_[2] = clamp(s); }
    void set_alpha(COLOF s) { col_[3] = clamp(s); }

    gym_color red  (COLOF s) const { return gym_color(clamp(s), col_[1], col_[2], col_[3]); }
    gym_color green(COLOF s) const { return gym_color(col_[0], clamp(s), col_[2], col_[3]); }
    gym_color blue (COLOF s) const { return gym_color(col_[0], col_[1], clamp(s), col_[3]); }
    gym_color alpha(COLOF s) const { return gym_color(col_[0], col_[1], col_[2], clamp(s)); }

    gym_color match_r(gym_color c) const { return gym_color(c.col_[0], col_[1], col_[2], col_[3]); }
    gym_color match_g(gym_color c) const { return gym_color(col_[0], c.col_[1], col_[2], col_[3]); }
    gym_color match_b(gym_color c) const { return gym_color(col_[0], col_[1], c.col_[2], col_[3]); }
    gym_color match_a(gym_color c) const { return gym_color(col_[0], col_[1], col_[2], c.col_[3]); }
    
#pragma mark -
    
    bool  visible()      const { return ( alpha_() > 0.0078125 ); }
    bool  invisible()    const { return ( alpha_() == 0 ); }
    bool  opaque()       const { return ( alpha_() == 255 ); }
    bool  transparent()  const { return ( alpha_() < 255 ); }
    
    COLOF transparency() const { return col_[3]; }
    COLOF normSqr()      const { return col_[0]*col_[0] + col_[1]*col_[1] + col_[2]*col_[2]; }

    /// Perceived brightness (Luminance), based on the ITU-R BT.601 standard
    COLOF brightness()   const { return (0.299f*col_[0] + 0.587f*col_[1] + 0.114f*col_[2]) * col_[3]; }

    COLOF difference(const gym_color& col) const
    {
        COLOF x = col_[0] - col.col_[0];
        COLOF y = col_[1] - col.col_[1];
        COLOF z = col_[2] - col.col_[2];
        return x*x + y*y + z*z;
    }

#pragma mark -

    gym_color darken(COLOF s) const
    {
        COLOF x = clamp(s);
        return gym_color(x*col_[0], x*col_[1], x*col_[2], col_[3]);
    }
    
    gym_color lighten(COLOF s) const
    {
        return gym_color(s*col_[0], s*col_[1], s*col_[2], col_[3]);
    }
    
    gym_color alpha_scaled(COLOF s) const
    {
        return gym_color(col_[0], col_[1], col_[2], clamp(s*col_[3]));
    }
    
    gym_color blend(gym_color C) const
    {
        COLOF s = a() + C.a();
        COLOF h = a()   / s;
        COLOF g = C.a() / s;
        return gym_color(h*col_[0]+g*C[0], h*col_[1]+g*C[1], h*col_[2]+g*C[2], 0.5*(h+g));
    }
    
    gym_color blend(COLOF g, gym_color B, COLOF h)
    {
        return gym_color(g*col_[0]+h*B[0], g*col_[1]+h*B[1], g*col_[2]+h*B[2], g*col_[3]+h*B[3]);
    }

    gym_color inverted() const
    {
        return gym_color(1.f-col_[0], 1.f-col_[1], 1.f-col_[2], col_[3]);
    }
    
    gym_color mix(gym_color C) const
    {
        constexpr COLOF A = 0.5, B = 0.5;
        return gym_color(col_[0]*A+C[0]*B, col_[1]*A+C[1]*B, col_[2]*A+C[2]*B, col_[3]);
    }

    gym_color tweak(uint32_t arg) const
    {
        return mix(gym_color(arg));
    }
    
#pragma mark -
    
    /// conversion function from RGB to HSV color space
    static void RGB2HSV(COLOF r, COLOF g, COLOF b, COLOF* h, COLOF* s, COLOF* v);
    
    /// conversion functions from HSV to RGB color space
    static void HSV2RGB(COLOF h, COLOF s, COLOF v, COLOF* r, COLOF* g, COLOF* b);
    
    
    /// set a RGB color from a factor in [-PI, PI], continuously varying through all colors
    static void set_hue_components(COLOF& r, COLOF& g, COLOF& b, COLOF h);

    /// return new saturated color with given Hue value `h` in [-PI, PI]
    static gym_color hue_color(COLOF h, COLOF alpha = 1.f);
    
    /// return new saturated color with Hue value `atan2(y, x)`
    static gym_color radial_colorXY(COLOF x, COLOF y, COLOF alpha);

    /// return color build from a normalized 3D vector {x, y, z}
    static gym_color radial_color(COLOF x, COLOF y, COLOF z, COLOF alpha);

    /// return new jet color for h in [0, 5] with specified alpha component
    static gym_color bgr_color(COLOF h, COLOF alpha = 1.f);

    /// return new jet color for h in [0, 5] with specified alpha component
    static gym_color jet_color(COLOF h, COLOF alpha = 1.f);
    
    /// return new jet color extended
    static gym_color dark_jet_color(COLOF h, COLOF alpha = 1.f);
    
    /// return new jet color extended
    static gym_color jet_color_alpha(COLOF h);

    /// print color in hexadecimal format (str must be of size 12)
    void hexadecimal(char* str) const;

    /// conversion of color to hexadecimal format
    std::string to_string() const;
    
};

/// input operator
std::istream& operator >> (std::istream&, gym_color&);

/// output operator
std::ostream& operator << (std::ostream&, const gym_color&);

    

#endif
