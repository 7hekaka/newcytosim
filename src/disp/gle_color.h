// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#ifndef GLE_COLOR_H
#define GLE_COLOR_H

#include "opengl.h"
#include <string>
#include <iostream>

/**
 gle_color implements colors with 4-components:
 - Red
 - Green
 - Blue
 - Alpha = transparency
 .
 
 This class implements the `RGBA` format using an 'unsigned integer'
 and an array of 4 floats.
 
 F. Nedelec -- Merged two older color classes on 23.08.2015
 */
/// Color with 4 components: red, green, blue, alpha (RGBA)
class gle_color
{
#pragma mark - Static methods
public:

    /// type used to quantify color components
    typedef GLfloat COLOF;

    
    static size_t stride()
    {
        return sizeof(gle_color) - 4 * sizeof(COLOF);
    }

private:

    /// concatenate 4 bytes into an int
    static uint32_t combine(uint32_t R, uint32_t G, uint32_t B, uint32_t A)
    {
        constexpr uint32_t K = 0xFF;
        return (R&K) << 24 | (G&K) << 16 | (B&K) << 8 | (A&K);
    }
    
    /// concatenate 4 bytes into an int
    static uint32_t combine(COLOF R, COLOF G, COLOF B, COLOF A)
    {
        return combine(uint32_t(255*R), uint32_t(255*G), uint32_t(255*B), uint32_t(255*A));
    }
    
    /// return value clamped to [0, 1]
    static COLOF clamp(COLOF s) { return std::max(COLOF(0), std::min(s, COLOF(1))); }
    
    /// array of 4 COLOF components, matching the `rgba_` integer
    COLOF col_[4];
    
    /// 32-bits integer containing 4 one-byte components: red, green, blue, alpha
    uint32_t rgba_;

#pragma mark - Private methods

    /// update 'rgba_' to match values in 'col_'
    void update_rgba()
    {
        rgba_ = combine(col_[0], col_[1], col_[2], col_[3]);
    }
    
    /// update 'col_' to match values in 'rgba_'
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
        rgba_ = 0xFFFFFFFF;
        col_[0] = 1;
        col_[1] = 1;
        col_[2] = 1;
        col_[3] = 1;
    }
    
    /// set to black
    void set_black()
    {
        rgba_ = 0x000000FF;
        col_[0] = 0;
        col_[1] = 0;
        col_[2] = 0;
        col_[3] = 1;
    }

    /// specify floating point components
    void set_float(COLOF r, COLOF g, COLOF b, COLOF a)
    {
        col_[0] = clamp(r);
        col_[1] = clamp(g);
        col_[2] = clamp(b);
        col_[3] = clamp(a);
        update_rgba();
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
        rgba_ = combine(uint32_t(r), uint32_t(g), uint32_t(b), uint32_t(a));
        update_float(rgba_);
    }

    /// export components as bytes
    void put_bytes(uint8_t& r, uint8_t& g, uint8_t& b, uint8_t& a) const
    {
        r = 0xFF & (uint8_t)( rgba_ >> 24 );
        g = 0xFF & (uint8_t)( rgba_ >> 16 );
        b = 0xFF & (uint8_t)( rgba_ >> 8 );
        a = 0xFF & (uint8_t)( rgba_ );
    }
    
#pragma mark - Constructors

    /// default constructor
    gle_color() : rgba_(0)
    {
        col_[0] = 0;
        col_[1] = 0;
        col_[2] = 0;
        col_[3] = 0;
    }
    
    /// constructor
    gle_color(const uint32_t& u)
    {
        rgba_ = u;
        update_float(u);
    }
    
    /// constructor from RGB values, with Alpha component = 1.0
    gle_color(const COLOF& r, const COLOF& g, const COLOF& b)
    {
        set_float(r,g,b,1.0f);
    }

    /// constructor from RGBA components
    gle_color(const COLOF& r, const COLOF& g, const COLOF& b, const COLOF& a)
    {
        set_float(r,g,b,a);
    }
    
#pragma mark - Public methods

    void operator = (const uint32_t& arg)
    {
        rgba_ = arg;
        update_float(arg);
    }
    
    bool operator ==(const gle_color col) const { return rgba_ == col.rgba_; }
    bool operator !=(const gle_color col) const { return rgba_ != col.rgba_; }
    
    COLOF const* colors() const { return col_; }

    /// access to float components
    COLOF& operator [] (int i) { return col_[i]; }

    COLOF r() const { return col_[0]; }
    COLOF g() const { return col_[1]; }
    COLOF b() const { return col_[2]; }
    COLOF a() const { return col_[3]; }
    
    COLOF red()   const { return col_[0]; }
    COLOF green() const { return col_[1]; }
    COLOF blue()  const { return col_[2]; }
    COLOF alpha() const { return col_[3]; }

    void set_red  (COLOF s) { col_[0] = clamp(s); update_rgba(); }
    void set_green(COLOF s) { col_[1] = clamp(s); update_rgba(); }
    void set_blue (COLOF s) { col_[2] = clamp(s); update_rgba(); }
    void set_alpha(COLOF s) { col_[3] = clamp(s); update_rgba(); }

    gle_color red  (COLOF s) const { return gle_color(clamp(s), col_[1], col_[2], col_[3]); }
    gle_color green(COLOF s) const { return gle_color(col_[0], clamp(s), col_[2], col_[3]); }
    gle_color blue (COLOF s) const { return gle_color(col_[0], col_[1], clamp(s), col_[3]); }
    gle_color alpha(COLOF s) const { return gle_color(col_[0], col_[1], col_[2], clamp(s)); }

    gle_color match_r(gle_color c) const { return gle_color(c.col_[0], col_[1], col_[2], col_[3]); }
    gle_color match_g(gle_color c) const { return gle_color(col_[0], c.col_[1], col_[2], col_[3]); }
    gle_color match_b(gle_color c) const { return gle_color(col_[0], col_[1], c.col_[2], col_[3]); }
    gle_color match_a(gle_color c) const { return gle_color(col_[0], col_[1], col_[2], c.col_[3]); }
    
#pragma mark -

    bool  visible()      const { return ( rgba_ & 0xFF ); }
    bool  invisible()    const { return ( rgba_ & 0xFF ) == 0; }
    COLOF normSqr()      const { return col_[0]*col_[0] + col_[1]*col_[1] + col_[2]*col_[2]; }
    COLOF brightness()   const { return normSqr() * col_[3]; }
    
    bool  opaque()       const { return ( (rgba_ & 0xFF) == 0xFF ); }
    bool  transparent()  const { return ( (rgba_ & 0xFF) != 0xFF ); }
    COLOF transparency() const { return col_[3]; }
    
    COLOF difference(gle_color back) const
    {
        COLOF x = col_[0] - back.col_[0];
        COLOF y = col_[1] - back.col_[1];
        COLOF z = col_[2] - back.col_[2];
        return x*x + y*y + z*z;
    }

#pragma mark -

    gle_color darken(COLOF s) const
    {
        COLOF x = clamp(s);
        return gle_color(x*col_[0], x*col_[1], x*col_[2], col_[3]);
    }
    
    gle_color lighten(COLOF s) const
    {
        return gle_color(s*col_[0], s*col_[1], s*col_[2], col_[3]);
    }
    
    gle_color alpha_scaled(COLOF s) const
    {
        return gle_color(col_[0], col_[1], col_[2], clamp(s*col_[3]));
    }
    
    gle_color blend(gle_color C) const
    {
        COLOF s = a() + C.a();
        COLOF h = a()   / s;
        COLOF g = C.a() / s;
        return gle_color(h*col_[0]+g*C[0], h*col_[1]+g*C[1], h*col_[2]+g*C[2], 0.5*(h+g));
    }
    
    gle_color blend(COLOF g, gle_color B, COLOF h)
    {
        return gle_color(g*col_[0]+h*B[0], g*col_[1]+h*B[1], g*col_[2]+h*B[2], g*col_[3]+h*B[3]);
    }

    gle_color inverted() const
    {
        return gle_color(1-col_[0], 1-col_[1], 1-col_[2], col_[3]);
    }
    
    gle_color tweak(uint32_t arg) const
    {
        gle_color C(arg);
        constexpr COLOF A = 0.5, B = 0.5;
        return gle_color(col_[0]*A+C[0]*B, col_[1]*A+C[1]*B, col_[2]*A+C[2]*B, col_[3]);
    }
    
#pragma mark -
    
    /// set current OpenGL color by calling glColor
    void load() const
    {
        //std::clog << "OpenGL load " << col_[0] << " " << col_[1] << " " << col_[2] << " " << col_[3] << "\n";
        glColor4fv(col_);
        //glColor4f(col_[0], col_[1], col_[2], col_[3]);
    }
    
    /// set current OpenGL color, but with `a` as alpha component
    void load(COLOF a) const
    {
        glColor4f(col_[0], col_[1], col_[2], clamp(a));
    }
    
    void load_clear() const
    {
        glClearColor(col_[0], col_[1], col_[2], col_[3]);
    }
    
    static void no_emission(GLenum face)
    {
        COLOF blk[4] = { 0, 0, 0, 1 };
        glMaterialfv(face, GL_EMISSION, blk);
    }

    /// set FRONT material property for lighting
    void load_front() const
    {
        glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, col_);
        no_emission(GL_FRONT);
    }
    
    /// set FRONT material property for lighting, and current color
    void load_load() const
    {
        glColor4fv(col_);
        glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, col_);
        no_emission(GL_FRONT);
    }
    
    /// set FRONT material property for lighting, and current color with given alpha-component
    void load_load(COLOF a) const
    {
        COLOF mat[4] = { col_[0], col_[1], col_[2], clamp(a) };
        glColor4fv(mat);
        glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, mat);
        no_emission(GL_FRONT);
    }

    /// set front OpenGL color, with `a` as alpha component
    void load_front(COLOF a) const
    {
        COLOF mat[4] = { col_[0], col_[1], col_[2], clamp(a) };
        glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, mat);
        no_emission(GL_FRONT);
    }

    /// set BACK material property for lighting
    void load_back() const
    {
        //glMaterialfv(GL_BACK, GL_AMBIENT_AND_DIFFUSE, col_);
        COLOF blk[4] = { 0, 0, 0, 1 };
        glMaterialfv(GL_BACK, GL_AMBIENT, col_);
        glMaterialfv(GL_BACK, GL_DIFFUSE, blk);
        glMaterialfv(GL_BACK, GL_EMISSION, blk);
    }
    
    /// set BACK material property for lighting, but with `a` as alpha component
    void load_back(COLOF a) const
    {
        COLOF mat[4] = { col_[0], col_[1], col_[2], clamp(a) };
        glMaterialfv(GL_BACK, GL_AMBIENT_AND_DIFFUSE, mat);
        no_emission(GL_BACK);
    }
    
    /// set FRONT and BACK material property for lighting
    void load_both() const
    {
#if 0
        COLOF blk[4] = { 0, 0, 0, 0 };
        glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, col_);
        glMaterialfv(GL_BACK, GL_AMBIENT, col_);
        glMaterialfv(GL_BACK, GL_DIFFUSE, blk);
#else
        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, col_);
#endif
        no_emission(GL_FRONT_AND_BACK);
    }
    
    /// set FRONT and BACK material property for lighting
    void load_both(COLOF a) const
    {
        COLOF blk[4] = { 0, 0, 0, 1 };
        COLOF mat[4] = { col_[0], col_[1], col_[2], clamp(a) };
        glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, mat);
        glMaterialfv(GL_BACK, GL_AMBIENT, mat);
        glMaterialfv(GL_BACK, GL_DIFFUSE, blk);
        no_emission(GL_FRONT_AND_BACK);
    }
    
    /// set FRONT and BACK material property for lighting
    void load_emission() const
    {
        COLOF blk[4] = { 0, 0, 0, 1 };
        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, blk);
        glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, col_);
    }
    
#pragma mark -
    
    /// conversion function from RGB to HSV color space
    static void RGB2HSV(COLOF r, COLOF g, COLOF b, COLOF* h, COLOF* s, COLOF* v);
    
    /// conversion functions from HSV to RGB color space
    static void HSV2RGB(COLOF h, COLOF s, COLOF v, COLOF* r, COLOF* g, COLOF* b);
    
    
    /// set a RGB color from a factor in [-PI, PI], continuously varying through all colors
    static void set_hue_components(COLOF& r, COLOF& g, COLOF& b, COLOF h);

    /// return new saturated color with given Hue value `h` in [-PI, PI]
    static gle_color hue_color(COLOF h, COLOF alpha = 1.0f);
    
    /// return new saturated color with Hue value `atan2(y, x)`
    static gle_color radial_color(COLOF x, COLOF y, COLOF alpha);

    /// return color build from a normalized 3D vector {x, y, z}
    static gle_color radial_color(COLOF x, COLOF y, COLOF z, COLOF alpha);

    /// return new jet color for h in [0, 5] with specified alpha component
    static gle_color jet_color(COLOF h, COLOF alpha = 1.0f);
    
    /// return new jet color extended
    static gle_color jet_color_dark(COLOF h, COLOF alpha = 1.0f);
    
    /// return new jet color extended
    static gle_color jet_color_alpha(COLOF h);

    /// print color in hexadecimal format (str must be of size 12)
    void hexadecimal(char* str) const;

    /// conversion of color to hexadecimal format
    std::string to_string() const;
    
};


/// input operator:
std::istream& operator >> (std::istream&, gle_color&);

/// output operator:
std::ostream& operator << (std::ostream&, const gle_color&);


#endif
