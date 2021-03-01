// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.

#ifndef FLUTE_H
#define FLUTE_H

/// accessory class to pack vertex data
struct flute2
{
    float xy[2];
    flute2() : xy{0, 0} {}
    flute2(float x, float y) : xy{x, y} {}
    flute2(Vector1 const& v) : xy{float(v.XX), 0} {}
    flute2(Vector1 const& v, float y) : xy{float(v.XX), y} {}
    flute2(Vector1 const& v, double y) : xy{float(v.XX), float(y)} {}
    flute2(Vector2 const& v) : xy{float(v.XX), float(v.YY)} {}
    static flute2 cast(double x, double y) { return flute2(float(x), float(y)); }
};


/// accessory class to pack vertex data
struct flute3
{
    float xyz[3];
    flute3() : xyz{0, 0, 0} {}
    flute3(float x, float y, float z) : xyz{x, y, z} {}
    //flute3(float* ptr) : xyz{ptr[0], ptr[1], ptr[2]} {}
    flute3(Vector1 const& v, float y, float z) : xyz{float(v.XX), y, z} {}
    flute3(Vector2 const& v, float z) : xyz{float(v.XX), float(v.YY), z} {}
    flute3(Vector3 const& v) : xyz{float(v.XX), float(v.YY), float(v.ZZ)} {}
    static flute3 cast(double x, double y, double z) { return flute3(float(x), float(y), float(z)); }
    float operator[](size_t i) const { return xyz[i]; }
    friend flute3 operator +(flute3 const& a, flute3 const& b) { return flute3{a[0]+b[0], a[1]+b[1], a[2]+b[2]}; }
    friend flute3 operator -(flute3 const& a, flute3 const& b) { return flute3{a[0]-b[0], a[1]-b[1], a[2]-b[2]}; }
    friend flute3 operator *(float const& a, flute3 const& b) { return flute3{a*b[0], a*b[1], a*b[2]}; }
    friend float normSqr(flute3 const& b) { return b[0]*b[0] + b[1]*b[1] + b[2]*b[2]; }
    friend flute3 normalize(flute3 const& b) { return (1/sqrt(normSqr(b))) * b; }
    //{ float a = 0.5; return flute3{a*b[0], a*b[1], a*b[2]}; }
};

/// accessory class to pack vertex or color data
struct flute4
{
    float xyz[4];
    flute4() : xyz{0, 0, 0, 0} {}
    flute4(float x, float y, float z, float t) : xyz{x, y, z, t} {}
    flute4(const float c[]) : xyz{c[0], c[1], c[2], c[3]} {}
    flute4(Vector3 const& v) : xyz{float(v.XX), float(v.YY), float(v.ZZ), 0} {}
    flute4(gle_color const& c) : xyz{c.r(), c.g(), c.b(), c.a()} {}
    static flute4 cast(double x, double y, double z, double t)
    { return flute4(float(x), float(y), float(z), float(t)); }
};


/// accessory class to pack vertex and color data together
struct flute6
{
    float xyz[6];
    flute6() : xyz{0, 0, 0, 0, 0, 0} {}
    flute6(float x, float y, float z, float a, float b, float c) : xyz{x, y, z, a, b, c} {}
    flute6(Vector1 const& v, gle_color const& c) : xyz{float(v.XX), 0, c.r(), c.g(), c.b(), c.a()} {}
    flute6(Vector2 const& v, gle_color const& c) : xyz{float(v.XX), float(v.YY), c.r(), c.g(), c.b(), c.a()} {}
    flute6(flute2 const& v, gle_color const& c) : xyz{v.xy[0], v.xy[1], c.r(), c.g(), c.b(), c.a()} {}
    flute6(Vector3 const& v, Vector3 const& w) : xyz{float(v.XX), float(v.YY), float(v.ZZ), float(w.XX), float(w.YY), float(w.ZZ)} {}
    static flute6 cast(double x, double y, double z, double a, double b, double c)
    { return flute6(float(x), float(y), float(z), float(a), float(b), float(c)); }
};


/// accessory class to pack vertex and color data together
struct flute8
{
    float xyz[8];
    flute8() : xyz{0, 0, 0, 0, 0, 0, 0, 0} {}
    flute8(float x, float y, float z, float t, float r, float g, float b, float a) : xyz{x, y, z, t, r, g, b, a} {}
    flute8(Vector1 const& v, gle_color const& c) : xyz{float(v.XX), 0, 0, 1, c.r(), c.g(), c.b(), c.a()} {}
    flute8(Vector2 const& v, gle_color const& c) : xyz{float(v.XX), float(v.YY), 0, 1, c.r(), c.g(), c.b(), c.a()} {}
    flute8(Vector3 const& v, gle_color const& c) : xyz{float(v.XX), float(v.YY), float(v.ZZ), 1, c.r(), c.g(), c.b(), c.a()} {}
    flute8(flute3 const& v, gle_color const& c) : xyz{v.xyz[0], v.xyz[1], v.xyz[2], 1, c.r(), c.g(), c.b(), c.a()} {}
    flute8(gle_color const& c, gle_color const& d) : xyz{c.r(), c.g(), c.b(), c.a(), d.r(), d.g(), d.b(), d.a()} {}
};

#endif /* FLUTE_H */
