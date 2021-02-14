// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University
#ifndef VECTOR_FLOAT_H
#define VECTOR_FLOAT_H


/// accessory class to facilitates using Vector2 in OpenGL
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


/// accessory class to facilitates using Vector3 in OpenGL
struct flute3
{
    float xyz[3];
    flute3() : xyz{0, 0, 0} {}
    flute3(float x, float y, float z) : xyz{x, y, z} {}
    flute3(Vector1 const& v, float y, float z) : xyz{float(v.XX), y, z} {}
    flute3(Vector2 const& v, float z) : xyz{float(v.XX), float(v.YY), z} {}
    flute3(Vector3 const& v) : xyz{float(v.XX), float(v.YY), float(v.ZZ)} {}
    static flute3 cast(double x, double y, double z) { return flute3(float(x), float(y), float(z)); }
};

/// accessory class to set RGBA colors in OpenGL
struct flute4
{
    float xyz[4];
    flute4() : xyz{0, 0, 0, 0} {}
    flute4(float x, float y, float z, float t) : xyz{x, y, z, t} {}
    flute4(const float c[]) : xyz{c[0], c[1], c[2], c[3]} {}
    flute4(gle_color const& c) : xyz{c.r(), c.g(), c.b(), c.a()} {}
    static flute4 cast(double x, double y, double z, double t)
    { return flute4(float(x), float(y), float(z), float(t)); }
};


/// accessory class to facilitates using Vector3 in OpenGL
struct flute6
{
    float xyz[6];
    flute6() : xyz{0, 0, 0, 0, 0, 0} {}
    flute6(float x, float y, float z, float a, float b, float c) : xyz{x, y, z, a, b, c} {}
    flute6(Vector3 const& v, Vector3 const& w) :
    xyz{float(v.XX), float(v.YY), float(v.ZZ), float(w.XX), float(w.YY), float(w.ZZ)} {}
    static flute6 cast(double x, double y, double z, double a, double b, double c)
    { return flute6(float(x), float(y), float(z), float(a), float(b), float(c)); }
};


/// accessory class to facilitates using RGBA colors in OpenGL
struct flute8
{
    float xyz[8];
    flute8() : xyz{0, 0, 0, 0, 0, 0, 0, 0} {}
    flute8(float x, float y, float z, float a, float b, float c, float d, float e)
    : xyz{x, y, z, a, b, c, d, e} {}
    flute8(gle_color const& c, gle_color const& d) :
    xyz{c.r(), c.g(), c.b(), c.a(), d.r(), d.g(), d.b(), d.a()} {}
};

#endif /* VECTOR_FLOAT_H */
