// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University
#ifndef VECTOR_FLOAT_H
#define VECTOR_FLOAT_H


/// accessory class to facilitates using Vector2 in OpenGL
struct float2
{
    float xy[2];
    float2() : xy{0, 0} {}
    float2(float x, float y) : xy{x, y} {}
    float2(Vector1 const& v, float y) : xy{float(v.XX), y} {}
    float2(Vector1 const& v, double y) : xy{float(v.XX), float(y)} {}
    float2(Vector2 const& v) : xy{float(v.XX), float(v.YY)} {}
    static float2 cast(double x, double y) { return float2(float(x), float(y)); }
};


/// accessory class to facilitates using Vector3 in OpenGL
struct float3
{
    float xyz[3];
    float3() : xyz{0, 0, 0} {}
    float3(float x, float y, float z) : xyz{x, y, z} {}
    float3(Vector1 const& v, float y, float z) : xyz{float(v.XX), y, z} {}
    float3(Vector2 const& v, float z) : xyz{float(v.XX), float(v.YY), z} {}
    float3(Vector3 const& v) : xyz{float(v.XX), float(v.YY), float(v.ZZ)} {}
    static float3 cast(double x, double y, double z) { return float3(float(x), float(y), float(z)); }
};

/// accessory class to set RGBA colors in OpenGL
struct float4
{
    float xyz[4];
    float4() : xyz{0, 0, 0, 0} {}
    float4(float x, float y, float z, float t) : xyz{x, y, z, t} {}
    float4(const float c[]) : xyz{c[0], c[1], c[2], c[3]} {}
    float4(gle_color const& c) : xyz{c.r(), c.g(), c.b(), c.a()} {}
    static float4 cast(double x, double y, double z, double t)
    { return float4(float(x), float(y), float(z), float(t)); }
};


/// accessory class to facilitates using Vector3 in OpenGL
struct float6
{
    float xyz[6];
    float6() : xyz{0, 0, 0, 0, 0, 0} {}
    float6(float x, float y, float z, float a, float b, float c) : xyz{x, y, z, a, b, c} {}
    float6(Vector3 const& v, Vector3 const& w) :
    xyz{float(v.XX), float(v.YY), float(v.ZZ), float(w.XX), float(w.YY), float(w.ZZ)} {}
    static float6 cast(double x, double y, double z, double a, double b, double c)
    { return float6(float(x), float(y), float(z), float(a), float(b), float(c)); }
};


/// accessory class to facilitates using RGBA colors in OpenGL
struct float8
{
    float xyz[8];
    float8() : xyz{0, 0, 0, 0, 0, 0, 0, 0} {}
    float8(float x, float y, float z, float a, float b, float c, float d, float e)
    : xyz{x, y, z, a, b, c, d, e} {}
    float8(gle_color const& c, gle_color const& d) :
    xyz{c.r(), c.g(), c.b(), c.a(), d.r(), d.g(), d.b(), d.a()} {}
};

#endif /* VECTOR_FLOAT_H */
