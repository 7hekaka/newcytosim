// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University.

#ifndef FLUTE_H
#define FLUTE_H

#include <cmath>


//#include "vector1.h"
//#include "vector2.h"
//#include "vector3.h"


/// accessory class to pack 1D vertex data
struct flute1
{
    float xy[2];
    flute1() : xy{0, 0} {}
    flute1(float x) : xy{x, 0} {}
    void operator = (const float v[1]) { xy[0] = v[0]; xy[1] = 0.f; }
    void operator = (const double v[1]) { xy[0] = float(v[0]); xy[1] = 0.f; }
};

/// accessory class to pack 2D vertex data
struct flute2
{
    float xy[2];
    flute2() : xy{0, 0} {}
    flute2(float x, float y) : xy{x, y} {}
    //flute2(double x, double y) : xy{float(x), float(y)} {}

    flute2(const float v[2]) : xy{v[0], v[1]} {}
    flute2(const double v[2]) : xy{float(v[0]), float(v[1])} {}
    void operator = (const float v[2]) { xy[0] = v[0]; xy[1] = v[1]; }
    void operator = (const double v[2]) { xy[0] = float(v[0]); xy[1] = float(v[1]); }

    void set(float x, float y) { xy[0] = x; xy[1] = y; }
    void set(double x, float y) { xy[0] = float(x); xy[1] = y; }
    void set(double x, double y) { xy[0] = float(x); xy[1] = float(y); }
#if 0
    void operator = (Vector1 const& v) { xy[0] = float(v.XX); xy[1] = 0.f; }
    flute2(Vector2 const& v) : xy{float(v.XX), float(v.YY)} {}
    void operator = (Vector2 const& v) { xy[0] = float(v.XX); xy[1] = float(v.YY); }
#endif
};


/// accessory class to pack 3D vertex data
struct flute3
{
    float xyz[3];
    flute3() : xyz{0, 0, 0} {}
    flute3(float x, float y, float z) : xyz{x, y, z} {}
    //flute3(double x, double y, double z) : xyz{float(x), float(y), float(z)} {}

    //flute3(float* ptr) : xyz{ptr[0], ptr[1], ptr[2]} {}
    void set(float x, float y, float z) { xyz[0] = x; xyz[1] = y; xyz[2] = z; }
    void set(double x, double y, double z) { xyz[0] = float(x); xyz[1] = float(y); xyz[2] = float(z); }
    
    flute3(void*);
    flute3(const float v[3]) : xyz{v[0], v[1], v[2]} {}
    flute3(const double v[3]) : xyz{float(v[0]), float(v[1]), float(v[2])} {}
    void operator = (const float v[3]) { xyz[0] = v[0]; xyz[1] = v[1]; xyz[2] = v[2]; }
    void operator = (const double v[3]) { xyz[0] = float(v[0]); xyz[1] = float(v[1]); xyz[2] = float(v[2]); }
    operator float const*() const { return xyz; }
    operator float*() { return xyz; }

#if 0
    flute3(Vector3 const& v) : xyz{float(v.XX), float(v.YY), float(v.ZZ)} {}
    //void operator = (Vector1 const& v) { xyz[0] = float(v.XX); xyz[1] = 0.f; xyz[2] = 0.f; }
    //void operator = (Vector2 const& v) { xyz[0] = float(v.XX); xyz[1] = float(v.YY); xyz[2] = 0.f; }
    void operator = (Vector3 const& v) { xyz[0] = float(v.XX); xyz[1] = float(v.YY); xyz[2] = float(v.ZZ); }
#endif
    /// elementary vector operations:
    float operator[](size_t i) const { return xyz[i]; }
    friend flute3 operator +(flute3 const& a, flute3 const& b) { return flute3{a[0]+b[0], a[1]+b[1], a[2]+b[2]}; }
    friend flute3 operator -(flute3 const& a, flute3 const& b) { return flute3{a[0]-b[0], a[1]-b[1], a[2]-b[2]}; }
    friend flute3 operator -(flute3 const& a) { return flute3{-a[0], -a[1], -a[2]}; }
    friend flute3 operator *(float const& a, flute3 const& b) { return flute3{a*b[0], a*b[1], a*b[2]}; }
    void operator *=(float const& a) { xyz[0]*=a; xyz[1]*=a; xyz[2]*=a; }
    friend flute3 operator *(flute3 const& b, float const& a) { return flute3{a*b[0], a*b[1], a*b[2]}; }
    friend float normSqr(flute3 const& b) { return b[0]*b[0] + b[1]*b[1] + b[2]*b[2]; }
    friend float norm(flute3 const& b) { return std::sqrt(normSqr(b)); }
    friend flute3 normalize(flute3 const& b) { return (1/sqrt(normSqr(b))) * b; }
    friend flute3 cross(flute3 const& a, flute3 const& b) { return flute3{a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]}; }
    friend float dot(flute3 const& a, flute3 const& b) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
    //{ float a = 0.5; return flute3{a*b[0], a*b[1], a*b[2]}; }
};


/// accessory class to pack vertex or color data
struct flute4
{
    float xyz[4];
    flute4() : xyz{0, 0, 0, 0} {}
    flute4(float x, float y, float z, float t) : xyz{x, y, z, t} {}
    flute4(const float c[4]) : xyz{c[0], c[1], c[2], c[3]} {}
    //flute4(double x, double y, double z, double t) : xyz{float(x), float(y), float(z), float(t)} {}
};


/// accessory class to pack 1D vertex and color data together
struct flute5
{
    float xyz[6];
    flute5() : xyz{0, 0, 0, 0, 0, 0} {}
    flute5(float const c[4], const float d[1]) : xyz{c[0], c[1], c[2], c[3], d[0], 0.f} {}
    flute5(float const c[4], const double d[1]) : xyz{c[0], c[1], c[2], c[3], float(d[0]), 0.f} {}
    void setY(float y) { xyz[5] = y; }
    void scale(float S) { xyz[4] *= S; }
};


/// accessory class to pack 2D vertex and color data together
struct flute6
{
    float xyz[6];
    flute6() : xyz{0, 0, 0, 0, 0, 0} {}
    flute6(float x, float y, float z, float a, float b, float c) : xyz{x, y, z, a, b, c} {}
    //flute6(double x, double y, double z, double a, double b, double c) : xyz{float(x), float(y), float(z), float(a), float(b), float(c)} {}

    flute6(float const c[4], float x, float y) : xyz{c[0], c[1], c[2], c[3], x, y} {}
    flute6(float const c[4], const float d[2]) : xyz{c[0], c[1], c[2], c[3], d[0], d[1]} {}
    flute6(float const c[4], const double d[2]) : xyz{c[0], c[1], c[2], c[3], float(d[0]), float(d[1])} {}
    flute6(const float c[4], flute2 const& v) : xyz{c[0], c[1], c[2], c[3], v.xy[0], v.xy[1]} {}
    flute6(flute3 const& a, flute3 const& b) : xyz{a.xyz[0], a.xyz[1], a.xyz[2], b.xyz[0], b.xyz[1], b.xyz[2]} {}
    flute6(flute3 const& a, float x, float y, float z) : xyz{a.xyz[0], a.xyz[1], a.xyz[2], x, y, z} {}

#if 0
    flute6(const float c[4], Vector1 const& v) : xyz{c[0], c[1], c[2], c[3], float(v.XX), 0.f, } {}
    flute6(const float c[4], Vector2 const& v) : xyz{c[0], c[1], c[2], c[3], float(v.XX), float(v.YY)} {}
    flute6(Vector3 const& v, Vector3 const& w) : xyz{float(v.XX), float(v.YY), float(v.ZZ), float(w.XX), float(w.YY), float(w.ZZ)} {}
#endif
    void setY(float y) { xyz[5] = y; }
    void scale(float S) { xyz[4] *= S; xyz[5] *= S; }
};


/// accessory class to pack 3D vertex and color data together
struct flute8
{
    float xyz[8];
    flute8() : xyz{0, 0, 0, 0, 0, 0, 0, 0} {}
    flute8(float r, float g, float b, float a, float x, float y, float z, float t) : xyz{r, g, b, a, x, y, z, t} {}
    flute8(const float c[4], float x, float y, float z) : xyz{c[0], c[1], c[2], c[3], x, y, z, 1.f} {}
    flute8(const float c[4], flute3 const& v) : xyz{c[0], c[1], c[2], c[3], v.xyz[0], v.xyz[1], v.xyz[2], 1.f} {}
    flute8(const float c[4], const float d[4]) : xyz{c[0], c[1], c[2], c[3], d[0], d[1], d[2], 1.0f} {}
    flute8(const float c[4], const double d[4]) : xyz{c[0], c[1], c[2], c[3], float(d[0]), float(d[1]), float(d[2]), 1.0f} {}
#if 0
    flute8(const float c[4], Vector1 const& v) : xyz{c[0], c[1], c[2], c[3], float(v.XX), 0.f, 0.f, 1.f} {}
    flute8(const float c[4], Vector2 const& v) : xyz{c[0], c[1], c[2], c[3], float(v.XX), float(v.YY), 0.f, 1.f} {}
    flute8(const float c[4], Vector3 const& v) : xyz{c[0], c[1], c[2], c[3], float(v.XX), float(v.YY), float(v.ZZ), 1.f} {}
#endif
    void setY(float y) { xyz[5] = y; }
    void scale(float S) { xyz[4] *= S; xyz[5] *= S; xyz[6] *= S; }
};


#endif /* FLUTE_H */
