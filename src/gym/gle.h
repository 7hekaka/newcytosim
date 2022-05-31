// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University

#ifndef GLE_H
#define GLE_H

#include "real.h"
#include "gle_color.h"

class Vector1;
class Vector2;
class Vector3;

/// Simple geometrical objects drawn with OpenGL
/**
 @todo namespace gle should be a class -> we use GL.vertex(v)
 Problem: the gle prefix is already used by standard OpenGL elements
 */
namespace gle
{
    /// this defines the number of triangles used to draw shapes
    /** Higher finesse improves the rendering:
     4 is okay, 8 is nice and 16 is very nice */
    constexpr size_t finesse = 8;
    
    /// number of circle points stored in buffer
    constexpr size_t pi_5half = finesse * 15;
    constexpr size_t pi_twice = finesse * 12;
    constexpr size_t pi_3half = finesse * 9;
    constexpr size_t pi_once = finesse * 6;
    constexpr size_t pi_half = finesse * 3;

    /// values of cosine, sine over a full circle
    extern float circle_[4*pi_twice+8];

    /// access to precomputed cosine
    inline float cos_(size_t n) { return circle_[2*n]; }
    
    /// access to precomputed sine
    inline float sin_(size_t n) { return circle_[1+2*n]; }

    /// initialize the arrays
    void initialize();
    
    /// release requested memory
    void quit();

    /// calculate sine and cosine
    void set_arc(size_t cnt, float CS[], double rad, double start, double delta, float cx, float cy);

    /// calculate sine and cosine for a circular arc
    void compute_arc(size_t cnt, float CS[], double rad, double start, double angle, float cx, float cy);

    /// initialize more buffer objects
    void setBuffers();

    //------------------------------------------------------------------------------
#pragma mark -

    /// draw 2D circle of radius 1 in XY plane, with +Z as normal
    void circle(float line_width);
    /// draw 2D circle of radius 1 in XY plane, with +Z as normal, fewer points
    void circle2(float line_width);
    /// draw 2D circle of radius 1 in XY plane, with +Z as normal, dotted
    void circle_dotted(float line_width);
    /// draw 2D oval within -1 to 1
    void capsule(float line_width);
    /// paint 2D within -1 to 1
    void capsule();

    /// draw 2D disc of radius 1 in XY plane, with +Z as normal
    void disc();
    /// draw 2D disc of radius 1 in XY plane, with +Z as normal
    void disc1();
    /// nicer 2D disc of radius 1 in XY plane, with +Z as normal
    void disc2();
    /// draw 2D disc of radius 1 at Z=1, with +Z as normal
    void discTop1();
    /// draw 2D disc of radius 1 at Z=1, with +Z as normal
    void discTop2();
    /// draw 2D disc of radius 1 in XY plane, with -Z as normal
    void discBottom1();
    /// nicer 2D disc of radius 1 in XY plane, with -Z as normal
    void discBottom2();

    /// paint a disc in XY plane, covering all points at distance to origin [ R0, R1 ]
    void paint_halo(float R0, float R1);
    /// paint spherocylinder in 2D, using current color
    void paint_capsule(float left, float right, float rad, size_t inc=1);
    /// draw spherocylinder contour in 2D
    void stroke_capsule(float left, float right, float rad, float line_width, size_t inc=1);
    /// paint two spherocylinder in 2D, joined by a cylinder of size tube x clos
    void paint_bicapsule(float left, float right, float rad, float clos, float tube, size_t inc=1);
    /// draw two spherocylinder contours in 2D, joined by a cylinder of size tube x clos
    void stroke_bicapsule(float left, float right, float rad, float clos, float tube, float line_width, size_t inc=1);

    /// draw a tetrahedron of radius 1
    void tetrahedron();
    /// upside-down tetrahedron of radius 1
    void upsideTetra();
    /// draw a octahedron of radius 1
    void octahedron();
    /// draw a icosahedron of radius 1
    void icosahedron();
    /// draw a icosahedron of radius 1
    void ICOSAHEDRON();

    /// draw a roughly spherical shape made of few triangles
    void blob();
    /// draw a icosahedron without normals
    void icoid();
    /// draw a centered blob of radius 1 with a cone extending up in Z
    void needle();
    /// draw a Cube of side 2
    void cube();
    /// draw a Cube of side 2
    void wireCube(float);
    /// draw a Cube of side 1
    void cuboid();
    /// draw a stellated octahedron
    void star();
    
    /// returns tetrahedron or octahedron
    inline void (*hedron(bool x))() { return x ? tetrahedron : cube; }

    /// display 3 arrow fins aligned with the Z axis, or radius 1, lenth 2, Z=[-0.5, 1.5]
    void arrowTail();

    /// draw an open tube along Z, of diameter 1 and length 1
    void tube1();
    /// draw an open tube along Z, of diameter 1 and length 1
    void tube2();
    /// draw an open tube along Z, of radius 1 covering Z [0, 1+epsilon]
    void tubeS();
    /// draw an open tube along Z, of radius 1 covering Z [-epsilon, 1+epsilon]
    void tubeM();
    /// draw an open tube along Z, of radius 1 covering Z [-epsilon, 1]
    void tubeE();
    /// draw an open tube along Z, of radius 1 covering Z [-epsilon, 1+epsilon]
    void doubleTubeE();
    /// draw an open tube along Z, of radius 1 covering Z [0, 1+epsilon]
    void doubleTubeF();
    /// draw a nice open tube along Z, of radius 1 and length 1
    void tube4();
    /// draw a nicer open tube along Z, of radius 1 and length 1
    void tube8();
    /// draw a tube along Z, of diameter 1 and length 1.5, Z=[-4, 256]
    void longTube1();
    /// draw a nicer tube along Z, of diameter 1 and length 1.5, Z=[-4, 256]
    void longTube2();
    /// draw a nicer tube along Z, of diameter 1 and length 1.5, Z=[-4, 256]
    void longTube4();
    /// draw a tube along Z, of diameter 1 with Z=[-256, 0]
    void halfTube1();
    /// draw a nicer tube along Z, of diameter 1 with Z=[-256, 0]
    void halfTube2();
    /// draw a nicer tube along Z, of diameter 1 with Z=[-256, 0]
    void halfTube4();
    /// draw a cylinder along Z, of hexagonal crosssection with Z=[0, 1]
    void hexTube();
    /// draw a cylinder along Z, of hexagonal crosssection with Z=[0, 1]
    void thinTube();
    /// draw a cylinder along Z, of hexagonal crosssection with Z=[0, 256]
    void thinLongTube();
    /// display a cone of axis Z, radius 1 at Z=0, summit at Z=1
    void cone1();
    /// display a nicer cone of axis Z, radius 1 at Z=0, summit at Z=1
    void cone2();
    /// display a long cone of axis Z, radius 1 at Z=0, summit at Z=3
    void cone3();

    /// display a cylinder of axis Z, radius 1 in Z=[0, 1]
    void cylinder1();
    /// display a cylinder of axis Z, radius 1 in Z=[-0.5, 0.5]
    void cylinder2();

    /// display a cone of axis Z, radius 1 at Z=0, summit at Z=1
    void cone();
    /// display a closed cone directed along Z, of radius 1 in Z=[-1, +2]
    void longCone();
    /// display a closed cone directed along Z, of radius 1.5 in Z=[0.7, +1.4]
    void shortCone();
    /// display an open cone directed along Z, of radius 1 at Z=0
    void truncatedCone();
    /// draw a 3-portion cylinder with a larger central section
    void barrel();
    /// draw a 3-portion cylinder with a larger central section
    void dualPassBarrel();
    /// display a dumbbell aligned with the Z axis, or radius 1/3, lenth 1
    void dumbbell();
    /// draw Torus of radius `rad` and thickness `thick`
    void torusZ(float rad, float thick, size_t inc = 1);

    /// draw a circular band composed of little triangles
    void arrowStrip(float width, size_t inc);
    /// draw 3 Arrowed Bands defining 8 quadrants on the sphere of radius 1
    void threeArrowStrip(float width, size_t inc);
 
    //------------------------------------------------------------------------------
    
    /// draw something
    void thing();
    /// do not draw
    void nothing();
    
    /// draw a sphere of radius 1 at the origin
    void sphere1();
    /// draw a nice sphere of radius 1 at the origin
    void sphere2();
    /// draw a very nice sphere of radius 1 at the origin
    void sphere4();
    /// draw a refined sphere of radius 1 at the origin
    void sphere8();
    
    void dualPassSphere1();
    void dualPassSphere2();
    void dualPassSphere4();
    void dualPassSphere8();

    /// draw half a sphere of radius 1 in Z < 0
    void hemisphere1();
    /// draw a nice half-sphere of radius 1 in Z < 0
    void hemisphere2();
    /// draw a very nice half-sphere of radius 1 in Z < 0
    void hemisphere4();

#if 1
    /// primitives to draw the MINUS ends of fibers:
    inline void capedTube1() { halfTube1(); hemisphere1(); }
    inline void capedTube2() { halfTube2(); hemisphere2(); }
    inline void capedTube4() { halfTube4(); hemisphere4(); }
#else
    /// primitives to draw the MINUS ends of fibers:
    inline void capedTube1() { halfTube1(); discBottom1(); }
    inline void capedTube2() { halfTube2(); discBottom2(); }
    inline void capedTube4() { halfTube4(); discBottom2(); }
#endif
    /// primitives to draw the PLUS ends of fibers:
    inline void endedTube1() { halfTube1(); discBottom1(); }
    inline void endedTube2() { halfTube2(); discBottom2(); }
    inline void endedTube4() { halfTube4(); discBottom2(); }

    //------------------------------------------------------------------------------
    #pragma mark -

    /// draw 'obj' with its ends at [a,b], of specified radius
    void stretchTube(Vector1 const& A, Vector1 const& B, float rad, void (*obj)());
    void stretchTube(Vector2 const& A, Vector2 const& B, float rad, void (*obj)());
    void stretchTube(Vector3 const& A, Vector3 const& B, float rad, void (*obj)());

    void drawTube(Vector1 const& A, float rad, Vector1 const& B, void (*obj)());
    void drawTube(Vector2 const& A, float rad, Vector2 const& B, void (*obj)());
    void drawTube(Vector3 const& A, float rad, Vector3 const& B, void (*obj)());
    
    /// draw a band from A to B, with specified radius
    void drawBand(Vector2 const& A, Vector2 const& B, real);
    void drawBand(Vector3 const& A, Vector3 const& B, real);

    /// draw a band from A to B, with specified radius in A and B
    void drawBand(Vector1 const& a, float, Vector1 const& b, float);
    void drawBand(Vector2 const& a, float, Vector2 const& b, float);
    
    /// draw a band from A to B, with specified radius and colors in A and B
    void drawBand(Vector1 const& a, float, gle_color, Vector1 const& b, float, gle_color);
    void drawBand(Vector2 const& a, float, gle_color, Vector2 const& b, float, gle_color);

    /// draw symbol linking A to B
    void drawHourglass(Vector2 const& a, Vector2 const&, Vector2 const& b, Vector2 const&);
    void drawHourglass(Vector2 const&, Vector2 const&, gle_color,
                Vector2 const&, Vector2 const&, gle_color);
    /// draw symbol linking A to B
    void drawCross(Vector2 const& a, Vector2 const&, Vector2 const& b, Vector2 const&, real);
    void drawBar(Vector3 const& a, Vector3 const& da, Vector3 const& b, Vector3 const& db, real);
    
    /// draw two discs in A and B, connected with a line
    void drawDumbbell(Vector2 const& A, Vector2 const& B, float diameter);

    /// display cone, dir should be normalized
    void drawCone(Vector1 const& center, Vector1 const& dir, float rad);
    /// display arrow-head, dir should be normalized
    void drawCone(Vector2 const& center, Vector2 const& dir, float rad);
    /// display arrow-head, dir should be normalized
    void drawCone(Vector3 const& center, Vector3 const& dir, float rad);
    
    /// display arrow-head, dir should be normalized
    void drawCylinder(Vector1 const& center, Vector1 const& dir, float rad);
    /// display arrow-head, dir should be normalized
    void drawCylinder(Vector2 const& center, Vector2 const& dir, float rad);
    /// display arrow-head, dir should be normalized
    void drawCylinder(Vector3 const& center, Vector3 const& dir, float rad);

    /// display arrow-head, dir should be normalized
    void drawArrowTail(Vector1 const& center, Vector1 const& dir, float rad);
    /// display arrow-head, dir should be normalized
    void drawArrowTail(Vector2 const& center, Vector2 const& dir, float rad);
    /// display arrow-head, dir should be normalized
    void drawArrowTail(Vector3 const& center, Vector3 const& dir, float rad);

    /// draw an arrow with ends [a,b], of specified radius
    void drawArrow(Vector1 const& A, Vector1 const& B, float rad);
    void drawArrow(Vector2 const& A, Vector2 const& B, float rad);
    void drawArrow(Vector3 const& A, Vector3 const& B, float rad);

    /// draw a set of 2 or 3 axes, depending on `dim`
    void drawAxes(float size, int dim);

    void drawCuboid(Vector3 const& A, Vector3 const& B, float w);
}


#endif
