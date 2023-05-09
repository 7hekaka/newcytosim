// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University

#ifndef GLE_H
#define GLE_H

#include "real.h"
#include "gym_color.h"

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
    /// `finesse` affects the number of triangles used to draw shapes such as cylinders
    /** Higher values are better: 2 is okay, 4 is good, 8 is nice and 16 is very nice */
    constexpr size_t finesse = 4;
    
    /// number of circle points stored in buffer
    /** We use multiples of 5 to match circles and the icosahedron */
    constexpr size_t pi_5half = finesse * 25;
    constexpr size_t pi_twice = finesse * 20;
    constexpr size_t pi_3half = finesse * 15;
    constexpr size_t pi_once = finesse * 10;
    constexpr size_t pi_half = finesse * 5;

    /// values of cosine, sine over two full circonvolutions
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
    void circle1(float stroke_width);
    /// draw 2D circle of radius 1 in XY plane, with +Z as normal, fewer points
    void circle2(float stroke_width);
    /// draw 2D circle of radius 1 in XY plane, with +Z as normal, dotted
    void dottedCircle(float point_size);
    /// draw 2D circle of radius 1 in XY plane, with +Z as normal
    void circle(float radius, float stroke_width);

    /// draw 2D oval within -1 to 1
    void stroke_capsule(float stroke_width);
    /// 2D oval within -1 to 1
    void paint_capsule();
    /// 2D cross within -1.5 to 1.5
    void paint_cross();
    /// 2D cross within -1.5 to 1.5
    void stroke_cross();

    /// draw nice 2D disc of radius 1 in XY plane, with +Z as normal
    void disc1();
    /// draw 2D disc of radius 1 in XY plane, with +Z as normal
    void disc2();
    /// draw nice 2D disc of radius 1 at Z=1, with +Z as normal
    void discTop1();
    /// draw 2D disc of radius 1 at Z=1, with +Z as normal
    void discTop2();
    /// draw nice 2D disc of radius 1 in XY plane, with -Z as normal
    void discBottom1();
    /// draw 2D disc of radius 1 in XY plane, with -Z as normal
    void discBottom2();
    /// draw 2D disc of radius 1 in XY plane, at Z = 0.5
    void discMid2();
    /// draw 2D disc covering radius [1, 1.4142] in XY plane, at Z = 0
    void ring();
    /// draw 2D disc covering radius [1, 1.2] in XY plane, at Z = 0
    void thin_ring();

    /// draw nice 2D disc of radius 1 in XY plane, with +Z as normal
    inline void disc() { disc1(); }

    /// paint a disc in XY plane, covering all points at distance to origin [ R0, R1 ]
    void paint_halo(float R0, float R1);
    /// paint spherocylinder in 2D, using current color
    void paint_capsule(float left, float right, float rad, size_t inc=1);
    /// draw spherocylinder contour in 2D
    void stroke_capsule(float left, float right, float rad, float stroke_width, size_t inc=1);
    /// paint two spherocylinder in 2D, joined by a cylinder of size tube x clos
    void paint_bicapsule(float left, float right, float rad, float clos, float tube, size_t inc=1);
    /// draw two spherocylinder contours in 2D, joined by a cylinder of size tube x clos
    void stroke_bicapsule(float left, float right, float rad, float clos, float tube, float stroke_width, size_t inc=1);

    /// draw a tetrahedron of radius 1
    void tetrahedron();
    /// upside-down tetrahedron of radius 1
    void upsideTetra();
    /// draw a octahedron of radius 1
    void octahedron();
    /// draw squared-based pyramid of radius 1
    void pyramid();
    /// draw squared-based pyramid of radius 1 pointing down
    void invPyramid();
    /// draw a icosahedron of radius 1
    void icosahedron();
    /// draw a icosahedron of radius 1
    void ICOSAHEDRON();

    /// draw a roughly spherical shape made of few triangles
    void blob();
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
    
    /// draw a icosahedron without normals, as the first approximation of a sphere
    void icoid();
    /// draw a icosahedron without normals, as the first approximation of a sphere
    void icoidF();

    /// returns tetrahedron or octahedron
    inline void (*hedron(bool x))() { return x ? tetrahedron : cube; }

    /// display 3 arrow fins aligned with the Z axis, or radius 1, lenth 2, Z=[-0.5, 1.5]
    void arrowTail();

    /// draw an very nice open tube along Z, of diameter 1 and length 1
    void tube1();
    /// draw a nice open tube along Z, of diameter 1 and length 1
    void tube2();
    /// draw an open tube along Z, of radius 1 and length 1
    void tube4();
    /// draw a rough open tube along Z, of radius 1 and length 1
    void tube8();

    /// draw an open tube along Z, of radius 1 covering Z [0, 1+epsilon]
    void tubeS();
    /// draw an open tube along Z, of radius 1 covering Z [-epsilon, 1+epsilon]
    void tubeM();
    /// draw an open tube along Z, of radius 1 covering Z [-epsilon, 1]
    void tubeE();
    
    /// draw a tube along Z, of diameter 1 and length 1.5, Z=[-4, 256]
    void longTube1();
    /// draw a nicer tube along Z, of diameter 1 and length 1.5, Z=[-4, 256]
    void longTube2();
    /// draw a nicer tube along Z, of diameter 1 and length 1.5, Z=[-4, 256]
    void longTube4();
    /// draw a nice tube along Z, of diameter 1 with Z=[-256, 0]
    void halfTube1();
    /// draw a tube along Z, of diameter 1 with Z=[-256, 0]
    void halfTube2();
    /// draw a rough tube along Z, of diameter 1 with Z=[-256, 0]
    void halfTube4();
    /// draw a rough tube along Z, of diameter 1 with Z=[-256, 0], closed at Z=0
    void endedTube4();
    /// draw a closed cylinder along Z, of hexagonal crosssection with Z=[0, 1]
    void hexTube();
    /// draw a closed cylinder along Z, of hexagonal crosssection with Z=[0, 1]
    void thinTube();
    /// draw a closed cylinder along Z, of hexagonal crosssection with Z=[0, 256]
    void thinLongTube();
    /// draw tube with grooves
    void flutedTube();

    /// display a super nice cone of axis Z, radius 1 at Z=0, summit at Z=1
    void cone1();
    /// display a nicer cone of axis Z, radius 1 at Z=0, summit at Z=1
    void cone2();
    /// display a rough cone of axis Z, radius 1 at Z=0, summit at Z=1
    void cone4();
    /// display a closed cone directed along Z, of radius 1 in Z=[-1, +2]
    void longCone();
    /// display an open cone directed along Z, of radius 1 at Z=0
    void truncatedCone();

    /// display a cylinder of axis Z, radius 1 in Z=[0, 1]
    inline void cylinder1() { tube2(); discBottom2(); discTop2(); }
    /// display a cylinder of axis Z, radius 1 in Z=[-0.5, 0.5]
    void cylinderT();
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
    
    /// draw a rough sphere of radius 1 at the origin
    void sphere8();
    /// draw a sphere of radius 1 at the origin
    void sphere4();
    /// draw a nice sphere of radius 1 at the origin
    void sphere2();
    /// draw a very nice sphere of radius 1 at the origin
    void sphere1();
    
    void dualPassSphere1();
    void dualPassSphere2();
    void dualPassSphere4();
    void dualPassSphere8();

    /// draw a very nice half-sphere of radius 1 in Z < 0
    void hemisphere1();
    /// draw a nice half-sphere of radius 1 in Z < 0
    void hemisphere2();
    /// draw half a sphere of radius 1 in Z < 0
    void hemisphere4();
    /// draw a falt hemisphere of radius 1 in Z < 0
    void hemisphereF();

    /// draw a blob with a pointy ends up in Z
    void droplet();
    /// draw nicest sphere available
    inline void sphere() { sphere1(); }
    /// draw nicest hemisphere available
    inline void hemisphere() { hemisphere1(); }
    /// draw 12 pentagons distributed on a sphere
    void footballPentagons();
    /// draw a line on the sphere
    void baseballSeamCurve(float R, float W);
    /// draw a line on the sphere
    void tennisballSeamCurve(float R, float W);
    /// draw a sphere decorated with 12 pentagons
    void football();
    /// draw a white baseball
    void baseball();
    /// draw a yellow tennisball
    void tennisball();
    
    /// draw lines on the surface of a tube
    void stripedTube(float width);
    /// primitive used to draw the central segments of fibers
    inline void innerTube() { longTube4(); }
    /// primitive used to draw the MINUS ends of fibers
    inline void capedTube() { halfTube4(); hemisphereF(); }
    //inline void capedTube() { halfTube2(); discBottom2(); }
    /// primitive used to draw the PLUS ends of fibers
    inline void endedTube() { endedTube4(); }
    //inline void endedTube() { halfTube4(); discBottom2(); }

    //------------------------------------------------------------------------------
    #pragma mark -
    
    /// draw a band from A to B, with specified radius
    void drawBand(Vector2 const& A, Vector2 const& B, real);
    void drawSpikyBand(Vector2 const& A, Vector2 const& B, real);

    /// draw a band from A to B, with specified radius in A and B
    void drawBand(Vector1 const& a, float, Vector1 const& b, float);
    void drawBand(Vector2 const& a, float, Vector2 const& b, float);
    
    /// draw a band from A to B, with specified radius and colors in A and B
    void drawBand(Vector1 const& a, float, gym_color, Vector1 const& b, float, gym_color);
    void drawBand(Vector2 const& a, float, gym_color, Vector2 const& b, float, gym_color);

    /// draw symbol linking A to B
    void drawHourglass(Vector2 const& a, Vector2 const&, Vector2 const& b, Vector2 const&);
    void drawHourglass(Vector2 const&, Vector2 const&, gym_color,
                Vector2 const&, Vector2 const&, gym_color);
    /// draw symbol linking A to B
    void drawCross(Vector2 const& a, Vector2 const&, Vector2 const& b, Vector2 const&, real);
    void drawBar(Vector3 const& a, Vector3 const& da, Vector3 const& b, Vector3 const& db, real);
    
    /// draw two discs in A and B, connected with a line
    void drawDumbbell(Vector2 const& A, Vector2 const& B, float diameter);
    
    /// draw Dome built on a
    void drawTipi(real*, int, real);
    
    /// display cone, dir should be normalized
    void drawCone(Vector1 const& center, Vector1 const& dir, float rad);
    /// display arrow-head, dir should be normalized
    void drawCone(Vector2 const& center, Vector2 const& dir, float rad);
    /// display arrow-head, dir should be normalized
    void drawCone(Vector3 const& center, Vector3 const& dir, float rad);
    
    /// display cylinder, dir should be normalized
    void drawCylinder(Vector1 const& center, Vector1 const& dir, float rad);
    /// display cylinder, dir should be normalized
    void drawCylinder(Vector2 const& center, Vector2 const& dir, float rad);
    /// display cylinder, dir should be normalized
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

    void strokeCuboid(Vector3 const& A, Vector3 const& B, float width);
    void paintCuboid(Vector3 const& A, Vector3 const& B, float rad);
    void paintTetrahedron(Vector3 const&, Vector3 const&, Vector3 const&, Vector3 const&);
    void paintPrism(Vector3 const&, Vector3 const&, Vector3 const&, Vector3 const&, Vector3 const&, Vector3 const&, Vector3 const&, Vector3 const&);
    void paintSpikyPrism(Vector3 const&, Vector3 const&, Vector3 const&, Vector3 const&, Vector3 const&, Vector3 const&, Vector3 const&, Vector3 const&);
}


#endif
