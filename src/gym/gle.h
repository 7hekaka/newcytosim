// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University

#ifndef GLE_H
#define GLE_H

#include "real.h"
#include "gym_color.h"
#include "gym_flute.h"

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
    /// `finesse` sets the number of triangles used to draw shapes such as cylinders
    /** Higher values are better: 2 is okay, 4 is good, 8 is nice and 16 is very nice */
    constexpr unsigned finesse = 6;
    
    /// number of circle points stored in buffer
    /** We use multiples of 5 to match circles and the icosahedron */
    constexpr unsigned pi_6half = finesse * 30;
    constexpr unsigned pi_5half = finesse * 25;
    constexpr unsigned pi_4half = finesse * 20;
    constexpr unsigned pi_3half = finesse * 15;
    constexpr unsigned pi_2half = finesse * 10;
    constexpr unsigned pi_1half = finesse * 5;

    /// values of cosine, sine over two full circonvolutions
    extern float circle_[4*pi_4half+8];

    /// access to precomputed cosine
    inline float cos_(unsigned n) { return circle_[2*n]; }
    
    /// access to precomputed sine
    inline float sin_(unsigned n) { return circle_[1+2*n]; }

    /// calculate sine and cosine
    void set_arc(size_t cnt, float CS[], double radius, double start, double delta, double cx, double cy);

    /// calculate sine and cosine for a circular arc
    void compute_arc(size_t cnt, float CS[], double radius, double start, double angle, double cx, double cy);

    /// initialize buffer object
    size_t setExplodedCube(flute6*, float, float, float, float, float, float);

    /// initialize the arrays
    void initialize();
    
    /// release requested memory
    void quit();

    //------------------------------------------------------------------------------
#pragma mark -

    /// draw 2D square ranging [-1, 1]
    void square1(float stroke_width);
    /// draw 2D circle of radius 1 in XY plane, with +Z as normal
    void circle1(float stroke_width);
    /// draw 2D circle of radius 1 in XY plane, with +Z as normal, fewer points
    void circle2(float stroke_width);
    /// draw 2D circle of radius 1 in XY plane, with +Z as normal, dotted
    void dottedCircle(float point_size);
    /// draw 2D circle of radius 1 in XY plane, with +Z as normal
    void circle(float radius, float stroke_width);

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
    void thinRing();

    /// draw nice 2D disc of radius 1 in XY plane, with +Z as normal
    inline void disc() { disc1(); }

    /// draw 2D oval within -1 to 1
    void strokeCapsule(float stroke_width);
    /// 2D oval within -1 to 1
    void paintCapsule();
    /// 2D cross within -1.5 to 1.5
    void paintCross();
    /// 2D cross within -1.5 to 1.5
    void strokeCross();

    /// paint a disc in XY plane, covering all points at distance to origin [ R0, R1 ]
    void paintHalo(float R0, float R1);
    /// paint spherocylinder in 2D, using current color
    void paintCapsule(float left, float right, float rad, unsigned inc=1);
    /// draw spherocylinder contour in 2D
    void strokeCapsule(float left, float right, float rad, float stroke_width, unsigned inc=1);
    /// paint two spherocylinder in 2D, joined by a cylinder of size tube x clos
    void paintBicapsule(float left, float right, float rad, float clos, float tube, unsigned inc=1);
    /// draw two spherocylinder contours in 2D, joined by a cylinder of size tube x clos
    void strokeBicapsule(float left, float right, float rad, float clos, float tube, float stroke_width, unsigned inc=1);

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
    /// draw a smaller cube
    void smallCube();
    /// draw top and bottom faces of a cube
    void cubeFaces();
    /// draw a Cube of side 2
    void cubeEdges(float);
    /// draw vertical edges of a cube
    void cubeVerticalEdges(float);

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
    
    /// draw a very nice tube along Z, of diameter 1 and length 1.5, Z=[-4, 256]
    void longTube1();
    /// draw a nice tube along Z, of diameter 1 and length 1.5, Z=[-4, 256]
    void longTube2();
    /// draw a tube along Z, of diameter 1 and length 1.5, Z=[-4, 256]
    void longTube4();
    /// draw a nice tube along Z, of diameter 1 with Z=[-256, 0]
    void halfTube1();
    /// draw a tube along Z, of diameter 1 with Z=[-256, 0]
    void halfTube2();
    /// draw a rough tube along Z, of diameter 1 with Z=[-256, 0]
    void halfTube4();
    /// draw a tube along Z, of diameter 1 with Z=[-0.5, 0.5], closed at Z=0
    void shutTubeC();
    /// draw a tube along Z, of diameter 1 with Z=[-256, 0], closed at Z=0
    void shutTube2();
    /// draw a rough tube along Z, of diameter 1 with Z=[-256, 0], closed at Z=0
    void shutLongTube2();
    /// draw a closed cylinder along Z, of hexagonal crosssection with Z=[0, 1]
    void hexTube();
    /// draw a closed cylinder along Z, of hexagonal crosssection with Z=[0, 1]
    void thinTube();
    /// draw a closed cylinder along Z, of hexagonal crosssection with Z=[0, 256]
    void thinLongTube();
    /// draw lines on the surface of a tube
    void stripedTube(float line_width, gym_color);
    /// draw helix on the surface of a tube
    void helix();

    /// display a super nice cone of axis Z, radius 1 at Z=0, summit at Z=1
    void cone1();
    /// display a nicer cone of axis Z, radius 1 at Z=0, summit at Z=1
    void cone2();
    /// display a cone of axis Z, radius 1 at Z=-1, summit at Z=1
    void coneC();
    /// display a closed cone directed along Z, of radius 1 in Z=[-1, +2]
    void longCone();
    /// display an open cone directed along Z, of radius 1 at Z=0
    void cutCone();

    /// display cylinder of axis Z, radius 1 covering Z=[0, 1]
    void cylinder1();
    /// display cylinder of axis Z, radius 0.5 covering Z=[-1, 1]
    void cylinderC();

    /// draw a 3-portion cylinder with a larger central section
    void barrel();
    /// draw a 3-portion cylinder with a larger central section
    void dualPassBarrel();
    /// display a dumbbell aligned with the Z axis, of radius 1/3, lenth 1
    void dumbbell();
    /// draw Torus of radius `rad` and thickness `thick`
    void torusZ(float rad, float thick, unsigned inc = 1);

    /// draw a circular band composed of little triangles
    void arrowStrip(float width, unsigned inc);
    /// draw 3 Arrowed Bands defining 8 quadrants on the sphere of radius 1
    void threeArrowStrip(float width, unsigned inc);
 
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
    
    // dual pass is used to draw transparent spheres:
    void dualPassSphere1();
    void dualPassSphere2();
    void dualPassSphere4();
    void dualPassSphere8();

    /// draw a very nice half-sphere of radius 1 in Z < 0
    void hemisphere1();
    /// draw a nice half-sphere of radius 1 in Z < 0
    void hemisphere2();
    /// draw a flattened dome of radius 1 in XY with Z in [-0.5, 0]
    void dome();
    /// draw a very flat dome of radius 1 in XY with Z in [-0.25, 0]
    void roof();
    /// draw a blob with a pointy ends up in Z
    void pin();
    /// draw a sphere stretched in Z
    void droplet();

    /// draw nicest sphere available
    inline void sphere() { sphere1(); }
    /// draw nicest hemisphere available
    inline void hemisphere() { hemisphere1(); }

    /// draw a sphere decorated with 12 black pentagons
    void football();
    /// draw a fine sphere decorated with 12 pentagons
    void football1(gym_color);
    /// draw a fine sphere decorated with 12 outlined pentagons
    void footballT(gym_color);
    /// draw a nice sphere decorated with 12 pentagons
    void football2(gym_color);
    /// draw a sphere decorated with 12 pentagons
    void football4(gym_color);
    /// draw a coarse sphere decorated with 12 pentagons
    void football8(gym_color);

    /// draw a line on the sphere
    void baseballSeamCurve(float R, float W);
    /// draw a line on the sphere
    void tennisballSeamCurve(float R, float W);
    /// draw a white baseball
    void baseball();
    /// draw a yellow tennisball
    void tennisball();
    
    /// draw multiple lines forming a tube
    void stripes(float line_width);
    /// primitive used to draw the central segments of fibers
    inline void innerTube() { longTube2(); }
    /// primitive used to draw the minus ends of fibers
    inline void capedTube() { halfTube2(); dome(); }
    //inline void capedTube() { halfTube2(); discBottom2(); }
    /// primitive used to draw the plus ends of fibers
    inline void endedTube() { shutLongTube2(); }
    //inline void endedTube() { halfTube2(); hemisphere(); }

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
