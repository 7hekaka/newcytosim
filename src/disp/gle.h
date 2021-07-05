// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University

#ifndef GLE_H
#define GLE_H

#include "real.h"
#include "opengl.h"
#include "gle_color.h"
#include "vector.h"
#include "flute.h"

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
    constexpr size_t pi_twice = finesse * 12;
    constexpr size_t pi_once = finesse * 6;
    constexpr size_t pi_half = finesse * 3;
    
    /// values of cosine, sine over a full circle
    extern const float* circle_;

    /// access to precomputed cosine
    inline float cos_(size_t n) { return circle_[4+2*n]; }
    
    /// access to precomputed sine
    inline float sin_(size_t n) { return circle_[5+2*n]; }

    /// initialize the arrays
    void initialize();
    
    /// release requested memory
    void release();

    /// calculate sine and cosine
    void compute_circle(size_t cnt, float CS[], double rad, double start = 0);

    /// calculate sine and cosine for a circular arc
    void compute_arc(size_t cnt, float CS[], double rad, double start, double angle, float cx, float cy);
    
    /// inverse square root
#if defined(__SSE3__)
    inline float invsqrt(float x) { return _mm_cvtss_f32(_mm_rsqrt_ss(_mm_set_ss(x))); }
#else
    inline float invsqrt(float x) { return x / sqrtf(x); }
#endif
   
    /// OpenGL buffers objects for streaming
    GLuint nextStream();
    
    /// initialize the Vertex Buffer Objects
    flute6* setTubeBuffers(flute6*, flute6* const);
    
    /// initialize more buffer objects
    flute6* setCubeBuffers(flute6*, flute6* const);

    /// initialize more buffer objects
    flute3* setBlobBuffers(flute3*, flute3* const);

    /// initialize more buffer objects
    void setBuffers();
    
    void bindBuffer();
    
    void unbindBuffer();
    
#pragma mark -
    
    inline void scale(float x)  { glScalef(x,x,x); }
    inline void scale(double x) { glScaled(x,x,x); }
   
    inline void translate(float x, float y, float z)    { glTranslatef(x, y, z); }
    inline void translate(double x, double y, double z) { glTranslated(x, y, z); }

    inline void transScale(float x, float y, float z, float s) { glTranslatef(x, y, z); glScalef(s,s,s); }

#if REAL_IS_DOUBLE
   
    inline void translate(Vector1 const& v) { glTranslated(v.XX, 0, 0); }
    inline void translate(Vector2 const& v) { glTranslated(v.XX, v.YY, 0); }
    inline void translate(Vector3 const& v) { glTranslated(v.XX, v.YY, v.ZZ); }

    inline void transScale(Vector1 const& v, double s) { glTranslated(v.XX, 0, 0); glScaled(s,s,s); }
    inline void transScale(Vector2 const& v, double s) { glTranslated(v.XX, v.YY, 0); glScaled(s,s,s); }
    inline void transScale(Vector3 const& v, double s) { glTranslated(v.XX, v.YY, v.ZZ); glScaled(s,s,s); }

#else
    
    inline void translate(Vector1 const& v) { glTranslatef(v.XX, 0, 0); }
    inline void translate(Vector2 const& v) { glTranslatef(v.XX, v.YY, 0); }
    inline void translate(Vector3 const& v) { glTranslatef(v.XX, v.YY, v.ZZ); }
    
    inline void transScale(Vector1 const& v, float s) { glTranslatef(v.XX, 0, 0); glScalef(s,s,s); }
    inline void transScale(Vector2 const& v, float s) { glTranslatef(v.XX, v.YY, 0); glScalef(s,s,s); }
    inline void transScale(Vector3 const& v, float s) { glTranslatef(v.XX, v.YY, v.ZZ); glScalef(s,s,s); }

#endif

    // colors that vary with the direction of a vector:
    inline gle_color radial_color(const Vector3& d) { return gle_color::radial_color((GLfloat)d.XX, (GLfloat)d.YY, (GLfloat)d.ZZ, 1.0f); }
    inline gle_color radial_color(const Vector2& d) { return gle_color::radial_color((GLfloat)d.XX, (GLfloat)d.YY, 1.0f); }
    inline gle_color radial_color(const Vector1& d) { if ( d.XX > 0 ) return gle_color(1,1,1); else return gle_color(0,1,0); }

    //------------------------------------------------------------------------------
#pragma mark -
    
    /// translate by T; rotate to align X with A, Y with B and Z with C
    void transRotate(Vector3 const& T, Vector3 const& A, Vector3 const& B, Vector3 const& C);

    /// translate by A; rotate to align Z with AB, Z replacing X. Scale in Z to put B at (0,0,1) Scale XY plane by `rad'
    void stretchAlignZ(Vector1 const& A, Vector1 const& B, float rad);
    /// translate by A; rotate to align Z with AB, Z replacing X. Scale in Z to put B at (0,0,1) Scale XY plane by `rad'
    void stretchAlignZ(Vector2 const& A, Vector2 const& B, float rad);
    /// translate by A; rotate to align Z with AB, Z replacing X. Scale XY plane by `rad'
    void stretchAlignZ(Vector3 const& A, Vector3 const& B, float rad);
    
    /// translate by pos; rotate to align Z with dir, scale XY plane by rad
    void transAlignZ(Vector1 const& pos, float rad, Vector1 const& dir);
    void transAlignZ(Vector2 const& pos, float rad, Vector2 const& dir);
    void transAlignZ(Vector3 const& pos, float rad, Vector3 const& dir);

    /// translate by pos; rotate to align Z with dir, given norm(dir)=1, scale XY by rad and Z by fac
    void stretchAlignZ1(Vector1 const& pos, float rad, Vector1 const& dir, float fac);
    void stretchAlignZ1(Vector2 const& pos, float rad, Vector2 const& dir, float fac);
    void stretchAlignZ1(Vector3 const& pos, float rad, Vector3 const& dir, float fac);

    /// translate by pos; rotate to align X to Z, scale by rad
    void transAlignZX(float pos, float rad, float dir);
    /// translate by A; rotate to align X to Z, scale XY by rad and Z by B-A
    void stretchAlignZX(float A, float B, float rad);

    void setClipPlane(GLenum, Vector1 const& dir, Vector1 const& pos);
    void setClipPlane(GLenum, Vector2 const& dir, Vector2 const& pos);
    void setClipPlane(GLenum, Vector3 const& dir, Vector3 const& pos);
    
    /// display back faces first, and then front faces
    void dualPass(void primitive());

    //------------------------------------------------------------------------------
#pragma mark -

    /// draw 2D circle of radius 1 in XY plane, with +Z as normal
    void circle();
    /// draw 2D circle of radius 1 in XY plane, with +Z as normal, fewer points
    void circle2();
    /// draw 2D circle of radius 1 in XY plane, with +Z as normal, dotted
    void circle_dotted();
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

    /// draw a tetrahedron of side 2 in 3D
    void tetrahedron();
    /// draw a octahedron of radius 1
    void octahedron();
    /// draw a icosahedron of radius 1
    void icosahedron();
    /// returns tetrahedron or octahedron
    inline void (*hedron(bool x))() { return x ? octahedron : tetrahedron; }

    /// draw a roughly spherical shape made of few triangles
    void blob();
    /// draw a centered blob of radius 1 with a cone extending up in Z
    void needle();
    /// draw a Cube of side 2
    void cube();
    /// draw a Cube of side 1
    void smallCube();
    /// draw a stellated octahedron
    void star();
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
    /// display a dumbbell aligned with the Z axis, or radius 1/3, lenth 1
    void dumbbell();
    /// draw Torus of radius `rad` and thickness `thick`
    void torusZ(float rad, float thick, size_t inc = 1);

    /// draw a circular band composed of little triangles
    void arrowStrip(float width, size_t inc);
    /// draw 3 Arrowed Bands defining 8 quadrants on the sphere of radius 1
    void threeArrowStrip(float width, size_t inc);
    
    /// a rectangle ( rect = [ left, bottom, right, top ] )
    void drawRectangle(const int rect[4]);
    void drawRectangle(float L, float B, float R, float T, float Z=0);

    //------------------------------------------------------------------------------
    
    /// draw something
    void thing();
    /// do not draw
    void nothing();
    
    /// draw a sphere at the origin of radius 1
    void sphere1();
    /// draw a nice sphere at the origin of radius 1
    void sphere2();
    /// draw a very nice sphere at the origin of radius 1
    void sphere4();
    /// draw a refined sphere at the origin of radius 1
    void sphere8();
    
    void dualPassSphere1();
    void dualPassSphere2();
    void dualPassSphere4();
    void dualPassSphere8();

    /// draw half a sphere in Z < 0
    void hemisphere1();
    /// draw a nice half-sphere in Z < 0
    void hemisphere2();
    /// draw a very nice half-sphere in Z < 0
    void hemisphere4();

#if 1
    /// primitives to draw the minus ends of fibers:
    inline void capedTube1() { halfTube1(); hemisphere1(); }
    inline void capedTube2() { halfTube2(); hemisphere2(); }
    inline void capedTube4() { halfTube4(); hemisphere4(); }
#else
    /// primitives to draw the minus ends of fibers:
    inline void capedTube1() { halfTube1(); discBottom1(); }
    inline void capedTube2() { halfTube2(); discBottom2(); }
    inline void capedTube4() { halfTube4(); discBottom2(); }
#endif
    /// primitives to draw the plus ends of fibers:
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
    void drawBand(Vector1 const& a, GLfloat, Vector1 const& b, GLfloat);
    void drawBand(Vector2 const& a, GLfloat, Vector2 const& b, GLfloat);
    
    /// draw a band from A to B, with specified radius and colors in A and B
    void drawBand(Vector1 const& a, GLfloat, gle_color, Vector1 const& b, GLfloat, gle_color);
    void drawBand(Vector2 const& a, GLfloat, gle_color, Vector2 const& b, GLfloat, gle_color);

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

    //------------------------------------------------------------------------------
#pragma mark -

    /// return height in pixel of GLUT font
    int fontHeight(void* font);

    /// compute size of text
    int maxTextWidth(const char text[], void* font, int& lines);
    
    /// display text on a rectangle of color `bcol`, in a corner of the center of the display window
    void drawText(const char text[], void* font, gle_color bcol, int position, int width, int height);
    
    /// draw text at the current OpenGL raster position and raster color
    void bitmapString(const char text[], void* font = nullptr, GLfloat vshift = 0);
    
    /// draw `text` at position `pos`
    void drawText(Vector1 const& pos, const char text[], void* font, float dx=0);
    /// draw `text` at position `pos`
    void drawText(Vector2 const& pos, const char text[], void* font, float dx=0);
    /// draw `text` at position `pos`
    void drawText(Vector3 const& pos, const char text[], void* font, float dx=0);
                        
    //------------------------------------------------------------------------------
#pragma mark -
    
    /// display rectangle specified in pixel-coordinates
    void drawRectangle(const int rect[4], int window_width, int window_height);
    
    /// draw a rectangle to indicate the GLUT window-resize handle
    void drawResizeBox(int window_width, int window_height);
    
    /// draw plane with squares of alternating colors
    void drawTiledFloor(int R, float T, float Z, gle_color col1, gle_color col2);

    /// draw a set of 2 or 3 axes, depending on `dim`
    void drawAxes(float size, int dim);

#pragma mark -

    /// convert OpenGL error code to string
    const char* errorString(GLenum code);

    /// check and print OpenGL error(s)
    void reportErrors(FILE*, const char* msg);
 
    /// print some info for debugging purpose
    void dump();
    
    /// print current color properties of OpenGL context
    void print_color_materials(std::ostream& os);
}

#ifdef NDEBUG
#  define CHECK_GL_ERROR(ARG) ((void) 0)
#else
#  define CHECK_GL_ERROR(ARG) gle::reportErrors(stderr, ARG)
#endif


#endif
