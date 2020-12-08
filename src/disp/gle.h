// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University

#ifndef GLE_H
#define GLE_H

#include "real.h"
#include "opengl.h"
#include "gle_color.h"
#include "vector.h"


/// Simple geometrical objects drawn with OpenGL
/**
 @todo namespace gle should be a class -> we use GL.vertex(v)
 Problem: the gle prefix is already used by standard OpenGL elements
 */
namespace gle
{
    /// this defines the number of triangles used to draw shapes
    /** higher finesse improves the rendering: 8 is good, 16 is nice and 32 is very nice */
    constexpr size_t finesse = 12;
    
    /// number of circle points stored in buffer
    constexpr size_t ncircle = finesse * 8;

    /// initialize the arrays
    void initialize();
    
    /// release requested memory
    void release();
    
    /// initialize the Vertex Buffer Objects
    void initTubeBuffers();

    /// initialize the Vertex Buffer Objects
    void initSphereBuffers();
   
    /// initialize more buffer objects
    void initBuffers();
    
    /// calculate sinus and cosinus for a circle
    void circle(size_t cnt, GLfloat C[], GLfloat S[], double rad, double start = 0);
    
    /// calculate sinus and cosinus
    void circle(size_t cnt, GLfloat CS[], double rad, double start = 0);

    /// calculate sinus and cosinus for a circular arc
    void arc(size_t cnt, GLfloat C[], GLfloat S[], double rad, double start, double end, GLfloat cenx, GLfloat ceny);

#pragma mark -
    
    inline void scale(float x)                          { glScalef(x,x,x); }
    inline void scale(double x)                         { glScaled(x,x,x); }

    inline void scale(float x, float y, float z)        { glScalef(x,y,z); }
    inline void scale(double x, double y, double z)     { glScaled(x,y,z); }
   
    inline void translate(float x, float y, float z)    { glTranslatef(x, y, z); }
    inline void translate(double x, double y, double z) { glTranslated(x, y, z); }

#if REAL_IS_DOUBLE
    
    inline void gleVertex(Vector1 const& v)       { glVertex2d(v.XX, 0); }
    inline void gleVertex(Vector2 const& v)       { glVertex2d(v.XX, v.YY); }
    inline void gleVertex(Vector3 const& v)       { glVertex3d(v.XX, v.YY, v.ZZ); }
    inline void gleVertex(real x, real y)         { glVertex2d(x, y); }
    inline void gleVertex(real x, real y, real z) { glVertex3d(x, y, z); }
    
    inline void gleNormal(Vector1 const& v)       { glNormal3d(v.XX, 0, 0); }
    inline void gleNormal(Vector2 const& v)       { glNormal3d(v.XX, v.YY, 0); }
    inline void gleNormal(Vector3 const& v)       { glNormal3d(v.XX, v.YY, v.ZZ); }
    
    inline void translate(Vector1 const& v)       { glTranslated(v.XX, 0, 0); }
    inline void translate(Vector2 const& v)       { glTranslated(v.XX, v.YY, 0); }
    inline void translate(Vector3 const& v)       { glTranslated(v.XX, v.YY, v.ZZ); }

    inline void rasterPos(Vector1 const& v)       { glRasterPos2d(v.XX, 0); }
    inline void rasterPos(Vector2 const& v)       { glRasterPos2d(v.XX, v.YY); }
    inline void rasterPos(Vector3 const& v)       { glRasterPos3d(v.XX, v.YY, v.ZZ); }

#else

    inline void gleVertex(Vector1 const& v)       { glVertex2f(v.XX, 0); }
    inline void gleVertex(Vector2 const& v)       { glVertex2f(v.XX, v.YY); }
    inline void gleVertex(Vector3 const& v)       { glVertex3f(v.XX, v.YY, v.ZZ); }
    inline void gleVertex(real x, real y)         { glVertex2f(x, y); }
    inline void gleVertex(real x, real y, real z) { glVertex3f(x, y, z); }

    inline void gleNormal(Vector1 const& v)       { glNormal3f(v.XX, 0, 0); }
    inline void gleNormal(Vector2 const& v)       { glNormal3f(v.XX, v.YY, 0); }
    inline void gleNormal(Vector3 const& v)       { glNormal3f(v.XX, v.YY, v.ZZ); }
    
    inline void translate(Vector1 const& v)       { glTranslatef(v.XX, 0, 0); }
    inline void translate(Vector2 const& v)       { glTranslatef(v.XX, v.YY, 0); }
    inline void translate(Vector3 const& v)       { glTranslatef(v.XX, v.YY, v.ZZ); }

    inline void rasterPos(Vector1 const& v)       { glRasterPos2f(v.XX, 0); }
    inline void rasterPos(Vector2 const& v)       { glRasterPos2f(v.XX, v.YY); }
    inline void rasterPos(Vector3 const& v)       { glRasterPos3f(v.XX, v.YY, v.ZZ); }
 
#endif

    // colors that vary with the direction of a vector:
    inline gle_color radial_color(const Vector3& d) { return gle_color::radial_color((GLfloat)d.XX, (GLfloat)d.YY, (GLfloat)d.ZZ, 1.0f); }
    inline gle_color radial_color(const Vector2& d) { return gle_color::radial_color((GLfloat)d.XX, (GLfloat)d.YY, 1.0f); }
    inline gle_color radial_color(const Vector1& d) { if ( d.XX > 0 ) return gle_color(1,1,1); else return gle_color(0,1,0); }

    //------------------------------------------------------------------------------
#pragma mark -
        
    /// translate by T, then rotate to align X with A, Y with B and Z with C
    void transRotate(Vector3 const& T, Vector3 const& A, Vector3 const& B, Vector3 const& C);

    /// translate by A; rotate to align Z with AB, Z replaces X. The X-Y plane is scaled by R
    void stretchAlignZ(Vector2 const& A, Vector2 const& B, float rad);
    /// translate by T, then rotate to align Z with dir
    void stretchAlignZ(Vector3 const& A, Vector3 const& B, float rad);
    
    /// translate by T, then rotate to align Z with dir, scaling X and Y by radis
    void transAlignZ(Vector1 const& pos, float rad, Vector1 const& dir);
    void transAlignZ(Vector2 const& pos, float rad, Vector2 const& dir);
    void transAlignZ(Vector3 const& pos, float rad, Vector3 const& dir);

    void setClipPlane(GLenum, Vector1 const& dir, Vector1 const& pos);
    void setClipPlane(GLenum, Vector2 const& dir, Vector2 const& pos);
    void setClipPlane(GLenum, Vector3 const& dir, Vector3 const& pos);

    //------------------------------------------------------------------------------
#pragma mark -

    /// draw 2D circle of radius 1 in XY plane, with +Z as normal
    void circle();
    /// draw 2D disc of radius 1 in XY plane, with +Z as normal
    void discUp();
    void discDown();
    /// draw nicer 2D disc of radius 1 in XY plane, with +Z as normal
    void disc2();

    /// draw a tetrahedron of side 2 in 3D
    void tetrahedron();
    /// draw a octahedron of radius 1
    void octahedron();
    /// draw a icosahedron of radius 1
    void icosahedron();
    /// returns tetrahedron or octahedron
    inline void (*hedron(bool x))() { return x ? octahedron : tetrahedron; }
    
    /// draw a Cube of side 2
    void cube();
    /// display 3 arrow fins aligned with the Z axis, or radius 1, lenth 2, Z=[-0.5, 1.5]
    void arrowTail();

    /// draw a sphere of radius 1 at origin
    void sphere1U();
    /// draw a nice sphere of radius 1 at origin
    void sphere2U();
    /// draw a very nice sphere of radius 1 at origin
    void sphere4U();
    /// draw a very nice sphere of radius 1 at origin
    void sphere8U();
    /// draw Torus of radius `rad` and thickness `thick`
    void torusZ(GLfloat rad, GLfloat thick, size_t inc = 1);
    
    /// draw an open tube from B to T along Z, of diameter 1
    void tubeZ(GLfloat B, GLfloat T, int inc);
    /// draw an open tube from B to T along Z, of diameter 1
    void tubeZ(GLfloat B, GLfloat RB, GLfloat T, GLfloat RT, int inc);
    /// draw an open tube along Z, of diameter 1 and length 1, Z=[0, 1]
    void hexTubeZ(GLfloat Zmin, GLfloat Zmax);
    
    /// draw an open tube along Z, of diameter 1 and length 1
    void tube1();
    /// draw an open tube along Z, of diameter 1 and length 1
    void tube2();
    /// draw a nice open tube along Z, of diameter 1 and length 1
    void tube4();
    /// draw a nicer open tube along Z, of diameter 1 and length 1
    void tube8();
    /// draw a tube along Z, of diameter 1 and length 1.5, Z=[-0.25, 1.25]
    void longTube1();
    /// draw a nicer tube along Z, of diameter 1 and length 1.5, Z=[-0.25, 1.25]
    void longTube2();
    /// draw a nicer tube along Z, of diameter 1 and length 1.5, Z=[-0.25, 1.25]
    void longTube4();
    /// draw a tube along Z, of diameter 1 with Z=[-16, 0]
    void halfTube1();
    /// draw a nicer tube along Z, of diameter 1 with Z=[-16, 0]
    void halfTube2();
    /// draw a nicer tube along Z, of diameter 1 with Z=[-16, 0]
    void halfTube4();
    /// draw a nicer tube along Z, of diameter 1 with Z=[-16, 0]
    void hexTube();

    
    /// draw a closed tube along Z, or diameter 1 and length 1
    void cylinder();
    /// spherocylinder of length L, radius R, centered and aligned with axis Z
    void capsuleZ(GLfloat L, GLfloat R);

    /// draw a 3-portion cylinder with a larger central section
    void barrel();
    /// display a cone directed along Z, of radius R at Z=B, and 0 at Z=T
    void coneZ(GLfloat R, GLfloat B, GLfloat T, bool closed);
    /// display an open cone directed along Z, of radius 1 at Z=0
    inline void cone() { tubeZ(0, 1, 1, 0.25, 4); }
    /// display a closed cone directed along Z, of radius 1 in Z=[-1, +2]
    inline void longCone() { coneZ(1, -1, 2, true); }
    /// display a cylindrical box, directed along Z, of length 1, radius 1 in Z=[-0.5, +0.5]
    void cylinderZ();
    /// display a dumbbell aligned with the Z axis, or radius 1/3, lenth 1
    void dumbbell();

    /// draw a circular band composed of little triangles
    void drawArrowedBand(size_t nb_triangles, float width);
    /// draw 3 Arrowed Bands defining 8 quadrants on the sphere of radius 1
    void drawThreeBands(size_t nb_triangles);
    
    /// a rectangle ( rect = [ left, bottom, right, top ] )
    void drawRectangle(const int rect[4]);
    
    /// a rectangle with cut corners
    void drawNiceRectangle(const int rect[4], int);

    //------------------------------------------------------------------------------
    
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
    /// primitives to draw the ends of spherocylinders:
    inline void capedTube1() { halfTube1(); hemisphere1(); }
    inline void capedTube2() { halfTube2(); hemisphere2(); }
    inline void capedTube4() { halfTube4(); hemisphere4(); }
#else
    /// primitives to draw the ends of spherocylinders:
    inline void capedTube1() { halfTube1(); discDown(); }
    inline void capedTube2() { halfTube2(); discDown(); }
    inline void capedTube4() { halfTube4(); discDown(); }
#endif
    /// primitives to draw the ends of spherocylinders:
    inline void endedTube1() { halfTube1(); discDown(); }
    inline void endedTube2() { halfTube2(); discDown(); }
    inline void endedTube4() { halfTube4(); discDown(); }
    
    inline void arrowHead() { tube1(); glTranslatef(0, 0, 1); glScalef(3.f, 3.f, 3.f); gle::longCone(); }

    //------------------------------------------------------------------------------
#pragma mark -
    
    /// display back first, and then front
    void dualPass(void primitive());
 
    /// draw 'obj' scaled by radius at position 'x', oriented along 'd'
    void gleObject(Vector1 const& pos, Vector1 const& dir, float rad, void (*obj)());
    void gleObject(Vector2 const& pos, Vector2 const& dir, float rad, void (*obj)());
    void gleObject(Vector3 const& pos, Vector3 const& dir, float rad, void (*obj)());
    
    //------------------------------------------------------------------------------

    /// draw 'obj' with its ends at [a,b], of specified radius
    void gleTube(Vector1 const& A, Vector1 const& B, float rad, void (*obj)());
    void gleTube(Vector2 const& A, Vector2 const& B, float rad, void (*obj)());
    void gleTube(Vector3 const& A, Vector3 const& B, float rad, void (*obj)());

    void drawTube(Vector1 const& A, float rad, Vector1 const& B, void (*obj)());
    void drawTube(Vector2 const& A, float rad, Vector2 const& B, void (*obj)());
    void drawTube(Vector3 const& A, float rad, Vector3 const& B, void (*obj)());
    
    /// draw a band from A to B, with specified radius
    void drawBand(Vector2 const& A, Vector2 const& B, real);
    void drawBand(Vector3 const& A, Vector3 const& B, real);

    /// draw a band from A to B, with specified radius in A and B
    void drawBand(Vector1 const& a, real, Vector1 const& b, real);
    void drawBand(Vector2 const& a, real, Vector2 const& b, real);
    
    /// draw a band from A to B, with specified radius and colors in A and B
    void drawBand(Vector1 const& a, real, gle_color, Vector1 const& b, real, gle_color);
    void drawBand(Vector2 const& a, real, gle_color, Vector2 const& b, real, gle_color);

    /// draw symbol linking A to B
    void drawHourglass(Vector2 const& a, Vector2 const&, Vector2 const& b, Vector2 const&);
    void drawHourglass(Vector2 const&, Vector2 const&, gle_color,
                Vector2 const&, Vector2 const&, gle_color);
    /// draw symbol linking A to B
    void drawCross(Vector2 const& a, Vector2 const&, Vector2 const& b, Vector2 const&, real);
    void drawBar(Vector3 const& a, Vector3 const& da, Vector3 const& b, Vector3 const& db, real);
    
    /// draw two discs in A and B, connected with a line
    void drawDumbbell(Vector2 const& A, Vector2 const& B, float diameter);
    /// draw two spheres in A and B, connected with a cylinder
    void drawDumbbell(Vector3 const& A, Vector3 const& B, float diameter);

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
    void bitmapText(const char text[], void* font = nullptr, GLfloat vshift = 0);
    
    /// draw `text` at position `pos`
    void drawText(Vector1 const& pos, const char text[], void* font, float dx=0);
    /// draw `text` at position `pos`
    void drawText(Vector2 const& pos, const char text[], void* font, float dx=0);
    /// draw `text` at position `pos`
    void drawText(Vector3 const& pos, const char text[], void* font, float dx=0);
                        
    //------------------------------------------------------------------------------
#pragma mark -
    
    /// draw pixel array `rgba` containing 4 bytes per pixels
    void drawPixels(int width, int height, int nbc, GLubyte rgba[], Vector2 pos, Vector2 dx, Vector2 dy);
    
    /// display rectangle specified in pixel-coordinates
    void drawRectangle(const int rect[4], int window_width, int window_height);
    
    /// draw a rectangle to indicate the GLUT window-resize handle
    void drawResizeBox(int window_width, int window_height);
    
    /// draw plane with squares of alternating colors
    void drawTiledFloor(int R, float T, float Z, gle_color col1, gle_color col2);

    /// draw a set of 2 or 3 axes, depending on `dim`
    void drawAxes(real size, int dim);
    
    /// convert OpenGL error code to string
    const char* errorString(GLenum code);

    /// check and print OpenGL error(s)
    void reportErrors(FILE*, const char* msg);
 
    /// print some info for debugging purpose
    void dump();
}

//#define CHECK_GL_ERROR(ARG) ((void) 0)
#define CHECK_GL_ERROR(ARG) gle::reportErrors(stderr, ARG)



#endif
