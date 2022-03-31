
#include "gym_check.h"
#include "gle_color.h"
#include "opengl.h"

namespace gym
{

char const* errorString(unsigned code)
{
    switch ( code )
    {
        case GL_NO_ERROR:          return "GL_NO_ERROR";
        case GL_INVALID_ENUM:      return "GL_INVALID_ENUM";
        case GL_INVALID_VALUE:     return "GL_INVALID_VALUE";
        case GL_INVALID_OPERATION: return "GL_INVALID_OPERATION";
        case GL_STACK_OVERFLOW:    return "GL_STACK_OVERFLOW";
        case GL_STACK_UNDERFLOW:   return "GL_STACK_UNDERFLOW";
        case GL_OUT_OF_MEMORY:     return "GL_OUT_OF_MEMORY";
        case GL_TABLE_TOO_LARGE:   return "GL_TABLE_TOO_LARGE";
        default:                   return "GL_UNKNOWN_ERROR";
    }
}

/**
 This is similar to glutReportError,
 but the additional argument can provide useful feedback for debugging
 */
void reportErrors(FILE * out, const char* msg)
{
    GLenum e = glGetError();
    while ( e != GL_NO_ERROR )
    {
        fprintf(out, "OpenGL error `%s' %s\n", errorString(e), msg);
        e = glGetError();
    }
}


//--------------------------------------------------------------------------

static void print_rgba(FILE* f, const char* str, GLfloat rgb[4])
{
    fprintf(f, "%s ( %4.2f %4.2f %4.2f %4.2f )", str, rgb[0], rgb[1], rgb[2], rgb[3]);
}

static void print_cap(GLenum cap, const char * str)
{
    GLint i = glIsEnabled(cap);
    fprintf(stderr, "%s %i ", str, i);
}

void print_caps()
{
    GLfloat c[4] = { 0 };
    glGetFloatv(GL_CURRENT_COLOR, c);
    print_rgba(stderr, "color", c);
    
    print_cap(GL_LIGHTING, "light");
    print_cap(GL_BLEND, "blend");
    print_cap(GL_FOG, "fog");
    print_cap(GL_DEPTH_TEST, "depth");
    print_cap(GL_ALPHA_TEST, "alpha");
    print_cap(GL_STENCIL_TEST, "stencil");
    print_cap(GL_CULL_FACE, "cull");
    print_cap(GL_COLOR_LOGIC_OP, "logic");
    print_cap(GL_COLOR_ARRAY, "array");
    print_cap(GL_COLOR_MATERIAL, "material");
    print_cap(GL_LINE_STIPPLE, "stipple");
    
    std::clog << '\n';
    
#if ( 0 )
    GLint vp[4] = { 0 };
    glGetIntegerv(GL_VIEWPORT, vp);
    fprintf(stderr, "viewport ( %4i %4i %4i %4i )", vp[0], vp[1], vp[2], vp[3]);
#endif
}

    
/// print current color properties of OpenGL context
void print_color_materials(FILE * out)
{
    GLfloat mat[4] = { 0 };
    glGetMaterialfv(GL_FRONT, GL_AMBIENT, mat);
    print_rgba(out, "front  amb", mat);
    glGetMaterialfv(GL_FRONT, GL_DIFFUSE, mat);
    print_rgba(out, "front  dif", mat);
    glGetMaterialfv(GL_FRONT, GL_EMISSION, mat);
    print_rgba(out, "front  emi", mat);

    glGetMaterialfv(GL_BACK, GL_AMBIENT, mat);
    print_rgba(out, "back  amb", mat);
    glGetMaterialfv(GL_BACK, GL_DIFFUSE, mat);
    print_rgba(out, "back  dif", mat);
    glGetMaterialfv(GL_BACK, GL_EMISSION, mat);
    print_rgba(out, "back  emi", mat);
}

}
