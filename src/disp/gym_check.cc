
#include "gym_check.h"

namespace gym
{
//--------------------------------------------------------------------------
#pragma mark - Debugging

char const* errorString(GLenum code)
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

void print_cap(GLenum cap, const char * str)
{
    GLint i = glIsEnabled(cap);
    std::clog << str << " " << i << "   ";
}

void dump_cap()
{
    GLfloat c[4] = { 0 };
    glGetFloatv(GL_CURRENT_COLOR, c);
    std::clog << "color = " << c[0] << " " << c[1] << " " << c[2] << " " << c[3] << '\n';
    
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
    std::clog << "viewport = " << vp[0] << " " << vp[1] << " " << vp[2] << " " << vp[3] << '\n';
#endif
}

/// print current color properties of OpenGL context
void print_color_materials(std::ostream& os)
{
    GLfloat mat[4] = { 0 };
    glGetMaterialfv(GL_FRONT, GL_AMBIENT, mat);
    os << "front  amb" << gle_color::components(mat) << '\n';
    glGetMaterialfv(GL_FRONT, GL_DIFFUSE, mat);
    os << "front  dif" << gle_color::components(mat) << '\n';
    glGetMaterialfv(GL_FRONT, GL_EMISSION, mat);
    os << "front  emi" << gle_color::components(mat) << '\n';
    os << '\n';

    glGetMaterialfv(GL_BACK, GL_AMBIENT, mat);
    os << "back  amb" << gle_color::components(mat) << '\n';
    glGetMaterialfv(GL_BACK, GL_DIFFUSE, mat);
    os << "back  dif" << gle_color::components(mat) << '\n';
    glGetMaterialfv(GL_BACK, GL_EMISSION, mat);
    os << "back  emi" << gle_color::components(mat) << '\n';
    os << '\n';
}

}
