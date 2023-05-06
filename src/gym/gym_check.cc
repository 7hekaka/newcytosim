
#include "gym_check.h"
#include "gym_color.h"
#include "opengl.h"


char const* gym::errorString(unsigned code)
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
 This check for OpenGL errors,
 the argument 'msg' can be useful to provide feedback for debugging
 */
void gym::reportErrors(FILE * out, const char* msg)
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
    fprintf(f, "%s( %4.2f %4.2f %4.2f %4.2f )", str, rgb[0], rgb[1], rgb[2], rgb[3]);
}

static void print_cap(GLenum cap, const char * str)
{
    GLint i = glIsEnabled(cap);
    fprintf(stderr, "%s %i ", str, i);
}

static void print_get(GLenum cap, const char * str)
{
    GLint i[4];
    glGetIntegerv(cap, i);
    fprintf(stderr, "%s %i ", str, i[0]);
}

void gym::printCaps(const char str[])
{
    fprintf(stderr, "%s :", str);
    print_cap(GL_NORMALIZE, "normalize");
    print_cap(GL_BLEND, "blend");
    print_cap(GL_ALPHA_TEST, "alpha");
    print_cap(GL_DEPTH_TEST, "depth");
    print_get(GL_DEPTH_WRITEMASK, "depthmask");
    print_cap(GL_LIGHTING, "light");
    print_cap(GL_LIGHT0, "");
    print_cap(GL_LIGHT1, "");
    print_cap(GL_LIGHT2, "");
    print_cap(GL_FOG, "fog");
    print_cap(GL_STENCIL_TEST, "stencil");
    print_cap(GL_CULL_FACE, "cull");
    print_cap(GL_SCISSOR_TEST, "scissor");
    print_cap(GL_COLOR_LOGIC_OP, "logic");
    print_cap(GL_COLOR_ARRAY, "array");
    print_cap(GL_COLOR_MATERIAL, "material");
    print_cap(GL_LINE_STIPPLE, "stipple");
    print_cap(GL_TEXTURE_2D, "texture");
    fprintf(stderr, "\n");
}

void gym::printMatrices(FILE * f)
{
    GLint V[4] = { 0 };
    glGetIntegerv(GL_VIEWPORT, V);
    fprintf(f, "viewport ( %4i %4i %4i %4i )\n", V[0], V[1], V[2], V[3]);
    GLfloat M[16] = { 0 };
    glGetFloatv(GL_PROJECTION_MATRIX, M);
    for ( int i = 0; i < 4; ++i )
        fprintf(f, "projection ( %8.2f %8.2f %8.2f %8.2f )\n", M[i], M[4+i], M[8+i], M[12+i]);
    glGetFloatv(GL_MODELVIEW_MATRIX, M);
    for ( int i = 0; i < 4; ++i )
        fprintf(f, "model_view ( %8.2f %8.2f %8.2f %8.2f )\n", M[i], M[4+i], M[8+i], M[12+i]);
}

    
/// print current color properties of OpenGL context
void gym::printColors(FILE * out)
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


void gym::printDepthRange(int X, int Y, int W, int H)
{
    size_t S = W * H;
    float * pixels = new float[S];
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glReadPixels(X, Y, W, H, GL_DEPTH_COMPONENT, GL_FLOAT, pixels);
    
    if ( glGetError() == GL_NO_ERROR )
    {
        float i = 2, s = -1;
        for ( size_t u = 0; u < S; ++u )
        {
            i = std::min(i, pixels[u]);
            if ( pixels[u] < 1 ) s = std::max(s, pixels[u]);
        }
        printf("OpenGL depth buffer range [ %8.6f %8.6f ]\n", i, s);
    }
    delete[] pixels;
}
