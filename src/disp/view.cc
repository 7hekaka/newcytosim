// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "assert_macro.h"
#include "view.h"
#include "gle.h"
#include "vector_float.h"
#include "offscreen.h"
#include "glu_unproject.cc"
#include "glut.h"
#include "tictoc.h"

using namespace gle;

//------------------------------------------------------------------------------

View::View(const std::string& n)
: ViewProp(n)
{
    mWindowId = 0;
    displayCallback = nullptr;
    
    visRegion[0] = view_size;
    visRegion[1] = view_size;
    visRegion[2] = view_size;
    
    eyePosition[0] = 0;
    eyePosition[1] = 0;
    eyePosition[2] = -0.5f * view_size;
    
    hasROI = false;
    mROI[0].reset();
    mROI[1].reset();
    
    hasMatrices = false;
}


View::~View()
{
}


//------------------------------------------------------------------------------
#pragma mark -

void View::initGL()
{
    // let GL normalize the normals:
    glEnable(GL_NORMALIZE);
    glEnableClientState(GL_VERTEX_ARRAY);
    glDisable(GL_STENCIL_TEST);
    glDisable(GL_DITHER);
    
    if ( 1 )
    {
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    }
    else
    {
        //glDisable(GL_BLEND);
    }
    
    if ( multisample > 1 )
    {
        glEnable(GL_MULTISAMPLE);
        /*
         GLint s = 0;
         glGetIntegerv(GL_MAX_SAMPLES, &s);
         std::clog << "OpenGL samples = " << samples << "  max = " << s << '\n';
         */
    }
    else
    {
        glDisable(GL_MULTISAMPLE);
        if ( 1 )
        {
            glEnable(GL_POINT_SMOOTH);
            glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
            glEnable(GL_LINE_SMOOTH);
            glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
        }
        /*
         Do not enable POLYGON_SMOOTH, which destroys joints of triangulated surfaces
         glEnable(GL_POLYGON_SMOOTH);
         glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
         */
    }
    
    if ( depth_clamp )
        glEnable(GL_DEPTH_CLAMP);
    else
        glDisable(GL_DEPTH_CLAMP);

    if ( depth_test )
    {
        glEnable(GL_DEPTH_TEST);
        //glDepthFunc(GL_LESS);
        glDepthFunc(GL_LEQUAL);
        // enable Alpha Test to discard transparent pixels:
        glEnable(GL_ALPHA_TEST);
        glAlphaFunc(GL_GREATER, 0.05f);
    }
    else
    {
        //std::clog << "no depth-test" << '\n';
        glDisable(GL_DEPTH_TEST);
        glDisable(GL_ALPHA_TEST);
        glAlphaFunc(GL_ALWAYS, 0);
    }
}


void View::openDisplay()
{
    load();
    back_color.load_clear();
    setFog(fog_type, fog_param, fog_color);
    setLights();
    setClipping();
}


/**
 Displays Axes and Scalebar
 */
void View::closeDisplay() const
{
    endClipping();
    
    if ( axes )
        gle::drawAxes(axes_size, axes);
    
    if ( scalebar )
    {
        scalebar_color.load();
        drawScaleBar(scalebar, scalebar_length);
    }
    
    if ( label != "off" && label != "none" )
    {
        // display only first line of text:
        std::string msg = full_label.substr(0, full_label.find('\n'));
        front_color.load();
        drawText(msg, nullptr, 0x0, 0);
    }
}


/**
 add over-the-window features for the interactive display
*/
void View::drawInteractiveFeatures() const
{
    if ( hasROI )
        drawROI();

    if ( draw_memo && memo.size() )
    {
        glColor3f(1,1,1);
        drawText(memo, GLUT_BITMAP_8_BY_13, 0x000000CC, 4);
    }

    if ( top_message.size() )
    {
        front_color.load();
        drawText(top_message, nullptr, back_color.alpha(0.5), 3);
    }
    
    if ( label != "off" )
    {
        front_color.load();
        drawText(full_label, nullptr, 0x0, 0);
    }
    
#if ( 0 )
    // display FPS = frames per seconds
    static char buf[16];
    static size_t cnt = 0;
    static double sec = TicToc::seconds_today();
    ++cnt;
    double now = TicToc::seconds_today();
    if ( now > sec + 1.0 )
    {
        double fps = cnt / ( now - sec );
        snprintf(buf, sizeof(buf), "%6.2f fps", fps);
        sec = now;
        cnt = 0;
    }
    front_color.load();
    drawText(buf, nullptr, 0x0, 1);
#endif
}


/**
 Set two light sources
 */
void View::setLights(bool local) const
{
    glMatrixMode(GL_MODELVIEW);
    if ( local )
    {
        glPushMatrix();
        glLoadIdentity();
    }
    
    glShadeModel(GL_SMOOTH);
    
    GLfloat matWhite[]  = { 1.0f, 1.0f, 1.0f, 1.0f };
    GLfloat matGray[]   = { 0.2f, 0.2f, 0.2f, 1.0f };
    GLfloat matBlack[]  = { 0.0f, 0.0f, 0.0f, 1.0f };
    //GLfloat matBlue[]   = { 0.f, 0.f, 1.f, 1.f };
    
    glMaterialfv(GL_FRONT, GL_AMBIENT,   matBlack);
    glMaterialfv(GL_FRONT, GL_DIFFUSE,   matBlack);
    glMaterialfv(GL_FRONT, GL_SPECULAR,  matWhite);
    glMateriali (GL_FRONT, GL_SHININESS, 32);

    // set a gray color for the back-side of everything
    glMaterialfv(GL_BACK, GL_AMBIENT,  matGray);
    glMaterialfv(GL_BACK, GL_DIFFUSE,  matBlack);
    glMaterialfv(GL_BACK, GL_SPECULAR, matBlack);
    glMateriali (GL_BACK, GL_SHININESS, 8);
    
    GLfloat lightDiffuse[]  = { 0.8f, 0.8f, 0.8f, 1.0f };
    GLfloat lightSpecular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
    GLfloat lModelAmbient[] = { 0.4f, 0.4f, 0.4f, 1.0f };
    
    GLfloat light0[] = { 0.577f, -0.577f, 0.577f, 0.0f };
    glLightfv(GL_LIGHT0, GL_POSITION, light0);
    glLightfv(GL_LIGHT0, GL_DIFFUSE,  lightDiffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, lightSpecular);
    glEnable(GL_LIGHT0);
    
    GLfloat light1[] = {-0.7f, 0.0f, -0.7f, 0.0f };
    glLightfv(GL_LIGHT1, GL_POSITION, light1);
    glLightfv(GL_LIGHT1, GL_DIFFUSE,  lightDiffuse);
    glLightfv(GL_LIGHT1, GL_SPECULAR, lightSpecular);
    glEnable(GL_LIGHT1);
    /*
    GLfloat light2[] = { 0.0f, 0.0f, -1.0f, 0.0f };
    glLightfv(GL_LIGHT2, GL_POSITION, light2);
    glLightfv(GL_LIGHT2, GL_DIFFUSE,  lightDiffuse);
    glLightfv(GL_LIGHT2, GL_SPECULAR, lightSpecular);
    glEnable(GL_LIGHT2);
     */
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lModelAmbient);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
    glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_FALSE);

    if ( local )
        glPopMatrix();
}


//------------------------------------------------------------------------------
#pragma mark -


void View::reshape(int W, int H)
{
    //std::clog << "View::reshaped " << W << " " << H << '\n';
    window_size[0] = W;
    window_size[1] = H;
    glViewport(0, 0, W, H);
    load();
    /*
    if ( W > H )
    {
        GLfloat ratio = H / GLfloat(W);
        glOrtho(-1.0, 1.0, -ratio, ratio, 0, 1);
    }
    else
    {
        GLfloat ratio = W / GLfloat(H);
        glOrtho(-ratio, ratio, -1.0, 1.0, 0, 1);
    }
    */
}


void View::adjust()
{
    //calculate the visible region in the 3 directions:
    if ( window_size[0] > window_size[1] )
    {
        GLfloat R = GLfloat(window_size[1]) / GLfloat(window_size[0]);
        visRegion[0] = view_size;
        visRegion[1] = view_size * R;
    }
    else
    {
        GLfloat R = GLfloat(window_size[0]) / GLfloat(window_size[1]);
        visRegion[0] = view_size * R;
        visRegion[1] = view_size;
    }
    visRegion[2] = view_size;
    
    //std::clog << "View::adjust  " << visRegion[0] << " " << visRegion[1] << " " << visRegion[2] << "\n";
    //std::clog << " pixel_size = " << zoom * visRegion[0]/window_size[0] << '\n';
}


void View::setProjection()
{
    GLfloat X = visRegion[0] * 0.5f;
    GLfloat Y = visRegion[1] * 0.5f;
    GLfloat Z = visRegion[2];
    GLfloat S = view_size;

    //std::clog << "View::setProjection  " << mWindowId << "\n";
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    if ( perspective == 3 )
    {
        // this creates a stronger perspective:
        eyePosition[2] = -1.5f * S;
        glFrustum(-X, X, -Y, Y, Z, 5.0f*Z);
    }
    else if ( perspective == 2 )
    {
        // this creates a strong perspective:
        eyePosition[2] = -2.0f * S;
        glFrustum(-X, X, -Y, Y, Z, 6.0f*Z);
    }
    else if ( perspective )
    {
        // this creates a perspective:
        eyePosition[2] = -2.0f * S;
        glFrustum(-X, X, -Y, Y, Z, 11.0f*Z);
    }
    else
    {
        // The back-plane is set behind to avoid clipping
        eyePosition[2] = -0.5f * S;
        glOrtho(-X, X, -Y, Y, 0, Z);
    }
    
    glMatrixMode(GL_MODELVIEW);
}


void View::load()
{
    //std::clog << "View::load() win " << window() << "\n";
    adjust();
    setProjection();
    setModelView();
}


void View::setModelView() const
{
    //std::cerr<<"setModelView win " << window() << '\n';

    // setup the OpenGL transformation matrix
    glMatrixMode(GL_MODELVIEW);
#if ( 1 )
    GLfloat mat[16];
    rotation.setOpenGLMatrix(mat, eyePosition);
    glLoadMatrixf(mat);
#else
    glLoadIdentity();
    real axs[3];
    real ang = rotation.getAngle(axs) * 180 / M_PI;
    glTranslatef(eyePosition[0], eyePosition[1], eyePosition[2]);
    glRotatef(ang, axs[0], axs[1], axs[2]);
#endif
    gle::scale(zoom);
    
    // point-of-focus:
    gle::translate(-(focus+focus_shift));
    
#if ( 0 )
    std::clog << "View::setModelView eye      " << eyePosition << "\n";
    std::clog << "View::setModelView rotation " << rotation << "\n";
    std::clog << "View::setModelView zoom     " << zoom << "\n";
    std::clog << "View::setModelView focus    " << focus+focus_shift << "\n";
#endif
#if ( 0 )
    GLint vp[4] = { 0 };
    glGetIntegerv(GL_VIEWPORT, vp);
    std::clog << "viewport = " << vp[0] << " " << vp[1] << " " << vp[2] << " " << vp[3] << '\n';
#endif
}


/**
 This will change what is visible in the Z direction near and far
 from the observer, using clipping planes and fog.
 A a function of `mode`:
 - 0 : disabled
 - 1 : show ( Z > 0 ) with fog
 - 2 : show ( Z < 0 ) with fog
 - 3 : show slice ( -a < Z < a ) where a = 5% of view_size
 .
 
 */
void View::sliceView(int mode) const
{
    real off = view_size * 0.5;
    
    if ( mode == 1 )
    {
        setClipPlaneEye(GL_CLIP_PLANE1, Vector3(0,0,+1), off);
        if ( !depth_clamp )
            setFog(1, 0, fog_color);
    }
    else if ( mode == 2 )
    {
        setClipPlaneEye(GL_CLIP_PLANE1, Vector3(0,0,-1), -off);
        setFog(1, 1, fog_color);
    }
    else if ( mode == 3 )
    {
        real thk = view_size * 0.05;
        setClipPlaneEye(GL_CLIP_PLANE1, Vector3(0,0,-1), thk-off);
        setClipPlaneEye(GL_CLIP_PLANE2, Vector3(0,0,+1), thk+off);
    }
}


//------------------------------------------------------------------------------
#pragma mark -

void View::reset()
{
    zoom = 0.933033;
    auto_scale = 1;
    focus.reset();
    focus_shift.reset();
    rotation.set(1,0,0,0);
    setModelView();
}


void View::zoom_to(GLfloat z)
{
    //std::clog << "zoom_to " << z << " " << this << '\n';
    zoom = z;
    setModelView();
}


void View::move_to(const Vector3& d)
{
    focus = d;
    setModelView();
}


void View::move_shift(const Vector3& d)
{
    focus_shift = d;
    setModelView();
}


void View::rotate_to(const Quaternion<real> & q)
{
    rotation = normalize(q);
    setModelView();
}


/**
 This assumes that vector `dir` is normalized
 */
void View::align_with(const Vector3 & dir)
{
    // axis is obtained by vector product: axis = cross((1, 0, 0), a)
    Vector3 axis( 0, dir.ZZ, -dir.YY );
    // cosine is scalar product, sine is norm of vector-product:
    real cs = dir.XX, si = axis.norm();
    if ( si > REAL_EPSILON )
    {
        rotation.setFromAxis(axis, std::atan2(si, cs));
        setModelView();
    }
}

//------------------------------------------------------------------------------
#pragma mark -

void View::adjustROI(real Z)
{
    getMatrices();
    setROI(unproject(0, 0, Z), unproject(width(), height(), Z));
}


void View::matchROI()
{
    if ( hasROI )
    {
        focus = 0.5 * ( mROI[0] + mROI[1] );
        real r = 0.5 * ( mROI[0] - mROI[1] ).norm_inf();
        
        // zoom only if region is 7 pixels wide:
        if ( r > 7 * pixelSize() )
            zoom = view_size / (GLfloat)r;
        
        setModelView();
    }
}


/** Only check X and Y components */
bool View::insideROI(Vector3 pos) const
{
    bool inX = ((mROI[0].XX < pos.XX) & (pos.XX < mROI[1].XX));
    bool inY = ((mROI[0].YY < pos.YY) & (pos.YY < mROI[1].YY));
    return ( inX & inY );
}


void View::setROI(Vector3 a, Vector3 b)
{
    hasROI = true;
    mROI[0].set(std::min(a.XX, b.XX), std::min(a.YY, b.YY), std::min(a.ZZ, b.ZZ));
    mROI[1].set(std::max(a.XX, b.XX), std::max(a.YY, b.YY), std::max(a.ZZ, b.ZZ));
}


//------------------------------------------------------------------------------
#pragma mark -

void View::getMatrices()
{
    //get the transformation matrices, to be used for mouse control
    glGetIntegerv(GL_VIEWPORT,         mViewport);
    glGetDoublev(GL_PROJECTION_MATRIX, mProjection);
    glGetDoublev(GL_MODELVIEW_MATRIX,  mModelview);
    hasMatrices = true;
}

/**
 This set a matrix like glOrtho()
 */
void View::setOrthoMat(GLdouble * mat)
{
    GLdouble L =-0.5*visRegion[0];
    GLdouble R = 0.5*visRegion[0];
    GLdouble B =-0.5*visRegion[1];
    GLdouble T = 0.5*visRegion[1];
    GLdouble N = 0;
    GLdouble F = visRegion[2];

    for ( int i = 0; i < 16; ++i )
        mat[i] = 0;
    
    mat[ 0] =  2.0 / ( R - L );
    mat[ 5] =  2.0 / ( T - B );
    mat[10] = -2.0 / ( F - N );
    
    mat[12] = -( R + L ) / ( R - L );
    mat[13] = -( T + B ) / ( T - B );
    mat[14] = -( F + N ) / ( F - N );
    mat[15] = 1.0;
}


/**
 Transforms the given window coordinates into user coordinates.
 
 It uses the matrices obtained at the last call of getMatrices(),
 or the current matrices if get_matrices == true.
 
 For more info, try `man gluUnProject`
 */
Vector3 View::unproject(GLdouble x, GLdouble y, GLdouble z, bool get_matrices)
{
    GLdouble ux = 0, uy = 0, uz = 0;
    if ( get_matrices )
    {
        GLint      vp[4];
        GLdouble   mv[16];
        GLdouble   pj[16];
        
        glGetIntegerv(GL_VIEWPORT,         vp);
        glGetDoublev(GL_PROJECTION_MATRIX, pj);
        glGetDoublev(GL_MODELVIEW_MATRIX,  mv);

        setOrthoMat(pj);
        myUnproject(x, y, z, mv, pj, vp, &ux, &uy, &uz);
    }
    else if ( hasMatrices )
        myUnproject(x, y, z, mModelview, mProjection, mViewport, &ux, &uy, &uz);
    else
        std::cerr << "warning: View::unproject called without matrices\n";
    
    //printf("unproject( %.2f, %.2f, %.2f ) = ( %.2f, %.2f, %.2f )\n", x, y, z, ux, uy, uz);
    return Vector3(ux, uy, uz);
}


//------------------------------------------------------------------------------
#pragma mark -

void View::setFog(GLint type, GLfloat param, gle_color color) const
{
    GLint gl_type = 0;
    switch( type )
    {
        case 1: gl_type = GL_LINEAR; break;
        case 2: gl_type = GL_EXP;    break;
        case 3: gl_type = GL_EXP2;   break;
        default: glDisable(GL_FOG); return;
    }
   
    glEnable(GL_FOG);
    glFogi(GL_FOG_MODE, gl_type);
    
    if ( gl_type == GL_LINEAR )
    {
        glFogf(GL_FOG_START, param*visRegion[2]);
        glFogf(GL_FOG_END, (param*2+1)*visRegion[2]);
    }
    else
    {
        glFogf(GL_FOG_DENSITY, param/visRegion[2]);
    }
    
    glFogfv(GL_FOG_COLOR, color.colors());
}

void View::enableFog(const GLint type, const GLfloat param, gle_color color)
{
    fog_type = type;
    fog_param = param;
    fog_color = color;
}

void View::enableFog(const GLint type, const GLfloat param)
{
    fog_type = type;
    fog_param = param;
}

//------------------------------------------------------------------------------
#pragma mark -

/**
 The plane equations is relative to the model
 */
void View::setClipPlane(GLenum glp, Vector3 dir, real sca) const
{
    GLdouble eq[] = {dir.XX, dir.YY, dir.ZZ, sca};
    glClipPlane(glp, eq);
    glEnable(glp);
}

/**
 The plane equation is relative to the camera
 */
void View::setClipPlaneEye(GLenum glp, Vector3 dir, real sca) const
{
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    setClipPlane(glp, dir, sca);
    glPopMatrix();
}


void View::setClipping() const
{
    for ( int ix = 0; ix < NB_CLIP_PLANES; ++ix )
    {
        if ( clip_plane_mode[ix] == 1 )
            setClipPlane(GL_CLIP_PLANE0+ix, clip_plane_vector[ix], clip_plane_scalar[ix]);
        else if ( clip_plane_mode[ix] == 2 )
            setClipPlaneEye(GL_CLIP_PLANE0+ix, clip_plane_vector[ix], clip_plane_scalar[ix]);
        else
            glDisable(GL_CLIP_PLANE0+ix);
    }
    
    if ( slice )
        sliceView(slice);
}

/*
 Disable all clip planes, including the one set by sliceView()
 */
void View::endClipping() const
{
    for ( int ix = 0; ix < NB_CLIP_PLANES; ++ix )
        glDisable(GL_CLIP_PLANE0+ix);
}

void View::enableClipPlane(int ix, Vector3 dir, real sc, bool model)
{
    if ( ix < NB_CLIP_PLANES )
    {
        clip_plane_mode[ix]   = ( model ? 1 : 2 );
        clip_plane_vector[ix] = dir;
        clip_plane_scalar[ix] = sc;
    }
}

void View::disableClipPlane(int ix)
{
    if ( ix < NB_CLIP_PLANES )
    {
        clip_plane_mode[ix] = 0;
        glDisable(GL_CLIP_PLANE0+ix);
    }
}

int View::hasClipPlane(int ix) const
{
    if ( ix < NB_CLIP_PLANES )
        return clip_plane_mode[ix];
    return false;
}

//------------------------------------------------------------------------------
#pragma mark -


void View::displayMagnifier(GLint Z, Vector3 foc, GLint mX, GLint mY) const
{
#if ( 1 )
    int W = width();
    int H = height();
    int P = ( W > H ? H : W ) / 4;
    int M = 2 * P;

/*
    GLint readbuf = 0, drawbuf = 0;
    glGetIntegerv(GL_READ_BUFFER, &readbuf);
    glGetIntegerv(GL_DRAW_BUFFER, &drawbuf);
    printf("normal buffers: read %i draw %i\n", readbuf, drawbuf);
*/
    
    // create off-screen buffer
    if ( OffScreen::createBuffer(M, M, 0) )
    {
        // operate with a copy of the current view:
        View view = *this;
        glMatrixMode(GL_PROJECTION);
        glPushMatrix();
        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();
        view.initGL();
        view.view_size = M * pixelSize() / Z;
        view.reshape(M, M);
        view.move_to(foc);
        view.zoom_to(1);
        setLights();
        setClipping();
        displayCallback(view, 1+Z/3);
        endClipping();
        glMatrixMode(GL_PROJECTION);
        glPopMatrix();
        glMatrixMode(GL_MODELVIEW);
        glPopMatrix();
        glViewport(0, 0, W, H);
/*
        GLubyte * tmp = (GLubyte*)malloc(3*M*M*sizeof(GLubyte));
        glReadPixels(0, 0, M, M, GL_RGB, GL_UNSIGNED_BYTE, tmp);
        SaveImage::savePixels("pixels.png", "png", tmp, W, H, 0);
        free(tmp);
*/
        //print_cap("read", GL_READ_BUFFER);
        //print_cap("draw", GL_DRAW_BUFFER);
        
        // restore writing destination:
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
        //glDrawBuffer(drawbuf);
        
        glBlitFramebuffer(0, 0, M, M, mX-P, mY-P, mX+P, mY+P, GL_COLOR_BUFFER_BIT, GL_LINEAR);
        //checkError("glBlitFramebuffer()");
        
        OffScreen::releaseBuffer();
        //glReadBuffer(readbuf);
    }
#endif
}


void View::drawText(std::string const& str, void* font, gle_color col, int pos) const
{
    gle::drawText(str.c_str(), font, col, pos, width(), height());
}

//------------------------------------------------------------------------------
#pragma mark -

void View::drawCuboid(Vector3 const& A, Vector3 const& B)
{
    GLfloat AX = A.XX, AY = A.YY, AZ = A.ZZ;
    GLfloat BX = B.XX, BY = B.YY, BZ = B.ZZ;
    
    GLfloat pts[24] = {AX, AY, AZ, BX, AY, AZ, BX, BY, AZ, AX, BY, AZ,
                       AX, AY, BZ, BX, AY, BZ, BX, BY, BZ, AX, BY, BZ};
    
    glVertexPointer(3, GL_FLOAT, 0, pts);
    glDrawArrays(GL_LINE_LOOP, 0, 4);
    glDrawArrays(GL_LINE_LOOP, 4, 4);
}


void View::drawROI() const
{
    glLineWidth(1);
    front_color.load();
    glDisable(GL_LIGHTING);
    drawCuboid(mROI[0], mROI[1]);
}


/// init vertical ticks over ] -cnt*d, +cnt*d [
size_t setTicksH(flute4* pts, int cnt, float d, float a, float b)
{
    size_t n = 0;
    for ( int i = 1; i < cnt; ++i )
    {
        pts[n++] = flute4{-i*d, a, -i*d, b};
        pts[n++] = flute4{ i*d, a,  i*d, b};
    }
    pts[n++] = flute4{0, a, 0, b};
    return n*2;
}

/// draw horizontal ticks over ] -cnt*d, +cnt*d [
size_t setTicksV(flute4* pts, int cnt, float d, float a, float b)
{
    size_t n = 0;
    for ( int i = 1; i < cnt; ++i )
    {
        pts[n++] = flute4{a, -i*d, b, -i*d};
        pts[n++] = flute4{a,  i*d, b,  i*d};
    }
    pts[n++] = flute4{a, 0, b, 0};
    return n*2;
}


/**
 This will draw:
 - a horizontal box of length scale, bounded with Y=a and Y=b
 - lines every scale/10, of width (b-a)/5
 - lines every scale/100, of width (b-a)/25
 - lines every scale/1000, of width (b-a)/125
 .
 */
void View::drawScaleH(GLfloat s, GLfloat a, GLfloat b) const
{
    glLineWidth(1);
    gle::drawRectangle(-s/2, a, s/2, b, 0);
    flute4 pts[24];
    glVertexPointer(2, GL_FLOAT, 0, pts);

    // draw bars
    s /= 10;
    GLfloat w(2);
    glLineWidth(w);
    setTicksH(pts, 5, s, a, b);
    glDrawArrays(GL_LINES, 0, 18);
    
    // draw tick marks
    char str[16] = {0};
    do {
        w -= 0.5;
        s /= 10;
        a /= 10;
        b /= 10;
        if ( s > 4 * pixelSize() )
        {
            glLineWidth(w);
            setTicksH(pts, 10, s, a, b);
            glDrawArrays(GL_LINES, 0, 36);
            glRasterPos2f(s-6*pixelSize(), b-12*pixelSize());
            snprintf(str, sizeof(str), "%g", s);
            gle::bitmapText(str);
        }
    } while ( w >= 0.5 );
}


/**
 This will draw:
 - a vertical box of length scale, bounded at X=a and X=b
 - lines every scale/10, of width (b-a)/5
 - lines every scale/100, of width (b-a)/25
 - lines every scale/1000, of width (b-a)/125
 .
 */
void View::drawScaleV(GLfloat s, GLfloat a, GLfloat b) const
{
    glLineWidth(1);
    gle::drawRectangle(a, -s/2, b, s/2, 0);
    flute4 pts[24];
    glVertexPointer(2, GL_FLOAT, 0, pts);

    // draw bars
    s /= 10;
    GLfloat w(2);
    glLineWidth(w);
    setTicksV(pts, 5, s, a, b);
    glDrawArrays(GL_LINES, 0, 18);
    
    // draw tick marks
    char str[16] = {0};
    do {
        w -= 0.5;
        s /= 10;
        a /= 10;
        b /= 10;
        if ( s > 4 * pixelSize() )
        {
            glLineWidth(w);
            setTicksV(pts, 10, s, a, b);
            glDrawArrays(GL_LINES, 0, 36);
            glRasterPos2f(b+pixelSize(), s-4*pixelSize());
            snprintf(str, sizeof(str), "%g", s);
            gle::bitmapText(str);
        }
    } while ( w >= 0.5 );
}


/**
 This will draw a centered cross with :
 - lines every scale/10, of width 1
 - small lines every scale/100, of width 0.5
 - tiny lines every scale/1000, of width 0.25
 .
 */
void View::drawScaleX(GLfloat scale) const
{
    GLfloat s(scale);
    GLfloat a( scale/20);
    GLfloat b(-scale/20);
    GLfloat w(2);

    flute4 pts[24] = {
        {-s, a,-s, b}, {s, a, s, b},
        { a,-s, b,-s}, {a, s, b, s},
        {-s, 0, s, 0}, {0,-s, 0, s}};

    glLineWidth(w);
    glVertexPointer(2, GL_FLOAT, 0, pts);
    glDrawArrays(GL_LINES, 0, 12);

    do {
        w -= 0.5;
        s /= 10;
        if ( s > 2 * pixelSize() )
        {
            glLineWidth(w);
            setTicksV(pts, 10, s, a, b);
            glDrawArrays(GL_LINES, 0, 36);
            setTicksH(pts, 10, s, a, b);
            glDrawArrays(GL_LINES, 0, 36);
        }
        a /= 10;
        b /= 10;
    } while ( w >= 0.5 );
}


/**
 */
void View::drawScaleBar(int mode, const real scale) const
{
    glPushAttrib(GL_ENABLE_BIT);
    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_BLEND);
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    
    GLfloat shift(32 * pixelSize() * zoom);
    GLfloat S(scale);
    
    switch( mode )
    {
        case 0:
            break;
        case 1:
            gle::translate(0, shift-0.5*visRegion[1], 0);
            gle::scale(zoom);
            drawScaleH(S, S/10, 0);
            break;
        case 2:
            gle::translate(0.5*visRegion[0]-shift, 0, 0);
            gle::scale(zoom);
            drawScaleV(S, -S/10, 0);
            break;
        case 3: {
            gle::scale(zoom);
            drawScaleX(S);
        } break;
    }
    
    glPopMatrix();
    glPopAttrib();
}
