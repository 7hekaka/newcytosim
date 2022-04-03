// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University

#include "assert_macro.h"
#include "view.h"
#include "gle.h"
#include "flute.h"
#include "offscreen.h"
#include "gym_check.h"
#include "gym_matrix.h"
#include "glu_unproject.cc"
#include "glut.h"
#include "time_date.h"


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
    
    viewport_[0] = 0;
    viewport_[1] = 0;
    viewport_[2] = 800;
    viewport_[3] = 800;

    hasROI = false;
    mROI[0].reset();
    mROI[1].reset();
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


void View::toggleDepthClamp()
{
    if ( glIsEnabled(GL_DEPTH_CLAMP) )
    {
        depth_clamp = false;
        glDisable(GL_DEPTH_CLAMP);
    }
    else
    {
        depth_clamp = true;
        glEnable(GL_DEPTH_CLAMP);
    }
}


void View::openDisplay()
{
    load();
    glDepthMask(GL_TRUE);
    back_color.load_clear();
    glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
    setFog(fog_type, fog_param, fog_color);
    setLights();
    setClipping();
    
    if ( floor_radius > 1 )
    {
        GLboolean depth = glIsEnabled(GL_DEPTH_TEST);
        glDisable(GL_DEPTH_TEST);
        glDisable(GL_LIGHTING);
        glDepthMask(GL_FALSE);
        gle::drawTiledFloor(floor_radius, floor_tile, floor_height, floor_color, back_color);
        glDepthMask(GL_TRUE);
        if ( depth ) glEnable(GL_DEPTH_TEST);
    }
}


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
        drawText(BITMAP_9_BY_15, msg, 0x0, 0);
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
        drawText(BITMAP_8_BY_13, memo, 0x000000CC, 4);
    }

    if ( top_message.size() )
    {
        front_color.load();
        drawText(BITMAP_9_BY_15, top_message, back_color.alpha(0.5), 3);
    }
    
    if ( label != "off" )
    {
        front_color.load();
        drawText(BITMAP_9_BY_15, full_label, 0x0, 0);
    }
    
#if ( 0 )
    // display FPS = frames per seconds
    static char buf[16];
    static size_t cnt = 0;
    static double sec = TimeDate::seconds_today();
    ++cnt;
    double now = TimeDate::seconds_today();
    if ( now > sec + 1.0 )
    {
        double fps = cnt / ( now - sec );
        snprintf(buf, sizeof(buf), "%6.2f fps", fps);
        sec = now;
        cnt = 0;
    }
    front_color.load();
    drawText(BITMAP_9_BY_15, buf, 0x0, 1);
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
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
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
    viewport_[2] = W;
    viewport_[3] = H;
    load();
}


/// calculate the visible region in the 3 directions:
void View::adjust()
{
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
    float X = visRegion[0] * 0.5f;
    float Y = visRegion[1] * 0.5f;
    float Z = visRegion[2];
    float S = view_size;

    if ( perspective == 3 )
    {
        // this creates a stronger perspective:
        eyePosition[2] = -1.5f * S;
        gym::mat_frustum(projection_, -X, X, -Y, Y, Z, 5.0f*Z);
    }
    else if ( perspective == 2 )
    {
        // this creates a strong perspective:
        eyePosition[2] = -2.0f * S;
        gym::mat_frustum(projection_, -X, X, -Y, Y, Z, 6.0f*Z);
    }
    else if ( perspective )
    {
        // this creates a perspective:
        eyePosition[2] = -2.0f * S;
        gym::mat_frustum(projection_, -X, X, -Y, Y, Z, 11.0f*Z);
    }
    else
    {
        // The back-plane is set behind to avoid clipping
        eyePosition[2] = -0.5f * S;
        gym::mat_ortho(projection_,-X, X, -Y, Y, 0, Z);
    }

    //std::clog << "View::setProjection  " << mWindowId << "\n";
    glMatrixMode(GL_PROJECTION);
    glLoadMatrixf(projection_);
    glMatrixMode(GL_MODELVIEW);
}


void View::load()
{
    //std::clog << "View::load() win " << window() << "\n";
    adjust();
    setProjection();
    setModelView();
}


void View::setModelView()
{
    rotation.setOpenGLMatrix(modelview_, zoom, eyePosition);
    Vector3 V = - ( focus + focus_shift );
    gym::mat_translate(modelview_, V.XX, V.YY, V.ZZ);
    glMatrixMode(GL_MODELVIEW);
    glLoadMatrixf(modelview_);

    //std::cerr<<"setModelView win " << window() << '\n';
#if ( 0 )
    std::clog << "View:eye      " << eyePosition << "\n";
    std::clog << "View:rotation " << rotation << "\n";
    std::clog << "View:zoom     " << zoom << "\n";
    std::clog << "View:focus    " << focus+focus_shift << "\n";
#endif
    //gym::print_matrices(stderr);
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
    switch ( mode )
    {
        case 1: {
            setClipPlaneEye(GL_CLIP_PLANE1, Vector3(0,0,+1), off);
            if ( !depth_clamp )
                setFog(1, 0, fog_color);
        } break;
        case 2: {
            setClipPlaneEye(GL_CLIP_PLANE1, Vector3(0,0,-1), -off);
            setFog(1, 1, fog_color);
        } break;
        case 3: {
            real thk = view_size * 0.05;
            setClipPlaneEye(GL_CLIP_PLANE1, Vector3(0,0,-1), thk-off);
            setClipPlaneEye(GL_CLIP_PLANE2, Vector3(0,0,+1), thk+off);
        } break;
        case 4: {
            setClipPlane(GL_CLIP_PLANE1, -depthAxis(), 0);
        } break;
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
    real C = dir.XX, S = axis.norm();
    if ( S > REAL_EPSILON )
    {
        rotation.setFromAxis(axis, std::atan2(S, C));
        setModelView();
    }
}

//------------------------------------------------------------------------------
#pragma mark -

void View::adjustROI(real Z)
{
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

/**
 return axis orthogonal to the display plane, and corresponding to depth
 obtained from the current modelview transformation
 */
Vector3 View::depthAxis() const
{
    return normalize(Vector3(modelview_[2], modelview_[6], modelview_[10]));
}


/**
 Transforms the given window coordinates into user coordinates.
 */
Vector3 View::unproject(GLfloat x, GLfloat y, GLfloat z)
{
    float un[4] = { 0 };
    myUnproject(x, y, z, modelview_, projection_, viewport_, un);
    return Vector3(un[0], un[1], un[2]);
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
        uint8_t * tmp = (uint8_t*)malloc(3*M*M*sizeof(uint8_t));
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


void View::drawText(FontType font, std::string const& str, gle_color col, int pos) const
{
    gym::drawText(font, str.c_str(), col, pos, width(), height());
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
        pts[n++] = {-i*d, a, -i*d, b};
        pts[n++] = { i*d, a,  i*d, b};
    }
    pts[n++] = {0, a, 0, b};
    return n*2;
}

/// draw horizontal ticks over ] -cnt*d, +cnt*d [
size_t setTicksV(flute4* pts, int cnt, float d, float a, float b)
{
    size_t n = 0;
    for ( int i = 1; i < cnt; ++i )
    {
        pts[n++] = {a, -i*d, b, -i*d};
        pts[n++] = {a,  i*d, b,  i*d};
    }
    pts[n++] = {a, 0, b, 0};
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
void View::drawScaleH(float s, float a, float b) const
{
    glLineWidth(1);
    gle::drawRectangle(-s/2, a, s/2, b, 0);
    flute4 pts[24];
    GLfloat w(2);

    // draw bars
    s /= 10;
    setTicksH(pts, 5, s, a, b);
    glLineWidth(w);
    glVertexPointer(2, GL_FLOAT, 0, pts);
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
            gym::bitmapString(BITMAP_HELVETICA_12, str, 15);
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
void View::drawScaleV(float s, float a, float b) const
{
    glLineWidth(1);
    gle::drawRectangle(a, -s/2, b, s/2, 0);
    flute4 pts[24];
    GLfloat w(2);

    // draw bars
    s /= 10;
    setTicksV(pts, 5, s, a, b);
    glLineWidth(w);
    glVertexPointer(2, GL_FLOAT, 0, pts);
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
            gym::bitmapString(BITMAP_HELVETICA_12, str, 15);
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
void View::drawScaleX(float scale) const
{
    float s(scale);
    float a( scale/20);
    float b(-scale/20);
    float w(2);

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
void View::drawScaleBar(int mode, const float S) const
{
    GLboolean depth = glIsEnabled(GL_DEPTH_TEST);
    GLboolean light = glIsEnabled(GL_LIGHTING);
    GLboolean blend = glIsEnabled(GL_BLEND);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_LIGHTING);
    glDisable(GL_BLEND);

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    
    GLfloat shift(32 * pixelSize() * zoom);
    
    switch( mode )
    {
        case 0:
            break;
        case 1:
            gle::transScale(0, shift-0.5*visRegion[1], 0, zoom);
            drawScaleH(S, S/10, 0);
            break;
        case 2:
            gle::transScale(0.5*visRegion[0]-shift, 0, 0, zoom);
            drawScaleV(S, -S/10, 0);
            break;
        case 3: {
            gle::transScale(0, 0, 0, zoom);
            drawScaleX(S);
        } break;
    }
    
    glPopMatrix();
    if ( depth ) glEnable(GL_DEPTH_TEST);
    if ( light ) glEnable(GL_LIGHTING);
    if ( blend ) glEnable(GL_BLEND);
}
