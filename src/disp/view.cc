// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University

#include "assert_macro.h"
#include "view.h"
#include "gle.h"
#include "flute.h"
#include "gym_flute.h"
#include "gym_flat.h"
#include "gym_check.h"
#include "gym_matrix.h"
#include "gym_view.h"
#include "gym_draw.h"
#include "gym_cap.h"
#include "gym_vect.h"
#include "fg_font.h"
#include "fg_stroke.h"
#include "time_date.h"


//------------------------------------------------------------------------------

View::View(const std::string& n, int depth)
: ViewProp(n)
{
    window_ = 0;
    drawCallback = nullptr;
    drawMagFunc = nullptr;
    depth_test = depth;
    
    visRange[0] = view_scale;
    visRange[1] = view_scale;
    eyeDistance = -0.25f * view_scale;

    viewport_[0] = 0;
    viewport_[1] = 0;
    viewport_[2] = 800;
    viewport_[3] = 800;

    gym::mat_diagonal(modelview_, 1);
    gym::mat_diagonal(projection_, 1);
    
    hasROI = false;
    mROI[0].reset();
    mROI[1].reset();
}


View::~View()
{
}


//------------------------------------------------------------------------------
#pragma mark -

void View::initGL() const
{
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnable(GL_NORMALIZE);
    glDisable(GL_STENCIL_TEST);
    glDisable(GL_DITHER);
    
    gym::initializeBlending();
    
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
    }
    
#ifdef GL_VERSION_2_1
    if ( depth_clamp )
        glEnable(GL_DEPTH_CLAMP);
    else
        glDisable(GL_DEPTH_CLAMP);
#endif
    
    if ( depth_test )
    {
        glClearDepth(1);
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
#ifdef GL_VERSION_2_1
    if ( depth_clamp )
    {
        depth_clamp = false;
        glDisable(GL_DEPTH_CLAMP);
    }
    else
    {
        depth_clamp = true;
        glEnable(GL_DEPTH_CLAMP);
    }
#endif
}


void View::clearPixels()
{
    gym::clearPixels(back_color);
}


void View::openDisplay() const
{
    adjust();
    loadView();
    gym::openDepthMask();
    gym::clearPixels(back_color);
    setFog(fog_type, fog_param, fog_color);
    setLights();
    setClipping();
    
    if ( floor_radius > 1 )
    {
        gym::disableDepthTest();
        gym::disableLighting();
        gym::closeDepthMask();
        gym::color(floor_color);
        gym::drawTiledFloor(floor_radius, floor_tile, floor_height);
        gym::openDepthMask();
        gym::restoreDepthTest();
    }
    
    if ( axes )
        gle::drawAxes(axes_size, axes);
}


void View::closeDisplay() const
{
    endClipping();
    
    if ( scalebar )
        drawScaleBar(scalebar, scalebar_length, scalebar_color);
    
#if 0
    if ( label != "off" )
    {
        // display only first line of text:
        char str[256];
        size_t n = full_label.find('\n');
        full_label.copy(str, sizeof(str), n);
        str[n] = '\0';
        strokeString(str);
    }
#endif
}

//------------------------------------------------------------------------------
#pragma mark - Text


void View::strokeString(const char str[], float width) const
{
    gym::color(front_color);
    gym::one_view(viewport_[2], viewport_[3]);
    fgStrokeString(8, 4, 1.0, 1, str, width, 0, -0.25);
    loadView();
}


void View::strokeString(float X, float Y, float Z, const char str[]) const
{
    gym::translate_ref(X, Y, Z);
    fgStrokeString(0, 0, pixelSize(), 1, str, 1);
}


/// width = text_width; height = text_heigth, (W, H) = window_size
float View::textPosition(float& px, float& py, int width, int height, int lines,
                        int W, int H, const int position) const
{
    assert_true( W > 0 );
    assert_true( H > 0 );
    
    switch( position )
    {
        case 0:
            //bottom-left, text going up
            px = height/2;
            py = height/2;
            return height;
        case 1:
            //bottom-right, text going up
            px = W - width - height/2;
            if ( px < 0 ) px = 0;
            py = height/2;
            return height;
        case 2:
            //top-right, text going down
            px = W - width - height/2;
            if ( px < 0 ) px = 0;
            py = H - height;
            return -height;
        default:
        case 3:
            //top-left, text going down
            px = height/2;
            py = H - height;
            return -height;
        case 4:
            //center, text going down
            px = ( W - width ) / 2;
            if ( px < 0 ) px = 0;
            py = ( H + lines*height ) / 2;
            return -height;
    }
    return height;
}

    
/**
 The text is displayed in the current color.
 A background rectangle is displayed only if `bcol` is visible.
 
 Possible values for `position`:
 - 0: bottom-left, text going up
 - 1: bottom-right, text going up
 - 2: top-right, text going down
 - 3: top-left, text going down
 - 4: center, text going down
 .
 
 Note: width and height are the current size of the viewport (window)
 */
void View::placeText(int position, FontType font, const float color[4],
                     const char text[], const float back[4], int W, int H) const
{
    int lines = 1;
    int height = fgFontHeight(font);
    int width = fgTextWidth(font, text, lines);
    
    float X, Y;
    float vshift = textPosition(X, Y, width, height, lines, W, H, position);
    
    if ( back && back[3] > 0 )
    {
        float E = height;
        float T = Y + lines * vshift;
        float B = std::min(Y, T) - E/4;
        T = std::max(Y, T) + E + E/4;
        float R = X + width + E;
        gym::paintOctagon(X-E, B, R, T, back, 5);
        if ( position == 4 )
            gym::drawOctagon(X-E, B, R, T, color, 5, 1);
    }
    
    fgBitmapString(X, Y, 1.f, font, color, text, vshift);
    //fgStrokeString(X, Y, 1.f, 0, text, 1, 0, 0);
}

/**
 draw text at position `vec`
 */
void View::drawText(Vector3 const& vec, const float color[4], const char str[], const float offset, FontType) const
{
    gym::disableLighting();
    gym::disableAlphaTest();
    gym::disableDepthTest();
    gym::cancelRotation();
    gym::translate_ref(vec.XX, vec.YY, vec.ZZ);
#if 0
    int H = fgFontHeight(font);
    fgBitmapString(offset, -H/3, pixelSize(), font, color, str, H);
#else
    gym::color(color);
    fgStrokeString(0, 0, pixelSize(), 1, str, 1);
#endif
    gym::restoreDepthTest();
    gym::restoreAlphaTest();
    gym::restoreLighting();
}

//------------------------------------------------------------------------------
#pragma mark -

// display FPS = frames per seconds
void View::drawFPS(size_t arg) const
{
    static char str[16];
    static size_t cnt = 0;
    static double sec = TimeDate::seconds_today();
    cnt += arg;
    double now = TimeDate::seconds_today();
    if ( now > sec + 1.0 )
    {
        double fps = cnt / ( now - sec );
        snprintf(str, sizeof(str), "%6.2f", fps);
        sec = now;
        cnt = 0;
    }
    strokeString(str, 2);
}


void View::drawLabel() const
{
    int W = width(), H = height();
    //set pixel coordinate system:
    gym::disableLighting();
    gym::disableAlphaTest();
    gym::disableDepthTest();
    gym::one_view(W, H);
    placeText(0, BITMAP_9_BY_15, front_color, full_label.c_str(), nullptr, W, H);
    gym::restoreDepthTest();
    gym::restoreAlphaTest();
    gym::restoreLighting();
}


/**
 add over-the-window features for the interactive display
*/
void View::drawInteractiveFeatures() const
{
    int W = width(), H = height();
    gym::disableLighting();

    if ( hasROI )
    {
        gym::color(front_color);
        gle::strokeCuboid(mROI[0], mROI[1], 1);
    }
    //set pixel coordinate system:
    gym::disableAlphaTest();
    gym::disableDepthTest();
    gym::one_view(W, H);

    if ( top_message.size() )
    {
        placeText(3, BITMAP_9_BY_15, front_color, top_message.c_str(), nullptr, W, H);
    }
    
    if ( label != "off" )
    {
        placeText(0, BITMAP_9_BY_15, front_color, full_label.c_str(), nullptr, W, H);
    }
    
    if ( draw_memo && memo.size() )
    {
        float white[4] = {1,1,1,1};
        float black[4] = {0,0,0,0.9};
        placeText(4, BITMAP_8_BY_13, white, memo.c_str(), black, W, H);
    }

    if ( flash.size() )
    {
        float yellow[4] = { 0.6f, 0.6f, 1.f, 1.f };
        placeText(2, BITMAP_9_BY_15, yellow, flash.c_str(), nullptr, W, H);
        if ( TimeDate::seconds_since_1970() > flash_end )
            flash = "";
    }
    
    gym::restoreDepthTest();
    gym::restoreAlphaTest();
    gym::restoreLighting();
}


/**
 Set two light sources
 */
void View::setLights() const
{
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
}

void View::setLightsEye() const
{
    gym::eye_view(eyeDistance, zoom);
    setLights();
}

//------------------------------------------------------------------------------
#pragma mark -

void View::loadViewport() const
{
    glViewport(viewport_[0], viewport_[1], viewport_[2], viewport_[3]);
}

void View::setViewport(int x, int y, size_t w, size_t h) const
{
    glViewport(x, y, w, h);
}


void View::resize(int mag)
{
    //std::clog << "View::reshaped " << W << " " << H << '\n';
    viewport_[2] = mag * window_size[0];
    viewport_[3] = mag * window_size[1];
    adjust(viewport_[2], viewport_[3]);
}


void View::reshape(int W, int H)
{
    window_size[0] = W;
    window_size[1] = H;
    resize(1);
    loadView();
    loadViewport();
}


/// calculate the window size in real units
void View::adjust(int W, int H) const
{
    if ( W > H )
    {
        visRange[0] = view_scale;
        visRange[1] = H * view_scale / W;
    }
    else
    {
        visRange[0] = W * view_scale / H;
        visRange[1] = view_scale;
    }
}


void View::setProjection() const
{
    float X = visRange[0] * 0.5f;
    float Y = visRange[1] * 0.5f;
    float S = view_scale;

    if ( perspective == 3 )
    {
        // this creates a stronger perspective:
        eyeDistance = -1.5f * S;
        gym::mat_frustum(projection_, -X, X, -Y, Y, S, 5.0f*S);
    }
    else if ( perspective == 2 )
    {
        // this creates a strong perspective:
        eyeDistance = -2.0f * S;
        gym::mat_frustum(projection_, -X, X, -Y, Y, S, 6.0f*S);
    }
    else if ( perspective )
    {
        // this creates a perspective:
        eyeDistance = -2.0f * S;
        gym::mat_frustum(projection_, -X, X, -Y, Y, S, 11.0f*S);
    }
    else
    {
        // The back-plane is set behind to avoid clipping
        eyeDistance = -1.25f * S;
        gym::mat_ortho(projection_,-X, X, -Y, Y, 0.25f*S, 2.25f*S);
    }

    //std::clog << "View::setProjection  " << window_ << "\n";
    gym::set_projection(projection_);
}


void View::setPickProjection(float X, float Y, float W, float H) const
{
    float mat[16];
    gym::mat_pick(mat, X, Y, W, H, viewport_);
    gym::mat_multiply(mat, projection_);
    gym::set_projection(mat);
}


void View::setModelView() const
{
    float T[4] = { 0, 0, eyeDistance, 1 };
    rotation.setOpenGLMatrix(modelview_, zoom, T);
    Vector3 V = focus + focus_shift;
    gym::mat_translate(modelview_, -V.XX, -V.YY, -V.ZZ);
    gym::set_view(modelview_);
    
    //std::clog << "setModelView win " << window() << '\n';
#if ( 0 )
    std::clog << "View:eye      " << eyeDistance << "\n";
    std::clog << "View:rotation " << rotation << "\n";
    std::clog << "View:zoom     " << zoom << "\n";
    std::clog << "View:focus    " << focus+focus_shift << "\n";
#endif
    //gym::printMatrices(stderr);
}

float View::pixelSize() const
{
    //float a = view_scale / ( zoom * std::max(width(), height()) );
    float b = visRange[0] / ( zoom * viewport_[2] );
    //float c = visRange[1] / ( zoom * viewport_[3] );
    assert_small(b - visRange[1] / ( zoom * viewport_[3] ));
    return b;
}

void View::loadView() const
{
    //std::clog << "View::loadView() win " << window() << "\n";
    setProjection();
    setModelView();
}


/**
 This will change what is visible in the Z direction near and far
 from the observer, using clipping planes and fog.
 A a function of `mode`:
 - 0 : disabled
 - 1 : show ( Z > 0 )
 - 2 : show ( Z > 0 ) with fog
 - 3 : show ( Z < 0 ) with fog
 - 4 : show slice ( -a < Z < a ) where a = 5% of view_scale
 - 5 : show ( H < 0 ), where H is depth relative to camera
 - 6 : show ( H > 0 ), where H is depth relative to camera
 .
 */
void View::sliceView(int mode) const
{
    if ( 0 == ( mode & 7 ))
        return;
    real off = 0.5;
    gym::ref_view();
    Vector3 V = depthAxis();
    switch ( mode )
    {
        case 1:
            gym::enableClipPlane(2, V.XX, V.YY, V.ZZ, off);
        break;
        case 2:
            gym::enableClipPlane(2,-V.XX,-V.YY,-V.ZZ,-off);
            setFog(1, 1, fog_color);
        break;
        case 3:
            gym::enableClipPlane(2, V.XX, V.YY, V.ZZ, off);
            if ( !depth_clamp )
                setFog(3, 2, fog_color);
        break;
        case 4: {
            real thk = view_scale * 0.1;
            gym::enableClipPlane(2, V.XX, V.YY, V.ZZ, thk);
            gym::enableClipPlane(1,-V.XX,-V.YY,-V.ZZ, thk);
        } break;
        case 5: {
            real thk = view_scale * 0.02;
            gym::enableClipPlane(2, V.XX, V.YY, V.ZZ, thk);
            gym::enableClipPlane(1,-V.XX,-V.YY,-V.ZZ, thk);
        } break;
        case 6:
            gym::enableClipPlane(2, V.XX, V.YY, V.ZZ, 0);
        break;
        case 7:
            gym::enableClipPlane(2,-V.XX,-V.YY,-V.ZZ, 0);
        break;
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


void View::set_scale(float s)
{
    //std::clog << "set_scale " << s << '\n';
    view_scale = s;
    zoom_in(0.933033);
    --auto_scale;
    adjust();
}


void View::zoom_to(float z)
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


void View::rotate_to(const Quaternion<real> & q)
{
    rotation = normalize(q);
    setModelView();
}


void View::align_with(const real dir[3])
{
    rotation.setRotationToVector(dir);
    setModelView();
}

//------------------------------------------------------------------------------
#pragma mark -

void View::adjustROI(float Z)
{
    setROI(unproject(0, 0, Z), unproject(width(), height(), Z));
}


void View::matchROI()
{
    if ( hasROI )
    {
        focus = 0.5 * ( mROI[0] + mROI[1] );
        float R = 0.5 * ( mROI[0] - mROI[1] ).norm_inf();
        
        // zoom only if region is 7 pixels wide:
        if ( R > 7 * pixelSize() )
            zoom = view_scale / R;
        
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


void View::project(float& H, float& V, const real XYZ[3]) const
{
    float vec[4] = { float(XYZ[0]), float(XYZ[1]), float(XYZ[2]), 1.0f };
    float out[4] = { 0 };
    gym::mat_mulvec(out, modelview_, vec);
    gym::mat_mulvec(vec, projection_, out);
    H = vec[0];
    V = vec[1];
}

/**
 Transforms the given window coordinates into user coordinates.
 */
Vector3 View::unproject(float x, float y, float z) const
{
    float un[4] = { 0 };
    gym::unproject(x, y, z, modelview_, projection_, viewport_, un);
    return Vector3(un[0], un[1], un[2]);
}


//------------------------------------------------------------------------------
#pragma mark -

void View::setFog(GLint type, float param, gym_color color) const
{
    GLint gl_type = 0;
    switch( type )
    {
        case 1: gl_type = GL_LINEAR; break;
        case 2: gl_type = GL_LINEAR; break;
        case 3: gl_type = GL_EXP;    break;
        case 4: gl_type = GL_EXP2;   break;
        default: glDisable(GL_FOG); return;
    }
   
    glEnable(GL_FOG);
    glFogi(GL_FOG_MODE, gl_type);
    
    if ( gl_type == GL_LINEAR )
    {
        glFogf(GL_FOG_START, (param)*view_scale);
        glFogf(GL_FOG_END, (param*2+1)*view_scale);
    }
    else
    {
        glFogf(GL_FOG_DENSITY, param/view_scale);
    }
    
    glFogfv(GL_FOG_COLOR, color.colors());
}

void View::enableFog(const GLint type, const float param, gym_color color)
{
    fog_type = type;
    fog_param = param;
    fog_color = color;
}

void View::enableFog(const GLint type, const float param)
{
    fog_type = type;
    fog_param = param;
}

//------------------------------------------------------------------------------
#pragma mark -


void View::setClipping() const
{
    for ( int i = 0; i < NB_CLIP_PLANES; ++i )
    {
        if ( clip_plane_mode[i] )
        {
            // can make the plane relative the viewing 'eye'
            if ( clip_plane_mode[i] == 2 )
                gym::eye_view(eyeDistance, zoom);
            Vector4 const& V = clip_plane[i];
            gym::setClipPlane(i, V.XX, V.YY, V.ZZ, V.TT);
        }
        else
            gym::disableClipPlane(i);
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
        gym::disableClipPlane(ix);
}

void View::enableClipPlane(int ix, real X, real Y, real Z, real S, bool mode)
{
    if ( ix < NB_CLIP_PLANES )
    {
        clip_plane_mode[ix] = ( mode ? 1 : 2 );
        clip_plane[ix].set(X, Y, Z, S);
    }
}

void View::disableClipPlane(int ix)
{
    if ( ix < NB_CLIP_PLANES )
    {
        clip_plane_mode[ix] = 0;
        gym::disableClipPlane(ix);
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

void View::drawMagnifier(float mag, Vector3 foc, Vector3 cen, int mX, int mY, int R) const
{
    if ( drawMagFunc )
    {
        // operate with a copy of the current view:
        View view = *this;
        view.magnify = mag;
        view.zoom *= mag;
        view.focus = foc;
        view.focus_shift = ( cen - foc ) / mag;
        view.focus.ZZ = 0;
        view.focus_shift.ZZ = 0;
        glEnable(GL_SCISSOR_TEST);
        glScissor(mX-R, mY-R, 2*R, 2*R);
        drawMagFunc(view);
        glDisable(GL_SCISSOR_TEST);
    }
}


void View::flashText(std::string const& str)
{
    if ( str != flash )
    {
        //std::clog << " flashText " << str << "\n";
        flash = str;
        flash_end = TimeDate::seconds_since_1970() + 3.0;
    }
}


//------------------------------------------------------------------------------
#pragma mark -


/// set horizontal lines over [ -cnt*d, +cnt*d ]
void setLadderH(float* pts, int cnt, float d, float a, float b)
{
    flute4* flu = (flute4*)pts;
    flu[0] = {0, a, 0, b};
    for ( int i = 1; i <= cnt; ++i )
    {
        flu[2*i-1] = {-i*d, a, -i*d, b};
        flu[2*i  ] = { i*d, b,  i*d, a};
    }
    // duplicate point to be able to draw a rectangle
    flu[2*cnt+1] = flu[2*cnt-1];
}

/// set vertical lines over [ -cnt*d, +cnt*d ]
void setLadderV(float* pts, int cnt, float d, float a, float b)
{
    flute4* flu = (flute4*)pts;
    flu[0] = {a, 0, b, 0};
    for ( int i = 1; i <= cnt; ++i )
    {
        flu[2*i-1] = {a, -i*d, b, -i*d};
        flu[2*i  ] = {b,  i*d, a,  i*d};
    }
    // duplicate point to be able to draw a rectangle
    flu[2*cnt+1] = flu[2*cnt-1];
}


/**
 This will draw:
 - a horizontal box of length `s`, bounded within Y=a and Y=b
 - lines every scale/10, of width (b-a)/5
 - lines every scale/100, of width (b-a)/25
 - lines every scale/1000, of width (b-a)/125
 .
 */
void View::drawScaleHV(float S, float a, float b, void (*func)(float*, int cnt, float, float, float), float border) const
{
    float W(2);
    S /= 10;
    
    float* flt = (float*)gym::mapBufferV2(24);
    func(flt, 5, S, a, b);
    gym::unmapBufferV2();
    gym::drawLineStrip(W+border, 18, 5);
    gym::drawLines(W+border, 0, 18);

    // draw tick marks
    char str[16] = {0};
    do {
        W -= 0.5;
        S /= 10;
        a /= 10;
        b /= 10;
        float Z = 10 * pixelSize();
        if ( S > Z )
        {
            flt = (float*)gym::mapBufferV2(44);
            func(flt, 10, S, a, b);
            gym::unmapBufferV2();
            gym::drawLines(W+border, 2, 36);
            snprintf(str, sizeof(str), "%g", S);
            fgStrokeString(S-Z, b+S+Z, pixelSize(), 1, str, 1);
        }
    } while ( W >= 0.5 );
}

/**
 This will draw a centered cross with :
 - lines every scale/10, of width 1
 - small lines every scale/100, of width 0.5
 - tiny lines every scale/1000, of width 0.25
 .
 */
void View::drawScaleX(float scale, float border) const
{
    float s(scale);
    float a( scale/20);
    float b(-scale/20);
    float W(2);

    flute4* flu = (flute4*)gym::mapBufferV2(12);
    flu[0] = {-s, a,-s, b};
    flu[1] = { s, a, s, b};
    flu[2] = { a,-s, b,-s};
    flu[3] = { a, s, b, s};
    flu[4] = {-s, 0, s, 0};
    flu[5] = { 0,-s, 0, s};
    gym::unmapBufferV2();
    gym::drawLines(W+border, 0, 12);

    do {
        W -= 0.5;
        s /= 10;
        if ( s > 2 * pixelSize() )
        {
            float* flt = (float*)gym::mapBufferV2(44);
            setLadderV(flt, 10, s, a, b);
            gym::unmapBufferV2();
            gym::drawLines(W+border, 2, 36);
            flt = (float*)gym::mapBufferV2(44);
            setLadderH(flt, 10, s, a, b);
            gym::unmapBufferV2();
            gym::drawLines(W+border, 2, 36);
        }
        a /= 10;
        b /= 10;
    } while ( W >= 0.5 );
}


/**
 */
void View::drawScaleBar(int mode, const float S, const float color[4]) const
{
    gym::disableDepthTest();
    gym::disableLighting();

    float shift(32 * pixelSize() * zoom);
    
    switch( mode )
    {
        case 0:
            break;
        case 1:
            gym::eye_view(0, shift-0.5*visRange[1], eyeDistance, zoom);
            gym::color(color);
            drawScaleHV(S, S/10, 0, setLadderH);
            break;
        case 2:
            gym::eye_view(0.5*visRange[0]-shift, 0, eyeDistance, zoom);
            gym::color(color);
            drawScaleHV(S, -S/10, 0, setLadderV);
            break;
        case 3: {
            gym::eye_view(0, 0, eyeDistance, zoom);
            gym::color(color);
            drawScaleX(S);
        } break;
    }

    gym::restoreLighting();
    gym::restoreDepthTest();
}
