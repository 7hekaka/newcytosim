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
    
    viewScale[0] = view_scale;
    viewScale[1] = view_scale;
    eyeDistance = 0.25f * view_scale;

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
    gym::clearPixels(back_color);
    setFog(fog_type, fog_param, fog_color);
    setLights();
    setClipping();
    
    if ( floor_radius > 1 )
    {
        gym::disableLighting();
        gym::color(floor_color);
        gym::drawTiledFloor(floor_radius, floor_tile, floor_height);
    }
    
    gym::openDepthMask();

    if ( axes )
    {
        gym::enableLighting();
        gym::enableCullFace(GL_BACK);
        gle::drawAxes(axes_size, axes);
        gym::disableCullFace();
    }
}


void View::closeDisplay() const
{
    endClipping();
    
    //drawFPS();
    
#if 0
    if ( label != "off" )
    {
        // display only first line of text:
        size_t n = full_label.find('\n');
        std::string str = full_label.substr(0, n);
        strokeString(str.c_str(), 1.f);
    }
#endif
}

//------------------------------------------------------------------------------
#pragma mark - Text


void View::strokeString(const char str[], float size) const
{
    int W = width(), H = height();
    gym::disableLighting();
    gym::disableAlphaTest();
    gym::disableDepthTest();
    gym::one_view(W, H);
    gym::color(back_color);
    fgStrokeString(0, H-18, 1, 1, str, 4+size, 4+size, 0);
    gym::color(front_color);
    fgStrokeString(0, H-18, 1, 1, str, size, size, 0);
    gym::restoreDepthTest();
    gym::restoreAlphaTest();
    gym::restoreLighting();
    loadView();
}


void View::strokeString(float X, float Y, float Z, const char str[], float scale) const
{
    gym::translate_ref(X, Y, Z);
    fgStrokeString(0, 0, scale*pixelSize(), 1, str, 2*scale, 2*scale);
}


/// width = text_width; height = text_heigth, (W, H) = window size
static float textPosition(float& X, float& Y, int width, int height, int lines,
                          int W, int H, const int position)
{
    float B = 12;
    assert_true( W > 0 );
    assert_true( H > 0 );
    
    switch( position )
    {
        case 0: //bottom-left, text going up
            X = B;
            Y = B;
            return height;
        case 1: //bottom-right, text going up
            X = std::max(0, W - width - height/2);
            Y = B;
            return height;
        case 2: //top-right, text going down
            X = std::max(0, W - width - height/2);
            Y = H - height;
            return -height;
        default:
        case 3: //top-left, text going down
            X = B;
            Y = H - height;
            return -height;
        case 4: //centered, text going down
            X = std::max(0, ( W - width ) / 2);
            Y = ( H + lines*height ) / 2;
            return -height;
        case 5: //centered, 1/4 up from the bottom
            X = std::max(0, ( W - width ) / 2);
            Y = ( H + (lines-1)*height ) / 8;
            return -height;
    }
    return 0;
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
 
 Note: W and H are the current size of the viewport (i.e. the window)
 */
void View::frameText(int position, FontType font, const float color[4],
                     const char text[], const float back[4], int W, int H) const
{
    float mag = 1.f;
    int lines = 1;
    int height = fgFontHeight(font);
    int width = fgTextWidth(font, text, lines);
    
    float X = 0, Y = 0;
    float dY = textPosition(X, Y, width, height, lines, W, H, position);
    
    if ( back && back[3] > 0 )
    {
        float E = height;
        float T = Y + lines * dY;
        float B = std::min(Y, T) + E/4;
        T = std::max(Y, T) + E + E/4;
        float R = X + width + E;
        gym::paintOctagon(X-E, B, R, T, back, 6);
        if ( position == 4 )
            gym::drawOctagon(X-E, B, R, T, color, 6, 1);
    }
    fgBitmapText(X, Y, mag, font, color, text, dY);
    //fgStrokeString(X, Y, mag, 0, text, 1, 0, 0);
}

/**
 draw text at position `{X, Y, Z}`
 */
void View::drawText(float X, float Y, float Z, const float color[4], const char str[], const float offset, FontType) const
{
    gym::disableLighting();
    gym::disableAlphaTest();
    gym::disableDepthTest();
    gym::cancelRotation();
    gym::translate_ref(X, Y, Z);
#if 0
    int H = fgFontHeight(font);
    fgBitmapString(offset, -H/3, pixelSize(), font, color, str, H);
#else
    gym::color(color);
    fgStrokeString(offset, 0, pixelSize(), 1, str, 1);
#endif
    gym::restoreDepthTest();
    gym::restoreAlphaTest();
    gym::restoreLighting();
}

//------------------------------------------------------------------------------
#pragma mark -

/// display frames per seconds
void View::drawFPS() const
{
    static char str[16] = { 0 };
    static size_t cnt = 0;
    static double sec = TimeDate::seconds_today();
    cnt++;
    double now = TimeDate::seconds_today();
    if ( now > sec + 1.0 )
    {
        double fps = cnt / ( now - sec );
        snprintf(str, sizeof(str), "%6.2f fps", fps);
        //printf("%s\n", str);
        sec = now;
        cnt = 0;
    }
    strokeString(str, 1);
}


void View::drawTextAndScale(std::string const& text) const
{
    int W = width(), H = height();
    //set pixel coordinate system:
    gym::disableDepthTest();
    gym::disableLighting();

    if ( scalebar )
        drawScaleBar(scalebar, scalebar_length, scalebar_color);
    
    gym::disableAlphaTest();
    gym::one_view(W, H);

    if ( label != "off" )
        frameText(0, BITMAP_9_BY_15, front_color, full_label.c_str(), nullptr, W, H);
    
    if ( text.size() )
    {
        float white[4] = {1,1,1,1};
        float black[4] = {0,0,0,0.9};
        frameText(5, BITMAP_TIMES_ROMAN_24, white, text.c_str(), black, W, H);
    }

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
        gym::ref_view();
        gym::color(front_color);
        gle::strokeCuboid(mROI[0], mROI[1], 1);
    }
    //set pixel coordinate system:
    gym::disableAlphaTest();
    gym::disableDepthTest();
    gym::one_view(W, H);

    if ( top_message.size() )
    {
        frameText(3, BITMAP_9_BY_15, front_color, top_message.c_str(), nullptr, W, H);
    }
    
    if ( subtitle.size() )
    {
        float white[4] = {1,1,1,1};
        float black[4] = {0,0,0,0.9};
        frameText(5, BITMAP_TIMES_ROMAN_24, white, subtitle.c_str(), black, W, H);
    }

    if ( draw_memo && memo.size() )
    {
        float white[4] = {1,1,1,1};
        float black[4] = {0,0,0,0.9};
        frameText(4, BITMAP_8_BY_13, white, memo.c_str(), black, W, H);
    }

    if ( flash.size() )
    {
        float yellow[4] = { 0.6f, 0.6f, 1.f, 1.f };
        frameText(2, BITMAP_9_BY_15, yellow, flash.c_str(), nullptr, W, H);
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

void View::setViewport(int x, int y, unsigned w, unsigned h) const
{
    glViewport(x, y, w, h);
}


void View::resize(int mag)
{
    viewport_[2] = mag * window_size[0];
    viewport_[3] = mag * window_size[1];
    adjust(viewport_[2], viewport_[3]);
    //printf(" View::resize(%6i, %6i)\n", viewport_[2], viewport_[3]);
}


void View::reshape(int W, int H)
{
    window_size[0] = W;
    window_size[1] = H;
    adjust(W, H);
    resize(1);
    loadView();
    loadViewport();
}


/// calculate the window size in real units
void View::adjust(int W, int H) const
{
    if ( W > H )
    {
        viewScale[0] = view_scale;
        viewScale[1] = H * view_scale / W;
    }
    else
    {
        viewScale[0] = W * view_scale / H;
        viewScale[1] = view_scale;
    }
    //printf(" View::viewScale(%6.3f, %6.3f)\n", viewScale[0], viewScale[1]);
}


void View::setProjection() const
{
    float X = viewScale[0] * 0.5f;
    float Y = viewScale[1] * 0.5f;
    float S = view_scale;

    if ( perspective == 3 )
    {
        // this creates a stronger perspective:
        eyeDistance = 1.5f * S;
        gym::mat_frustum(projection_, -X, X, -Y, Y, S, 5.0f*S);
    }
    else if ( perspective == 2 )
    {
        // this creates a strong perspective:
        eyeDistance = 2.0f * S;
        gym::mat_frustum(projection_, -X, X, -Y, Y, S, 6.0f*S);
    }
    else if ( perspective )
    {
        // this creates a perspective:
        eyeDistance = 2.0f * S;
        gym::mat_frustum(projection_, -X, X, -Y, Y, S, 11.0f*S);
    }
    else
    {
        // The back-plane is set behind to avoid clipping
        eyeDistance = 1.25f * S;
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
    float T[4] = { 0, 0, -eyeDistance, 1 };
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
    float b = viewScale[0] / ( zoom * viewport_[2] );
    //float c = viewScale[1] / ( zoom * viewport_[3] );
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
    Vector3 V = depthAxis();
    gym::ref_view();
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
                setFog(1, 1, fog_color);
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
    if ( s > 0 )
    {
        //std::clog << "set_scale " << s << '\n';
        view_scale = s;
        zoom_in(0.933033);
        --auto_scale;
        adjust();
    }
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


/// adjust the current rotation to align with `dir`:
void View::align_with(const real dir[3], bool can_flip, real scale)
{
    real vec[3] = { 0, 0, 0 };
    rotation.rotateVector(vec, dir);
    if ( can_flip )
    {
        real flip = std::copysign(1.0, vec[0]);
        vec[0] *= flip;
        vec[1] *= flip;
        vec[2] *= flip;
    }
    //printf("\ndir: %+9.3f %+9.3f %+9.3f", dir[0], dir[1], dir[2]);
    //printf(" : %+9.3f %+9.3f %+9.3f : ", vec[0], vec[1], vec[2]);
    Quaternion<real> R;
    R.setRotationToVector(vec, scale);
    rotation = ( R * rotation ).normalized();
    //R.print(stdout, true); //rotation.print(stdout, true);
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
 return first 'X' axis in the display plane, in the current modelview transformation
 */
Vector3 View::firstAxis() const
{
    const float * ptr = modelview_;
    return normalize(Vector3(ptr[0], ptr[4], ptr[8]));
}

/**
 return axis orthogonal to the display plane, corresponding to depth axis
 in the current modelview transformation
 */
Vector3 View::depthAxis() const
{
    const float * ptr = modelview_;
    return normalize(Vector3(ptr[2], ptr[6], ptr[10]));
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
    float pt[4] = { x, y, z, 1 };
    gym::unproject(pt, modelview_, projection_, viewport_);
    return Vector3(pt[0], pt[1], pt[2]) - focus_shift;
}


//------------------------------------------------------------------------------
#pragma mark - Fog

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
        // fog to start at the edge of Simulation's volume
        float start = eyeDistance - 0.5*view_scale;
        glFogf(GL_FOG_START, start);
        glFogf(GL_FOG_END, start+param*view_scale);
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
#pragma mark - Clipping planes


void View::setClipping() const
{
    for ( int i = 0; i < NB_CLIP_PLANES; ++i )
    {
        if ( clip_plane_mode[i] )
        {
            gym::ref_view();
            // can make the plane relative the viewing 'eye'
            if ( clip_plane_mode[i] == 2 )
                gym::eye_view(eyeDistance, zoom);
            Vector4 const& V = clip_plane[i];
            gym::enableClipPlane(i, V.XX, V.YY, V.ZZ, V.TT);
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
        view.focus.ZZ = 0;
        view.focus_shift = ( cen - foc ) / mag;
        view.focus_shift.ZZ = 0;
        glEnable(GL_SCISSOR_TEST);
        glScissor(mX-R, mY-R, 2*R, 2*R);
        drawMagFunc(view);
        glDisable(GL_SCISSOR_TEST);
    }
}


bool View::flashText(std::string const& str)
{
    if ( str != flash )
    {
        //std::clog << " flashText " << str << "\n";
        flash = str;
        flash_end = TimeDate::seconds_since_1970() + 3.0;
        return window() == 1;
    }
    return 0;
}


//------------------------------------------------------------------------------
#pragma mark - scale bar


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
    float shift(32 * pixelSize() * zoom);
    
    switch( mode )
    {
        case 0:
            break;
        case 1:
            gym::eye_view(0, shift-0.5*viewScale[1], -eyeDistance, zoom);
            gym::color(color);
            drawScaleHV(S, S/10, 0, setLadderH);
            break;
        case 2:
            gym::eye_view(0.5*viewScale[0]-shift, 0, -eyeDistance, zoom);
            gym::color(color);
            drawScaleHV(S, -S/10, 0, setLadderV);
            break;
        case 3: {
            gym::eye_view(0, 0, -eyeDistance, zoom);
            gym::color(color);
            drawScaleX(S);
        } break;
    }
}
