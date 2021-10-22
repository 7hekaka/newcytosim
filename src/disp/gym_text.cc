// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University

#include "glut.h"
#include "gym_text.h"
#include "gle_color.h"
#include "vector.h"

namespace gym
{
    int fontHeight(void* font)
    {
        if ( font == GLUT_BITMAP_8_BY_13 )        return 13;
        if ( font == GLUT_BITMAP_9_BY_15 )        return 15;
        if ( font == GLUT_BITMAP_TIMES_ROMAN_10 ) return 11;
        if ( font == GLUT_BITMAP_TIMES_ROMAN_24 ) return 26;
        if ( font == GLUT_BITMAP_HELVETICA_10 )   return 11;
        if ( font == GLUT_BITMAP_HELVETICA_12 )   return 15;
        if ( font == GLUT_BITMAP_HELVETICA_18 )   return 22;
        return 13;
    }
    
    
    /**
     Compute the max width of all the lines in the given text
     This uses GLUT, which should be initialized.
     */
    int maxTextWidth(const char text[], void* font, int& lines)
    {
        int res = 0;
        lines = 0;
        int w = 0;
        for (const char* c = text; *c != '\0' ; ++c)
        {
            if ( *c == '\n' )
            {
                res = std::max(res, w);
                ++lines;
                w = 0;
            }
            else if ( isspace(*c))
            {
                w += glutBitmapWidth(font, ' ');
            }
            else if ( isprint(*c))
            {
                w += glutBitmapWidth(font, *c);
            }
        }
        res = std::max(res, w);
        if ( res > 0 )
            lines = std::max(1, lines);
        return res;
    }
    

    /**
     draw the string character per character using:
     glutBitmapCharacter()
     */
    void bitmapString(const char text[], void* font, GLfloat vshift)
    {
        if ( !font )
        {
            font = GLUT_BITMAP_HELVETICA_12;
            vshift = sign_real(vshift) * fontHeight(font);
        }
        if ( vshift == 0 )
            vshift = -fontHeight(font);
        
        GLfloat ori[4], pos[4];
        glGetFloatv(GL_CURRENT_RASTER_POSITION, ori);
        
        for (const char* p = text; *p; ++p)
        {
            if ( *p == '\n' )
            {
                glGetFloatv(GL_CURRENT_RASTER_POSITION, pos);
                glBitmap(0, 0, 0, 0, ori[0]-pos[0], vshift, nullptr);
            }
            else if ( isspace(*p) )
            {
                glutBitmapCharacter(font, ' ');
            }
            else if ( isprint(*p) )
            {
                glutBitmapCharacter(font, *p);
            }
        }
    }
    
    
    /**
     draw text at position `vec`, if this corresponds to a valid raster position
     */
    void drawText(Vector3 const& vec, const char text[], void* font, float dx)
    {
        GLboolean valid = false;
        glRasterPos3f(vec.x(), vec.y(), vec.z());
        glGetBooleanv(GL_CURRENT_RASTER_POSITION_VALID, &valid);

        if ( valid == GL_TRUE )
        {
            GLboolean depth = glIsEnabled(GL_DEPTH_TEST);
            GLboolean alpha = glIsEnabled(GL_ALPHA_TEST);
            GLboolean light = glIsEnabled(GL_LIGHTING);
            glDisable(GL_DEPTH_TEST);
            glDisable(GL_ALPHA_TEST);
            glDisable(GL_LIGHTING);
            int L = 1;
            int H = fontHeight(font);
            int W = maxTextWidth(text, font, L);
            // center text:
            glBitmap(0,0,0,0,-W*dx,-H/3,nullptr);
            bitmapString(text, font, H);
            if ( depth ) glEnable(GL_DEPTH_TEST);
            if ( alpha ) glEnable(GL_ALPHA_TEST);
            if ( light ) glEnable(GL_LIGHTING);
        }
    }
    
    void drawText(Vector2 const& w, const char text[], void* font, float dx)
    {
        drawText(Vector3(w.XX, w.YY, 0), text, font, dx);
    }
    
    void drawText(Vector1 const& w, const char text[], void* font, float dx)
    {
        drawText(Vector3(w.XX, 0, 0), text, font, dx);
    }

    
    void paintOctogon(const int rec[4], const int rad)
    {
        GLfloat L(rec[0]), B(rec[1]), R(rec[2]), T(rec[3]);
        GLfloat D(rad);
        GLfloat pts[16] = {R,B+D, R,T-D, R-D,B, R-D,T, L+D,B, L+D,T, L,B+D, L,T-D};
        glVertexPointer(2, GL_FLOAT, 0, pts);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 8);
    }
    
    void drawOctogon(const int rec[4], const int rad)
    {
        GLfloat L(rec[0]), B(rec[1]), R(rec[2]), T(rec[3]);
        GLfloat D(rad);
        GLfloat pts[18] = {L,B+D, L+D,B, R-D,B, R,B+D, R,T-D, R-D,T, L+D,T, L,T-D, L,B+D};
        glVertexPointer(2, GL_FLOAT, 0, pts);
        glDrawArrays(GL_LINE_STRIP, 0, 9);
    }

    /**
     The text is displayed in the current color.
     A background rectangle is displayed only if `bcol` is visible.
     
         glColor3f(1,1,1);
         drawText(fKeyString, GLUT_BITMAP_8_BY_13, 0x0, 1);
     
     Possible values for `position`:
     - 0: bottom-left, text going up
     - 1: bottom-right, text going up
     - 2: top-right, text going down
     - 3: top-left, text going down
     - 4: center, text going down
     .
     
     Note: width and height are the current size of the viewport (window)
     */
    void drawText(const char text[], void* font, const gle_color back,
                  const int position, int width, int height)
    {
        assert_true( width > 0 );
        assert_true( height > 0 );
        
        if ( !font )
            font = GLUT_BITMAP_9_BY_15;
        
        int n_lines = 1;
        int lineHeight = fontHeight(font);
        int textWidth = maxTextWidth(text, font, n_lines);
        
        GLint px, py;
        switch( position )
        {
            case 0: {
                //bottom-left, text going up
                px = lineHeight/2;
                py = lineHeight/2;
            } break;
            case 1: {
                //bottom-right, text going up
                px = width - textWidth - lineHeight/2;
                if ( px < 0 ) px = 0;
                py = lineHeight/2;
            } break;
            case 2: {
                //top-right, text going down
                px = width - textWidth - lineHeight/2;
                if ( px < 0 ) px = 0;
                py = height - lineHeight;
                lineHeight *= -1;
            } break;
            default:
            case 3: {
                //top-left, text going down
                px = lineHeight/2;
                py = height - lineHeight;
                lineHeight *= -1;
            } break;
            case 4: {
                //center, text going down
                px = ( width - textWidth ) / 2;
                if ( px < 0 ) px = 0;
                py = ( height + n_lines*lineHeight ) / 2;
                lineHeight *= -1;
            } break;
        }
        
        //set pixel coordinate system:
        GLboolean depth = glIsEnabled(GL_DEPTH_TEST);
        GLboolean alpha = glIsEnabled(GL_ALPHA_TEST);
        GLboolean light = glIsEnabled(GL_LIGHTING);
        glDisable(GL_DEPTH_TEST);
        glDisable(GL_ALPHA_TEST);
        glDisable(GL_LIGHTING);

        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();
        glLoadIdentity();
        
        glMatrixMode(GL_PROJECTION);
        glPushMatrix();
        glLoadIdentity();
        glOrtho(0, width, 0, height, 0, 1);
        
        glRasterPos2i(0, 0);
        glBitmap(0, 0, 0, 0, px, py, nullptr);
        
        if ( back.visible() )
        {
            float col[4];
            glGetFloatv(GL_CURRENT_COLOR, col);
            int R = abs(lineHeight);
            int B = std::min(py, py + n_lines * lineHeight);
            int T = std::max(py, py + n_lines * lineHeight);
            int rec[4] = { px-R, B, px+textWidth+R, T+R+R/2+R/4 };
            back.load();
            paintOctogon(rec, 3);
            if ( position == 4 )
            {
                glLineWidth(0.5);
                drawOctogon(rec, 3);
            }
            glColor4fv(col);
        }
        
        bitmapString(text, font, lineHeight);
        
        glMatrixMode(GL_PROJECTION);
        glPopMatrix();
        glMatrixMode(GL_MODELVIEW);
        glPopMatrix();
        if ( depth ) glEnable(GL_DEPTH_TEST);
        if ( alpha ) glEnable(GL_ALPHA_TEST);
        if ( light ) glEnable(GL_LIGHTING);
    }
    
    
}
