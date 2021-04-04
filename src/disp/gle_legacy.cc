

    /// display a cone directed along Z, of radius R at Z=B, and 0 at Z=T
    void coneZ(GLfloat R, GLfloat B, GLfloat T);
    /// display a cone directed along Z, of radius R at Z=B, and 0 at Z=T
    void discZ(GLfloat R, GLfloat Z, GLfloat N);
    /// draw an open tube from B to T along Z, of diameter 1
    void tubeZ(GLfloat B, GLfloat T, int inc);
    /// draw an open tube from B to T along Z, of diameter 1
    void tubeZ(GLfloat B, GLfloat rB, GLfloat T, GLfloat rT, int inc);
    /// draw an open tube along Z, of diameter 1 and length 1, Z=[0, 1]
    void hexTubeZ(GLfloat Zmin, GLfloat Zmax);
    /// draw Torus of radius `rad` and thickness `thick`
    void torusZ(GLfloat rad, GLfloat thick, size_t inc = 1);
    /// spherocylinder of length L, radius R, centered and aligned with axis Z
    void capsuleZ(GLfloat B, GLfloat T, GLfloat R);
    /// draw ellipse
    void ellipseZ(GLfloat rX, GLfloat rY, GLfloat rZ);
    /// draw circle on the ellipse
    void ellipse_circleZ(GLfloat rX, GLfloat rY, GLfloat rZ, GLfloat u);


    //-----------------------------------------------------------------------
    #pragma mark - Legacy 3D objects
    
    void coneZ(GLfloat R, GLfloat B, GLfloat T)
    {
        glBegin(GL_TRIANGLE_FAN);
        glNormal3f( 0, 0, 1 );
        glVertex3f( 0, 0, T );
        const GLfloat L(T-B);
        const GLfloat Y = 1.f/sqrtf(L*L+1);
        const GLfloat X = Y * L;
        for ( size_t n = 0; n <= pi_twice; ++n )
        {
            GLfloat S = sin_(n), C = cos_(n);
            glNormal3f(X*C, X*S, Y);
            glVertex3f(R*C, R*S, B);
        }
        glEnd();
    }

    void tubeZ(GLfloat B, GLfloat T, int inc)
    {
        assert_true( B <= T );
        glBegin(GL_TRIANGLE_STRIP);
        for( size_t n = 0; n <= pi_twice; n += inc )
        {
            GLfloat S = sin_(n), C = cos_(n);
            glNormal3f(C, S, 0);
            glVertex3f(C, S, T);
            glVertex3f(C, S, B);
        }
        glEnd();
    }

    void tubeZ(GLfloat B, GLfloat rB, GLfloat T, GLfloat rT, int inc)
    {
        assert_true( B <= T );
        const GLfloat H(T-B);
        const GLfloat N = 1.f/sqrtf(H*H+(rT-rB)*(rT-rB));
        const GLfloat tC = N * (T-B);
        const GLfloat tS = N * (rB-rT);
        glBegin(GL_TRIANGLE_STRIP);
        for( size_t n = 0; n <= pi_twice; n += inc )
        {
            GLfloat S = sin_(n), C = cos_(n);
            glNormal3f(tC*C, tC*S, tS);
            glVertex3f(rT*C, rT*S, T);
            glVertex3f(rB*C, rB*S, B);
        }
        glEnd();
    }

    void tubeZ(GLfloat B, GLfloat R, gle_color col, GLfloat T, GLfloat D, gle_color loc)
    {
        glBegin(GL_TRIANGLE_STRIP);
        for( size_t n = 0; n <= pi_twice; ++n )
        {
            GLfloat S = sin_(n), C = cos_(n);
            col.load_load();
            glNormal3f(  C,   S, 0);
            glVertex3f(D*C, D*S, B);
            loc.load_load();
            glNormal3f(  C,   S, 0);
            glVertex3f(R*C, R*S, T);
        }
        glEnd();
    }
    
    /// draw spherocylinder of radius R, of axis Z with Z in [B, T]
    void capsuleZ(GLfloat B, GLfloat T, GLfloat R)
    {
        const size_t fin = pi_twice >> 2;
        const size_t inc = 4;
        //display strips along the side of the volume:
        for ( size_t t = 0; t < 4*fin; t += inc )
        {
            //compute the transverse angles:
            GLfloat cb = cos_(t),     sb = sin_(t);
            GLfloat ca = cos_(t+inc), sa = sin_(t+inc);
            GLfloat cB = R * cb, sB = R * sb;
            GLfloat cA = R * ca, sA = R * sa;
            
            //draw one srip of the oval:
            glBegin(GL_TRIANGLE_STRIP);
            for ( size_t i=0; i <= fin; i += inc )
            {
                GLfloat x = cos_(i), y = sin_(i);
                glNormal3f(ca*y, sa*y, x);
                glVertex3f(cA*y, sA*y, T+R*x);
                glNormal3f(cb*y, sb*y, x);
                glVertex3f(cB*y, sB*y, T+R*x);
            }
            for ( int i=fin; i >= 0; i -= inc )
            {
                GLfloat x = -cos_(i), y = sin_(i);
                glNormal3f(ca*y, sa*y, x);
                glVertex3f(cA*y, sA*y, B+R*x);
                glNormal3f(cb*y, sb*y, x);
                glVertex3f(cB*y, sB*y, B+R*x);
            }
            glEnd();
        }
    }
    
    /// draw a Torus of radius R and a thickness 2*T
    void torusZ(GLfloat R, GLfloat T, size_t inc)
    {
        for ( size_t n = 0; n < pi_twice; n += inc )
        {
            GLfloat X0 = cos_(n    ), Y0 = sin_(n    );
            GLfloat X1 = cos_(n+inc), Y1 = sin_(n+inc);
            glBegin(GL_TRIANGLE_STRIP);
            for ( size_t p = 0; p <= pi_twice; p += 2*inc )
            {
                GLfloat S = sin_(p), C = cos_(p);
                glNormal3f(X0*C, Y0*C, S);
                glVertex3f(X0*(R+T*C), Y0*(R+T*C), T*S);
                glNormal3f(X1*C, Y1*C, S);
                glVertex3f(X1*(R+T*C), Y1*(R+T*C), T*S);
            }
            glEnd();
        }
    }
    
    void hexTubeZ(GLfloat A, GLfloat B)
    {
        /// draw hexagon that has the same surface as a disc of radius 1.
        constexpr GLfloat R = 1.0996361107912678f; //std::sqrt( 2 * M_PI / ( 3 * std::sqrt(3) ));
        constexpr GLfloat C = 0.8660254037844386f; //std::sqrt(3)/2;
        constexpr GLfloat S = 0.5f;
        constexpr GLfloat H = R * C, X = R * S;
        
        const GLfloat pts[] = {
             R,  0, B,  R,  0, A,
             X,  H, B,  X,  H, A,
            -X,  H, B, -X,  H, A,
            -R,  0, B, -R,  0, A,
            -X, -H, B, -X, -H, A,
             X, -H, B,  X, -H, A,
             R,  0, B,  R,  0, A };
        
        constexpr GLfloat dir[] = {
             1,  0, 0,  1,  0, 0,
             S,  C, 0,  S,  C, 0,
            -S,  C, 0, -S,  C, 0,
            -1,  0, 0, -1,  0, 0,
            -S, -C, 0, -S, -C, 0,
             S, -C, 0,  S, -C, 0,
             1,  0, 0,  1,  0, 0 };

        glEnableClientState(GL_VERTEX_ARRAY);
        glEnableClientState(GL_NORMAL_ARRAY);
        glVertexPointer(3, GL_FLOAT, 0, pts);
        glNormalPointer(GL_FLOAT, 0, dir);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 14);
        glDisableClientState(GL_VERTEX_ARRAY);
        glEnableClientState(GL_NORMAL_ARRAY);
    }

    void ellipseZ(GLfloat rX, GLfloat rY, GLfloat rZ)
    {
        GLfloat iX(1.f/rX), iY(1.f/rY), iZ(1.f/rZ);
        /*
         A vector orthogonal to the ellipse surface at position (X, Y, Z) is
         ( X / rX^2, Y / rY^2, Z / rZ^2 )
          */
        for ( size_t n = 0; n < pi_once; ++n )
        {
            GLfloat uC = cos_(n  ), uS = sin_(n  );
            GLfloat lC = cos_(n+1), lS = sin_(n+1);
            GLfloat uX = uS * rX, uY = uS * rY, uZ = uC * rZ;
            GLfloat lX = lS * rX, lY = lS * rY, lZ = lC * rZ;
            GLfloat xu = uS * iX, yu = uS * iY, zu = uC * iZ;
            GLfloat xl = lS * iX, yl = lS * iY, zl = lC * iZ;
            glBegin(GL_TRIANGLE_STRIP);
            for ( size_t p = 0; p <= pi_twice; ++p )
            {
                GLfloat S = sin_(p), C = cos_(p);
                glNormal3f(C*xu, S*yu, zu);
                glVertex3f(C*uX, S*uY, uZ);
                glNormal3f(C*xl, S*yl, zl);
                glVertex3f(C*lX, S*lY, lZ);
            }
            glEnd();
        }
    }

    void ellipse_circleZ(GLfloat rX, GLfloat rY, GLfloat rZ, GLfloat u)
    {
        GLfloat R(std::sqrt((1-u)*(1+u)));
        GLfloat iX(1.f/rX), iY(1.f/rY), iZ(u/rZ);
        GLfloat vX(R*rX), vY(R*rY), vZ(u*rZ);
        glBegin(GL_LINE_LOOP);
        for ( size_t n = 0; n <= pi_twice; ++n )
        {
            GLfloat S = sin_(n), C = cos_(n);
            glNormal3f(C*iX, S*iY, iZ);
            glVertex3f(C*vX, S*vY, vZ);
        }
        glEnd();
    }


