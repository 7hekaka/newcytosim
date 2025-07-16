// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University.

#include "primitives.h"
#include "assert_macro.h"
#include "exceptions.h"
#include "quaternion.h"
#include "iowrapper.h"
#include "tokenizer.h"
#include "evaluator.h"
#include "glossary.h"
#include "random.h"
#include "space.h"


int Cytosim::skip_space(std::string const& str, size_t& sci, bool eat_line)
{
    assert_true(sci <= str.length());
    int c = str[sci];
    while ( isspace(c) )
    {
        if (( c == '\n' ) & !eat_line )
            break;
        c = str[++sci];
    }
    return c;
}

/// return `true` if stream contains unread non-space character(s)
size_t Cytosim::has_trail(std::string const& str, size_t pos)
{
    while ( isspace(str[pos]) )
        ++pos;
    if ( str[pos] )
        return pos;
    return 0;
}

/** With C++11, the extracted value is zeroed even upon failure.
extract(), on the other hand, will preserve the value of 'var' if no read occurs.
Returns `true` if a value was set. This is used with T=real and T=Vector */
template < typename T >
static bool extract(std::string const& arg, size_t& sci, T& var)
{
    T backup = var;
    std::stringstream iss(arg.substr(sci));
    if ( iss >> var )
    {
        sci += iss.gcount();
        return true;
    }
    var = backup;
    return false;
}


/** Specialized version for floating point values: `var` is unchanged upon failure */
static bool extract(std::string const& arg, size_t& sci, real& var)
{
    const char* str = arg.c_str();
    const char* ptr = str + sci;
    while ( isspace(*ptr) )
        ++ptr;
    errno = 0;
    char * end = nullptr;
#if REAL_IS_DOUBLE
    real x = strtod(ptr, &end);
#else
    real x = strtof(ptr, &end);
#endif
    if ( errno == 0 && end > ptr )
    {
        var = x;
        sci = end - str;
        //std::clog << "extract(" << ptr << ") --float--> " << var << "\n";
        return true;
    }
    //std::clog << "extract(" << ptr << ") not-a-float\n";
    return false;
}


/** Specialized version for Vector: `vec` is unchanged upon failure. */
static bool extract(std::string const& arg, size_t& sci, Vector& vec)
{
    if ( extract(arg, sci, vec.XX) )
    {
        [[maybe_unused]] size_t pos = sci;
        real v;
        if ( extract(arg, sci, v) )
        {
#if ( DIM > 1 )
            vec.YY = v;
#else
            if ( v ) {
                sci = pos;
                return true;
            }
#endif
            pos = sci;
            if ( extract(arg, sci, v) )
            {
#if ( DIM > 2 )
                vec.ZZ = v;
#else
                if ( v ) {
                    sci = pos;
                    return true;
                }
#endif
            }
        }
        else
        {
#if ( DIM > 1 )
            vec.YY = 0;
#endif
#if ( DIM > 2 )
            vec.ZZ = 0;
#endif
        }
        //std::clog << "extract(" << arg << ") --vector--> " << vec << "\n";
        return true;
    }
    return false;
}


/// extract a rotation angle specified in Radian
static real get_angle(std::string const& arg, size_t& sci)
{
    real a = 0;
    if ( extract(arg, sci, a) )
        return a;
    // if no angle is specified, set it randomly:
    return RNG.sreal() * M_PI;
}


/// return 'X' 'Y' or 'Z' if this is the last character of the string
static char get_axis(std::string const& str, size_t pos)
{
    if ( str.length() > pos )
    {
        char c = str.at(pos);
        if ( 'X' <= c && c <= 'Z' )
            return c;
    }
    return 'Z';
}

//------------------------------------------------------------------------------
#pragma mark - Position

/**
 There are different ways to specify a position:
 
 keyword & parameters | Position (X, Y, Z)                                     |
 ---------------------|---------------------------------------------------------
 `A B C`              | The specified vector (A,B,C)
 `inside`             | A random position inside the current Space
 `edge E`             | At distance E from the edge of the current Space
 `surface E`          | On the surface of the current Space\n By projecting a point at distance E from the surface.
 `line L T`           | L: Length, T: thickness. Selected randomly with `-L/2 < X < L/2; norm(Y,Z) < T`
 `sphere R T`         | At distance R +/- T/2 from the origin\n `R-T/2 < norm(X,Y,Z) < R+T/2`
 `ball R`             | At distance R at most from the origin\n `norm(X,Y,Z) < R`
 `disc R T`           | in 2D, a disc in the XY-plane \n in 3D, a disc in the XY-plane of thickness T in Z
 `discXZ R T`         | Disc in the XZ-plane of radius R, thickness T
 `discYZ R T`         | Disc in the YZ-plane of radius R, thickness T
 `equator R T`        | On the sphere of radius R, and at distance T from the plane `Z=0`
 `cap R T`            | On the sphere of radius R, and at distance T from the plane `Z=R`
 `circle R T`         | At distance T from the circle of radius R, in the plane `Z=0`
 `cylinder L R`       | Cylinder of axis X, L: length in X, R: radius in YZ plane
 `cylinderZ L R`      | Cylinder of axis Z, L: length in Z, R: radius in XY plane
 `ring L R T`         | Surface of a cylinder of axis X, L=length in X, R=radius in YZ, T = thickness
 `ellipse A B C`      | Inside the ellipse or ellipsoid of main axes 2A, 2B and 2C
 `arc L Theta`        | A piece of circle of length L and covering an angle Theta
 `stripe L R`         | Random vector with `L < X < R`
 `square R`           | Random vector with `-R < X < R`; `-R < Y < R`; `-R < Z < R`;
 `rectangle A B C`    | Random vector with `-A < X < A`; `-B < Y < B`; `-C < Z < C`;
 `gradient S E`       | Linear density gradient along X, of value 0 at `X=S` and 1 at `X=E`
 `gradient S E R`     | Linear density gradient, contained inside a cylinder of radius R
 `exponential S L`    | Exponential density gradient of length scale L, starting at S
 `exponential S L R`  | Exponential density gradient, contained inside a cylinder of radius R
 `XY Z`               | randomly in the XY plane and within the Space, at specified Z (a value)

 Each primitive describes a certain area in Space, and in most cases the returned position is
 chosen randomly inside this area following a uniform probability.
 */

Vector Cytosim::readPositionPrimitive(std::string const& arg, size_t& sci, Space const* spc)
{
    int c = skip_space(arg, sci, false);

    if ( isalpha(c) )
    {
        std::string tok = Tokenizer::get_symbol(arg, sci);
        //StreamFunc::mark_line(std::cerr, is);

        // primitives are alphabetically sorted in this switch:
        switch ( c )
        {
            case 'a':
                
                if ( tok == "arc" )
                {
                    real L = -1, A = 1.57;
                    if ( !extract(arg, sci, L) || L < 0 )
                        throw InvalidParameter("length L must be >= 0 in `arc L`");
                    extract(arg, sci, A);
                    
                    real x = 0, y = 0;
                    if ( A == 0 ) {
                        x = 0;
                        y = L * RNG.shalf();
                    }
                    else {
                        real R = L / A;
                        real a = A * RNG.shalf();
                        x = R * std::cos(a) - R; // origin centered on arc
                        y = R * std::sin(a);
                    }
                    return Vector(x, y, 0);
                }
                break;

            case 'b':
                
                if ( tok == "ball" )
                {
                    real R = -1;
                    if ( !extract(arg, sci, R) || R < 0 )
                        throw InvalidParameter("radius R must be >= 0 in `ball R`");
                    return Vector::randB(R);
                }
                break;
                
            case 'c':
                
                if ( tok.compare(0, 3, "cap") == 0 )
                {
                    real R = 0, T = 0;
                    if ( !extract(arg, sci, R) || R < 0 )
                        throw InvalidParameter("radius R must be >= 0 in `cap R`");
                    if ( extract(arg, sci, T) && T < 0 )
                        throw InvalidParameter("thickness T must be >= 0 in `cap R T`");
                    real Z = std::max(R - T * RNG.preal(), -R);
                    real r = std::sqrt(R*R - Z*Z);
#if ( DIM >= 3 )
                    real C, S;
                    RNG.urand2(C, S, r);
                    switch ( get_axis(tok, 3) )
                    {
                        case 'X': return Vector(Z, C, S);
                        case 'Y': return Vector(C, Z, S);
                        default : return Vector(C, S, Z);
                    }
#else
                    return Vector(r*RNG.sflip(), Z, 0);
#endif
                }
                
                if ( tok.compare(0, 8, "cylinder") == 0 )
                {
                    real L = -1, R = -1;
                    if ( !extract(arg, sci, L) || L < 0 )
                        throw InvalidParameter("length L must be >= 0 in `cylinder L R`");
                    if ( !extract(arg, sci, R) || R < 0 )
                        throw InvalidParameter("radius R must be >= 0 in `cylinder L R`");
                    const Vector2 V = Vector2::randB(R);
                    switch ( get_axis(tok, 8) )
                    {
                        case 'X': return Vector(L*RNG.shalf(), V.XX, V.YY);
                        case 'Y': return Vector(V.XX, L*RNG.shalf(), V.YY);
                        default : return Vector(V.XX, V.YY, L*RNG.shalf());
                    }
                }
 
                if ( tok.compare(0, 6, "circle") == 0 )
                {
                    real R = -1, T = 0;
                    if ( !extract(arg, sci, R) || R < 0 )
                        throw InvalidParameter("radius R must be >= 0 in `circle R T`");
                    if ( extract(arg, sci, T) && T < 0 )
                        throw InvalidParameter("thickness T must be >= 0 in `circle R T`");
#if ( DIM >= 3 )
                    real C, S;
                    RNG.urand2(C, S, R);
                    Vector3 W = (0.5*T) * Vector3::randU();
                    switch( get_axis(tok, 6) )
                    {
                        case 'X': return Vector3(0, C, S) + W;
                        case 'Y': return Vector3(C, 0, S) + W;
                        default : return Vector3(C, S, 0) + W;
                    }
#endif
                    return Vector::randU(R) + Vector::randU(T*0.5);
                }
                
                if ( tok == "cube" )
                {
                    real R = 1;
                    extract(arg, sci, R);
                    return Vector::randH()*R;
                }
                
                if ( tok == "center" || tok == "origin" )
                    return Vector(0,0,0);

                break;
                
            case 'd':
                
                if ( tok.compare(0, 4, "disc") == 0 )
                {
                    real R = -1, T = 0;
                    if ( !extract(arg, sci, R) || R < 0 )
                        throw InvalidParameter("radius R must be >= 0 in `disc R`");
                    if ( extract(arg, sci, T) && T < 0 )
                        throw InvalidParameter("thickness T must be >= 0 in `disc R T`");
#if ( DIM >= 3 )
                    //in 3D, a disc in the XY-plane of thickness T in Z-direction
                    Vector2 V = Vector2::randB(R);
                    switch( get_axis(tok, 4) )
                    {
                        case 'X': return Vector(T*RNG.shalf(), V.XX, V.YY);
                        case 'Y': return Vector(V.XX, T*RNG.shalf(), V.YY);
                        default : return Vector(V.XX, V.YY, T*RNG.shalf());
                    }
#endif
                    //in 2D, a disc in the XY-plane
                    return Vector::randB(R);
                }
                break;
                
            case 'e':
                
                if ( tok == "ellipse" )
                {
                    Vector S(1, 1, 0);
                    extract(arg, sci, S);
                    return S.e_mul(Vector::randB());
                }
                
                if ( tok == "ellipse_surface" )
                {
                    Vector S(1, 1, 0);
                    extract(arg, sci, S);
                    return S.e_mul(Vector::randU());
                }
                 
                if ( tok == "equator" )
                {
                    real R = 1, T = 0;
                    if ( extract(arg, sci, R) && R < 0 )
                        throw InvalidParameter("radius R must be >= 0 in `equator R T`");
                    if ( extract(arg, sci, T) && T < 0 )
                        throw InvalidParameter("thickness T must be >= 0 in `equator R T`");
                    real C, S;
                    RNG.urand2(C, S, R);
                    return Vector(C, S, T*RNG.shalf());
                }

                if ( spc && tok == "edge" )
                {
                    real R = 0;
                    if ( !extract(arg, sci, R) || R < REAL_EPSILON )
                        throw InvalidParameter("distance R must be > 0 in `edge R`");
                    return spc->placeNearEdge(R);
                }
                
                if ( spc && tok == "exponential" )
                {
                    real B = -10, E = 1, R = 0;
                    extract(arg, sci, B);
                    extract(arg, sci, E);
                    extract(arg, sci, R);
                    if ( R == 0 )
                    {
                        Vector vec;
                        real p;
                        do {
                            vec = spc->place();
                            p = std::exp( ( B - vec.XX ) / E );
                        } while ( p < 0 || p > 1 || p < RNG.preal() );
                        return vec;
                    }
                    real x = RNG.exponential();
#if ( DIM < 3 )
                    return Vector(B+E*x, R*RNG.sreal(), 0);
#else
                    real C, S;
                    RNG.urand2(C, S, R);
                    return Vector(B+E*x, C, S);
#endif
                }
               break;
                
            case 'g':
                
                if ( spc && tok == "gradient" )
                {
                    real B = -10, E = 10, R = 0;
                    extract(arg, sci, B);
                    extract(arg, sci, E);
                    extract(arg, sci, R);
                    if ( R == 0 )
                    {
                        Vector vec;
                        real p;
                        do {
                            vec = spc->place();
                            p = ( vec.XX - B ) / ( E - B );
                        } while ( p < 0 || p > 1 || p < RNG.preal() );
                        return vec;
                    }
                    real x = std::sqrt(RNG.preal());
#if ( DIM < 3 )
                    return Vector(B+(E-B)*x, R*RNG.sreal(), 0);
#else
                    real C, S;
                    RNG.urand2(C, S, R);
                    return Vector(B+(E-B)*x, C, S);
#endif
                }
                break;

            case 'i':
                
                if ( spc && tok == "inside" )
                    return spc->place();
                break;
 
            case 'l':
                
                if ( tok == "line" )
                {
                    real L = -1, T = 0;
                    if ( !extract(arg, sci, L) || L < 0 )
                        throw InvalidParameter("length L must be >= 0 in `line L`");
                    if ( extract(arg, sci, T) && T < 0 )
                        throw InvalidParameter("thickness T must be >= 0 in `line L T`");
#if ( DIM >= 3 )
                    const Vector2 V = Vector2::randB(T);
                    return Vector(L*RNG.shalf(), V.XX, V.YY);
#endif
                    return Vector(L*RNG.shalf(), T*RNG.shalf(), 0);
                }
                break;

            case 'o':
                
                if ( spc && tok == "outside_sphere" )
                {
                    real R = 0;
                    extract(arg, sci, R);
                    if ( R < 0 )
                        throw InvalidParameter("distance R must be >= 0 in `outside_sphere R`");
                    Vector P;
                    do
                        P = spc->place();
                    while ( P.norm() < R );
                    return P;
                }
                break;
                
            case 'r':
                
                if ( spc && tok == "random" )
                    return spc->place();

                if ( tok == "rectangle" )
                {
                    Vector S(0, 0, 0);
                    extract(arg, sci, S);
                    return S.e_mul(Vector::randH());
                }

                if ( tok.compare(0, 4, "ring") == 0 )
                {
                    real L = -1, R = -1, T = 0;
                    if ( !extract(arg, sci, L) || L < 0 )
                        throw InvalidParameter("length L must be >= 0 in `ring L R`");
                    if ( !extract(arg, sci, R) || R < 0 )
                        throw InvalidParameter("radius R must be >= 0 in `ring L R`");
                    if ( extract(arg, sci, T) && T < 0 )
                        throw InvalidParameter("thickness T must be >= 0 in `ring L R T`");
                    real C, S;
                    RNG.urand2(C, S, R * ( 1.0 + RNG.shalf()*T ));
                    switch ( get_axis(tok, 4) )
                    {
                        case 'X': return Vector(L*RNG.shalf(), C, S);
                        case 'Y': return Vector(C, L*RNG.shalf(), S);
                        default : return Vector(C, S, L*RNG.shalf());
                    }
                }
                break;

            case 's':
                
                if ( spc && tok == "surface" )
                {
                    real R = 1;
                    extract(arg, sci, R);
                    if ( R < REAL_EPSILON )
                        throw InvalidParameter("distance R must be > 0 in `surface R`");
                    return spc->placeOnEdge(R);
                }

                if ( spc && tok == "stripe" )
                {
                    real S = -0.5, E = 0.5;
                    extract(arg, sci, S);
                    extract(arg, sci, E);
                    Vector inf, sup;
                    spc->boundaries(inf, sup);
                    Vector pos = inf + (sup-inf).e_mul(Vector::randP());
                    pos.XX = RNG.real_uniform(S, E);
                    return pos;
                }
                
                if ( tok == "sphere" )
                {
                    real R = -1, T = 0;
                    if ( !extract(arg, sci, R) || R < 0 )
                        throw InvalidParameter("radius R must be >= 0 in `sphere R`");
                    if ( extract(arg, sci, T) && T < 0 )
                        throw InvalidParameter("thickness T must be >= 0 in `sphere R T`");
                    return Vector::randU(R) + Vector::randU(T*0.5);
                }
                
                if ( tok == "square" )
                {
                    real R = 1;
                    extract(arg, sci, R);
                    return Vector(R*RNG.shalf(), R*RNG.shalf(), 0);
                }
#if ( 1 )
                /// A contribution from Beat Rupp
                if ( tok == "segment" )
                {
                    real B = 0, L = 0, T = 0, R = 0;
                    extract(arg, sci, B);
                    extract(arg, sci, L);
                    extract(arg, sci, T);
                    extract(arg, sci, R);
                    real x = T * RNG.shalf();
                    real y = L * RNG.preal();
                    if ( B > 0 ) {
                        real radius = L / (B * M_PI);
                        real inner = radius -T/2.0;
                        real theta = abs_real( L / radius );
                        real angle = RNG.preal() * theta;
                        // substract R to have the arc start from 0,0:
                        x = (inner + T * RNG.preal()) * std::cos(angle) - radius;
                        y = (inner + T * RNG.preal()) * std::sin(angle);
                    }
                    real C = std::cos(R);
                    real S = std::sin(R);
                    return Vector(C*x + S*y , -S*x + C*y, 0 );
                }
#endif
                break;
                
            case 'X':
                if ( spc && tok == "XY" )
                {
                    Vector V = spc->place();
                    real H = 0;
                    extract(arg, sci, H);
#if ( DIM > 2 )
                    V.ZZ = H;
#endif
                    return V;
                }
#if ( DIM > 1 )
                if ( spc && tok == "XZ" )
                {
                    real H = 0;
                    extract(arg, sci, H);
                    Vector V = spc->place();
                    V.YY = H;
                    return V;
                }
#endif
                break;
                
            case 'Y':
                if ( spc && tok == "YZ")
                {
                    real H = 0;
                    extract(arg, sci, H);
                    Vector V = spc->place();
                    V.XX = H;
                    return V;
                }
                break;
                
        }
        if ( spc )
            throw InvalidParameter("Unknown position specification `"+tok+"'");
        else
            throw InvalidParameter("Unknown position specification `"+tok+"' (with no space defined)");
    }

    // accept a Vector:
    Vector vec(0,0,0);
    
    if ( c == 0 || extract(arg, sci, vec) )
        return vec;

    throw InvalidParameter("cannot extract vector from `"+arg.substr(sci)+"`");
}


/**
 Modify a vector, according to further instructions from `is`.
 This applies `at`, `blur`, `if`, `or`...
*/
void Cytosim::modifyPosition(std::string const& arg, size_t& sci, Space const* spc, Vector& pos)
{
    while ( arg[sci] )
    {
        size_t isp = sci;
        std::string tok = Tokenizer::get_symbol(arg, sci);
       
        // Translation is specified with 'at' or 'move'
        if ( tok == "at" )
        {
            Vector vec(0,0,0);
            extract(arg, sci, vec);
            pos += vec;
        }
#if ( DIM == 3 )
        // translate in Z
        else if ( tok == "addZ" )
        {
            real Z = 0;
            extract(arg, sci, Z);
            pos.ZZ += Z;
        }
#endif
        // Convolve with shape
        else if ( tok == "add" )
        {
            Vector vec = readPositionPrimitive(arg, sci, spc);
            pos += vec;
        }
        // Alignment with a vector is specified with 'align'
        else if ( tok == "align" )
        {
            Vector vec = readDirection(arg, sci, pos, spc);
            Rotation rot = Rotation::randomRotationToVector(vec);
            pos = rot * pos;
        }
        // Rotation is specified with 'turn'
        else if ( tok == "turn" )
        {
            Rotation rot = readOrientation(arg, sci, pos, spc);
            pos = rot * pos;
        }
        // apply central symmetry with 'flip'
        else if ( tok == "flip" )
        {
            pos = -pos;
        }
        // Gaussian noise specified with 'blur'
        else if ( tok == "blur" )
        {
            real blur = 0;
            extract(arg, sci, blur);
            pos += Vector::randG(blur);
        }
        // extend along one of the main axis
        else if ( tok.compare(0, 6, "extend") == 0 )
        {
            real B = 0, T = 0;
            extract(arg, sci, B);
            extract(arg, sci, T);
#if ( DIM >= 3 )
            switch ( get_axis(tok, 6) )
            {
                case 'Z': pos.ZZ += B + ( T - B ) * RNG.preal(); break;
                case 'Y': pos.YY += B + ( T - B ) * RNG.preal(); break;
                default : pos.XX += B + ( T - B ) * RNG.preal(); break;

            }
#else
            pos.XX += B + ( T - B ) * RNG.preal();
#endif
        }
        // returns a random position between the two points specified
        else if ( tok == "to" )
        {
            Vector vec = readPositionPrimitive(arg, sci, spc);
            pos += ( vec - pos ) * RNG.preal();
            return;
        }
        // returns one of the two points specified
        else if ( tok == "or" )
        {
            Vector alt = readPositionPrimitive(arg, sci, spc);
            if ( RNG.flip() ) pos = alt;
        }
        else if ( tok == "if" )
        {
            int c = skip_space(arg, sci, false);
            if ( c == '(' )
                tok = Tokenizer::get_blocked_text(arg, ++sci, '(', ')');
            else
                throw InvalidParameter("expected condition block `(...)' after if");
            constexpr real alpha = 1.0 / ( DIM==3? M_SQRT3: M_SQRT2 );
            real S = pos.e_sum() * alpha;
            Evaluator evaluator{{"X", pos.x()}, {"Y", pos.y()}, {"Z", pos.z()}, {"S", S}};
            try {
                if ( 0 == evaluator.eval(tok) )
                {
                    pos.set(nan(""), nan(""), nan(""));
                    return;
                }
            }
            catch( Exception& e ) {
                e.message(e.message()+" in `"+tok+"'");
                throw;
            }
        }
        else
        {
            // unget last token:
            sci = isp;
            //throw InvalidParameter("unexpected `"+tok+"'");
            break;
        }
    }
}


//------------------------------------------------------------------------------
/**
 A position is defined with a SHAPE followed by a number of TRANSFORMATION.
 
 TRANSFORMATION         | Result                                               |
 -----------------------|-------------------------------------------------------
 `at X Y Z`             | Translate by specified vector (X,Y,Z)
 `add SHAPE`            | Translate by a vector chosen according to SHAPE
 `align VECTOR`         | Rotate to align parallel with specified vector
 `turn ROTATION`        | Apply specified rotation
 `extend B T`           | Extend along the Z axis, between Z=B and Z=T
 `blur REAL`            | Add centered Gaussian noise of variance REAL
 `to X Y Z`             | Interpolate with the previously specified position
 `or POSITION`          | flip randomly between two specified positions
 
 A vector is set according to SHAPE, and the transformations are applied one after
 the other, in the order in which they were given.\n

 Examples:
 
   position = 1 0 0
   position = circle 3 at 1 0
   position = square 3 align 1 1 0 at 1 1
 
 */
Vector Cytosim::readPosition(std::string const& arg, size_t& sci, Space const* spc)
{
    Vector vec = readPositionPrimitive(arg, sci, spc);
    assert_true( vec.valid() );
    if ( arg[sci] )
    {
        int c = skip_space(arg, sci, false);
        if ( isalpha(c) )
            modifyPosition(arg, sci, spc, vec);
    }
    return vec;
}


/// convert string to a position
Vector Cytosim::readPosition(std::string const& arg, Space const* spc)
{
    size_t sci = 0;
    Vector vec = Cytosim::readPosition(arg, sci, spc);
    size_t t = has_trail(arg, sci);
    if ( t )
        throw InvalidSyntax("unexpected trailing `"+arg.substr(t)+"' in position `"+arg+"'");
    return vec;
}


/// convert string to a position, bailing out after many trials
Vector Cytosim::findPosition(std::string const& arg, Space const* spc)
{
    long max_trials = 1 << 14;
    while ( --max_trials >= 0 )
    {
        size_t sci = 0;
        Vector vec = Cytosim::readPosition(arg, sci, spc);
        if ( vec.valid() )
        {
            size_t t = has_trail(arg, sci);
            if ( t )
                throw InvalidSyntax("unexpected trailing `"+arg.substr(t)+"' in position `"+arg+"'");
            return vec;
        }
    }
    throw InvalidParameter("failed to determine Vector from `"+arg+"'");
    return Vector(0,0,0);
}


//------------------------------------------------------------------------------
#pragma mark - Direction

/**
 Reads a direction which is a unit vector (norm = 1):
 
 Keyword                                     | Resulting Vector                                          |
 --------------------------------------------|------------------------------------------------------------
 `REAL REAL REAL`                            | the vector of norm 1 co-aligned with given vector
 `parallel REAL REAL REAL`                   | one of the two vectors of norm 1 parallel with given vector
 `orthogonal REAL REAL REAL`                 | a vector of norm 1 perpendicular to the given vector
 `horizontal` \n `parallel X`                | (+1,0,0) or (-1,0,0), randomly chosen with equal chance
 `vertical`\n `parallel Y`                   | (0,+1,0) or (0,-1,0), randomly chosen with equal chance
 `parallel Z`                                | (0,0,+1) or (0,0,-1), randomly chosen with equal chance
 `parallel XY`\n`parallel XZ`\n`parallel YZ` | A random vector in the specified plane
 `radial`                                    | directed from the origin to the current point
 `antiradial`                                | directed from the current point to the origin
 `circular`                                  | perpendicular to axis joining the current point to the origin
 `or DIRECTION`                              | flip randomly between two specified directions

 
 If a Space is defined, one may also use:
 
 Keyword         | Resulting Vector                       |
 ----------------|-----------------------------------------
 `tangent`       | parallel to the surface of the Space
 `normal`        | perpendicular to the surface
 `inward`        | normal to the surface, directed outward
 `outward`       | normal to the surface, directed inward


 Note: when the rotation is not uniquely determined in 3D (eg. `horizontal`), 
 cytosim will pick uniformly among all the possible rotations that fulfill the requirements.
 */

Vector Cytosim::readDirectionPrimitive(std::string const& arg, size_t& sci, Vector const& pos, Space const* spc)
{
    int c = skip_space(arg, sci, false);
    
    if ( c == 0 )
        return Vector::randU();

    if ( isalpha(c) )
    {
        const std::string tok = Tokenizer::get_symbol(arg, sci);
        
        if ( tok == "random" )
            return Vector::randU();

        if ( tok == "X" )
            return Vector(RNG.sflip(), 0, 0);
        if ( tok == "Y" )
            return Vector(0, RNG.sflip(), 0);
        if ( tok == "XY" )
        {
            real C, S;
            RNG.urand2(C, S);
            return Vector(C, S, 0);
        }
#if ( DIM >= 3 )
        if ( tok == "Z" )
            return Vector(0, 0, RNG.sflip());
        if ( tok == "XZ" )
        {
            real C, S;
            RNG.urand2(C, S);
            return Vector(C, 0, S);
        }
        if ( tok == "YZ" )
        {
            real C, S;
            RNG.urand2(C, S);
            return Vector(0, C, S);
        }
#endif

        if ( tok == "align" )
        {
            Vector vec;
            if ( extract(arg, sci, vec) )
                return vec.normalized(RNG.sflip());
            throw InvalidParameter("expected vector after `align`");
        }
        
        if ( tok == "parallel" )
        {
            Vector vec;
            if ( extract(arg, sci, vec) )
                return normalize(vec);
            throw InvalidParameter("expected vector after `parallel`");
        }

        if ( tok == "orthogonal" )
        {
            Vector vec;
            if ( extract(arg, sci, vec) )
                return normalize(vec).randOrthoU(1.0);
            throw InvalidParameter("expected vector after `orthogonal`");
        }

        if ( tok == "horizontal" )
        {
#if ( DIM >= 3 )
            real C, S;
            RNG.urand2(C, S);
            return Vector(C, S, 0);
#else
            return Vector(RNG.sflip(), 0, 0);
#endif
        }
        
        if ( tok == "vertical" )
        {
#if ( DIM >= 3 )
            return Vector(0, 0, RNG.sflip());
#else
            return Vector(0, RNG.sflip(), 0);
#endif
        }
        
        if ( tok == "radial" )
            return normalize(pos);
        
        if ( tok == "antiradial" )
            return -normalize(pos);

        if ( tok == "circular" )
            return normalize(pos).randOrthoU(1.0);
        
        if ( tok == "orthoradial" )
            return normalize(pos).randOrthoU(1.0);
        
        if ( tok == "inwardX" )
            return Vector(sign_real(pos.XX), 0, 0);
        
        if ( tok == "outwardX" )
            return Vector(-sign_real(pos.XX), 0, 0);

        if ( spc )
        {
            if ( tok == "tangent" )
            {
                size_t cnt = 0;
                Vector dir;
                do {
                    dir = Vector::randU();
                    dir = spc->project(pos+dir) - spc->project(pos);
                    if ( ++cnt > 128 )
                    {
                        printf("warning: tangent placement(%9.3f %9.3f %9.3f) failed\n", pos.x(), pos.y(), pos.z());
                        return Vector::randU();
                    }
                } while ( dir.normSqr() < REAL_EPSILON );
                return dir.normalized();
            }
            
#if ( DIM >= 3 )
            if ( tok == "clockwise" )
            {
                real ang = 0;
                extract(arg, sci, ang);
                Vector out = spc->normalToEdge(pos);
                Vector tan = cross(Vector(0,0,1), out);
                real C = std::cos(ang), S = std::sin(ang);
                Vector dir(C*tan.XX, C*tan.YY, S);
                real n = dir.norm();
                if ( n > REAL_EPSILON )
                    return dir.normalized();
                return out.randOrthoU(1.0);
            }
            
            if ( tok == "anticlockwise" )
            {
                real ang = 0;
                extract(arg, sci, ang);
                Vector out = spc->normalToEdge(pos);
                Vector tan = cross(Vector(0,0,-1), out);
                real C = std::cos(ang), S = std::sin(ang);
                Vector dir(C*tan.XX, C*tan.YY, S);
                real n = dir.norm();
                if ( n > REAL_EPSILON )
                    return dir.normalized();
                return out.randOrthoU(1.0);
            }
#elif ( DIM == 2 )
            if ( tok == "clockwise" )
                return cross(+1, spc->normalToEdge(pos));
            
            if ( tok == "anticlockwise" )
                return cross(-1, spc->normalToEdge(pos));
#endif
            
            if ( tok == "normal" )
                return RNG.sflip() * spc->normalToEdge(pos);
            
            if ( tok == "inward" )
                return -spc->normalToEdge(pos);
            
            if ( tok == "outward" )
                return spc->normalToEdge(pos);
        }
        
        throw InvalidParameter("Unknown direction specification `"+tok+"'");
    }
    
    // accept a Vector:
    Vector vec(0,0,0);
    if ( extract(arg, sci, vec) )
    {
        real n = vec.norm();
        if ( n < REAL_EPSILON )
            throw InvalidParameter("direction vector appears singular");
        return vec / n;
    }

    throw InvalidParameter("cannot extract direction from `"+arg+"`");
}


Vector Cytosim::readDirection(std::string const& arg, size_t& sci, Vector const& pos, Space const* spc)
{
    size_t ouf = 0, max_trials = 1<<14;
    std::string tok;
    Vector dir(1,0,0);
    size_t isp, start = sci;
    
restart:
    sci = start;
    dir = readDirectionPrimitive(arg, sci, pos, spc);
    
    while ( arg[sci] )
    {
        isp = sci;
        tok = Tokenizer::get_symbol(arg, sci);
        
        if ( tok.empty() )
            return dir;
        
        // Gaussian noise specified with 'blur'
        else if ( tok == "blur" )
        {
            real blur = 0;
            extract(arg, sci, blur);
#if ( DIM == 3 )
            Vector3 X, Y;
            dir.orthonormal(X, Y);
            X *= blur*RNG.gauss();
            Y *= blur*RNG.gauss();
            dir = normalize(dir+X+Y);
#elif ( DIM == 2 )
            dir = Rotation::rotation(blur*RNG.gauss()) * dir;
#endif
        }
        // returns one of the two points specified
        else if ( tok == "or" )
        {
            Vector alt = readDirectionPrimitive(arg, sci, pos, spc);
            if ( RNG.flip() ) dir = alt;
        }
        else if ( tok == "if" )
        {
            int c = skip_space(arg, sci, false);
            if ( c == '(' )
                tok = Tokenizer::get_blocked_text(arg, ++sci, '(', ')');
            else
                throw InvalidParameter("expected condition block `(...)' after if");
            Evaluator evaluator{{"X", dir.x()}, {"Y", dir.y()}, {"Z", dir.z()}};
            try {
                if ( 0 == evaluator.eval(tok) )
                {
                    if ( ++ouf < max_trials )
                        goto restart;
                    throw InvalidParameter("condition `"+tok+"' could not be fulfilled");
                    break;
                }
            }
            catch( Exception& e ) {
                e.message(e.message()+" in `"+tok+"'");
                throw;
            }
        }
        else
        {
            // unget last token
            sci = isp;
            break;
        }
    }
    return dir;
}


Vector Cytosim::readDirection(std::string const& arg, Vector const& pos, Space const* spc)
{
    size_t sci = 0;
    Vector vec(0, 0, 0);
    try {
        vec = Cytosim::readDirection(arg, sci, pos, spc);
    }
    catch ( Exception& e )
    {
        throw InvalidSyntax("could not determine direction from `"+arg+"'");
    }
    size_t t = has_trail(arg, sci);
    if ( t )
       throw InvalidSyntax("discarded `"+arg.substr(t)+"' in direction");
    return vec;
}

/**
 A rotation can be specified in 3D as follows:
 
 Keyword                 | Rotation / Result
 ------------------------|-----------------------------------------------------------
 `random`                | A rotation selected uniformly among all possible rotations
 `identity`, `off`       | The object is not rotated
 `X theta`               | A rotation around axis 'X' of angle `theta` in radians
 `Y theta`               | A rotation around axis 'Y' of angle `theta` in radians
 `Z theta`               | A rotation around axis 'Z' of angle `theta` in radians
 `angle A`               | Rotation around Z axis with angle of rotation in radians
 `angle A axis X Y Z`    | Rotation around given axis with angle of rotation in radians
 `axis X Y Z angle A`    | Rotation around given axis with angle of rotation in radians
 `degree A axis X Y Z`   | As specified by axis and angle of rotation in degrees
 `quat q0 q1 q2 q3`      | As specified by the Quaternion (q0, q1, q2, q3)
*/

Rotation Cytosim::readRotation(std::string const& arg, size_t& sci)
{
    std::string tok = Tokenizer::get_symbol(arg, sci);
    
    if ( tok == "random" )
        return Rotation::randomRotation();
    else if ( tok == "off" || tok == "identity" )
        return Rotation::one();
    else if ( tok == "align111" )
        return Rotation::align111();
    else if ( tok == "angle" )
    {
        real A = 0;
        extract(arg, sci, A);
        Vector3 dir(0,0,1);
        if ( Tokenizer::has_symbol(arg, sci, "axis") )
            extract(arg, sci, dir);
#if ( DIM >= 3 )
        return Rotation::rotationAroundAxis(normalize(dir), std::cos(A), std::sin(A));
#else
        return Rotation::rotation(std::cos(A), std::sin(A));
#endif
    }
    else if ( tok == "axis" )
    {
        Vector3 dir(0,0,1);
        extract(arg, sci, dir);
        real ang = 0;
        if ( Tokenizer::has_symbol(arg, sci, "angle") )
            extract(arg, sci, ang);
        else if ( Tokenizer::has_symbol(arg, sci, "degree") )
        {
            extract(arg, sci, ang);
            ang *= M_PI/180.0;
        }
#if ( DIM >= 3 )
        return Rotation::rotationAroundAxis(normalize(dir), std::cos(ang), std::sin(ang));
#else
        return Rotation::rotation(std::cos(ang), std::sin(ang));
#endif
    }
    else if ( tok == "degree" )
    {
        real ang = 0;
        extract(arg, sci, ang);
        ang *= M_PI/180.0;
        Vector3 dir(0,0,1);
        if ( Tokenizer::has_symbol(arg, sci, "axis") )
            extract(arg, sci, dir);
#if ( DIM >= 3 )
        return Rotation::rotationAroundAxis(normalize(dir), std::cos(ang), std::sin(ang));
#else
        return Rotation::rotation(std::cos(ang), std::sin(ang));
#endif
    }
#if ( DIM >= 3 )
    else if ( tok == "quat" )
    {
        Quaternion<real> quat;
        extract(arg, sci, quat);
        quat.normalize();
        Rotation rot;
        quat.setMatrix3(rot);
        return rot;
    }
    else if ( tok == "X" )
        return Rotation::rotationAroundX(get_angle(arg, sci));
    else if ( tok == "Y" )
        return Rotation::rotationAroundY(get_angle(arg, sci));
    else if ( tok == "Z" )
        return Rotation::rotationAroundZ(get_angle(arg, sci));
    else if ( tok == "Xflip" )
        return Rotation::rotationAroundX(M_PI);
    else if ( tok == "Yflip" )
        return Rotation::rotationAroundY(M_PI);
    else if ( tok == "Zflip" )
        return Rotation::rotationAroundZ(M_PI);
#else
    else if ( tok == "X" )
        return Rotation(0, RNG.sflip());
    else if ( tok == "Z" )
        return Rotation::rotation(get_angle(arg, sci));
    else if ( tok == "Xflip" )
        return Rotation(0, RNG.sflip());
#endif
    
    if ( tok.size() )
        throw InvalidSyntax("unexpected token `"+tok+"' in rotation");
    return Rotation(0, 1);
}


Rotation Cytosim::readRotation(std::string const& arg)
{
    size_t sci = 0;
    Rotation rot(0, 1);
    try {
        rot = Cytosim::readRotation(arg, sci);
        // can combine another rotation:
        if ( arg[sci] )
            rot = Cytosim::readRotation(arg, sci) * rot;
    }
    catch ( Exception& e )
    {
        throw InvalidSyntax("could not determine rotation from `"+arg+"'");
    }
    size_t t = has_trail(arg, sci);
    if ( t )
       throw InvalidSyntax("discarded `"+arg.substr(t)+"' in rotation");
    return rot;
}

/**
 The initial orientation of objects is defined by a rotation, but it is usually
 sufficient to specify a unit vector:
 
 Keyword                 | Rotation / Result
 ------------------------|------------------------------------------------------
 ROTATION                | see @ref Cytosim::readRotation()
 DIRECTION               | see @ref Cytosim::readDirection
 DIRECTION or DIRECTION  | flip randomly between two specified directions
 
 When a DIRECTION is specified, a rotation will be built that transforms (1, 0, 0)
 into the given vector (after normalization). In 3D, this does not define a rotation
 uniquely, and cytosim will randomly pick one of the possible rotations, with equal
 probability among all the possible rotation, by rotating around (1, 0, 0) beforehand.
*/

Rotation Cytosim::readOrientation(std::string const& arg, size_t& sci, Vector const& pos, Space const* spc)
{
    int c = skip_space(arg, sci, false);

    if ( isalpha(c) )
    {
        try {
            return readRotation(arg, sci);
        }
        catch ( Exception& e )
        {
            // Some keywords are handled by readDirection(), so we just warn here
            std::cerr << "Warning, " << e.message() << '\n';
        }
    }

    // normally a unit vector is specified:
    Vector vec = readDirection(arg, sci, pos, spc);
    
    /*
     A single Vector does not uniquely define a rotation in 3D:
     hence we return a random rotation that is picked uniformly
     among all possible rotations transforming (1,0,0) in vec.
     */
    return Rotation::randomRotationToVector(vec);
}

