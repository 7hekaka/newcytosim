// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University.

#include "primitives.h"
#include "assert_macro.h"
#include "stream_func.h"
#include "exceptions.h"
#include "quaternion.h"
#include "iowrapper.h"
#include "tokenizer.h"
#include "evaluator.h"
#include "glossary.h"
#include "random.h"
#include "space.h"


/* With C++11, the extracted value is zeroed even upon failure.
This on the other hand will preserve the value of 'var' if no read occurs.
Returns `true` if a value was set. This is used with T=real and T=Vector */
template < typename T >
static bool extract(std::istream& is, T& var)
{
    T backup = var;
    std::streampos isp = is.tellg();
    if ( is >> var )
        return true;
    is.clear();
    is.seekg(isp);
    var = backup;
    return false;
}


/// extract a rotation angle specified in Radian
static real get_angle(std::istream& is)
{
    real a = 0;
    if ( is >> a )
        return a;
    // if no angle is specified, set it randomly:
    return RNG.sreal() * M_PI;
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
 `line L T`           | L = Length, T = thickness. Selected randomly with -L/2 < X < L/2; norm(Y,Z) < T
 `sphere R T`         | At distance R +/- T/2 from the origin\n `R-T/2 < norm(X,Y,Z) < R+T/2`
 `ball R`             | At distance R at most from the origin\n `norm(X,Y,Z) < R`
 `disc R T`           | in 2D, a disc in the XY-plane \n in 3D, a disc in the XY-plane of thickness T in Z
 `discXZ R T`         | Disc in the XZ-plane of radius R, thickness T
 `discYZ R T`         | Disc in the YZ-plane of radius R, thickness T
 `equator R T`        | At distance R from the origin, and T from the XY plane:\n `norm(X,Y) < R` `norm(Z) < T`
 `circle R T`         | Circle of radius R and thickness T \n At distance T from the circle of radius R
 `cylinder L R`       | Cylinder of axis X, L=length in X, R=radius in YZ plane
 `cylinderZ L R`      | Cylinder of axis Z, L=length in Z, R=radius in XY plane
 `ring L R T`         | Surface of a cylinder of axis X, L=length in X, R=radius in YZ, T = thickness
 `ellipse A B C`      | Inside the ellipse or ellipsoid of main axes 2A, 2B and 2C
 `arc L Theta`        | A piece of circle of length L and covering an angle Theta
 `stripe L R`         | Random vector with L < X < R
 `square R`           | Random vector with -R < X < R; -R < Y < R; -R < Z < R;
 `rectangle A B C`    | Random vector with -A < X < A; -B < Y < B; -C < Z < C;
 `gradient S E`       | Linear density gradient along X, of value 0 at X=S and 1 at X=E
 `gradient S E R`     | Linear density gradient, contained inside a cylinder of radius R
 `exponential S L`    | Exponential density gradient of length scale L, starting at S
 `exponential S L R`  | Exponential density gradient, contained inside a cylinder of radius R
 `XY Z`               | randomly in the XY plane and within the Space, at specified Z (a value)

 Each primitive describes a certain area in Space, and in most cases the returned position is
 chosen randomly inside this area following a uniform probability.
 */

Vector Cytosim::readPositionPrimitive(std::istream& is, Space const* spc)
{
    int c = Tokenizer::skip_space(is, false);

    if ( c == EOF )
        return Vector(0,0,0);

    if ( isalpha(c) )
    {
        std::string tok = Tokenizer::get_symbol(is);
        //StreamFunc::mark_line(std::cerr, is);

        if ( spc )
        {
            if ( tok == "inside" || tok == "random" )
                return spc->place();
            
            if ( tok == "XY" )
            {
                Vector V = spc->place();
                real H = 0;
                extract(is, H);
#if ( DIM > 2 )
                V.ZZ = H;
#endif
                return V;
            }
            if ( tok == "YZ" || tok == "XZ" )
            {
                real H = 0;
                extract(is, H);
                Vector V = spc->place();
                if ( tok == "YZ" ) V.XX = H;
#if ( DIM > 1 )
                if ( tok == "XZ" ) V.YY = H;
#endif
                return V;
            }

            if ( tok == "edge" )
            {
                real R = 0;
                if ( !extract(is, R) || R < REAL_EPSILON )
                    throw InvalidParameter("distance R must be > 0 in `edge R`");
                return spc->placeNearEdge(R);
            }
            
            if ( tok == "surface" )
            {
                real R = 1;
                extract(is, R);
                if ( R < REAL_EPSILON )
                    throw InvalidParameter("distance R must be > 0 in `surface R`");
                return spc->placeOnEdge(R);
            }

            if ( tok == "outside_sphere" )
            {
                real R = 0;
                extract(is, R);
                if ( R < 0 )
                    throw InvalidParameter("distance R must be >= 0 in `outside_sphere R`");
                Vector P;
                do
                    P = spc->place();
                while ( P.norm() < R );
                return P;
            }
            
            if ( tok == "stripe" )
            {
                real S = -0.5, E = 0.5;
                extract(is, S);
                extract(is, E);
                Vector inf, sup;
                spc->boundaries(inf, sup);
                Vector pos = inf + (sup-inf).e_mul(Vector::randP());
                pos.XX = RNG.real_uniform(S, E);
                return pos;
            }
        
            if ( tok == "gradient" )
            {
                real S = -10, E = 10, R = 0;
                extract(is, S);
                extract(is, E);
                extract(is, R);
                if ( R == 0 )
                {
                    Vector vec;
                    real p;
                    do {
                        vec = spc->place();
                        p = ( vec.XX - S ) / ( E - S );
                    } while ( p < 0 || p > 1 || p < RNG.preal() );
                    return vec;
                }
                real x = std::sqrt(RNG.preal());
#if ( DIM < 3 )
                return Vector(S+(E-S)*x, R*RNG.sreal(), 0);
#else
                const Vector2 V = Vector2::randU();
                return Vector(S+(E-S)*x, R*V.XX, R*V.YY);
#endif
            }
        
            if ( tok == "exponential" )
            {
                real S = -10, E = 1, R = 0;
                extract(is, S);
                extract(is, E);
                extract(is, R);
                if ( R == 0 )
                {
                    Vector vec;
                    real p;
                    do {
                        vec = spc->place();
                        p = std::exp( ( S - vec.XX ) / E );
                    } while ( p < 0 || p > 1 || p < RNG.preal() );
                    return vec;
                }
                real x = RNG.exponential();
#if ( DIM < 3 )
                return Vector(S+E*x, R*RNG.sreal(), 0);
#else
                const Vector2 V = Vector2::randU();
                return Vector(S+E*x, R*V.XX, R*V.YY);
#endif
            }
        }
        
        if ( tok == "sphere" )
        {
            real R = -1, T = 0;
            if ( !extract(is, R) || R < 0 )
                throw InvalidParameter("radius R must be >= 0 in `sphere R`");
            if ( extract(is, T) && T < 0 )
                throw InvalidParameter("thickness T must be >= 0 in `sphere R T`");
            return Vector::randU(R) + Vector::randU(T*0.5);
        }
        
        if ( tok == "ball" )
        {
            real R = -1;
            if ( !extract(is, R) || R < 0 )
                throw InvalidParameter("radius R must be >= 0 in `ball R`");
            return Vector::randB(R);
        }
        
        if ( tok == "equator" )
        {
            real R = 0, T = 0;
            if ( extract(is, R) && R < 0 )
                throw InvalidParameter("radius R must be >= 0 in `equator R T`");
            if ( extract(is, T) && T < 0 )
                throw InvalidParameter("thickness T must be >= 0 in `equator R T`");
            const Vector2 V = Vector2::randU();
            return Vector(R*V.XX, R*V.YY, T*RNG.shalf());
        }
        
        if ( tok == "cap" )
        {
            real R = 0, T = 0;
            if ( !extract(is, R) || R < 0 )
                throw InvalidParameter("radius R must be >= 0 in `cap R`");
            if ( extract(is, T) && T < 0 )
                throw InvalidParameter("thickness T must be >= 0 in `cap R T`");
            real Z = std::max(R - T * RNG.preal(), -R);
            real r = std::sqrt(R*R - Z*Z);
#if ( DIM >= 3 )
            const Vector2 XY = Vector2::randU();
            return Vector(r*XY.XX, r*XY.YY, Z);
#else
            return Vector(r*RNG.sflip(), Z, 0);
#endif
        }

        if ( tok.compare(0, 8, "cylinder") == 0 )
        {
            real L = -1, R = -1;
            if ( !extract(is, L) || L < 0 )
                throw InvalidParameter("length L must be >= 0 in `cylinder L R`");
            if ( !extract(is, R) || R < 0 )
                throw InvalidParameter("radius R must be >= 0 in `cylinder L R`");
            const Vector2 V = Vector2::randB(R);
            std::string axis = tok.substr(8);
            if ( axis == "Z" )
                return Vector(V.XX, V.YY, L*RNG.shalf());
            else if ( axis == "Y" )
                return Vector(V.XX, L*RNG.shalf(), V.YY);
            return Vector(L*RNG.shalf(), V.XX, V.YY);
        }
        
        if ( tok.compare(0, 4, "ring") == 0 )
        {
            real L = -1, R = -1, T = 0;
            if ( !extract(is, L) || L < 0 )
                throw InvalidParameter("length L must be >= 0 in `ring L R`");
            if ( !extract(is, R) || R < 0 )
                throw InvalidParameter("radius R must be >= 0 in `ring L R`");
            if ( extract(is, T) && T < 0 )
                throw InvalidParameter("thickness T must be >= 0 in `ring L R T`");
            const Vector2 V = Vector2::randU(R) * ( 1.0 + RNG.shalf()*T );
            std::string axis = tok.substr(4);
            if ( axis == "Z" )
                return Vector(V.XX, V.YY, L*RNG.shalf());
            else if ( axis == "Y" )
                return Vector(V.XX, L*RNG.shalf(), V.YY);
            return Vector(L*RNG.shalf(), V.XX, V.YY);
        }

        if ( tok.compare(0, 6, "circle") == 0 )
        {
            real R = -1, T = 0;
            if ( !extract(is, R) || R < 0 )
                throw InvalidParameter("radius R must be >= 0 in `circle R T`");
            if ( extract(is, T) && T < 0 )
                throw InvalidParameter("thickness T must be >= 0 in `circle R T`");
#if ( DIM >= 3 )
            Vector2 V = Vector2::randU(R);
            Vector3 W = (0.5*T) * Vector3::randU();
            std::string axis = tok.substr(6);
            if ( axis == "X" )
                return Vector3(0, V.XX, V.YY) + W;
            else if ( axis == "Y" )
                return Vector3(V.XX, 0, V.YY) + W;
            else
                return Vector3(V.XX, V.YY, 0) + W;
#endif
            return Vector::randU(R) + Vector::randU(T*0.5);
        }
        
        if ( tok.compare(0, 4, "disc") == 0 )
        {
            real R = -1, T = 0;
            if ( !extract(is, R) || R < 0 )
                throw InvalidParameter("radius R must be >= 0 in `disc R`");
            if ( extract(is, T) && T < 0 )
                throw InvalidParameter("thickness T must be >= 0 in `disc R T`");
#if ( DIM >= 3 )
            //in 3D, a disc in the XY-plane of thickness T in Z-direction
            Vector2 V = Vector2::randB(R);
            std::string axis = tok.substr(4);
            if ( axis == "X" )
                return Vector(T*RNG.shalf(), V.XX, V.YY);
            else if ( axis == "Y" )
                return Vector(V.XX, T*RNG.shalf(), V.YY);
            else
                return Vector(V.XX, V.YY, T*RNG.shalf());
#endif
            //in 2D, a disc in the XY-plane
            return Vector::randB(R);
        }
        
        if ( tok == "ellipse" )
        {
            Vector S(1, 1, 0);
            extract(is, S);
            return S.e_mul(Vector::randB());
        }
        
        if ( tok == "ellipse_surface" )
        {
            Vector S(1, 1, 0);
            extract(is, S);
            return S.e_mul(Vector::randU());
        }
        
        if ( tok == "line" )
        {
            real L = -1, T = 0;
            if ( !extract(is, L) || L < 0 )
                throw InvalidParameter("length L must be >= 0 in `line L`");
            if ( extract(is, T) && T < 0 )
                throw InvalidParameter("thickness T must be >= 0 in `line L T`");
#if ( DIM >= 3 )
            const Vector2 V = Vector2::randB(T);
            return Vector(L*RNG.shalf(), V.XX, V.YY);
#endif
            return Vector(L*RNG.shalf(), T*RNG.shalf(), 0);
        }
        
        if ( tok == "arc" )
        {
            real L = -1, A = 1.57;
            if ( !extract(is, L) || L < 0 )
                throw InvalidParameter("length L must be >= 0 in `arc L`");
            extract(is, A);
            
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
        
        if ( tok == "square" )
        {
            real x = 1;
            extract(is, x);
            return Vector::randS(x);
        }
        
        if ( tok == "rectangle" )
        {
            Vector S(0, 0, 0);
            extract(is, S);
            return S.e_mul(Vector::randH());
        }

#if ( 1 )
        /// A contribution from Beat Rupp
        if ( tok == "segment" || tok == "newsegment" )
        {
            real B = 0, L = 0, T = 0, R = 0;
            extract(is, B);
            extract(is, L);
            extract(is, T);
            extract(is, R);
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
        if ( tok == "center" || tok == "origin" )
            return Vector(0,0,0);
        
        if ( spc )
            throw InvalidParameter("Unknown position specification `"+tok+"'");
        else
            throw InvalidParameter("Unknown position specification `"+tok+"' (with no space defined)");
    }
    
    // accept a Vector:
    Vector vec(0,0,0);
    if ( extract(is, vec) )
        return vec;
    throw InvalidParameter("expected a vector specifying a `position`");
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
Vector Cytosim::readPosition(std::istream& is, Space const* spc)
{
    Vector pos = readPositionPrimitive(is, spc);
    assert_true( pos.valid() );
    is.clear();
    int c = Tokenizer::skip_space(is, false);
    if ( isalpha(c) )
        return modifyPosition(is, spc, pos);
    return pos;
}


Vector Cytosim::modifyPosition(std::istream& is, Space const* spc, Vector pos)
{
    std::string tok;
    std::streampos isp;
    
    while ( !is.eof() )
    {
        isp = is.tellg();
        tok = Tokenizer::get_symbol(is);
        //StreamFunc::mark_line(std::cerr, is);

        if ( tok.empty() )
            return pos;
        
        // Translation is specified with 'at' or 'move'
        if ( tok == "at"  ||  tok == "move" )
        {
            Vector vec(0,0,0);
            extract(is, vec);
            pos = pos + vec;
        }
        // Convolve with shape
        else if ( tok == "add" )
        {
            Vector vec = readPositionPrimitive(is, spc);
            pos = pos + vec;
        }
        // Alignment with a vector is specified with 'align'
        else if ( tok == "align" )
        {
            Vector vec = readDirection(is, pos, spc);
            Rotation rot = Rotation::randomRotationToVector(vec);
            pos = rot * pos;
        }
        // Rotation is specified with 'turn'
        else if ( tok == "turn" )
        {
            Rotation rot = readOrientation(is, pos, spc);
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
            extract(is, blur);
            pos += Vector::randG(blur);
        }
        // extend long the X axis
        else if ( tok == "extend" )
        {
            real B = 0, T = 0;
            extract(is, B);
            extract(is, T);
            pos.XX += B + ( T - B ) * RNG.preal();
        }
#if ( DIM > 2 )
        // extend long the Z axis
        else if ( tok == "extendZ" )
        {
            real B = 0, T = 0;
            extract(is, B);
            extract(is, T);
            pos.ZZ += B + ( T - B ) * RNG.preal();
        }
#endif
        // returns a random position between the two points specified
        else if ( tok == "to" )
        {
            Vector vec = readPositionPrimitive(is, spc);
            return pos + ( vec - pos ) * RNG.preal();
        }
        // returns one of the two points specified
        else if ( tok == "or" )
        {
            Vector alt = readPositionPrimitive(is, spc);
            if ( RNG.flip() ) pos = alt;
        }
        else if ( tok == "if" )
        {
            tok = Tokenizer::get_token(is);
            real S = pos.e_sum()/std::sqrt(DIM);
            Evaluator evaluator{{"X", pos.x()}, {"Y", pos.y()}, {"Z", pos.z()}, {"S", S}};
            try {
                if ( 0 == evaluator.eval(tok) )
                    return Vector(nan(""), nan(""), nan(""));
            }
            catch( Exception& e ) {
                e.message(e.message()+" in `"+tok+"'");
                throw;
            }
        }
        else
        {
            // unget last token:
            is.clear();
            is.seekg(isp);
#if 0
            /*
             We need to work around a bug in the stream extraction operator,
             which eats extra characters ('a','n','e','E') if a double is read
             19.10.2015
             */
            is.seekg(-1, std::ios_base::cur);
            int c = is.peek();
            if ( c=='a' || c=='b' )
                continue;
            is.seekg(1, std::ios_base::cur);
#endif
            //throw InvalidParameter("unexpected `"+tok+"'");
            break;
        }
    }
    return pos;
}


/// convert string to a position
Vector Cytosim::readPosition(std::string const& arg, Space const* spc)
{
    std::istringstream iss(arg);
    Vector vec = Cytosim::readPosition(iss, spc);
    if ( StreamFunc::has_trail(iss) )
    {
        std::string str;
        std::getline(iss, str);
        throw InvalidSyntax("unexpected trailing `"+str+"' in position `"+arg+"'");
    }
    return vec;
}


/// convert string to a position, bailing out after many trials
Vector Cytosim::findPosition(std::string const& arg, Space const* spc)
{
    long max_trials = 1 << 14;
    while ( --max_trials >= 0 )
    {
        std::istringstream iss(arg);
        Vector vec = Cytosim::readPosition(iss, spc);
        if ( vec.valid() )
        {
            if ( StreamFunc::has_trail(iss) )
            {
                std::string str;
                std::getline(iss, str);
                throw InvalidSyntax("unexpected trailing `"+str+"' in position `"+arg+"'");
            }
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

Vector Cytosim::readDirectionPrimitive(std::istream& is, Vector const& pos, Space const* spc)
{
    int c = Tokenizer::skip_space(is, false);
    
    if ( c == EOF )
        return Vector::randU();

    if ( isalpha(c) )
    {
        const std::string tok = Tokenizer::get_symbol(is);
        
        if ( tok == "random" )
            return Vector::randU();

        if ( tok == "X" )
            return Vector(RNG.sflip(), 0, 0);
        if ( tok == "Y" )
            return Vector(0, RNG.sflip(), 0);
        if ( tok == "XY" )
        {
            const Vector2 V = Vector2::randU();
            return Vector(V.XX, V.YY, 0);
        }
#if ( DIM >= 3 )
        if ( tok == "Z" )
            return Vector(0, 0, RNG.sflip());
        if ( tok == "XZ" )
        {
            const Vector2 V = Vector2::randU();
            return Vector(V.XX, 0, V.YY);
        }
        if ( tok == "YZ" )
        {
            const Vector2 V = Vector2::randU();
            return Vector(0, V.XX, V.YY);
        }
#endif

        if ( tok == "align" )
        {
            Vector vec;
            if ( extract(is, vec) )
                return vec.normalized(RNG.sflip());
            throw InvalidParameter("expected vector after `align`");
        }
        
        if ( tok == "parallel" )
        {
            Vector vec;
            if ( extract(is, vec) )
                return normalize(vec);
            throw InvalidParameter("expected vector after `parallel`");
        }

        if ( tok == "orthogonal" )
        {
            Vector vec;
            if ( extract(is, vec) )
                return vec.randOrthoU(1);
            throw InvalidParameter("expected vector after `orthogonal`");
        }

        if ( tok == "horizontal" )
        {
#if ( DIM >= 3 )
            return Vector(Vector2::randU());
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
            return pos.randOrthoU(1.0);
        
        if ( tok == "orthoradial" )
            return pos.randOrthoU(1.0);
        
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
                extract(is, ang);
                Vector vec = spc->normalToEdge(pos);
                Vector tan = cross(Vector(0,0,1), vec);
                real C = std::cos(ang), S = std::sin(ang);
                Vector dir(C*tan.XX, C*tan.YY, S);
                real n = dir.norm();
                if ( n > REAL_EPSILON )
                    return dir.normalized();
                return vec.randOrthoU(1.0);
            }
            
            if ( tok == "anticlockwise" )
            {
                real ang = 0;
                extract(is, ang);
                Vector vec = spc->normalToEdge(pos);
                Vector tan = cross(Vector(0,0,-1), vec);
                real C = std::cos(ang), S = std::sin(ang);
                Vector dir(C*tan.XX, C*tan.YY, S);
                real n = dir.norm();
                if ( n > REAL_EPSILON )
                    return dir.normalized();
                return vec.randOrthoU(1.0);
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
    
    {
        // accept a Vector:
        Vector vec(0,0,0);
        if ( extract(is, vec) )
        {
            real n = vec.norm();
            if ( n < REAL_EPSILON )
                throw InvalidParameter("direction vector appears singular");
            return vec / n;
        }
    }
    throw InvalidParameter("expected a vector specifying a `direction`");
}


Vector Cytosim::readDirection(std::istream& is, Vector const& pos, Space const* spc)
{
    size_t ouf = 0, max_trials = 1<<14;
    std::string tok;
    Vector dir(1,0,0);
    std::streampos isp, start = is.tellg();
    
restart:
    if ( is.fail() )
        return dir;
    dir = readDirectionPrimitive(is, pos, spc);
    is.clear();
    
    while ( !is.eof() )
    {
        isp = is.tellg();
        tok = Tokenizer::get_symbol(is);
        
        if ( tok.empty() )
            return dir;
        
        // Gaussian noise specified with 'blur'
        else if ( tok == "blur" )
        {
            real blur = 0;
            extract(is, blur);
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
            Vector alt = readDirectionPrimitive(is, pos, spc);
            if ( RNG.flip() ) dir = alt;
        }
        else if ( tok == "if" )
        {
            tok = Tokenizer::get_token(is);
            Evaluator evaluator{{"X", dir.x()}, {"Y", dir.y()}, {"Z", dir.z()}};
            try {
                if ( 0 == evaluator.eval(tok) )
                {
                    if ( ++ouf < max_trials )
                    {
                        is.seekg(start);
                        goto restart;
                    }
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
            is.clear();
            is.seekg(isp);
            break;
        }
    }
    return dir;
}


Vector Cytosim::readDirection(std::string const& arg, Vector const& pos, Space const* spc)
{
    std::istringstream iss(arg);
    Vector vec(0, 0, 0);
    try {
        vec = Cytosim::readDirection(iss, pos, spc);
    }
    catch ( Exception& e )
    {
        throw InvalidSyntax("could not determine direction from `"+arg+"'");
    }
    if ( StreamFunc::has_trail(iss) )
    {
        std::string sub = arg.substr(iss.tellg());
        throw InvalidSyntax("discarded `"+sub+"' in direction");
    }
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

Rotation Cytosim::readRotation(std::istream& is)
{
    std::string tok = Tokenizer::get_symbol(is);
    
    if ( tok == "random" )
        return Rotation::randomRotation();
    else if ( tok == "identity" || tok == "off" || tok == "none" )
        return Rotation::one();
    else if ( tok == "align111" )
        return Rotation::align111();
    else if ( tok == "angle" )
    {
        real A = 0;
        extract(is, A);
        Vector3 dir(0,0,1);
        if ( Tokenizer::has_symbol(is, "axis") )
            extract(is, dir);
#if ( DIM >= 3 )
        return Rotation::rotationAroundAxis(normalize(dir), std::cos(A), std::sin(A));
#else
        return Rotation::rotation(std::cos(A), std::sin(A));
#endif
    }
    else if ( tok == "axis" )
    {
        Vector3 dir(0,0,1);
        extract(is, dir);
        real ang = 0;
        if ( Tokenizer::has_symbol(is, "angle") )
            extract(is, ang);
        else if ( Tokenizer::has_symbol(is, "degree") )
        {
            extract(is, ang);
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
        extract(is, ang);
        ang *= M_PI/180.0;
        Vector3 dir(0,0,1);
        if ( Tokenizer::has_symbol(is, "axis") )
            extract(is, dir);
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
        extract(is, quat);
        quat.normalize();
        Rotation rot;
        quat.setMatrix3(rot);
        return rot;
    }
    else if ( tok == "X" )
        return Rotation::rotationAroundX(get_angle(is));
    else if ( tok == "Y" )
        return Rotation::rotationAroundY(get_angle(is));
    else if ( tok == "Z" )
        return Rotation::rotationAroundZ(get_angle(is));
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
        return Rotation::rotation(get_angle(is));
    else if ( tok == "Xflip" )
        return Rotation(0, RNG.sflip());
#endif
    
    if ( tok.size() )
        throw InvalidSyntax("unexpected token in rotation");
    return Rotation(0, 1);
}


Rotation Cytosim::readRotation(std::string const& arg)
{
    std::istringstream iss(arg);
    Rotation rot(0, 1);
    try {
        rot = Cytosim::readRotation(iss);
        // can combine a second rotation:
        if ( iss.good() )
            rot = Cytosim::readRotation(iss) * rot;
    }
    catch ( Exception& e )
    {
        throw InvalidSyntax("could not determine rotation from `"+arg+"'");
    }
    if ( StreamFunc::has_trail(iss) )
    {
        std::string sub = arg.substr(iss.tellg());
        throw InvalidSyntax("discarded `"+sub+"' in rotation");
    }
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

Rotation Cytosim::readOrientation(std::istream& is, Vector const& pos, Space const* spc)
{
    int c = Tokenizer::skip_space(is, false);

    if ( isalpha(c) )
    {
        try {
            return readRotation(is);
        }
        catch ( Exception& e )
        {
            // Some keywords are handled by readDirection(), so we just warn here
            std::cerr << "Warning, " << e.message() << '\n';
        }
    }

    // normally a unit vector is specified:
    Vector vec = readDirection(is, pos, spc);
    
    /*
     A single Vector does not uniquely define a rotation in 3D:
     hence we return a random rotation that is picked uniformly
     among all possible rotations transforming (1,0,0) in vec.
     */
    return Rotation::randomRotationToVector(vec);
}

