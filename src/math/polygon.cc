// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University.
#include "polygon.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <errno.h>
#include <string>
#include <cmath>
#include "exceptions.h"


Polygon::Polygon()
: pts_(nullptr), npts_(0)
{
}


Polygon::~Polygon()
{
    delete[] pts_;
}


void Polygon::allocate(unsigned s)
{
    delete[] pts_;
    pts_  = new Point2D[s+2];
    maxp_ = s;
}


void Polygon::set(unsigned ord, real rad, real ang)
{
    allocate(ord);
    npts_ = ord;
    real a = 2 * M_PI / (real)ord;
    for ( unsigned i = 0; i < ord; ++i )
    {
        pts_[i].xx = rad * std::cos(ang);
        pts_[i].yy = rad * std::sin(ang);
        ang += a;
    }
    wrap();
}


void Polygon::setPoint(unsigned i, real x, real y, long c)
{
    if ( pts_ && i < npts_ )
    {
        pts_[i].xx = x;
        pts_[i].yy = y;
        pts_[i].spot = c;
    }
    else
        throw InvalidParameter("invalid index in Polygon::setPoint()");
}


/**
 Each point should be on its own line: X Y
 This will set coordinates in pts_[] within the limit given by `alc`.
 but it will return the number of points in the stream, even if this is greater than `alc`.
 Thus two calls to this function should be enough to:
 - determine the number of points in the file
 - allocate an appropriate array
 - read the coordinates
 */
size_t Polygon::read(std::istream& in, Point2D* pts, size_t pts_size)
{
    size_t i = 0;
    char str[2048];
    real x, y;
    long k;

    while ( in.good() )
    {
        in.getline(str, sizeof(str));

        if ( in.fail() && in.gcount() )
            throw InvalidIO("Could not read polygon coordinate file");
       
        char const* ptr = str;
        char * end = nullptr;
        
        x = strtod(ptr, &end);
        
        if ( end == ptr )
            continue;

        ptr = end;
        y = strtod(ptr, &end);
        
        if ( end == ptr )
        {
            y = 0;
            k = 0;
        }
        else
        {
            ptr = end;
            errno = 0;
            k = strtol(ptr, &end, 10);
            if ( errno || end == ptr )
                k = 0;
        }
        
        if ( i < pts_size )
        {
            //std::clog << "polygon["<<i<<"]: "<<str<<"| "<<x<<" "<<y<<" "<<k<<"\n";
            pts[i].xx = x;
            pts[i].yy = y;
            pts[i].spot = k;
        }
        ++i;
    }
    return i;
}


void Polygon::read(std::istream& in)
{
    size_t s = read(in, nullptr, 0);
    unsigned n = (unsigned)s;
    if ( s != n )
        throw InvalidIO("Error: too many points in polygon coordinate file");
    allocate(n);
    npts_ = n;
    in.clear();
    in.seekg(0);
    read(in, pts_, npts_);
    wrap();
}


void Polygon::read(std::string const& file)
{
    if ( file.empty() )
        throw InvalidParameter("a polygon file should be specified");
    
    std::ifstream in(file.c_str(), std::ifstream::in);
    
    if ( ! in.good() )
        throw InvalidParameter("polygon file `"+file+"' not found");
    
    //std::clog << "Reading polygon from " << file << '\n';
    
    read(in);
    in.close();

    if ( nbPoints() < 3 )
        throw InvalidParameter("polygon: too few points specified in `"+file+"'");
}


void Polygon::write(std::ostream& os) const
{
    os.precision(6);
    os << std::fixed;
    for ( unsigned i = 1; i < npts_; ++i )
    {
        os << std::setw(12) << pts_[i].xx << "  " << std::setw(12) << pts_[i].yy;
        if ( pts_[i].spot )
            os << " " << pts_[i].spot;
        std::endl(os);
    }
}


/**
 This makes a clockwise Polygon counter-clockwise and vice-versa
 */
void Polygon::flip()
{
    unsigned n = 1;
    unsigned p = npts_-1;
    
    while ( n < p )
    {
        Point2D X = pts_[p];
        pts_[p] = pts_[n];
        pts_[n] = X;
        ++n;
        --p;
    }
    wrap();
}


void Polygon::transform(real sx, real sy, real dx, real dy)
{
    for ( unsigned n = 0; n < npts_; ++n )
    {
        pts_[n].xx = sx * pts_[n].xx + dx;
        pts_[n].yy = sy * pts_[n].yy + dy;
    }
    wrap();
}


/**
 This will increase the area if the Polygon is counterclockwise
 */
void Polygon::inflate(real eps)
{
    complete(REAL_EPSILON);
    
    if ( npts_ < 3 )
        return;
    
    // tangent 'T' to previous segment
    real tx = pts_[npts_-1].dx;
    real ty = pts_[npts_-1].dy;
    
    // previous point shifted out by 'eps'
    real px = pts_[npts_-1].xx + eps * ty;
    real py = pts_[npts_-1].yy - eps * tx;

    for ( unsigned n = 0; n < npts_; ++n )
    {
        // normal 'N' to current segment
        real nx =  pts_[n].dy;
        real ny = -pts_[n].dx;

        /*
         Calculate interestion of the shifted lines
         First line is parametric:     X = P + a * T
         Second line is defined by:  ( X - Q ) . N = eps
         where P = previous point, Q = current point
         */
        real s = tx * nx + ty * ny;
        
        if ( abs_real(s) < REAL_EPSILON )
        {
            px = pts_[n].xx + eps * nx;
            py = pts_[n].yy + eps * ny;
        }
        else
        {
            real a = eps + ( pts_[n].xx - px ) * nx + ( pts_[n].yy - py ) * ny;
            px += tx * a / s;
            py += ty * a / s;
        }
        //std::clog << " n " << n << "  " << px << "  " << py << "\n";
        
        tx = pts_[n].dx;
        ty = pts_[n].dy;
    
        pts_[n].xx = px;
        pts_[n].yy = py;
    }
    wrap();
}


/**
 box[] = { xmin, xmax, ymin, ymax }
 
 result is undefined if ( npts == 0 ).
 */
void Polygon::find_extremes(real box[4]) const
{
    if ( npts_ > 0 )
    {
        box[0] = pts_[0].xx;
        box[1] = pts_[0].xx;
        box[2] = pts_[0].yy;
        box[3] = pts_[0].yy;
    }
    else
    {
        box[0] = 0;
        box[1] = 0;
        box[2] = 0;
        box[3] = 0;
    }
    
    for ( unsigned i = 1; i < npts_; ++i )
    {
        box[0] = std::min(box[0], pts_[i].xx);
        box[1] = std::max(box[1], pts_[i].xx);
        box[2] = std::min(box[2], pts_[i].yy);
        box[3] = std::max(box[3], pts_[i].yy);
    }
}


void Polygon::wrap()
{
    if ( npts_ > 1 ) pts_[npts_] = pts_[0];
    if ( npts_ > 2 ) pts_[npts_+1] = pts_[1];
}

/**
 pre-calculate offset of successive points,
 and length of segments, used in project for efficiency.
 Also duplicate first two points to simplify some later calculations:
 - point[L  ] <- point[0]
 - point[L+1] <- point[1]
 .
 
 The array should be allocated to hold (npts+2) Point2D
 */
int Polygon::complete(const real epsilon)
{
    //std::clog << "Polygon::complete(" << epsilon << ")\n";
    if ( npts_ > 1 )
    {
        unsigned i = 1, n = 0, p = 0;
        do {
            real dx = 0, dy = 0;
            // skip consecutive points that are too close from each other:
            do {
                if ( ++n > npts_ )
                    goto finish;
                dx = pts_[n].xx - pts_[p].xx;
                dy = pts_[n].yy - pts_[p].yy;
                //std::clog << "poly  " << n << "   " << dx << " " << dy << "\n";
            } while ( std::max(abs_real(dx), abs_real(dy)) < epsilon );
            //std::clog << i << " <-- " << n << "\n";
            pts_[i] = pts_[n];
            // normalize the vector:
            real d = std::sqrt( dx * dx + dy * dy );
            pts_[i].dx  = dx / d;
            pts_[i].dy  = dy / d;
            pts_[i].len = d;
            p = n;
            ++i;
        } while ( p <= npts_ );

finish:
        if ( i < npts_ )
        {
            std::clog << "Polygon had " << npts_-i << " degenerate vertices\n";
            npts_ = i;
            return 1;
        }
        wrap();
    }

    return 0;
}

    
/**
 calculate volume of polygon, using an algorithm that return a negative value
 for a polygon defined clockwise and a positive value for anti-clockwise.
 http://mathworld.wolfram.com/PolygonArea.html
 */
real Polygon::surface() const
{
    if ( npts_ < 3 )
        return 0;
    
    real S = pts_[npts_-1].xx * ( pts_[0].yy - pts_[npts_-2].yy );
    for ( unsigned i = 2; i < npts_; ++i )
        S += pts_[i-1].xx * ( pts_[i].yy - pts_[i-2].yy );
    
    return S * 0.5;
}


/**
 Count the number of time a ray from (xx, yy) to (infinity, yy) crosses the polygon
 The point is inside if the result is odd. 
 This method works for clockwise and anticlockwise polygons.
 
 @return
 0    : point is outside
 1    : point is inside
 edge : point is near the boundary, within distance `threshold`
 .
*/
int Polygon::inside(real xx, real yy, int edge, real threshold) const
{
    int cross = 0;
        
    // check sides of polygon against horizontal line from (xx, yy) to (+inf, yy)
    for ( unsigned i = 0; i < npts_; ++i )
    {
        const Point2D & p1 = pts_[i]; // not included in edge
        const Point2D & p2 = pts_[i+1];

        //p1.write(std::clog); p2.write(std::clog); std::clog << "\n";
        // check if edge cannot intersect with ray
        bool above = ( yy <= p1.yy) && ( yy < p2.yy );
        bool below = ( yy >= p1.yy) && ( yy > p2.yy );
        if ( above || below )
            continue;
        
        // ray may go through p2
        if ( yy == p2.yy )
        {
            // check for horizontal edge
            if ( p1.yy == p2.yy )
            {
                bool left = ( xx < p1.xx ) && ( xx < p2.xx );
                bool right = ( xx > p1.xx ) && ( xx > p2.xx );
                if ( left | right )
                    continue;
                return edge;
            }
            
            if ( p2.xx < xx )
                continue;
            
            if ( xx == p2.xx )
                return edge;
            
            ++cross;

            // next vertex
            const Point2D& p3 = pts_[i+1];
         
            // check that p2 is not a corner
            bool updown = ( p1.yy < yy ) && ( yy < p3.yy );
            bool downup = ( p3.yy < yy ) && ( yy < p1.yy );
            cross += updown | downup;
            continue;
        }
        
        // xx is left of edge
        if ((xx <= p1.xx) | (xx <= p2.xx))
        {
            // intersection of ray with edge
            real xi = ( yy - p1.yy ) * ( p2.xx - p1.xx ) / ( p2.yy - p1.yy ) + p1.xx;
            
            // overlies on an edge
            if ( abs_real( xx - xi ) < threshold )
                return edge;
                
            // xx left of intersection
            if ( xx < xi )
                ++cross;
        }
    }
    
    //std::clog << " polygon::inside " << cross << " for (" << xx << " " << yy << ")\n";
    return ( cross & 1 );
}

/**
 Find the closest point on the polygon to (xx, yy)
 
 @return
 1 : projects on an edge
 0 : projects on a vertex
 .
 The function will also set `hit` to be the index of the point,
 or the segment where the projection landed

 */
int Polygon::project(real xx, real yy, real& pX, real& pY, unsigned& hit) const
{
    if ( npts_ < 1 )
        throw InvalidParameter("cannot project on uninitialized polygon");
    int res = 0;
    
    //initialize with first point:
    hit = 0;
    pX = pts_[0].xx;
    pY = pts_[0].yy;
    
    real best = square(xx-pX) + square(yy-pY);
    
    for ( unsigned i = 0; i < npts_; ++i )
    {
        real x = xx - pts_[i].xx;
        real y = yy - pts_[i].yy;
        // distance to this polygon point:
        real d = x * x + y * y;
        // abscissa of projection on segment [i, i+1] of the polygon:
        real a = pts_[i].dx * x + pts_[i].dy * y;
        
        if ( a > 0 )
        {
            if ( a < pts_[i].len )
            {
                // distance from (xx, yy) to segment:
                real da = d - a * a;
                
                if ( da < best )
                {
                    best = da;
                    pX = pts_[i].xx + a * pts_[i].dx;
                    pY = pts_[i].yy + a * pts_[i].dy;
                    hit = i;
                    res = 1;
                }
            }
        }
        else if ( d < best )
        {
            best = d;
            pX = pts_[i].xx;
            pY = pts_[i].yy;
            hit = i;
            res = 0;
        }
    }
    return res;
}


void Polygon::dump(std::ostream& os) const
{
    const int W = 10;
    os << "polygon " << npts_ << "\n";
    for ( unsigned n = 0; n < npts_; ++n )
    {
        os << " " << std::setw(W) << pts_[n].xx << " " << std::setw(W) << pts_[n].yy << " " << pts_[n].spot;
        os << " " << std::setw(W) << pts_[n].dx << " " << std::setw(W) << pts_[n].dy << "\n";
    }
}


void Polygon::print(FILE * f) const
{
    fprintf(f, "polygon %u\n", npts_);
    for ( unsigned n = 0; n < npts_; ++n )
    {
        fprintf(f, "%10.2f %10.2f %4li", pts_[n].xx, pts_[n].yy, pts_[n].spot);
        fprintf(f, "  %10.2f %10.2f\n", pts_[n].dx, pts_[n].dy);
    }
}

