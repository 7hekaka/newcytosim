// Cytosim was created by Francois Nedelec. Copyright Cambridge University 2019-2020

#include <iostream>
#include <numeric>
#include <list>
#include <set>
#include "tokenizer.h"
#include "stream_func.h"
#include "organizer.h"
#include "matrix22.h"
#include "random_pcg.h"


/// width of columns in formatted output, in number of characters
static int column_width = 10;

/// use this macro at the beginning of a line of comment
#define COM "\n% " << std::setw(column_width-2)

/// use this macro at the beginning of a new line of data
#define LIN '\n' << std::setw(column_width)

/// use this macro to separate between values on a line of data
#define SEP ' ' << std::setw(column_width-1)

#include "accumulator.h"

/// pad string by adding white-space on the right up to size 'n * column_width - p'
static std::string ljust(std::string const& str, size_t n, size_t p = 0)
{
    size_t s = n * (size_t)column_width - p;
    if ( str.size() < s )
        return str + std::string(s-str.size(), ' ');
    else
        return str;
}

/// pad string by adding white-space on the left up to size 'n * column_width - p'
static std::string rjust(std::string const& str, size_t n, size_t p = 1)
{
    size_t s = n * (size_t)column_width - p;
    if ( str.size() < s )
        return std::string(s-str.size(), ' ') + str;
    else
        return str;
}

/// repeat string DIM times with X, Y, Z endings as appropriate
static std::string repeatXYZ(std::string const& str)
{
    std::string res(rjust(str+"X", 1, 1));
#if ( DIM > 1 )
    res += " " + rjust(str+"Y", 1, 1);
#endif
#if ( DIM > 2 )
    res += " " + rjust(str+"Z", 1, 1);
#endif
    return res;
}

/// remove any 's' at the end of the argument
static void remove_plural(std::string & str)
{
    if ( str.size() > 2  &&  str.at(str.size()-1) == 's' )
        str.resize(str.size()-1);
}


/**
 combines multiple report, if `what` has multiple instructions separated by ','
 for example:
 
     report fiber:force,fiber:length
 
 parameters can be specifed as in:
 
     report fiber:cluster{couple=1},fiber:length

 */
void Simul::poly_report(std::ostream& out, std::string what, Glossary& opt, int frm) const
{
    int ver = 2;
    opt.set(ver, "verbose");
    if (( ver & 2 ) && frm >= 0 )
        out << "% frame   " << frm << '\n';
    std::stringstream is(what);
    while ( is.good() )
    {
        std::string arg = Tokenizer::get_polysymbol(is, false);
        std::string blk = Tokenizer::get_block(is, '{');
        try {
            if ( blk.empty() )
            {
                //out << "\nSimul::report(" << arg << ")";
                mono_report(out, arg, opt, ver);
            }
            else
            {
                //out << "\nSimul::report(" << arg << ", " << blk << ")";
                Glossary glos(opt);
                glos.read(blk, 0);
                mono_report(out, arg, glos, ver);
                opt.add_counts(glos);
            }
        }
        catch( Exception & e )
        {
            out << '\n' << e.brief();
        }
        // another report instruction should be separated by a comma:
        if ( is.peek() != ',' )
            break;
        is.get();
    }
    if ( ver & 2 )
    {
        out << "% end report\n\n";
    }
    {
        // check for unused characters in instruction stream
        std::string str;
        std::getline(is, str);
        if ( str.size() > 0 )
            throw InvalidParameter("unexpected `" + str + "' in report string");
    }
    //opt.print_counts(std::cerr);
}


/**
 Surround the report with comments to identify start/end
 */
void Simul::mono_report(std::ostream& out, std::string const& arg, Glossary& opt, int ver) const
{
    std::streamsize p = 4;
    opt.peek(p, "precision");
    opt.peek(column_width, "column") || opt.peek(column_width, "width");
    
    // adjust floating-point notation:
    out.setf(std::ios_base::fixed, std::ios_base::floatfield);
    std::streamsize op = out.precision();
    out.precision(p);

    if ( ver & 1 )
    {
        //out << "% start\n";
        out << "% time " << std::to_string(time()) << '\n';
    }
    if ( ver > 0 )
    {
        out << "% report " << arg << " " << opt.to_string();
        report_one(out, arg, opt);
        out << "\n";
    }
    else
    {
        std::stringstream ss;
        report_one(ss, arg, opt);
        StreamFunc::skip_lines(out, ss, '%', false);
    }
    if ( ver & 1 )
        out << "% end\n";
    out.precision(op);
}


/**
 Split 'arg' into who:what and call report_one()
 Surround the report with comments to identify start/end
 */
void Simul::report_one(std::ostream& out, std::string const& arg, Glossary& opt) const
{
    std::string who = arg, what;
    
    // split argument string into who:what
    std::string::size_type pos = arg.find(':');
    if ( pos != std::string::npos )
    {
        who  = arg.substr(0, pos);
        what = arg.substr(pos+1);
    }
    
    // allow for 's' or not at the ends of words:
    remove_plural(who);
    remove_plural(what);
    
    //std::clog << "report("<< who << "|" << what << ")\n";
    if ( isCategory(who) )
    {
        int split = false;
        if ( opt.set(split, "split") && split )
        {
            int com = 1;
            opt.peek(com, "verbose");
            // generate a separate report for all classes in this category:
            PropertyList plist = properties.find_all(who);
            if ( com && plist.size() )
            {
                out << COM << "split:";
                for ( Property const* sel : plist )
                    out << SEP << sel->name();
            }
            for ( Property const* sel : plist )
            {
                if ( com )
                {
                    report_one(out, who, sel, what, opt);
                    com = 0;
                }
                else
                {
                    std::stringstream ss;
                    report_one(ss, who, sel, what, opt);
                    StreamFunc::skip_lines(out, ss, '%');
                }
            }
        }
        else
        {
            // pool all objects in this category:
            report_one(out, who, nullptr, what, opt);
        }
    }
    else
    {
        // check if name corresponds to a property:
        Property const* sel = nullptr;
        sel = properties.find(who);
        if ( sel )
            report_one(out, sel->category(), sel, what, opt);
        else
            report_one(out, who, nullptr, what, opt);
    }
}


/**
 WHAT            | Output
 ----------------|--------------------------------------------------------------
 `bead`          | Position of beads
 `couple`        | Summary with number of couples in each state
 `fiber`         | Length and position of the ends of fibers
 `single`        | Number and state of singles
 `solid`         | Position of center and first point of solids
 `sphere`        | Position of center and first point of spheres
 `organizer`     | Position of the center of asters and other organizers
 `field`         | Total quantity of substance in field and Lattices
 
 
 WHAT               | Output
 -------------------|----------------------------------------------------------
 `simul:time`       | Time
 `simul:inventory`  | summary list of objects
 `simul:property`   | All properties
 `simul:parameter`  | global parameters
 `simul:NAME`       | parameters for Property 'NAME'


 WHAT                    | Output
 ------------------------|------------------------------------------------------
 `fiber:position`        | Position and orientation of fibers
 `fiber:age`             | Average age of fibers
 `fiber:length`          | Average length and variance of lengths of fibers
 `fiber:distribution`    | length distribution of fiber lengths (option: `max` and `interval`)
 `fiber:dynamic`         | Number of fiber classified by PLUS_END Dynamic state
 `fiber:point`           | coordinates of vertices of all fibers
 `fiber:displacement`    | mean squared displacement of fibers since the last call
 `fiber:moments`         | standard deviation of vertices of all fibers
 `fiber:speckle`         | coordinates of points randomly distributed along all fibers (option: `interval`)
 `fiber:sample`          | coordinates of points newly distributed along all fibers (option: `interval`)
 `fiber:segment`         | information about lengths of segments, number of kinks
 `fiber:end`             | Positions and dynamic states of all fiber ends
 `fiber:force`           | Position of vertices and Forces acting on vertices
 `fiber:tension`         | Internal stress along fibers
 `fiber:energy`          | Fiber's elastic bending energy
 `fiber:confinement`     | Force applied by fibers on their confinement Space
 `fiber:binder`          | Positions of bridging hands along each fiber
 `fiber:lattice`         | Total quantity on fiber's lattices
 `fiber:mesh`            | Total quantity on fiber's meshes
 `fiber:intersection`    | Intersections point of fibers
 `fiber:hand`            | Position of hands attached to fibers
 `fiber:link`            | Positions of attached hands for which stiffness > 0
 `fiber:cluster`         | Clusters made of fibers connected by Couples


 WHAT                    | Output
 ------------------------|------------------------------------------------------
 `bead:all`              | Position of beads
 `bead:single`           | Number of Beads with no single attached, 1 single attached etc.
 `solid:hand`            | Number of hands and number of attached hands on Solids
 `solid:orientation`     | Solid's position and orientation vectors
 `spindle:indice`        | Two scalar indices that caracterize the organization of fibers
 `spindle:profile`       | Number of right- and left-pointing fiber as a function of position
 `single:all`            | Position and force of singles
 `single:force`          | Average and maximum force of singles
 `single:NAME`           | Position and force of singles of class NAME
 `NAME:position`         | Position and force of singles of class NAME
 `couple:state`          | Position and state of all couples
 `couple:NAME`           | Position and state of couples of class NAME
 `couple:link`           | detailed information on doubly-attached couples
 `couple:configuration`  | number of Couples in { X, P, A, V, T } states
 `couple:force`          | Average and maximum of tension in the couple links
 `couple:histogram`      | Histogram of tension in the couple links
 `couple:active`         | Position of active couples
 `couple:anatomy`        | Composition of couples
 `NAME:position`         | Position of couples of class NAME
 `couple:hands`          | Composition of couples
 
 */
void Simul::report_one(std::ostream& out, std::string const& who, Property const* sel,
                       std::string const& what, Glossary& opt) const
{
    if ( who == "fiber" )
    {
        if ( what.empty() || what == "position" )
            return reportFibers(out, sel);
        if ( what == "plus_end" )
            return reportFiberEnds(out, PLUS_END, sel);
        if ( what == "minus_end" )
            return reportFiberEnds(out, MINUS_END, sel);
        if ( what == "end" )
            return reportFiberEnds(out, BOTH_ENDS, sel);
        if ( what == "point" )
            return reportFiberPoints(out, sel);
        if ( what == "displacement" )
            return reportFiberDisplacement(out, sel);
        if ( what == "direction" )
            return reportFiberDirections(out, sel);
        if ( what == "moment" )
            return reportFiberMoments(out);
        if ( what == "speckle" )
            return reportFiberSpeckles(out, opt);
        if ( what == "sample" )
            return reportFiberSamples(out, opt);
        if ( what == "length" )
            return reportFiberLengths(out, sel);
        if ( what == "distribution" || what == "histogram" )
            return reportFiberLengthHistogram(out, opt);
        if ( what == "tension" )
            return reportFiberTension(out, opt);
        if ( what == "segment" )
            return reportFiberSegments(out);
        if ( what == "energy" )
            return reportFiberBendingEnergy(out);
        if ( what == "extension" )
            return reportFiberExtension(out);
        if ( what == "nematic" )
            return reportFiberNematic(out);
        if ( what == "end_state" || what == "dynamic" )
        {
            reportFiberEndState(out, PLUS_END, sel);
            std::stringstream ss;
            reportFiberEndState(ss, MINUS_END, sel);
            StreamFunc::skip_lines(out, ss, '%');
            return;
        }
        if ( what == "plus_state" )
            return reportFiberEndState(out, PLUS_END, sel);
        if ( what == "minus_state" )
            return reportFiberEndState(out, MINUS_END, sel);
        if ( what == "force" )
            return reportFiberForces(out);
        if ( what == "confine_force" )
            return reportFiberConfineForce(out);
        if ( what == "confinement" )
            { reportFiberConfinement(out); return; }
        if ( what == "cluster" )
            return reportClusters(out, opt);
        if ( what == "age" )
            return reportFiberAge(out);
        if ( what == "intersection" )
            return reportFiberIntersections(out, opt);
        if ( what == "hand" )
            return reportFiberHands(out);
        if ( what == "link" )
            return reportFiberLinks(out);
        if ( what == "lattice" )
            return reportFiberLattice(out, sel);
        if ( what == "mesh" )
            return reportFiberMeshAverage(out, false, sel);
        if ( what == "mesh_density" )
            return reportFiberMeshAverage(out, true, sel);
        if ( what == "mesh_all" )
            return reportFiberMesh(out, false, sel);
        if ( what == "connector" )
            return reportFiberConnectors(out, opt);
        if ( what == "collision" )
            return reportFiberCollision(out, sel, opt);

        throw InvalidSyntax("I can only report fiber: position, end, minus_end, plus_end, "\
                            "point, moment, speckle, sample, segment, dynamic, length, extension,"\
                            "distribution, tension, force, cluster, age, energy, hand, link");
    }
    if ( who == "bead" )
    {
        if ( what == "position" || what.empty() )
            return reportBeadPosition(out, sel);
        if ( what == "single" )
            return reportBeadSingles(out);
        throw InvalidSyntax("I can only report bead: position, single");
    }
    if ( who == "solid" )
    {
        if ( what == "hand" )
            return reportSolidHands(out, sel);
        if ( what == "orientation" )
            return reportSolidOrientation(out, sel);
        if ( what == "position" || what.empty() )
            return reportSolidPosition(out, sel);
        throw InvalidSyntax("I can only report solid: hand, position");
    }
    if ( who == "space" )
    {
        if ( what == "force" )
            return reportSpaceForce(out);
        if ( what.empty() )
            return reportSpace(out);
        throw InvalidSyntax("I can only report `space` and space:force");
    }
    if ( who == "sphere" )
    {
        if ( what == "position" || what.empty() )
            return reportSpherePosition(out, sel);
        throw InvalidSyntax("I can only report sphere:position");
    }
    if ( who == "single" )
    {
        if ( what.empty() )
            return reportSingle(out, sel);
        if ( what == "link" )
            return reportSingleLink(out, sel);
        if ( what == "state" )
            return reportSingleState(out, sel);
        if ( what == "force" )
            return reportSingleForce(out, sel);
        if ( what == "position" )
            return reportSinglePosition(out, sel);
        throw InvalidSyntax("I can only report single: link, state, force, position");
    }
    if ( who == "couple" )
    {
        if ( what.empty() )
            return reportCouple(out, sel);
        if ( what == "list" )
            return reportCoupleList(out, sel);
        if ( what == "state" )
            return reportCoupleState(out, sel);
        if ( what == "link" )
            return reportCoupleLink(out, sel);
        if ( what == "configuration" )
            return reportCoupleConfiguration(out, sel, opt);
        if ( what == "force" )
            return reportCoupleForce(out, sel);
        if ( what == "histogram" )
            return reportCoupleForceHistogram(out, opt);
        if ( what == "active" )
            return reportCoupleActive(out, sel);
        if ( what == "anatomy" )
            return reportCoupleAnatomy(out, sel);
        throw InvalidSyntax("I can only report couple: state, link, configuration, active, force, anatomy");
    }
    if ( who == "organizer" )
    {
        if ( what.empty() )
            return reportOrganizer(out);
        throw InvalidSyntax("I can only report `organizer'");
    }
    if ( who == "aster" )
    {
        if ( what.empty() )
            return reportAster(out);
        throw InvalidSyntax("I can only report `aster'");
    }
    if ( who == "field" )
    {
        return reportField(out);
    }
    if ( who == "simul" || who.empty() )
    {
        if ( what.empty() )
            return reportSimul(out);
        if ( what == "time" )
            return reportTime(out);
        if ( what == "inventory" )
            return reportInventory(out);
        if ( what == "property" || what == "parameter" )
            return writeProperties(out, false);
        throw InvalidSyntax("I can only report simul: time, inventory, property");
    }
    if ( who == "property" )
    {
        if ( what.empty() || what=="all" )
            return writeProperties(out, false);
        Property const* p = findProperty(what);
        if ( p )
            return p->write(out);
        throw InvalidSyntax("unknown property");
    }
    if ( who == "spindle" )
    {
        if ( what == "indice" )
            return reportIndices(out);
        if ( what == "profile" )
            return reportProfile(out);
        throw InvalidSyntax("I can only report spindle: indices, profile");
    }
    if ( who == "network" )
    {
        if ( what == "bridge" )
            return reportNetworkBridges(out, opt);
        if ( what == "size" )
            return reportNetworkSize(out);
    }
    if ( who == "ring" )
        return reportRing(out);
    if ( who == "platelet" )
        return reportPlatelet(out);
    if ( who == "ashbya" )
        return reportAshbya(out);
    if ( who == "something" )
        return reportSomething(out);

    throw InvalidSyntax("unknown report `"+who+":"+what+"'");
}

//------------------------------------------------------------------------------
#pragma mark - Fiber Aggregated Properties

/**
 Export average length and variance of lengths for each class of fiber
 */
void Simul::reportFiberAge(std::ostream& out) const
{
    out << COM << ljust("class", 2, 2) << SEP << "count" << SEP << "avg_birth";
    out << SEP << "var_birth" << SEP << "avg_age" << SEP << "min_age" << SEP << "max_age";
    
    const real now = time();

    for ( Property const* i : properties.find_all("fiber") )
    {
        FiberProp const* fp = static_cast<FiberProp const*>(i);
        ObjectList objs = fibers.collect(fp);
        size_t cnt = 0;
        real avg = 0, var = 0, mn = INFINITY, mx = -INFINITY;
        fibers.infoBirthtime(objs, cnt, avg, var, mn, mx);
        out << LIN << ljust(fp->name(), 2);
        out << SEP << cnt;
        out << SEP << avg;
        out << SEP << var;
        out << SEP << now-mx;
        out << SEP << now-avg;
        out << SEP << now-mn;
    }
}


/**
Export average length and variance of length for a class of fiber
*/
void Simul::reportFiberLengths(std::ostream& out, FiberProp const* sel) const
{
    size_t cnt = 0;
    real avg = 0, var = 0, mn = INFINITY, mx = -INFINITY;
    ObjectList objs = fibers.collect(sel);
    fibers.infoLength(objs, cnt, avg, var, mn, mx);
    
    out << LIN << ljust(sel->name(), 2);
    out << SEP << cnt;
    out.precision(3);
    out << SEP << std::fixed << avg;
    out << SEP << std::fixed << var;
    out << SEP << std::fixed << mn;
    out << SEP << std::fixed << mx;
    out.precision(1);
    out << SEP << std::fixed << avg*cnt;
}

/**
 Export average length and variance for each class of fiber
 */
void Simul::reportFiberLengths(std::ostream& out, Property const* sel) const
{
    out << COM << ljust("class", 2, 2) << SEP << "count" << SEP << "avg_len" << SEP << "var_len";
    out << SEP << "min_len" << SEP << "max_len" << SEP << "total";
    
    std::streamsize p = out.precision();
    if ( sel )
        reportFiberLengths(out, static_cast<FiberProp const*>(sel));
    else
    {
        for ( Property const* i : properties.find_all("fiber") )
            reportFiberLengths(out, static_cast<FiberProp const*>(i));
    }
    out.precision(p);
}


/**
 Export length histograms for each class of fiber
 */
void Simul::reportFiberLengthHistogram(std::ostream& out, Glossary & opt) const
{
    const size_t BMAX = 256;
    unsigned cnt[BMAX+1];

    real sup = 0;
    for ( Fiber const* fib = fibers.first(); fib; fib = fib->next() )
        sup = std::max(sup, fib->length());

    size_t nbin = BMAX;
    real delta = ( sup > 2 ) ? 1 : 0.1;
    opt.set(delta, "interval", "bins");
    if ( !opt.set(nbin, "interval", 1, "bins", 1) )
        nbin = std::ceil(sup/delta);
    nbin = std::min(nbin, BMAX);
    
    if ( 1 )
    {
        out << COM << "bin " << delta << " count " << fibers.size();
        out << LIN << ljust("scale", 2);
        std::streamsize p = out.precision();
        out.precision(2);
        for ( size_t u = 0; u <= nbin; ++u )
            out << " " << std::setw(5) << delta * ( u + 0.5 );
        out.precision(p);
    }
    
    for ( Property const* i : properties.find_all("fiber") )
    {
        FiberProp const* fp = static_cast<FiberProp const*>(i);
        
        for ( size_t u = 0; u <= nbin; ++u )
            cnt[u] = 0;
        
        for ( Fiber const* fib = fibers.first(); fib; fib = fib->next() )
        {
            if ( fib->prop == fp )
            {
                size_t u = (size_t)std::floor( fib->length() / delta );
                ++cnt[std::min(u, nbin)];
            }
        }

        out << LIN << ljust(fp->name(), 2);
        for ( size_t u = 0; u <= nbin; ++u )
            out << " " << std::setw(5) << cnt[u];
    }
}


/**
 Export number of fiber, classified according to dynamic state of one end
 */
void Simul::reportFiberEndState(std::ostream& out, FiberEnd end, Property const* sel) const
{
    std::string name;
    if ( sel )
        name = sel->name() + ":";
    name.append(end==PLUS_END ?"plus":"minus");
    
    out << COM << ljust("class", 2, 2) << SEP << "total" << SEP << "static";
    out << SEP << "green" << SEP << "yellow" << SEP << "orange" << SEP << "red";
    
    constexpr size_t MAX = 5;
    size_t cnt[MAX+1] = { 0 };
    size_t sum = 0;
    
    for ( Fiber const* fib = fibers.first(); fib; fib = fib->next() )
    {
        if ( !sel || sel == fib->prop )
        {
            ++sum;
            state_t x = std::min(fib->endState(end), (state_t)MAX);
            ++cnt[x];
        }
    }
    
    out << LIN << ljust(name, 2) << SEP << sum;
    for ( size_t i = 0; i < MAX; ++i )
        out << SEP << cnt[i];
}


//------------------------------------------------------------------------------
#pragma mark - Fiber conformation


void Simul::reportFiberSegments(std::ostream& out) const
{
    out << COM << ljust("class", 2, 2) << SEP << "fibers" << SEP << "joints";
    out << SEP << "kinks" << SEP << "min_seg" << SEP << "max_seg" << SEP << "err_seg";
    
    for ( Property const* i : properties.find_all("fiber") )
    {
        FiberProp const* fp = static_cast<FiberProp const*>(i);
        
        size_t cnt, points;
        real mn = 0, mx = 0, dv = 0;
        
        ObjectList objs = fibers.collect(fp);
        fibers.infoSegments(objs, cnt, points, mn, mx, dv);
        out << LIN << ljust(fp->name(), 2);
        out << SEP << cnt;
        out << SEP << points - 2 * cnt;
        out << SEP << fibers.nbKinks(objs);
        out << SEP << std::fixed << mn;
        out << SEP << std::fixed << mx;
        out << SEP << std::fixed << dv;
    }
}


/**
 Export fiber elastic bending energy
 */
void Simul::reportFiberBendingEnergy(std::ostream& out) const
{
    out << COM << ljust("bending_energy",2,2) << SEP << "count";
    out << SEP << "sum" << SEP << "avg" << SEP << "var" << SEP << "rigidity";
    
    for ( Property const* i : properties.find_all("fiber") )
    {
        FiberProp const* fp = static_cast<FiberProp const*>(i);
        ObjectList objs = fibers.collect(fp);
        size_t cnt = 0;
        real avg = 0, var = 0;
        fibers.infoBendingEnergy(objs, cnt, avg, var);
        if ( cnt > 0 )
        {
            out << LIN << ljust(fp->name(), 2);
            out << SEP << cnt;
            out << SEP << avg*cnt;
            out << SEP << avg;
            out << SEP << var;
            out << SEP << fp->rigidity;
        }
    }
}


void Simul::reportFiberExtension(std::ostream& out) const
{
    out << COM << ljust("end_to_end_dist",2,2) << SEP << "count";
    out << SEP << "avg" << SEP << "var" << SEP << "min" << SEP << "max";
    
    for ( Property const* i : properties.find_all("fiber") )
    {
        size_t cnt = 0;
        real avg = 0, var = 0, in = INFINITY, ax = -INFINITY;
        FiberProp const* p = static_cast<FiberProp const*>(i);
        for ( Object const* f : fibers.collect(p) )
        {
            Fiber const* fib = Fiber::toFiber(f);
            real x = ( fib->posEndP() - fib->posEndM() ).norm();
            in = std::min(in, x);
            ax = std::max(ax, x);
            avg += x;
            var += x * x;
            ++cnt;
        }
        if ( cnt > 0 )
        {
            avg /= cnt;
            var -= square(avg)*cnt;
            if ( cnt > 1 )
                var /= real(cnt-1);
            out << LIN << ljust(p->name(), 2);
            out << SEP << cnt;
            out << SEP << avg;
            out << SEP << var;
            out << SEP << in;
            out << SEP << ax;
        }
    }
}


void Simul::reportFiberNematic(std::ostream& out) const
{
    out << COM << ljust("class",2,2) << SEP << "count" << SEP << "order";
    out << SEP << "dirX" << SEP << "dirY" << SEP << "dirZ";
    
    for ( Property const* i : properties.find_all("fiber") )
    {
        FiberProp const* p = static_cast<FiberProp const*>(i);
        ObjectList objs = fibers.collect(p);
        if ( objs.size() > 0 )
        {
            real vec[9] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };
            real S = FiberSet::infoNematic(objs, vec);
            out << LIN << ljust(p->name(), 2);
            out << SEP << objs.size();
            out << SEP << S;
            out << SEP << vec[0] << SEP << vec[1] << SEP << vec[2];
        }
    }
}


//------------------------------------------------------------------------------
#pragma mark - Fiber Lattice & Mesh


/**
 Report quantity of substance in the fiber's Lattice
 */
void Simul::reportFiberLattice(std::ostream& out, Property const* sel) const
{
    out << COM << ljust("class", 2, 2) << SEP << "count" << SEP << "vacant";
    out << SEP << "avg" << SEP << "avg_if" << SEP << "min" << SEP << "max";
    
    size_t cnt = 0, vac = 0;
    real sum = 0, mn = INFINITY, mx = -INFINITY;
    
    for ( Fiber const* fib = fibers.first(); fib; fib = fib->next() )
    {
        if ( !sel || sel == fib->prop )
            fib->infoLattice(cnt, vac, sum, mn, mx);
    }
    
    std::streamsize p = out.precision();
    out << LIN << ljust("fiber:lattice", 2);
    out << SEP << cnt;
    out << SEP << vac;
    out << SEP << std::setprecision(4) << sum / (real)cnt;
    out << SEP << std::setprecision(4) << sum / (real)(cnt-vac);
    out << SEP << std::fixed << std::setprecision(6) << mn;
    out << SEP << std::fixed << std::setprecision(6) << mx;
    out.precision(p);
}


/**
 Report quantity of substance in the fiber's Mesh
 */
void Simul::reportFiberMeshAverage(std::ostream& out, bool density, Property const* sel) const
{
    out << COM << ljust("class", 2, 2) << SEP << "total";
    out << SEP << "avg" << SEP << "min" << SEP << "max" << SEP << "length";
    
    size_t cnt = 0;
    real len = 0, sum = 0, mn = INFINITY, mx = -INFINITY;
    
    for ( Fiber const* fib = fibers.first(); fib; fib = fib->next() )
    {
        if ( !sel || sel == fib->prop )
            fib->infoMesh(len, cnt, sum, mn, mx, density);
    }
    
    std::streamsize p = out.precision();
    out << LIN << ljust("fiber:mesh", 2);
    out << SEP << sum;
    out << SEP << std::setprecision(4) << sum / (real)cnt;
    out << SEP << std::fixed << std::setprecision(6) << mn;
    out << SEP << std::fixed << std::setprecision(6) << mx;
    out << SEP << std::setprecision(3) << len;
    out.precision(p);
}

/**
 Report quantity of substance in the fiber's Mesh
 */
void Simul::reportFiberMesh(std::ostream& out, bool density, Property const* sel) const
{
    out << COM << ljust("fiber", 2, 2) << SEP << "total";
    out << SEP << "avg" << SEP << "min" << SEP << "max" << SEP << "length";
    
    for ( Fiber const* fib = fibers.firstID(); fib; fib = fibers.nextID(fib) )
    {
        size_t cnt = 0;
        real len = 0, sum = 0, mn = INFINITY, mx = -INFINITY;
        if ( !sel || sel == fib->prop )
        {
            fib->infoMesh(len, cnt, sum, mn, mx, density);
            std::streamsize p = out.precision();
            out << LIN << ljust(fib->reference(), 2);
            out << SEP << sum;
            out << SEP << std::setprecision(4) << sum / (real)cnt;
            out << SEP << std::fixed << std::setprecision(6) << mn;
            out << SEP << std::fixed << std::setprecision(6) << mx;
            out << SEP << std::setprecision(3) << len;
            out.precision(p);
        }
    }
}


//------------------------------------------------------------------------------
#pragma mark - Fiber Individual Properties

/**
 Export length, position and directions at center of fibers
 */
void Simul::reportFiber(std::ostream& out, Fiber const* fib) const
{
    out << LIN << fib->prop->number();
    out << SEP << fib->identity();
    out << SEP << fib->length();
    out << SEP << fib->posEnd(CENTER);
    out << SEP << fib->dirEnd(CENTER);
    out << SEP << (fib->posEndM()-fib->posEndP()).norm();
    real C = dot(fib->dirEndM(), fib->dirEndP());
#if ( DIM >= 3 )
    real S = cross(fib->dirEndM(), fib->dirEndP()).norm();
#else
    real S = cross(fib->dirEndM(), fib->dirEndP());
#endif
    out << SEP << C << SEP << S;
    out << SEP << organizers.findOrganizerID(fib);
}

    
/// qsort function comparing length of Fibers
static int compareFibers(Object const* A, Object const* B)
{
    real a = static_cast<Fiber const*>(A)->length();
    real b = static_cast<Fiber const*>(B)->length();
    return ( a < b ) - ( a > b ); //sort in decreasing length
    //return ( a > b ) - ( a < b );
}
static int comp(void const* A, void const* B)
{
    real a = static_cast<Fiber const*>(A)->length();
    real b = static_cast<Fiber const*>(B)->length();
    return ( a < b ) - ( a > b ); //sort in decreasing length
    //return ( a > b ) - ( a < b );
}


/**
 Export length, position and directions at center of fibers
 */
void Simul::reportFibersSorted(std::ostream& out, Property const* sel) const
{
    // sort fibers by decreasing length:
    const_cast<FiberSet&>(fibers).pool_.blinksort(compareFibers);
    
    out << COM << "class" << SEP << "identity" << SEP << "length";
    out << SEP << repeatXYZ("pos") << SEP << repeatXYZ("dir") << SEP << "endToEnd";
    out << SEP << "cos" << SEP << "sin" << SEP << "organizer";

    for ( Fiber const* fib = fibers.first(); fib; fib = fib->next() )
    {
        if ( !sel || sel == fib->prop )
            reportFiber(out, fib);
    }
}


/**
 Export length, position and directions at center of fibers
 */
void Simul::reportFibers(std::ostream& out, Property const* sel) const
{
    out << COM << "class" << SEP << "identity" << SEP << "length";
    out << SEP << repeatXYZ("pos") << SEP << repeatXYZ("dir") << SEP << "endToEnd";
    out << SEP << "cos" << SEP << "sin" << SEP << "organizer";
    
    // list fibers in the order of the inventory:
    for ( Fiber const* fib = fibers.firstID(); fib; fib = fibers.nextID(fib) )
    {
        if ( !sel || sel == fib->prop )
            reportFiber(out, fib);
    }
}


/**
 Export dynamic state, positions and directions of fiber
 Argument `end` can be MINUS_END, PLUS_END or BOTH_ENDS
 */
void Simul::reportFiberEnds(std::ostream& out, FiberEnd end, Property const* sel) const
{
    out << COM << "class" << SEP << "identity" << SEP << "length";
    if ( end & PLUS_END )
        out << SEP << "stateP" << SEP << repeatXYZ("posP") << SEP << repeatXYZ("dirP");
    if ( end & MINUS_END )
        out << SEP << "stateM" << SEP << repeatXYZ("posM") << SEP << repeatXYZ("dirM");
    
    for ( Fiber const* fib = fibers.firstID(); fib; fib = fibers.nextID(fib) )
    {
        if ( sel && sel != fib->prop )
            continue;
        out << LIN << fib->prop->number();
        out << SEP << fib->identity();
        out << SEP << fib->length();
        if ( end & PLUS_END )
        {
            out << SEP << fib->endStateP();
            out << SEP << fib->posEndP();
            out << SEP << fib->dirEndP();
        }
        if ( end & MINUS_END )
        {
            out << SEP << fib->endStateM();
            out << SEP << fib->posEndM();
            out << SEP << fib->dirEndM();
        }
    }
}


/**
 Export Fiber-number, position of vertices
 */
void Simul::reportFiberPoints(std::ostream& out, Property const* sel) const
{
    out << COM << "identity" << SEP << repeatXYZ("pos") << SEP << "curvature";

    // list fibers in the order of the inventory:
    for ( Fiber const* fib = fibers.firstID(); fib; fib = fibers.nextID(fib) )
    {
        if ( sel && sel != fib->prop )
            continue;
        out << COM << "fiber " << fib->reference() << "  " << fib->segmentation();
        
        for ( size_t p = 0; p < fib->nbPoints(); ++p )
        {
            out << LIN << fib->identity();
            out << SEP << fib->posPoint(p);
            out << SEP << fib->curvature(p);
        }
    }
}


/**
 Export positions of points taken randomly along all fibers,
 but that remain static with respect to the lattice of each fiber,
 during the life-time of this fiber.
 
 This is meant to simulate the `speckle microscopy` that is obtained
 in microcscopy with a low amount of fluorescent-monomers.
 
 The distance between the speckles follows an exponential distribution
 with an average defined by the parameter `spread`.
 */
void Simul::reportFiberSpeckles(std::ostream& out, Glossary& opt) const
{
    real spread = 1;
    if ( opt.set(spread, "density") )
        spread = 1.0 / spread;
    else
        opt.set(spread, "interval");
    constexpr real TINY = 0x1p-32;

    Fiber const* fib = fibers.first();
    while ( fib )
    {
        out << COM << "fiber " << fib->reference();
        
        // generate speckles below the origin of abscissa
        if ( fib->abscissaM() < 0 )
        {
            uint64_t Z = pcg32_init(fib->signature());
            real a = spread * std::log(pcg32(Z)*TINY);
            while ( a > fib->abscissaP() )
            {
                a += spread * std::log(pcg32(Z)*TINY);
            }
            while ( a >= fib->abscissaM() )
            {
                out << '\n' << fib->pos(a);
                a += spread * std::log(pcg32(Z)*TINY);
            }
        }
        // generate speckles above the origin of abscissa
        if ( fib->abscissaP() > 0 )
        {
            uint64_t Z = pcg32_init(~fib->signature());
            real a = -spread * std::log(pcg32(Z)*TINY);
            while ( a < fib->abscissaM() )
                a -= spread * std::log(pcg32(Z)*TINY);
            while ( a <= fib->abscissaP() )
            {
                out << '\n' << fib->pos(a);
                a -= spread * std::log(pcg32(Z)*TINY);
            }
        }
        
        fib = fib->next();
    }
}


/**
 Export positions of points taken randomly along all fibers,
 changing the distribution at every time.
 */
void Simul::reportFiberSamples(std::ostream& out, Glossary& opt) const
{
    real spread = 1;
    if ( opt.set(spread, "density") )
        spread = 1.0 / spread;
    else
        opt.set(spread, "interval");
    
    Array<FiberSite> loc(1024, 1024);
    fibers.uniFiberSites(loc, spread);
    
    Fiber const* ofib = nullptr;
    for ( FiberSite & i : loc )
    {
        if ( ofib != i.fiber() )
        {
            out << COM << "fiber " << i.fiber()->reference();
            ofib = i.fiber();
        }
        
        out << LIN << i.pos();
    }
}


/**
 Export Mean Squared Displacement of fiber's MINUS_ENDs since the last call
 to this function.
 \todo: it would be simpler to add a 'Vector oldPos;' directly in class Fiber.
 */
void Simul::reportFiberDisplacement(std::ostream& out, Property const* sel) const
{
    typedef std::map <ObjectID, Vector> fiber_map;
    static fiber_map positions;
    static double past = 0;
    
    out << COM << "delta_time nb_fibers mean_squared_displacement";
    
    real sum = 0;
    size_t cnt = 0;
    for ( Fiber const* fib = fibers.first(); fib; fib = fib->next() )
    {
        if ( sel && sel != fib->prop )
            continue;
        /* using minus end to avoid effects of plus-end growth,
         but there is no check that the minus-end is not growing... */
        Vector pos = fib->posEndM();
        fiber_map::iterator i = positions.find(fib->identity());
        if ( i != positions.end() )
        {
            ++cnt;
            sum += distanceSqr(pos, i->second);
            i->second = pos;
        }
        else
        {
            positions[fib->identity()] = pos;
        }
    }
    
    if ( cnt > 0 )
        out << LIN << time() - past << SEP << cnt << SEP << sum / (real)cnt;
    else
        out << LIN << time() - past << SEP << 0 << SEP << 0;
    
    past = time();
}


/** This is a hack to maintain some credibility */
static bool isSymmetricAroundAxisZ(std::string const& shape)
{
    if ( shape == "sphere" ) return true;
    if ( shape == "cylinderZ" ) return true;
    if ( shape == "ellipse" ) return true;
    if ( shape == "torus" ) return true;
    if ( shape == "ring" ) return true;
    if ( shape == "disc" ) return true;
    return false;
}

void Simul::reportFiberDirections(std::ostream& out, Property const* sel) const
{
    Space const* spc = spaces.master();
    if ( sel )
        spc = static_cast<FiberProp const*>(sel)->confine_space_ptr;
    if ( !isSymmetricAroundAxisZ(spc->prop->shape) )
        throw InvalidParameter("reportFiberDirections() cannot handle non symmetric Space");

    real sum = 0;
    Vector eZ(0, 0, 1);
    Vector2 avg(0, 0);
    Matrix22 mat(0, 0);
#if ( DIM == 3 )
    for ( Fiber const* fib = fibers.first(); fib; fib = fib->next() )
    {
        if ( sel && sel != fib->prop )
            continue;
        for ( size_t p = 0; p < fib->nbSegments(); ++p )
        {
            Vector pos = fib->midPoint(p, 0.5);
            Vector dir = fib->dirSegment(p);
            Vector nor = spc->normalToEdge(pos);
            Vector tan = cross(eZ, nor);
            real n = normSqr(tan);
            if ( n > 0.5 )
            {
                Vector2 V(dot(dir, tan/sqrt(n)), dir.ZZ);
                V.normalize();
                mat(0,0) += V.XX * V.XX;
                mat(1,0) += V.XX * V.YY;
                mat(1,1) += V.YY * V.YY;
                avg += V;
                ++sum;
            }
        }
    }
#endif
    real X = 0, Y = 0, S = 0;
    if ( sum > 0 )
    {
        avg /= sum;
        mat *= 2.0 / sum;
        mat(0,1) = mat(1,0);
        X = std::sqrt(mat(0,0) * 0.5);
        Y = std::sqrt(mat(1,1) * 0.5);
        // subtract trace:
        mat(0,0) -= 1.0;
        mat(1,1) -= 1.0;
        //std::clog << mat << " ";
        // eigenvalue of a 2x2 traceless symmetric matrix:
        S = sqrt(square(mat(0,0))+square(mat(1,0)));
    }
    // polar order parameter:
    real M = norm(avg);
    out << COM << "n_seg" << SEP << "nematic" << SEP << "polar" << SEP << "avg_t" << SEP << "avg_z" << SEP << "mom_t" << SEP << "mom_z";
    out << LIN << sum << SEP << S << SEP << M << SEP << avg.XX << SEP << avg.YY << SEP << X << SEP << Y;
}


/**
 Export first and second-order moments of vertices for each class of Fiber
 */
void Simul::reportFiberMoments(std::ostream& out) const
{
    out << COM << ljust("class", 2, 2) << SEP << "sum";
    out << SEP << "avgX" << SEP << "avgY" << SEP << "avgZ";
    out << SEP << "varX" << SEP << "varY" << SEP << "varZ" << SEP << "var_sum";
    out << std::fixed;
    
    Accumulator acc;
    
    for ( Property const* i : properties.find_all("fiber") )
    {
        FiberProp const* fp = static_cast<FiberProp const*>(i);
        
        acc.reset();
        
        for ( Fiber const* fib = fibers.first(); fib; fib = fib->next() )
        {
            if ( fib->prop == fp )
            {
                const real w = fib->segmentation();
                acc.add(0.5*w, fib->posEndM());
                for ( size_t n = 1; n < fib->lastPoint(); ++n )
                    acc.add(w, fib->posPoint(n));
                acc.add(0.5*w, fib->posEndP());
            }
        }
        
        acc.subtract_mean();
        out << LIN << ljust(fp->name(), 2);
        acc.print(out, 0);
    }
}

/**
 Export first and second-order moments of vertices for each class of Fiber
 */
void Simul::reportNetworkSize(std::ostream& out) const
{
    out << COM << "polymer" << SEP << "surface";
    out << std::fixed;
    Accumulator acc;
    for ( Fiber const* fib = fibers.first(); fib; fib = fib->next() )
    {
        const real w = fib->segmentation();
        acc.add(0.5*w, fib->posEndM());
        for ( size_t n = 1; n < fib->lastPoint(); ++n )
            acc.add(w, fib->posPoint(n));
        acc.add(0.5*w, fib->posEndP());
    }
    acc.subtract_mean();
    real S = 2 * M_PI * acc.total_variance();
    out << LIN << acc.total_length() << SEP << S;
}

//------------------------------------------------------------------------------
#pragma mark - Fiber forces

/**
 Export Fiber-number, position of vertices and tension in each segment
 */
void Simul::reportFiberForces(std::ostream& out) const
{
    computeForces();

    out << COM << "identity" << SEP << repeatXYZ("pos") << SEP << repeatXYZ("force") << SEP << "tension";
    
    // list fibers in the order of the inventory:
    for ( Fiber const* fib = fibers.firstID(); fib; fib = fibers.nextID(fib) )
    {
        out << COM << "fiber " << fib->reference();
        
        for ( size_t p = 0; p < fib->nbPoints(); ++p )
        {
            out << LIN << fib->identity();
            out << SEP << fib->posPoint(p);
            out << SEP << fib->netForce(p);
            if ( p == fib->lastPoint() )
                out << SEP << 0.0;
            else
                out << SEP << fib->tension(p);
        }
    }
}


/**
 Sum of the internal tensions from fiber segments that intersect a plane
 specified in `opt`.
 The plane is defined by <em> n.pos + a = 0 </em>

     plane = NORMAL, SCALAR

 */
void Simul::reportFiberTension(std::ostream& out, Glossary& opt) const
{
    computeForces();
    
    out << COM << "count" << SEP << "sum_force" << SEP << "min_force" << SEP << "max_force";

    Vector n(1,0,0);
    real ten = 0, inf = 0, sup = 0;
    size_t cnt = 0;
    if ( opt.value_is("plane", 0, "all") )
    {
        // extending the comments:
        for ( int d = 1; d < DIM; ++d )
            out << SEP << "count" << SEP << "force";
        
        // plane orthogonal to X:
        fibers.infoTension(cnt, ten, inf, sup, Vector(1,0,0), 0);
        out << LIN << cnt << SEP << ten << SEP << inf << SEP << sup;
#if ( DIM > 1 )
        // plane orthogonal to Y:
        fibers.infoTension(cnt, ten, inf, sup, Vector(0,1,0), 0);
        out << SEP << cnt << SEP << ten << SEP << inf << SEP << sup;
#endif
#if ( DIM > 2 )
        // plane orthogonal to Z:
        fibers.infoTension(cnt, ten, inf, sup, Vector(0,0,1), 0);
        out << SEP << cnt << SEP << ten << SEP << inf << SEP << sup;
#endif
    }
    else if ( opt.set(n, "plane") )
    {
        real a = 0;
        opt.set(a, "plane", 1);
        out << COM << "fiber tension orthogonal to plane: (" << n << ").pos = " << -a;
        fibers.infoTension(cnt, ten, inf, sup, n, a);
        out << LIN << cnt << SEP << ten << SEP << inf << SEP << sup;
    }
    else
    {
        // if no plane is specified, sum all tension from all segments
        fibers.infoTension(cnt, ten, inf, sup);
        out << LIN << cnt << SEP << ten << SEP << inf << SEP << sup;
    }
}


/** Attention: This does not handle all cases */
static bool confinementApplies(Confinement mode, Space const* spc, Vector const& pos)
{
    switch ( mode )
    {
        case CONFINE_OFF: return false;
        case CONFINE_INSIDE: return spc->outside(pos);
        case CONFINE_OUTSIDE: return spc->inside(pos);
        case CONFINE_ON: return true;
        case CONFINE_PLUS_END: return true;
        case CONFINE_MINUS_END: return true;
        case CONFINE_BOTH_ENDS: return true;
        case CONFINE_PLUS_OUT: return spc->inside(pos);
        default: return true;
    }
    return true;
}

/** Attention: This does not handle all cases */
static bool vertexIsConfined(Confinement mode, Fiber const* fib, size_t inx)
{
    switch ( mode )
    {
        case CONFINE_OFF: return false;
        case CONFINE_INSIDE: return true;
        case CONFINE_OUTSIDE: return true;
        case CONFINE_ON: return true;
        case CONFINE_PLUS_END: return ( inx == fib->lastPoint() );
        case CONFINE_MINUS_END: return ( inx == 0 );
        case CONFINE_BOTH_ENDS: return ( inx == 0 || inx == fib->lastPoint() );
        case CONFINE_PLUS_OUT: return ( inx == fib->lastPoint() );
        default: return true;
    }
    return true;
}


/**
 Export total magnitude of force exerted by Fiber on the confinement
 */
void Simul::reportFiberConfineForce(std::ostream& out) const
{
    out << COM << "confinement forces";
    out << COM << "identity" << SEP << "vertex" << SEP << repeatXYZ("pos") << SEP << repeatXYZ("force");
     
     // list fibers in the order of the inventory:
     for ( Fiber const* fib = fibers.firstID(); fib; fib = fibers.nextID(fib) )
     {
         out << COM << "fiber " << fib->reference();
         Space const* spc = findSpace(fib->prop->confine_space);
         const real stiff = fib->prop->confine_stiff[0];
         const Confinement mode = fib->prop->confine;

         for ( size_t p = 0; p < fib->nbPoints(); ++p )
         {
             Vector w, pos = fib->posPoint(p);
             if ( vertexIsConfined(mode, fib, p) && confinementApplies(mode, spc, pos) )
             {
                 w = spc->project(pos);
                 out << LIN << fib->identity();
                 out << SEP << p << SEP << pos;
                 out << SEP << stiff * ( w - pos );
             }
         }
     }
}

/**
 Export total magnitude of force exerted by Fiber on the confinement.
 The radial components of the forces are summed up, which is only meaningful
 in very particular systems, when the geometry is circular around the Z-axis!
 */
real Simul::reportFiberConfinement(std::ostream& out) const
{
    out << COM << "count" << SEP << repeatXYZ("force") << SEP << "radial";
    size_t cnt = 0;
    Vector sum(0,0,0);
    real   rad = 0;
    
#if ( DIM > 1 )
    for ( Fiber const* fib = fibers.first(); fib; fib = fib->next() )
    {
        Space const* spc = findSpace(fib->prop->confine_space);
        const real stiff = fib->prop->confine_stiff[0];
        const Confinement mode = fib->prop->confine;

        if ( !isSymmetricAroundAxisZ(spc->prop->shape) )
            throw InvalidParameter("reportFiberConfinement() cannot handle non symmetric Space");
        
        for ( size_t p = 0; p < fib->nbPoints(); ++p )
        {
            Vector w, pos = fib->posPoint(p);
            if ( vertexIsConfined(mode, fib, p) && confinementApplies(mode, spc, pos) )
            {
                ++cnt;
                w = spc->project(pos);
                // assuming the goemetry is rotational-symmetric around the Z axis:
                Vector dir = normalize(Vector(pos.XX, pos.YY, 0));
                Vector vec = stiff * ( pos - w );
                sum += vec;
                rad += dot(vec, dir);
            }
        }
    }
#endif
    out << LIN << cnt << SEP << sum << SEP << rad;
    return rad;
}


//------------------------------------------------------------------------------
#pragma mark - bound Hands per Fiber


void Simul::reportFiberHands(std::ostream& out) const
{
    out << COM << "fib_type" << SEP << "fib_id" << SEP << "class" << SEP << "abs";
    for ( Fiber const* fib = fibers.firstID(); fib; fib = fibers.nextID(fib) )
    {
        if ( fib->nbHands() > 0 )
        {
            out << COM << "on fiber " << fib->reference();
            fib->sortHands();
            for ( Hand * ha = fib->firstHand(); ha; ha = ha->next() )
            {
                out << LIN << fib->prop->number();
                out << SEP << fib->identity();
                out << SEP << ha->prop->number();
                out << SEP << ha->abscissa();
            }
        }
    }
}


void Simul::reportFiberLinks(std::ostream& out) const
{
    out << COM << "fib_type" << SEP << "fib_id" << SEP << "class" << SEP << "abs" << SEP << "position";
    for ( Fiber const* fib = fibers.firstID(); fib; fib = fibers.nextID(fib) )
    {
        if ( fib->nbHands() > 0 )
        {
            out << COM << "on fiber " << fib->reference();
            fib->sortHands();
            for ( Hand const* ha = fib->firstHand(); ha; ha = ha->next() )
            {
                if ( ha->linkStiffness() > 0 )
                {
                    out << LIN << fib->prop->number();
                    out << SEP << fib->identity();
                    out << SEP << ha->prop->number();
                    out << SEP << ha->abscissa();
                    out << SEP << ha->linkFoot();
                    out << SEP << ha->linkStiffness();
                }
            }
        }
    }
}

//------------------------------------------------------------------------------
#pragma mark - Networks


void Simul::reportFiberIntersections(std::ostream& out, Glossary& opt) const
{
    int details = 2;
    real up = 0;
    opt.set(up, "threshold");
    opt.set(details, "details");
    
    const real sup = up * up;
    
    if ( details == 2 )
    {
        out << COM << "id1" << SEP << "abs1";
        out << SEP << "id2" << SEP << "abs2" << SEP << repeatXYZ("pos");
    }
    Accumulator acc;
    
    for ( Fiber const* fib = fibers.firstID(); fib; fib = fibers.nextID(fib) )
    {
        unsigned cnt = 0;
        for ( Fiber const* fox = fibers.nextID(fib); fox; fox = fibers.nextID(fox) )
        {
            for ( size_t ii = 0; ii < fib->nbSegments(); ++ii )
            {
                FiberSegment seg(fib, ii);
                for ( size_t jj = 0; jj < fox->nbSegments(); ++jj )
                {
                    FiberSegment soc(fox, jj);
                    real abs1, abs2;
                    real dis2 = seg.shortestDistanceSqr(soc, abs1, abs2);
                    if (( dis2 < sup ) & seg.within(abs1) & soc.within(abs2))
                    {
                        ++cnt;
                        Vector pos1 = seg.pos(abs1/seg.len());
                        //Vector pos2 = loc2.pos(abs2/loc2.len());
                        if ( details == 2 )
                        {
                            out << LIN << fib->identity();
                            out << SEP << abs1 + seg.abscissa1();
                            out << SEP << fox->identity();
                            out << SEP << abs2 + soc.abscissa1();
                            out << SEP << pos1;
                        }
                        acc.add(pos1);
                    }
                }
            }
        }
        if ( cnt && details >= 1 )
        {
            out << COM << "total";
            out << SEP << fib->identity();
            out << SEP << cnt;
        }
    }
    acc.subtract_mean();
    acc.print_doc(out);
    acc.print(out, 1);
}


/// accessory class to analyse the connections in a network of fibers
struct Connector
{
    real a;
    real s;
    long f;
    long g;
    long h;
    Connector(real as, long fs) { a = as; f = fs; g = -1; h = -1; s = 0; }
    Connector(real as, long fs, long gs) { a = as; f = fs; g = gs; h = -1; s = 0; }
    Connector(real as, long fs, long gs, long hs) { a = as; f = fs; g = gs; h = hs; s = 0; }
    Connector(real as, long fs, long gs, long hs, real ss) { a = as; f = fs; g = gs; h = hs; s = ss; }
    bool operator < (const Connector& r) const { return a < r.a; }
};


/**
 This is the older version
 */
void Simul::reportFiberConnectors(std::ostream& out, Glossary& opt) const
{
    int details = 2;
    opt.set(details, "details");

    if ( details > 1 )
    {
        out << COM << "class  identity      abs1    fiber1     hand1      dist";
        out << "      abs2    fiber2     hand2      dist ...";
    }
    else
    {
        out << COM << "fiber connectors";
    }
    
    // used to calculate the size of the network from the position of connectors
    Accumulator acc;
    typedef std::vector<Connector> clist_t;
    typedef std::map<ObjectID, clist_t> map_t;
    
    map_t map;

    for ( Fiber const* fib = fibers.firstID(); fib; fib = fibers.nextID(fib) )
    {
        map.clear();
        // check all connecting Hands and record abscissa, depending on the fiber that is linked
        for ( Hand const* ha = fib->firstHand(); ha; ha = ha->next() )
        {
            Hand const* oh = ha->monitor()->otherHand(ha);
            if ( oh && oh->attached() )
            {
                ObjectID f2 = oh->fiber()->identity();
                map[f2].push_back(Connector(ha->abscissa(), ha->prop->number()));
            }
        }
        if ( map.size() )
        {
            if ( details > 1 )
            {
                out << LIN << fib->prop->number() << SEP << fib->identity();
            }
            
            clist_t list;
            // average all the abscissa linking to the same fiber:
            for ( map_t::const_iterator mi = map.begin(); mi != map.end(); ++mi )
            {
                clist_t const& sec = mi->second;
                real a = 0.0;
                // number of connector of each type
                int c1 = 0, c2 = 0;
                for ( clist_t::const_iterator ci = sec.begin(); ci != sec.end(); ++ci )
                {
                    a += ci->a;
                    if ( ci->f == 1 ) ++c1;
                    if ( ci->f == 2 ) ++c2;
                 }
                a /= sec.size();
                list.push_back(Connector(a, mi->first, c1, c2, sec.size()-c1-c2));
            }
            // sort the list in increasing abscissa
            std::sort(list.begin(), list.end());
            
            clist_t::const_iterator p = list.begin();
            for ( clist_t::const_iterator c = list.begin(); c != list.end(); ++c )
            {
                if ( details > 1 )
                {
                    out << SEP << c->a;
                    out << SEP << c->f;
                    out << SEP << std::to_string(c->g)+"+"+std::to_string(c->h);
                    // calculate direct distance to previous point of intersection:
                    if ( c != p )
                        out << SEP << ( fib->pos(p->a) - fib->pos(c->a) ).norm();
                    else
                        out << SEP << 0;
                }
                p = c;
                acc.add(fib->pos(p->a));
            }
            if ( details > 0 )
            {
                out << COM << "total";
                out << SEP << fib->prop->number();
                out << SEP << fib->identity();
                out << SEP << list.size();
            }
        }
    }
    acc.subtract_mean();
    acc.print_doc(out);
    acc.print(out, 1);
}

/**
 F. Nedelec, 18/08/2017
 */
void Simul::reportNetworkBridges(std::ostream& out, Glossary& opt) const
{
    int details = 0;
    opt.set(details, "details");

    out << COM << "length" << SEP << "speed" << SEP << "type1" << SEP << "type2";

    typedef std::vector<Connector> clist_t;
    typedef std::map<ObjectID, clist_t> map_t;
    
    map_t map;
    
    HandProp const* hp1 = findProperty<HandProp>("hand", 1);
    HandProp const* hp2 = findProperty<HandProp>("hand", 2);
    
    const real speedh1 = hp1->motorSpeed();
    const real speedh2 = hp2->motorSpeed();

    for ( Fiber const* fib = fibers.firstID(); fib; fib = fibers.nextID(fib) )
    {
        map.clear();
        // check all connecting Hands and record abscissa, depending on the fiber that is linked
        for ( Hand * ha = fib->firstHand(); ha; ha = ha->next() )
        {
            Hand const* oh = ha->monitor()->otherHand(ha);
            if ( oh && oh->attached() )
            {
                ObjectID f2 = oh->fiber()->identity();
                if ( ha->prop == hp1 )
                    map[f2].push_back(Connector(ha->abscissa(), 1));
                else if ( ha->prop == hp2 )
                    map[f2].push_back(Connector(ha->abscissa(), 2));
                else
                    out << COM << "report network:bridge can only handle 2 hand types";
            }
        }
        if ( map.size() )
        {
            clist_t list;
            // average all the abscissa linking to the same fiber:
            for ( map_t::const_iterator mi = map.begin(); mi != map.end(); ++mi )
            {
                clist_t const& sec = mi->second;
                // average abscissa:
                real a = 0.0;
                // number of connector of each type
                int c1 = 0, c2 = 0;
                for ( clist_t::const_iterator ci = sec.begin(); ci != sec.end(); ++ci )
                {
                    a += ci->a;
                    if ( ci->f == 1 ) ++c1;
                    if ( ci->f == 2 ) ++c2;
                }
                a /= sec.size();
                real speed = 0;
                if ( c1 > 0 && c2 > 0 )
                    speed = std::min(speedh1, speedh2);
                else if ( c1 > 0 )
                    speed = speedh1;
                else if ( c2 > 0 )
                    speed = speedh2;
                list.push_back(Connector(a, mi->first, c1, c2, speed));
            }
            // sort the list in increasing abscissa
            std::sort(list.begin(), list.end());
            
            if ( details > 0 )
            {
                out << COM << "connectors on fiber f" << fib->identity() << ":";
                // print all connector attachment positions:
                out << COM << "abscissa" << SEP << "fiber_id" << SEP << "speed" << SEP << "type";
                for ( clist_t::const_iterator c = list.begin(); c != list.end(); ++c )
                {
                    out << LIN << c->a;
                    out << SEP << c->f;
                    out << SEP << c->s;
                    out << SEP << std::to_string(c->g)+"+"+std::to_string(c->h);
                }
            }
            if ( list.size() > 1 )
            {
                // print all bridges
                if ( details > 0 )
                    out << COM << "length" << SEP << "speed" << SEP << "type1" << SEP << "type2";
#if ( 1 )
                for ( clist_t::const_iterator p = list.begin(); p != list.end(); ++p )
                for ( clist_t::const_iterator c = p+1; c != list.end(); ++c )
#else
                for ( clist_t::const_iterator p = list.begin(), c = p+1; c != list.end(); ++p, ++c )
#endif
                {
                    out << LIN << c->a - p->a;
                    out << SEP << c->s - p->s;
                    out << SEP << std::to_string(p->g)+"+"+std::to_string(p->h);
                    out << SEP << std::to_string(c->g)+"+"+std::to_string(c->h);
                }
            }
        }
    }
}


//------------------------------------------------------------------------------
#pragma mark - Beads, Solid, Space


void Simul::reportTime(std::ostream& out) const
{
    out << LIN << time();
}


void Simul::reportInventory(std::ostream& out) const
{
    //out << COM << "properties:";
    //properties.write_names(out, "");
    //out << COM << "objects:";
    spaces.report(out);
    fields.report(out);
    fibers.report(out);
    spheres.report(out);
    beads.report(out);
    solids.report(out);
    singles.report(out);
    couples.report(out);
    organizers.report(out);
    tubules.report(out);
    events.report(out);
}

template < typename SET >
void reportSimulSet(std::ostream& out, SET& set, PropertyList const& properties)
{
    for ( Property const* i : properties.find_all(set.title()) )
    {
        ObjectID id = 0;
        size_t points = 0, sup = 0;
        ObjectList objs = set.collect(match_property, i);
        for ( Object * o : objs )
        {
            Mecable * mec = Simul::toMecable(o);
            if ( mec )
            {
                points += mec->nbPoints();
                sup = std::max(sup, mec->nbPoints());
                id = std::max(id, mec->identity());
            }
        }
        if ( points > 0 )
        {
            out << LIN << ljust(i->name(), 2);
            out << SEP << objs.size() << SEP << points;
            out << SEP << sup << SEP << id;
        }
    }
}

void Simul::reportSimul(std::ostream& out) const
{
    out << COM << ljust("class", 2, 2) << SEP << "objects" << SEP << "vertices";
    out << SEP << "largest" << SEP << "identity";
    reportSimulSet(out,  fibers, properties);
    reportSimulSet(out,  solids, properties);
    reportSimulSet(out, spheres, properties);
    reportSimulSet(out,   beads, properties);
}

/**
 Export position of all organizers
 */
void Simul::reportOrganizer(std::ostream& out) const
{
    out << COM << "class" << SEP << "identity" << SEP << repeatXYZ("pos");

    for ( Organizer const* obj=organizers.first(); obj; obj=obj->next() )
    {
        out << LIN << obj->property()->number();
        out << SEP << obj->identity();
        out << SEP << obj->position();
        out << SEP << obj->nbOrganized();
    }
}


/**
 Export position of Asters
 */
void Simul::reportAster(std::ostream& out) const
{
    out << COM << "class" << SEP << "identity" << SEP << repeatXYZ("pos");
    
    for ( Organizer const* obj=organizers.first(); obj; obj=obj->next() )
    {
        if ( obj->tag() == Organizer::TAG_ASTER )
        {
            out << LIN << obj->property()->number();
            out << SEP << obj->identity();
            out << SEP << obj->position();
        }
    }
}


/**
 Export position of Beads
 */
void Simul::reportBeadPosition(std::ostream& out, Property const* sel) const
{
    out << COM << "class" << SEP << "identity" << SEP << repeatXYZ("pos");
    
    for ( Bead const* obj=beads.first(); obj; obj=obj->next() )
    {
       if ( sel && sel != obj->prop )
            continue;
        out << LIN << obj->prop->number();
        out << SEP << obj->identity();
        out << SEP << obj->position();
    }
}


/**
 Export number of beads classified as a function of
 the number of grafted Single that are attached to Fibers
 */
void Simul::reportBeadSingles(std::ostream& out) const
{
    out << COM << "identity" << "amount(nb_attached_hands)";
    
    std::map<ObjectID, int> cnt;
    
    for ( Single const* i=singles.firstA(); i; i=i->next() )
    {
        Bead const* obj = Bead::toBead(i->base());
        if ( obj )
            ++cnt[ obj->identity() ];
    }

    const int max = 12;
    int nb[max] = { 0 };
    for ( Bead const* obj=beads.first(); obj; obj=obj->next() )
        ++nb[ cnt[obj->identity()] ];
    
    for ( int c = 0; c < max; ++c )
        out << " " << std::setw(3) << nb[c];
}


/**
 Export position of Solids
 */
void Simul::reportSolidPosition(std::ostream& out, Property const* sel) const
{
    out << COM << "class" << SEP << "identity" << SEP << repeatXYZ("cen");
    out << SEP << repeatXYZ("point0") << SEP << repeatXYZ("point1");
        
    for ( Solid const* obj=solids.first(); obj; obj=obj->next() )
    {
        if ( sel && sel != obj->prop )
            continue;
        out << LIN << obj->prop->number();
        out << SEP << obj->identity();
        out << SEP << obj->centroid();
        out << SEP << obj->posPoint(0);
        if ( obj->nbPoints() > 1 )
        out << SEP << obj->posPoint(1);
        
        if ( modulo )
        {
            Vector pos = obj->centroid();
            modulo->fold(pos);
            out << SEP << pos;
        }
    }
}

/**
 Export orientation of Solids
 */
void Simul::reportSolidOrientation(std::ostream& out, Property const* sel) const
{
    out << COM << "class" << SEP << "identity" << SEP << repeatXYZ("cen");
    out << SEP << repeatXYZ("dir");
        
    for ( Solid const* obj=solids.first(); obj; obj=obj->next() )
    {
        if ( sel && sel != obj->prop )
            continue;
        out << LIN << obj->prop->number();
        out << SEP << obj->identity();
        out << SEP << obj->centroid();
        out << SEP << obj->orientation();
    }
}

/**
 Export position of Solids with counts of Hands and attached Hands
 */
void Simul::reportSolidHands(std::ostream& out, Property const* sel) const
{
    out << COM << "class" << SEP << "identity" << SEP << repeatXYZ("pos");
    out << SEP << "nb_hand" << SEP << "nb_link";
        
    for ( Solid const* obj = solids.firstID(); obj; obj = solids.nextID(obj) )
    {
        if ( sel && sel != obj->prop )
            continue;
        out << LIN << obj->prop->number();
        out << SEP << obj->identity();
        Vector pos = obj->centroid();
        if ( modulo ) modulo->fold(pos);
        out << SEP << pos;
        SingleList anchored = singles.collectWrists(obj);
        int cnt = 0;
        for ( Single const* s : anchored )
            cnt += s->attached();
        out << SEP << anchored.size() << SEP << cnt;
    }
}


/**
 Report position of Sphere
 */
void Simul::reportSpherePosition(std::ostream& out, Property const* sel) const
{
    out << COM << "class" << SEP << "identity";
    out << SEP << repeatXYZ("point0") << SEP << repeatXYZ("point1");
        
    for ( Sphere const* obj=spheres.first(); obj; obj=obj->next() )
    {
        if ( sel && sel != obj->prop )
            continue;
        out << LIN << obj->prop->number();
        out << SEP << obj->identity();
        out << SEP << obj->posPoint(0);
        if ( obj->nbPoints() > 1 )
        out << SEP << obj->posPoint(1);
    }
}


/**
 Report something about Space (incomplete)
 */
void Simul::reportSpace(std::ostream& out) const
{
    out << COM << "class" << SEP << "identity";
    
    for ( Space const* obj=spaces.firstID(); obj; obj=spaces.nextID(obj) )
    {
        out << LIN << obj->prop->name();
        out << SEP << obj->identity();
        out << SEP << std::fixed << obj->prop->shape;
    }
}


/**
 Report force on Space
 */
void Simul::reportSpaceForce(std::ostream& out) const
{
    out << COM << "class" << SEP << "identity" << SEP << "shape";
    
    for ( Space const* obj=spaces.first(); obj; obj=obj->next() )
    {
        out << LIN << obj->prop->name();
        out << SEP << obj->identity();
        out << SEP << obj->prop->shape;
        obj->report(out);
    }
}


/**
 Report quantity of substance in Field
 */
void Simul::reportField(std::ostream& out) const
{
    if ( fields.size() == 0 )
        return;
    out << COM << ljust("class", 2, 2);
    out << SEP << "total" << SEP << "avg" << SEP << "min" << SEP << "max";
    
    FieldCell s, a, n, x;
    // report total substance in each Field
    for ( Field const* obj=fields.first(); obj; obj=obj->next() )
    {
        const real alpha = 1.0 / obj->cellVolume();
        obj->infoValues(s, a, n, x);
        out << LIN << ljust(obj->prop->name(), 2);
        out << SEP << s;
        out << SEP << a * alpha;
        out << SEP << n * alpha;
        out << SEP << x * alpha;
    }
}


//------------------------------------------------------------------------------
#pragma mark - Single

void write(std::ostream& out, Single const* obj, Simul const* simul)
{
    out << LIN << obj->prop->number();
    out << SEP << obj->identity();
    out << SEP << obj->position();
    Fiber const* fib = obj->fiber();
    if ( fib )
    {
        out << SEP << obj->force();
        out << SEP << fib->identity();
        out << SEP << obj->abscissa();
        out << SEP << simul->organizers.findOrganizerID(fib);
    }
    else
    {
        out << SEP << Vector(0,0,0);
        out << SEP << "0";
        out << SEP << "nan";
        out << SEP << "0";
    }
}


/**
 Export details of Singles, possiby selecting for a certain kind
 */
void Simul::reportSingleState(std::ostream& out, Property const* sel) const
{
    out << COM << "class" << SEP << "identity";
    out << SEP << repeatXYZ("pos") << SEP << repeatXYZ("force");
    out << SEP << "fiber" << SEP << "abscissa" << SEP << "aster";
    
    for ( Single const* obj = singles.firstID(); obj; obj = singles.nextID(obj) )
        if ( !sel || sel == obj->prop )
            write(out, obj, this);
}


/**
 Export details of attached Singles
 */
void Simul::reportSinglePosition(std::ostream& out, Property const* sel) const
{
    out << COM << "class" << SEP << "identity" << SEP << repeatXYZ("pos");
    out << SEP << "fiber" << SEP << "abscissa";
        
    for ( Single const* obj = singles.firstID(); obj; obj = singles.nextID(obj) )
    {
        if ( sel && sel != obj->prop )
            continue;
        out << LIN << obj->prop->number();
        out << SEP << obj->identity();
        out << SEP << obj->posFoot();
        if ( obj->attached() )
        {
            out << SEP << obj->fiber()->identity();
            out << SEP << obj->abscissa();
        }
        else
        {
            out << SEP << 0 << SEP << 0;
        }
    }
}

/**
 Export details of attached Singles
 */
void Simul::reportSingleLink(std::ostream& out, Property const* sel) const
{
    out << COM << "class" << SEP << "identity";
    out << SEP << repeatXYZ("pos") << SEP << repeatXYZ("force");
    out << SEP << "fiber" << SEP << "abscissa" << SEP << "aster";
        
    for ( Single const* obj = singles.firstID(); obj; obj = singles.nextID(obj) )
    {
        if ( obj->attached()  && ( !sel || sel == obj->prop ))
            write(out, obj, this);
    }
}


/**
 Export number of Single in each state
 */
void Simul::reportSingle(std::ostream& out, Property const* sel) const
{
    constexpr size_t SUP = 128;
    
    size_t free[SUP+1] = { 0 }, bound[SUP+1] = { 0 }, based[SUP+1] = { 0 };
    
    for ( Single const* i = singles.firstF(); i ; i = i->next() )
    {
        assert_true(!i->attached());
        size_t x = i->prop->number();
        if ( x < SUP )
        {
            ++free[x];
            based[x] += ( i->base() != nullptr );
        }
    }
    
    for ( Single const* i=singles.firstA(); i ; i=i->next() )
    {
        assert_true(i->attached());
        size_t x = i->prop->number();
        if ( x < SUP )
        {
            ++bound[x];
            based[x] += ( i->base() != nullptr );
        }
    }
    
    out << COM << ljust("single", 2, 2);
    out << SEP << "total";
    out << SEP << "free";
    out << SEP << "bound";
    out << SEP << "based";
    
    for ( Property const* i : properties.find_all("single") )
    {
        if ( sel && i != sel )
            continue;
        out << LIN << ljust(i->name(), 2);
        size_t x = i->number();
        if ( x < SUP )
        {
            out << SEP << free[x] + bound[x];
            out << SEP << free[x];
            out << SEP << bound[x];
            out << SEP << based[x];
        }
        else
            out << SEP << " out-of-range ";
    }
}


/**
 Export average properties of Couples forces
 */
void Simul::reportSingleForce(std::ostream& out, Property const* sel) const
{
    constexpr size_t MAX = 8;
    real cnt[MAX+1] = { 0 };
    real avg[MAX+1] = { 0 };
    real sup[MAX+1] = { 0 };
    real len[MAX+1] = { 0 };

    // accumulate counts:
    for ( Single const* i=singles.firstA(); i; i=i->next() )
    {
        if ( i->hasLink() && ( !sel || sel == i->prop ))
        {
            size_t x = i->prop->number();
            if ( x < MAX )
            {
                real f = i->force().norm();
                avg[x] += f;
                cnt[x] += 1;
                sup[x] = std::max(sup[x], f);
                len[x] = std::max(len[x], i->stretch().norm());
            }
        }
    }
    
    out << COM << ljust("single", 2, 2) << SEP << "avg_force" << SEP << "max_force" << SEP << "max_len";
    for ( size_t i = 0; i < MAX; ++i )
    {
        if ( cnt[i] > 0 )
        {
            Property const* p = properties.find_or_die("single", i);
            out << LIN << ljust(p->name(), 2);
            out << SEP << avg[i] / cnt[i];
            out << SEP << sup[i];
            out << SEP << len[i];
        }
    }
}


//------------------------------------------------------------------------------
#pragma mark - Couple

void Simul::reportCoupleList(std::ostream& out, Property const* sel) const
{
    out << COM << "couples" << SEP << "count";

    size_t cnt = couples.collect(match_property, sel).size();

    if ( sel )
        out << LIN << sel->name();
    else
        out << LIN << "all";

    out << SEP << cnt;
 }


void write(std::ostream& out, Couple const* obj)
{
    out << LIN << obj->prop->number();
    out << SEP << obj->identity();
    out << SEP << obj->active();
    out << SEP << obj->position();

    Fiber const* fib = obj->fiber1();
    if ( fib )
    {
        out << SEP << fib->identity();
        out << SEP << obj->abscissa1();
    }
    else
    {
        out << SEP << "0";
        out << SEP << "nan";
    }

    fib = obj->fiber2();
    if ( fib )
    {
        out << SEP << fib->identity();
        out << SEP << obj->abscissa2();
    }
    else
    {
        out << SEP << "0";
        out << SEP << "nan";
    }
}

        
/**
 Export position of Couples of a certain kind
 */
void Simul::reportCoupleState(std::ostream& out, Property const* sel) const
{
    out << COM << "class" << SEP << "identity" << SEP << "active" << SEP << repeatXYZ("pos");
    out << SEP << "fiber1" << SEP << "abscissa1" << SEP << "fiber2" << SEP << "abscissa2";
    
    for ( Couple const* obj=couples.firstFF(); obj ; obj=obj->next() )
        if ( !sel || sel == obj->prop )
            write(out, obj);
    
    for ( Couple const* obj=couples.firstAF(); obj ; obj=obj->next() )
        if ( !sel || sel == obj->prop )
            write(out, obj);
    
    for ( Couple const* obj=couples.firstFA(); obj ; obj=obj->next() )
        if ( !sel || sel == obj->prop )
            write(out, obj);
    
    for ( Couple const* obj=couples.firstAA(); obj ; obj=obj->next() )
        if ( !sel || sel == obj->prop )
            write(out, obj);
}


/**
 Export position of active Couples of a certain kind
 */
void Simul::reportCoupleActive(std::ostream& out, Property const* sel) const
{
    out << COM << "state" << SEP << repeatXYZ("pos");
    
    for ( Couple const* obj=couples.firstFF(); obj ; obj=obj->next() )
        if ( obj->active()  &&  obj->prop == sel )
            out << "\n 0 " << obj->position();
   
    for ( Couple const* obj=couples.firstAF(); obj ; obj=obj->next() )
        if ( obj->prop == sel )
            out << "\n 1 " << obj->position();
    
    for ( Couple const* obj=couples.firstFA(); obj ; obj=obj->next() )
        if ( obj->prop == sel )
            out << "\n 2 " << obj->position();
    
    for ( Couple const* obj=couples.firstAA(); obj ; obj=obj->next() )
        if ( obj->prop == sel )
            out << "\n 3 " << obj->position();
}


/**
 Export position and force of Couples that are bound to 2 filaments
 */
void Simul::reportCoupleLink(std::ostream& out, Property const* sel) const
{
    out << COM << "class" << SEP << "identity";
    out << SEP << "fiber1" << SEP << "abscissa1";// << SEP << repeatXYZ("pos1");
    out << SEP << "fiber2" << SEP << "abscissa2";// << SEP << repeatXYZ("pos2");
    out << SEP << "force" << SEP << "cos_angle";
        
    for ( Couple const* obj=couples.firstAA(); obj ; obj=obj->next() )
    {
        if ( sel && sel != obj->prop )
            continue;
        out << LIN << obj->prop->number();
        out << SEP << obj->identity();
        
        out << SEP << obj->fiber1()->identity();
        out << SEP << obj->abscissa1();
        //out << SEP << obj->posHand1();

        out << SEP << obj->fiber2()->identity();
        out << SEP << obj->abscissa2();
        //out << SEP << obj->posHand2();

        out << SEP << obj->force().norm();
        out << SEP << obj->cosAngle();
    }
}


/**
 Export configuration of bridging couple, as
 P: parallel side-side links
 A: antiparallel side-side links
 X: other side-side links
 T: side-end
 V: end-end
 
 T and V are defined with respect to the `end`, at distance `threshold`,
 both can be set as parameters.
 
 by Jamie Li Rickman for
 Determinants of Polar versus Nematic Organization in Networks of Dynamic Microtubules
 and Mitotic Motors, Cell 2018
 */
void Simul::reportCoupleConfiguration(std::ostream& out, Property const* sel,
                                      Glossary& opt) const
{
    real dis = 0.010;  // max distance to end to constitute T or V link
    opt.set(dis, "distance", "threshold");
    
    size_t T[8] = { 0 };
    for ( Couple const* obj=couples.firstAA(); obj ; obj=obj->next() )
    {
        if ( !sel || sel == obj->prop )
            ++T[obj->configuration(dis)];
    }
    size_t S = T[0]+T[1]+T[2]+T[3]+T[4]+T[5]+T[6];
    
    out << COM << "couples" << SEP << "Total" << SEP << "P" << SEP << "A";
    out << SEP << "X" << SEP << "T+" << SEP << "V+" << SEP << "T-" << SEP << "V-";

    if ( sel )
        out << LIN << sel->name();
    else
        out << LIN << "all";

    out << SEP << S << SEP << T[0] << SEP << T[1] << SEP << T[2] << SEP << T[3] << SEP << T[4] << SEP << T[5] << SEP << T[6];
 }



/**
 Export average properties of Couples forces
 */
void Simul::reportCoupleForce(std::ostream& out, Property const* sel) const
{
    constexpr size_t MAX = 8;
    real cnt[MAX+1] = { 0 };
    real avg[MAX+1] = { 0 };
    real sup[MAX+1] = { 0 };
    real len[MAX+1] = { 0 };

    // accumulate counts:
    for ( Couple const* i=couples.firstAA(); i ; i = i->next() )
    {
        if ( !sel || sel == i->prop )
        {
            size_t x = i->prop->number();
            if ( x < MAX )
            {
                real f = i->force().norm();
                avg[x] += f;
                cnt[x] += 1;
                sup[x] = std::max(sup[x], f);
                len[x] = std::max(len[x], i->stretch().norm());
            }
        }
    }
        
    out << COM << ljust("couple", 2) << SEP << "avg_force" << SEP << "max_force" << SEP << "max_len";
    for ( size_t i = 0; i < MAX; ++i )
    {
        if ( cnt[i] > 0 )
        {
            Property const* p = properties.find_or_die("couple", i);
            out << LIN << ljust(p->name(), 2);
            out << SEP << avg[i] / cnt[i];
            out << SEP << sup[i];
            out << SEP << len[i];
        }
    }
}


/**
 Export histogram of Couples forces
 */
void Simul::reportCoupleForceHistogram(std::ostream& out, Glossary& opt) const
{
    const size_t IMAX = 8;
    const size_t BMAX = 256;
    size_t cnt[IMAX][BMAX+1];

    real delta = 0.5;
    size_t nbin = 64;
    opt.set(delta, "interval");
    opt.set(nbin, "interval", 1);
    nbin = std::min(nbin, BMAX);

    // reset counts:
    for ( size_t ii = 0; ii <  IMAX; ++ii )
    for ( size_t jj = 0; jj <= nbin; ++jj )
        cnt[ii][jj] = 0;
    
    // accumulate counts:
    for ( Couple const* i=couples.firstAA(); i ; i = i->next() )
    {
        size_t x = i->prop->number();
        if ( x < IMAX )
        {
            unsigned f = (unsigned)( i->force().norm() / delta );
            if ( f < nbin )
                ++cnt[x][f];
            else
                ++cnt[x][nbin];
        }
    }
    
    if ( 1 )
    {
        std::streamsize p = out.precision();
        out.precision(3);
        out << COM << "force_distribution" << " (`scale` indicates the center of each bin)";
        out << LIN << ljust("scale", 2);
        for ( size_t u = 0; u <= nbin; ++u )
            out << " " << std::setw(5) << delta * ( u + 0.5 );
        out.precision(p);
    }
    
    for ( size_t ii = 0; ii < IMAX; ++ii )
    {
        size_t sum = 0;
        for ( size_t jj = 0; jj <= nbin; ++jj )
            sum += cnt[ii][jj];
        if ( sum )
        {
            Property const* p = properties.find_or_die("couple", ii);
            out << LIN << ljust(p->name(), 2);
            for ( size_t jj = 0; jj <= nbin; ++jj )
                out << ' ' << std::setw(5) << cnt[ii][jj];
        }
    }
}


/**
 Export number of Couples in each state
 */
void Simul::reportCouple(std::ostream& out, Property const* sel) const
{
    constexpr size_t SUP = 128;
    int act[SUP] = { 0 }, cnt[SUP][4];
    
    //reset counts:
    for ( size_t i = 0; i < SUP; ++i )
    {
        cnt[i][0] = 0;
        cnt[i][1] = 0;
        cnt[i][2] = 0;
        cnt[i][3] = 0;
    }
    
    for ( Couple const* i=couples.firstFF(); i ; i = i->next() )
    {
        assert_true(!i->attached1() && !i->attached2());
        size_t x = i->prop->number();
        if ( x < SUP )
        {
            if ( i->active() ) ++act[x];
            ++(cnt[x][0]);
        }
    }
    
    for ( Couple const* i=couples.firstAF(); i ; i = i->next() )
    {
        assert_true(i->attached1() && !i->attached2());
        size_t x = i->prop->number();
        if ( x < SUP )
        {
            if ( i->active() ) ++act[x];
            ++(cnt[x][1]);
        }
    }
    for ( Couple const* i=couples.firstFA(); i ; i = i->next() )
    {
        assert_true(!i->attached1() && i->attached2());
        size_t x = i->prop->number();
        if ( x < SUP )
        {
            if ( i->active() ) ++act[x];
            ++(cnt[x][2]);
        }
    }
    
    for ( Couple const* i=couples.firstAA(); i ; i = i->next() )
    {
        assert_true(i->attached1() && i->attached2());
        size_t x = i->prop->number();
        if ( x < SUP )
        {
            if ( i->active() ) ++act[x];
            ++(cnt[x][3]);
        }
    }
    
    out << COM << ljust("couple", 2, 2);
    out << SEP << "total";
    out << SEP << "active";
    out << SEP << "FF";
    out << SEP << "AF";
    out << SEP << "FA";
    out << SEP << "AA";
    
    for ( Property const* i : properties.find_all("couple") )
    {
        if ( sel && i != sel )
            continue;
        out << LIN << ljust(i->name(), 2);
        size_t x = i->number();
        if ( x < SUP )
        {
            out << SEP << cnt[x][0]+cnt[x][1]+cnt[x][2]+cnt[x][3];
            out << SEP << act[x];
            for ( size_t d = 0; d < 4; ++d )
                out << SEP << cnt[x][d];
        }
        else
            out << SEP << "out-of-range";
    }
}


/**
 Export composition of Couple classes
 */
void Simul::reportCoupleAnatomy(std::ostream& out, Property const* sel) const
{
    out << COM << "hand_id" << SEP << rjust("hand_name", 2, 1);
    
    for ( Property const* i : properties.find_all("hand") )
    {
        HandProp const* p = static_cast<HandProp const*>(i);
        out << LIN << p->number();
        out << SEP << rjust(p->name(), 2);
    }

    out << COM << "couple_id" << SEP << rjust("couple_name", 2, 1);
    out << SEP << rjust("hand1", 2) << SEP << rjust("hand2", 2);
    
    for ( Property const* i : properties.find_all("couple") )
    {
        if ( sel && i != sel )
            continue;
        CoupleProp const* p = static_cast<CoupleProp const*>(i);
        out << LIN << p->number();
        out << SEP << rjust(p->name(), 2);
        out << SEP << rjust(p->hand1_prop->name(), 2);
        out << SEP << rjust(p->hand2_prop->name(), 2);
    }
}

//------------------------------------------------------------------------------
#pragma mark - Cluster Analysis

/**
 equalize flag() if connected by a Couple
 */
void Simul::flagClustersCouples() const
{
    for ( Couple const* X=couples.firstAA(); X ; X=X->next() )
    {
        ObjectFlag f = X->fiber1()->flag();
        ObjectFlag g = X->fiber2()->flag();
        if ( f != g )
            changeFlags(std::max(f, g), std::min(f, g));
    }
}

/**
 equalize flag() if connected by a Couple of given type
 */
void Simul::flagClustersCouples(Property const* sel) const
{
    for ( Couple const* X=couples.firstAA(); X ; X=X->next() )
    {
        if ( X->prop == sel )
        {
            ObjectFlag f = X->fiber1()->flag();
            ObjectFlag g = X->fiber2()->flag();
            if ( f != g )
                changeFlags(std::max(f, g), std::min(f, g));
        }
    }
}

/**
 equalize flag() if connected through Single of class Wrist:
 */
void Simul::flagClustersSingles() const
{
    for ( Single * X=singles.firstA(); X; X=X->next() )
    {
        Mecable const* B = X->base();
        if ( B )
        {
            ObjectFlag f = X->fiber()->flag();
            ObjectFlag g = B->flag();
            if ( f != g )
                changeFlags(std::max(f, g), std::min(f, g));
        }
    }
}


void Simul::flagClusters(bool C, bool S, bool M) const
{
    if ( ! ( C | S | M ) )
        throw InvalidSyntax("you must specify a cluster type: couple=1 or solid=1 or meca=1");

    setUniqueFlags();
    if ( C ) flagClustersCouples();
    if ( S ) flagClustersSingles();
    if ( M ) flagClustersMeca();
    flagLargestCluster(1UL);
}


/// class to store info about a Cluster
struct Cluster
{
    size_t     cnt;
    ObjectFlag flg;

    Cluster(ObjectFlag f, size_t n) : cnt(n), flg(f) {}
        
    /// Compare function
    bool operator < (Cluster const&b) const
    {
        if ( cnt == b.cnt ) return flg < b.flg;
        return cnt > b.cnt;
    }
};


/**
Set Mecable::flag() to 'f' for Mecables in the largest cluster
*/
void Simul::flagLargestCluster(ObjectFlag f) const
{
    // for large systems, it would be better to use a std:unordered_map here:
    typedef std::map<ObjectFlag, size_t> map_t;
    map_t map;

    // collect number of fibers with same 'flag' value:
    for ( Fiber* fib = fibers.first(); fib; fib = fib->next() )
        ++map[fib->flag()];
    
    size_t size = 0;
    ObjectFlag largest = 0;
    // find largest clusters in 'map':
    for ( map_t::value_type const& i : map )
    {
        if ( i.second > size )
        {
            largest = i.first;
            size = i.second;
        }
    }
    
    if ( size > 2 )
    {
        // swap 'largest' and 'f':
        for ( Fiber* fib = fibers.first(); fib; fib = fib->next() )
        {
            if ( fib->flag() == f )
                fib->flag(largest);
            else if ( fib->flag() == largest )
                fib->flag(f);
        }
    }
}

/**
Order Clusters from 1 to max, in order of decreasing size,
and set fiber:flag() to the corresponding cluster index.
@return number of clusters
*/
size_t reportOrderedClusters(std::ostream& out, Fiber* first, size_t threshold, size_t details)
{
    typedef std::vector<Fiber*> list_t;
    typedef std::map<ObjectFlag, list_t> map_t;
    // the std::set will keep its elements ordered:
    typedef std::multiset<Cluster> set_t;
    map_t map;
    set_t clusters;

    // extract clusters in 'map' and reset fiber's flag:
    for ( Fiber* F = first; F; F = F->next() )
    {
        //std::clog << fib->reference() << " " << fib->flag() << "\n";
        map[F->flag()].push_back(F);
        F->flag(0);
    }
    
    // insert clusters with size information to get them sorted:
    for ( map_t::value_type const& i : map )
    {
        size_t s = i.second.size();
        if ( s >= threshold )
            clusters.emplace(i.first, s);
    }
    
    if ( details > 1 )
    {
        out << COM << "cluster_id" << SEP << "nb_fibers :" << SEP << "fiber_id";
        out << LIN << clusters.size() << " clusters:";
    }
    else if ( details > 0 )
    {
        out << COM << "cluster_id" << SEP << "nb_fibers";
        out << LIN << clusters.size() << " clusters:";
    }
    
    // consider clusters by decreasing size:
    size_t idx = 0;
    for ( set_t::value_type const& c : clusters )
    {
        ++idx;
        list_t const& list = map[c.flg];
        
        for ( Fiber * i : list )
            i->flag(idx);

        if ( details > 0 )
        {
            out << LIN << c.flg << SEP << c.cnt << ":";
            size_t cnt = std::min(list.size(), details);
            for ( size_t i = 0; i < cnt; ++i )
                out << " " << list[i]->identity();
        }
    }
    return idx;
}


/**
 Export clusters defined by Simul::flagClusters()
 Clusters are ordered in decreasing size.
 */
void Simul::reportClusters(std::ostream& out, Glossary& opt) const
{
    size_t details = 128;
    bool C = false, S = false, M = false;

    opt.set(details, "details");
    opt.set(C, "couples", "couple");
    opt.set(S, "singles", "single");
    opt.set(M, "meca");
    
    flagClusters(C, S, M);
    
    out << COM << "cluster by couples " << C << " solids " << S << " meca " << M;
    reportOrderedClusters(out, fibers.first(), 2, details);
}


//------------------------------------------------------------------------------
#pragma mark - Ring

/**
 Evaluates if the Fiber distribution makes a connected ring around the Z-axis
 @returns number of rings
 FJN 8.07.2017, 8.11.2018, for Blood Platelet project
 */
size_t Simul::flagRing() const
{
    flagClusters(1, 0, 0);

    typedef std::list<ObjectFlag> list_t;
    list_t ring;
    
    // include all cluster into 'ring'
    for ( Fiber * fib = fibers.first(); fib; fib=fib->next() )
    {
        ObjectFlag f = fib->flag();
        if ( f > 0 )
        {
            list_t::const_iterator i = std::find(ring.begin(), ring.end(), f);
            if ( i == ring.end() )
                ring.push_back(f);
        }
    }
    
    // rotate plane around the Z-axis and find intersecting fibers
    for ( unsigned a = 0; a < 360; ++a )
    {
        real ang = a * M_PI / 180.0;
        Vector nor( std::cos(ang), std::sin(ang), 0.0);
        Vector dir(-std::sin(ang), std::cos(ang), 0.0);
        
        list_t sec;
        for ( Fiber const* fib=fibers.first(); fib; fib=fib->next() )
        {
            for ( size_t s = 0; s < fib->nbSegments(); ++s )
            {
                // check that fiber intersect with plane:
                real abs = fib->planarIntersect(s, nor, 0);
                if ( 0 <= abs  &&  abs < 1 )
                {
                    // check that intersection is located on 'dir' side of Z-axis:
                    Vector pos = fib->midPoint(s, abs);
                    if ( dot(pos, dir) > 0 )
                    {
                        // transfer cluster if it was already present before:
                        ObjectFlag f = fib->flag();
                        list_t::iterator i = std::find(ring.begin(), ring.end(), f);
                        if ( i != ring.end() )
                        {
                            ring.erase(i);
                            sec.push_back(f);
                        }
                    }
                }
            }
        }
        ring = sec;
    }
    
    if ( ring.size() == 1 )
    {
        ObjectFlag f = *ring.begin();
        for ( Fiber * fib = fibers.first(); fib; fib=fib->next() )
        {
            if ( fib->flag() == f )
                fib->flag(1);
            else
                fib->flag(0);
        }
    }
    else if ( ring.size() > 0 )
    {
        // unflag all fibers that are not part of a ring:
        for ( Fiber * fib = fibers.first(); fib; fib=fib->next() )
        {
            if ( std::find(ring.begin(), ring.end(), fib->flag()) == ring.end() )
                fib->flag(0);
        }
    }
    else
    {
        // unflag all
        for ( Fiber * fib = fibers.first(); fib; fib=fib->next() )
            fib->flag(0);
    }
    return ring.size();
}


/**
 Calculate the length of the ring and its mean radius
 FJN 15.09.2018, for Blood Platelet project
*/
void Simul::analyzeRing(ObjectFlag flg, real& length, real& radius) const
{
    length = 0.0;
    radius = 0.0;

    Vector cen_old;
    unsigned rad_cnt = 0;
    
    // rotate plane around the Z-axis and find intersecting fibers
    for ( unsigned a = 0; a < 360; ++a )
    {
        real ang = a * M_PI / 180.0;
        Vector nor( std::cos(ang), std::sin(ang), 0);
        Vector dir(-std::sin(ang), std::cos(ang), 0);
        
        unsigned cen_cnt = 0;
        Vector cen(0,0,0);
        for ( Fiber const* fib=fibers.first(); fib; fib=fib->next() )
        {
            // only consider fiber that are part of the ring:
            if ( fib->flag() == flg )
            {
                for ( size_t s = 0; s < fib->nbSegments(); ++s )
                {
                    // check that fiber intersect with plane:
                    real abs = fib->planarIntersect(s, nor, 0);
                    if ( 0 <= abs  &&  abs < 1 )
                    {
                        Vector pos = fib->midPoint(s, abs);
                        // check that intersection is located on 'dir' side of Z-axis:
                        real H = dot(pos, dir);
                        if ( H > 0 )
                        {
                            radius += H;
                            ++rad_cnt;
                            cen += pos;
                            ++cen_cnt;
                        }
                    }
                }
            }
        }
        if ( cen_cnt > 0 )
        {
            cen /= cen_cnt;
            if ( a > 0 )
                length += ( cen - cen_old ).norm();
            cen_old = cen;
        }
    }
    radius /= rad_cnt;
}

/**
Calculate the length of the ring and its radius
*/
void Simul::reportRing(std::ostream& out) const
{
    out << COM << "nb_rings" << SEP << "length" << SEP << "radius";
    size_t nring = flagRing();
    if ( nring == 1 )
    {
        real len, rad;
        analyzeRing(1, len, rad);
        out << LIN << 1 << SEP << std::fixed << len<< SEP << rad ;
    }
    else
    {
        out << LIN << nring << SEP << 0.0 << SEP << 0.0;
    }
}


/**
 Evaluates if the Fiber distribution makes a connected ring around the Z-axis
 FJN 8.07.2017, for Blood Platelet project
 */
void Simul::reportPlatelet(std::ostream& out) const
{
    size_t nfib;
    real pol, var, mn, mx;
    fibers.infoLength(fibers.collect(), nfib, pol, var, mn, mx);
    pol *= nfib;
    
    if ( nfib > 1024 )
    {
        // the calculation can take too much time with lots of fibers, so we cut it here:
        out << COM << "nb_fiber" << SEP << "polymer";
        out << LIN << nfib << SEP << std::fixed << pol;
        return;
    }
    
    computeForces();

    real t = 0, ten = 0, inf = 0, sup = 0;
    size_t c, cnt = 0;

    // rotate plane around the Z-axis and find intersecting fibers
    for ( unsigned a = 0; a < 180; a += 10 )
    {
        real ang = a * M_PI / 180.0;
        Vector dir(std::cos(ang), std::sin(ang), 0);
        fibers.infoTension(c, t, inf, sup, dir, 0);
        cnt += 2;  // every plane should intersect the ring twice
        ten += t;
    }
    ten /= (real)cnt;
    
    std::ofstream nos("/dev/null");
    real force = reportFiberConfinement(nos);

    real len = 0.0, rad = 0.0;
    if ( flagRing() == 1 )
        analyzeRing(1, len, rad);

    std::streamsize p = out.precision();
    out.precision(std::min(p,(std::streamsize)3));
    out << COM << "nb_fiber" << SEP << "polymer" << SEP << "tension" << SEP << "force" << SEP << "length" << SEP << "radius";
    out << LIN << nfib << SEP << std::fixed << pol << SEP << ten << SEP << force << SEP << len << SEP << rad;
    out.precision(p);
}


//------------------------------------------------------------------------------
#pragma mark - Misc

/**
 Export indices calculated by FiberSet::infoSpindle
 */
void Simul::reportIndices(std::ostream& out) const
{
    out << COM << "amount" << SEP << "radial" << SEP << "polar";
    real ixa, ixp;
    fibers.infoSpindle(ixa, ixp, Vector(1,0,0), 0, 10, 0.5);
    out << LIN << fibers.size();
    out << SEP << ixa;
    out << SEP << ixp;
}


/**
 Export number of Fibers pointing left and right,
 that intersect a plane parallel to YZ.
 The planes are distributed regularly every 0.5 um along the X-axis.
 */
void Simul::reportProfile(std::ostream& out) const
{
    out << COM << "position" << SEP << "leftward" << SEP << "rightward" << SEP << "right-left";
    Vector n(1,0,0);
    real m = 10, dm = 0.5;
    int R, L;
    for ( real x = -m ; x <= m ; x += dm )
    {
        fibers.infoPlane(R, L, n, -x);
        out << LIN << x;
        out << SEP << L;
        out << SEP << R;
        out << SEP << R-L;
    }
}


/**
 l'angle entre un vecteur 1 (centre du noyau --> SPB)
 et un vecteur 2 (axe de l'hyphe; gauche --> droite = sens du flow).
 */
void Simul::reportAshbya(std::ostream& out) const
{
    out << COM << "class id point_0, vector_1, angle";
    for ( Solid const* obj=solids.first(); obj; obj=obj->next() )
    {
        out << LIN << obj->prop->number();
        out << SEP << obj->identity();
        out << SEP << obj->posPoint(0);
        if ( obj->nbPoints() > 1 )
        {
            Vector vec = normalize(obj->diffPoints(0));
            Vector dir(1,0,0);
            out << SEP << vec;
            out << SEP << std::acos(dot(vec, dir));
        }
    }
}


/**
 categorize the configuration of two microtubules, with respect to collisions
 This calculated 3 boolean values: K = catastrophe, X = crossing, Z = zippered
 1.10.2021 -- 11.2022
 */
void Simul::reportFiberCollision(std::ostream& out, Property const* sel, Glossary& opt) const
{
#if ( DIM == 1 )
	throw InvalidParameter("fiber:collision meaningless in 1D");
#endif
	if ( fibers.size() > 2 )
		throw InvalidParameter("fiber:collision can only handle two fibers");
    
    bool print = 0;
    opt.set(print, "print");
    static real abs = -77; // abscissa of contact point
    static real ang = 777; // angle at first contact
    static real dis = INFINITY; // minimum distance reached
    static char cat = 'U'; // category

    Fiber const * fib = nullptr, *fox = nullptr;
    for ( Fiber const* f = fibers.first(); f; f = f->next() )
    {
        if ( f->prop == sel )
            fib = f;
        else
            fox = f;
    }
    bool K = 1, X = 0, T = 0, Z = 0;

    if ( fib && fox )
    {
        const real sup = 3 * fib->prop->steric_radius;
        
        // K = state at plus-end
        K = ( fib->endStateP() == 4 );
        
        // check if tip of 'fib' is close to 'fox':
        Vector tip = fib->posEndP();
        real d, aaa = fox->projectPoint(tip, d);
        Vector dir = fox->dir(aaa);
        dis = std::min(dis, d);
        
        // plus-tip of 'fil' is in contact with 'fox':
        bool contact = ( d < sup*sup );
        if ( contact )
        {
            // check direction of fib's tip to fox at closest point:
            real C = dot(fib->dirEndP(), dir);
            real A = std::acos(C);
            // the angle and abscissa are set at first contact:
            if ( ang > 700 )
            {
                ang = A;
                abs = aaa;
            }
            // 'zippering' implies 'being tangent' and 'moving along' the obstacle
            T = ( abs_real(C) > 0.94 );   // tangent within 20 degrees
            // distanced zipped is measured along the obstacle, with abscissa:
            Z = ( abs_real(aaa-abs) > 2 ); // distance zipped sufficient
        }
        else if ( fib->length() > 1 )
        {
            // the plus-tip may have crossed the other filament if it is not in contact
            // consider a point 1um back, and check if it is on opposite side of 'fox'
            Vector bak = fib->posFrom(1, PLUS_END);
            real ddd, bbb = fox->projectPoint(bak, ddd);
            bbb = 0.5 * ( aaa + bbb );
            Vector mid = fox->pos(bbb);
            Vector axs = fox->dir(bbb);
            Torque TP = cross(axs, tip-mid);
            Torque TM = cross(axs, bak-mid);
#if ( DIM == 3 )
            X = ( dot(TP, TM) < 0 );
#else
            X = ( TP * TM < 0 );
#endif
        }

        if ( cat == 'U' )
        {
            if ( K )
            {
                // recorded catastrophes must be at contact
                if ( contact )
                    cat = 'K';
                // but in any case, we can stop the simulation
                abort_time();
            }
            else if ( Z ) cat = 'Z';
        }
        // The 'X' may superseed the Z and U category
        if ( X && !K ) cat = 'X';
                
        // check if MT has reached its max_length:
        if ( fib->length() >= fib->prop->max_length - 0.01 )
        {
            cat = 'U';
            abort_time();
        }

        // since these states are final, we can terminate the simulation
        if ( cat == 'K' || cat == 'X' || cat == 'Z' )
            abort_time();
    }
    else
        abort_time();
    
    if ( print )
    {
        out << LIN << std::fixed << std::setprecision(5) << ang;
        out << SEP << std::fixed << std::setprecision(5) << std::sqrt(dis);
        out << SEP << K << SEP << X << SEP << Z << SEP << T << SEP << cat;
        // reset static variables for next round:
        abs = -77;
        ang = 777;
        dis = INFINITY;
        cat = 'U';
    }
}

/**
 Export end-to-end distance of Fiber
 */
void Simul::reportSomething(std::ostream& out) const
{
    for ( Fiber const* fib = fibers.firstID(); fib; fib = fibers.nextID(fib) )
    {
        Vector ee = fib->posEndP() - fib->posEndM();
        out << ee.norm() << " ";
    }
}
