// Cytosim was created by Francois Nedelec. Copyright Cambridge University 2019-2020

#include <iostream>
#include <numeric>
#include <list>
#include <set>
#include "tokenizer.h"
#include "stream_func.h"

/// width of columns in formatted output, in number of characters
int column_width = 10;

/// use this macro at the beginning of a line of comment
#define COM "\n% " << std::setw(column_width-2)

/// use this macro at the beginning of a new line of data
#define LIN '\n' << std::setw(column_width)

/// use this macro to separate between values on a line of data
#define SEP ' ' << std::setw(column_width-1)

#include "accumulator.h"

/// pad string by adding white-space on the right up to size 'n * column_width - p'
std::string ljust(std::string const& str, size_t n, size_t p = 0)
{
    size_t s = n * (size_t)column_width - p;
    if ( str.size() < s )
        return str + std::string(s-str.size(), ' ');
    else
        return str;
}

/// pad string by adding white-space on the left up to size 'n * column_width - p'
std::string rjust(std::string const& str, size_t n, size_t p = 1)
{
    size_t s = n * (size_t)column_width - p;
    if ( str.size() < s )
        return std::string(s-str.size(), ' ') + str;
    else
        return str;
}

/// repeat string DIM times with X, Y, Z endings as appropriate
std::string repeatXYZ(std::string const& str)
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
void remove_plural(std::string & str)
{
    if ( str.size() > 2  &&  str.at(str.size()-1) == 's' )
        str.resize(str.size()-1);
}



/**
 Surround the report with comments to identify start/end
 */
void Simul::report_wrap(std::ostream& out, std::string const& arg, Glossary& opt) const
{
    out << "% time " << std::fixed << std::setprecision(3) << prop->time;
    out << "\n% report " << arg;
    report(out, arg, opt);
    out << "% end\n";
}


/**
 calls report_one() multiple times, if `arg` contains multiple instructions,
 separated by a comma, for example:
 
     report fiber:force,fiber:length
 
 parameters can be includes as in:
 
     report fiber:cluster {couple=1},fiber:length

 */
void Simul::report(std::ostream& out, std::string what, Glossary& opt) const
{
    std::streamsize p = 4;
    opt.peek(p, "precision");
    opt.peek(column_width, "column") || opt.peek(column_width, "width");
    
    // use fixed notation:
    out.setf(std::ios_base::fixed, std::ios_base::floatfield);
    out.precision(p);

    std::stringstream is(what);
    while ( is.good() )
    {
        std::string arg = Tokenizer::get_polysymbol(is, false);
        std::string blok = Tokenizer::get_block(is, '{');
        try {
            if ( blok.empty() )
            {
                //out << "\nSimul::report(" << arg << ")";
                report_one(out, arg, opt);
            }
            else
            {
                //out << "\nSimul::report(" << arg << ", " << blok << ")";
                Glossary glos(opt);
                glos.read(blok, 0);
                report_one(out, arg, glos);
                opt.add_counts(glos);
            }
        }
        catch( Exception & e )
        {
            out << '\n' << e.brief();
        }
        char c = Tokenizer::get_character(is, true, false);
        if ( c == EOF )
            break;
        if ( c != ',' )
        {
            out << '\n';
            char rest[256] = { 0 };
            rest[0] = c;
            is.getline(rest+1, sizeof(rest)-1);
            throw InvalidParameter("unexpected `" + std::string(rest) + "' in report string");
        }
    }
    //opt.write_counts(std::cerr);
    out << '\n';
}


/**
 Split 'arg' into who:what and call report_one()
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
    
    // allow for approximate spelling (missing 's'):
    remove_plural(who);
    remove_plural(what);

    int com = 1, split = false;
    opt.peek(com, "verbose");
    opt.set(split, "split");
    
    //std::clog << "report("<< what << "|" << who << ")\n";
    if ( isCategory(who) )
    {
        if ( split )
        {
            // process every class in the category separatly:
            PropertyList plist = properties.find_all(who);
            for ( Property const* sel : plist )
            {
                out << COM << "class " + std::to_string(sel->number()) + " is `" + sel->name() << "'";
                report_one(out, who, sel, what, com, opt);
                com = 0;
            }
        }
        else
        {
            // process all objects irrespective of their class:
            report_one(out, who, nullptr, what, com, opt);
        }
    }
    else
    {
        // check if name corresponds to a property:
        Property const* sel = nullptr;
        sel = properties.find(who);
        if ( sel )
            report_one(out, sel->category(), sel, what, com, opt);
        else
            report_one(out, who, nullptr, what, com, opt);
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
 
 
 WHAT                | Output
 --------------------|----------------------------------------------------------
 `system:time`       | Time
 `system:inventory`  | summary list of objects
 `system:property`   | All properties
 `system:parameter`  | global parameters
 `system:NAME`       | parameters for Property 'NAME'


 WHAT                    | Output
 ------------------------|------------------------------------------------------
 `fiber:position`        | Position and orientation of fibers
 `fiber:age`             | Average age of fibers
 `fiber:length`          | Average length and standard deviation of fibers
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
                       std::string const& what, bool com, Glossary& opt) const
{
    if ( who == "fiber" )
    {
        if ( what.empty() || what == "position" )
            return reportFibers(out, sel, com);
        if ( what == "plus_end" )
            return reportFiberEnds(out, PLUS_END, sel, com);
        if ( what == "minus_end" )
            return reportFiberEnds(out, MINUS_END, sel, com);
        if ( what == "end" )
            return reportFiberEnds(out, BOTH_ENDS, sel, com);
        if ( what == "point" )
            return reportFiberPoints(out, sel, com);
        if ( what == "displacement" )
            return reportFiberDisplacement(out, sel, com);
        if ( what == "moment" )
            return reportFiberMoments(out);
        if ( what == "speckle" )
            return reportFiberSpeckles(out, opt);
        if ( what == "sample" )
            return reportFiberSamples(out, opt);
        if ( what == "segment" )
            return reportFiberSegments(out);
        if ( what == "length" )
            return reportFiberLengths(out);
        if ( what == "distribution" || what == "histogram" )
            return reportFiberLengthHistogram(out, opt);
        if ( what == "tension" )
            return reportFiberTension(out, opt);
        if ( what == "energy" )
            return reportFiberBendingEnergy(out);
        if ( what == "dynamic" )
        {
            reportFiberEndState(out, PLUS_END, sel, true);
            reportFiberEndState(out, MINUS_END, sel, false);
            return;
        }
        if ( what == "plus_state" )
            return reportFiberEndState(out, PLUS_END, sel, true);
        if ( what == "minus_state" )
            return reportFiberEndState(out, MINUS_END, sel, true);
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
            return reportFiberLattice(out);
        if ( what == "mesh" )
            return reportFiberMesh(out, false);
        if ( what == "mesh_density" )
            return reportFiberMesh(out, true);
        if ( what == "connector" )
            return reportFiberConnectors(out, opt);

        throw InvalidSyntax("I can only report fiber: position, end, minus_end, plus_end, "\
                            "point, moment, speckle, sample, segment, dynamic, length, "\
                            "distribution, tension, force, cluster, age, energy, hand, link");
    }
    if ( who == "bead" )
    {
        if ( what == "position" || what.empty() )
            return reportBeadPosition(out, sel, com);
        if ( what == "single" )
            return reportBeadSingles(out);
        throw InvalidSyntax("I can only report bead: position, single");
    }
    if ( who == "solid" )
    {
        if ( what == "hand" )
            return reportSolidHands(out);
        if ( what == "position" || what.empty() )
            return reportSolidPosition(out, sel, com);
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
            return reportSpherePosition(out, sel, com);
        throw InvalidSyntax("I can only report sphere:position");
    }
    if ( who == "single" )
    {
        if ( what.empty() )
            return reportSingle(out);
        if ( what == "link" )
            return reportSingleLink(out, sel, com);
        if ( what == "state" )
            return reportSingleState(out, sel, com);
        if ( what == "force" )
            return reportSingleForce(out, sel, com);
        if ( what == "position" )
            return reportSinglePosition(out, sel, com);
        throw InvalidSyntax("I can only report single: link, state, force, position");
    }
    if ( who == "couple" )
    {
        if ( what.empty() )
            return reportCouple(out);
        if ( what == "state" )
            return reportCoupleState(out, sel, com);
        if ( what == "link" )
            return reportCoupleLink(out, sel, com);
        if ( what == "configuration" )
            return reportCoupleConfiguration(out, sel, com, opt);
        if ( what == "force" )
            return reportCoupleForce(out, sel, com);
        if ( what == "histogram" )
            return reportCoupleForceHistogram(out, opt);
        if ( what == "active" )
            return reportCoupleActive(out, sel, com);
        if ( what == "anatomy" )
            return reportCoupleAnatomy(out);
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
    if ( who == "system" )
    {
        if ( what.empty() )
            return reportSystem(out);
        if ( what == "time" )
            return reportTime(out);
        if ( what == "inventory" )
            return reportInventory(out);
        if ( what == "property" || what == "parameter" )
            return writeProperties(out, false);
        throw InvalidSyntax("I can only report system: time, inventory, property");
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
    if ( who == "network" && what == "bridge" )
        return reportNetworkBridges(out, opt);
    if ( who == "ring" )
        return reportRing(out);
    if ( who == "platelet" )
        return reportPlatelet(out);
    if ( who == "ashbya" )
        return reportAshbya(out);
    if ( who == "custom" )
        return reportCustom(out);

    throw InvalidSyntax("unknown report `"+who+"'");
}

//------------------------------------------------------------------------------
#pragma mark - Fiber Aggregated Properties

/**
 Export average length and standard-deviation for each class of fiber
 */
void Simul::reportFiberAge(std::ostream& out) const
{
    out << COM << ljust("class", 2, 2) << SEP << "count" << SEP << "avg_birth";
    out << SEP << "dev_birth" << SEP << "avg_age" << SEP << "min_age" << SEP << "max_age";
    
    size_t cnt;
    real avg, dev, mn, mx;
    const real now = prop->time;

    for ( Property const* i : properties.find_all("fiber") )
    {
        FiberProp const* fp = static_cast<FiberProp const*>(i);
        ObjectList objs = fibers.collect(fp);
        fibers.infoBirthtime(objs, cnt, avg, dev, mn, mx);
        out << LIN << ljust(fp->name(), 2);
        out << SEP << cnt;
        out << SEP << avg;
        out << SEP << dev;
        out << SEP << now-mx;
        out << SEP << now-avg;
        out << SEP << now-mn;
    }
}

/**
 Export average length and standard-deviation for each class of fiber
 */
void Simul::reportFiberLengths(std::ostream& out) const
{
    out << COM << ljust("class", 2, 2) << SEP << "count" << SEP << "avg_len" << SEP << "std_dev";
    out << SEP << "min_len" << SEP << "max_len" << SEP << "total";

    size_t cnt;
    real avg, dev, mn, mx;
    
    std::streamsize p = out.precision();
    for ( Property const* i : properties.find_all("fiber") )
    {
        FiberProp const* fp = static_cast<FiberProp const*>(i);
        
        ObjectList objs = fibers.collect(fp);
        fibers.infoLength(objs, cnt, avg, dev, mn, mx);
        
        out << LIN << ljust(fp->name(), 2);
        out << SEP << cnt;
        out.precision(3);
        out << SEP << std::fixed << avg;
        out << SEP << std::fixed << dev;
        out << SEP << std::fixed << mn;
        out << SEP << std::fixed << mx;
        out.precision(1);
        out << SEP << std::fixed << avg*cnt;
    }
    out.precision(p);
}


/**
 Export average length and standard-deviation for each class of fiber
 */
void Simul::reportFiberLengthHistogram(std::ostream& out, Glossary & opt) const
{
    const size_t BMAX = 256;
    unsigned cnt[BMAX+1];

    real sup = 0;
    for ( Fiber const* fib = fibers.firstID(); fib; fib = fibers.nextID(fib) )
        sup = std::max(sup, fib->length());

    real delta = ( sup > 2 ) ? 1 : 0.1;
    size_t nbin = std::ceil(sup/delta);
    opt.set(delta, "interval");
    opt.set(nbin, "interval", 1);
    nbin = std::min(nbin, BMAX);
    
    if ( 1 )
    {
        out << COM << "fiber length histogram (bin size " << delta <<")";
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
        
        for ( Fiber const* fib = fibers.firstID(); fib; fib = fibers.nextID(fib) )
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
void Simul::reportFiberEndState(std::ostream& out, FiberEnd end, Property const* sel, bool com) const
{
    std::string name;
    if ( sel )
        name = sel->name() + ":";
    name.append(end==PLUS_END ?"plus_end":"minus_end");
    
    if ( com )
    {
        out << COM << ljust("class", 2, 2) << SEP << "total" << SEP << "static";
        out << SEP << "green" << SEP << "yellow" << SEP << "orange" << SEP << "red";
    }
    
    constexpr size_t MAX = 5;
    size_t cnt[MAX+1] = { 0 };
    size_t sum = 0;
    
    for ( Fiber const* fib = fibers.firstID(); fib; fib = fibers.nextID(fib) )
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
                    out << SEP << ha->linkBase();
                    out << SEP << ha->linkStiffness();
                }
            }
        }
    }
}


/**
 Report quantity of substance in the fiber's Lattice
 */
void Simul::reportFiberLattice(std::ostream& out) const
{
    out << COM << ljust("class", 2, 2);
    out << SEP << "total" << SEP << "avg" << SEP << "min" << SEP << "max" << SEP << "length";
    
    size_t cnt = 0;
    real len = 0, sm = 0, mn = INFINITY, mx = -INFINITY;
    
    for ( Fiber const* fib = fibers.firstID(); fib; fib = fibers.nextID(fib) )
        fib->infoLattice(len, cnt, sm, mn, mx);

    out << LIN << ljust("fiber:lattice", 2);
    out << SEP << sm;
    out << SEP << std::setprecision(4) << sm / (real)cnt;
    out << SEP << std::fixed << std::setprecision(6) << mn;
    out << SEP << std::fixed << std::setprecision(6) << mx;
    out << SEP << std::setprecision(3) << len;
}


/**
 Report quantity of substance in the fiber's Lattice
 */
void Simul::reportFiberMesh(std::ostream& out, bool density) const
{
    out << COM << ljust("class", 2, 2);
    out << SEP << "total" << SEP << "avg" << SEP << "min" << SEP << "max" << SEP << "length";
    
    size_t cnt = 0;
    real len = 0, sm = 0, mn = INFINITY, mx = -INFINITY;
    
    for ( Fiber const* fib = fibers.firstID(); fib; fib = fibers.nextID(fib) )
        fib->infoMesh(len, cnt, sm, mn, mx, density);

    out << LIN << ljust("fiber:mesh", 2);
    out << SEP << sm;
    out << SEP << std::setprecision(4) << sm / (real)cnt;
    out << SEP << std::fixed << std::setprecision(6) << mn;
    out << SEP << std::fixed << std::setprecision(6) << mx;
    out << SEP << std::setprecision(3) << len;
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
    out << SEP << dot(fib->dirEndM(), fib->dirEndP());
    out << SEP << organizers.findOrganizerID(fib);
}

    
/// to sort in ascending order: return -1 if ( a < b ) and +1 if ( a > b )
int compareFibers(Object const* A, Object const* B)
{
    real a = static_cast<Fiber const*>(A)->length();
    real b = static_cast<Fiber const*>(B)->length();
    return ( a < b ) - ( a > b );
    //return ( a > b ) - ( b > a );
}


/**
 Export length, position and directions at center of fibers
 */
void Simul::reportFibersSorted(std::ostream& out, Property const* sel, bool com) const
{
    // sort fibers in the ObjectPool:
    const_cast<Simul*>(this)->fibers.pool.blinksort(compareFibers);

    if ( com )
    {
        out << COM << "class" << SEP << "identity" << SEP << "length";
        out << SEP << repeatXYZ("pos") << SEP << repeatXYZ("dir");
        out << SEP << "endToEnd" << SEP << "cosinus" << SEP << "organizer";
    }

    for ( Fiber const* fib = fibers.first(); fib; fib = fib->next() )
    {
        if ( !sel || sel == fib->prop )
            reportFiber(out, fib);
    }
}


/**
 Export length, position and directions at center of fibers
 */
void Simul::reportFibers(std::ostream& out, Property const* sel, bool com) const
{
    if ( com )
    {
        out << COM << "class" << SEP << "identity" << SEP << "length";
        out << SEP << repeatXYZ("pos") << SEP << repeatXYZ("dir");
        out << SEP << "endToEnd" << SEP << "cosinus" << SEP << "organizer";
    }

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
void Simul::reportFiberEnds(std::ostream& out, FiberEnd end, Property const* sel, bool com) const
{
    if ( com )
    {
        out << COM << "class" << SEP << "identity" << SEP << "length";
        if ( end & PLUS_END )
            out << SEP << "stateP" << SEP << repeatXYZ("posP") << SEP << repeatXYZ("dirP");
        if ( end & MINUS_END )
            out << SEP << "stateM" << SEP << repeatXYZ("posM") << SEP << repeatXYZ("dirM");
    }
    
    for ( Fiber const* fib = fibers.firstID(); fib; fib = fibers.nextID(fib) )
    {
        if ( !sel || sel == fib->prop )
        {
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
}


/**
 Export Fiber-number, position of vertices
 */
void Simul::reportFiberPoints(std::ostream& out, Property const* sel, bool com) const
{
    if ( com )
        out << COM << "identity" << SEP << repeatXYZ("pos") << SEP << "curvature";

    // list fibers in the order of the inventory:
    for ( Fiber const* fib = fibers.firstID(); fib; fib = fibers.nextID(fib) )
    {
        if ( !sel || sel == fib->prop )
        {
            out << COM << "fiber " << fib->reference() << "  " << fib->segmentation();
            
            for ( size_t p = 0; p < fib->nbPoints(); ++p )
            {
                out << LIN << fib->identity();
                out << SEP << fib->posP(p);
                out << SEP << fib->curvature(p);
            }
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
            uint32_t z = fib->signature();
            real a = spread * std::log(z*TINY);
            while ( a > fib->abscissaP() )
            {
                z = lcrng2(z);
                a += spread * std::log(z*TINY);
            }
            while ( a >= fib->abscissaM() )
            {
                out << '\n' << fib->pos(a);
                z = lcrng2(z);
                a += spread * std::log(z*TINY);
            }
        }
        // generate speckles above the origin of abscissa
        if ( fib->abscissaP() > 0 )
        {
            uint32_t z = ~fib->signature();
            real a = -spread * std::log(z*TINY);
            while ( a < fib->abscissaM() )
            {
                z = lcrng1(z);
                a -= spread * std::log(z*TINY);
            }
            while ( a <= fib->abscissaP() )
            {
                out << '\n' << fib->pos(a);
                z = lcrng1(z);
                a -= spread * std::log(z*TINY);
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
 Export Mean Squared displacement of fibers since the last call
 to this function.
 */
void Simul::reportFiberDisplacement(std::ostream& out, Property const* sel, bool com) const
{
    typedef std::map <ObjectID, Vector> fiber_map;
    static fiber_map positions;
    static double past = 0;
    
    if ( com )
        out << COM << "delta_time nb_fibers mean_squared_displacement";
    
    real sum = 0;
    size_t cnt = 0;
    for ( Fiber const* fib = fibers.firstID(); fib; fib = fibers.nextID(fib) )
    {
        if ( !sel || sel == fib->prop )
        {
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
    }
    
    if ( cnt > 0 )
        out << LIN << time() - past << SEP << cnt << SEP << sum / (real)cnt;
    else
        out << LIN << time() - past << SEP << 0 << SEP << 0;
    
    past = time();
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
    
    Accumulator accum;
    
    for ( Property const* i : properties.find_all("fiber") )
    {
        FiberProp const* fp = static_cast<FiberProp const*>(i);
        
        accum.reset();
        
        for ( Fiber const* fib = fibers.firstID(); fib; fib = fibers.nextID(fib) )
        {
            if ( fib->prop == fp )
            {
                const real w = fib->segmentation();
                accum.add(0.5*w, fib->posEndM());
                for ( size_t n = 1; n < fib->lastPoint(); ++n )
                    accum.add(w, fib->posP(n));
                accum.add(0.5*w, fib->posEndP());
            }
        }
        
        accum.subtract_mean();
        out << LIN << ljust(fp->name(), 2);
        accum.print(out, 0);
    }
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
            out << SEP << fib->posP(p);
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


/**
 Export fiber elastic bending energy
 */
void Simul::reportFiberBendingEnergy(std::ostream& out) const
{
    out << COM << ljust("bending_energy",2,2) << SEP << "count";
    out << SEP << "sum" << SEP << "avg" << SEP << "dev" << SEP << "rigidity";
    
    size_t cnt;
    real avg, dev;
    
    for ( Property const* i : properties.find_all("fiber") )
    {
        FiberProp const* fp = static_cast<FiberProp const*>(i);
        ObjectList objs = fibers.collect(fp);
        fibers.infoBendingEnergy(objs, cnt, avg, dev);
        if ( cnt > 0 )
        {
            out << LIN << ljust(fp->name(), 2);
            out << SEP << cnt;
            out << SEP << std::setprecision(3) << avg*cnt;
            out << SEP << std::setprecision(3) << avg;
            out << SEP << std::setprecision(3) << dev;
            out << SEP << std::setprecision(3) << fp->rigidity;
        }
    }
}


/**
 Export total magnitude of force exerted by Fiber on the confinement
 */
void Simul::reportFiberConfineForce(std::ostream& out) const
{
    out << COM << "confinement forces";
    out << COM << "identity" << SEP << repeatXYZ("pos") << SEP << repeatXYZ("force");
     
     // list fibers in the order of the inventory:
     for ( Fiber const* fib = fibers.firstID(); fib; fib = fibers.nextID(fib) )
     {
         out << COM << "fiber " << fib->reference();
         Space const* spc = findSpace(fib->prop->confine_space);
         const real stiff = fib->prop->confine_stiffness;

         for ( size_t p = 0; p < fib->nbPoints(); ++p )
         {
             out << LIN << fib->identity();
             Vector w, pos = fib->posP(p);
             out << SEP << pos;
             if ( spc->outside(pos) )
             {
                 w = spc->project(pos);
                 out << SEP << stiff * ( w - pos );
             }
             else
             {
                 out << SEP << Vector(0,0,0);
             }
         }
     }
}


/**
 Export total magnitude of force exerted by Fiber on the confinement
 */
real Simul::reportFiberConfinement(std::ostream& out) const
{
    out << COM << "count" << SEP << repeatXYZ("force") << SEP << "radial";
    size_t cnt = 0;
    Vector sum(0,0,0);
    real   rad = 0;
    
#if ( DIM > 1 )
    for ( Fiber const* fib = fibers.firstID(); fib; fib = fibers.nextID(fib) )
    {
        Space const* spc = findSpace(fib->prop->confine_space);
        const real stiff = fib->prop->confine_stiffness;
        
        for ( size_t p = 0; p < fib->nbPoints(); ++p )
        {
            Vector w, pos = fib->posP(p);
            if ( spc->outside(pos) )
            {
                ++cnt;
                w = spc->project(pos);
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
    Accumulator accum;
    
    for ( Fiber const* fib = fibers.firstID(); fib; fib = fibers.nextID(fib) )
    {
        unsigned cnt = 0;
        for ( Fiber const* fob = fibers.nextID(fib); fob; fob = fibers.nextID(fob) )
        {
            for ( size_t ii = 0; ii < fib->nbSegments(); ++ii )
            {
                FiberSegment seg1(fib, ii);
                for ( size_t jj = 0; jj < fob->nbSegments(); ++jj )
                {
                    real abs1, abs2;
                    FiberSegment seg2(fob, jj);
                    if ( seg1.shortestDistance(seg2, abs1, abs2) <= sup )
                    {
                        if ( seg1.within(abs1) & seg2.within(abs2) )
                        {
                            ++cnt;
                            Vector pos1 = seg1.pos(abs1/seg1.len());
                            //Vector pos2 = loc2.pos(abs2/loc2.len());
                            if ( details == 2 )
                            {
                                out << LIN << fib->identity();
                                out << SEP << abs1 + seg1.abscissa1();
                                out << SEP << fob->identity();
                                out << SEP << abs2 + seg2.abscissa1();
                                out << SEP << pos1;
                            }
                            accum.add(pos1);
                        }
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
    accum.subtract_mean();
    accum.print_doc(out);
    accum.print(out, 1);
}


/// accessory class to analyse the connections in a network of fibers
struct Connector
{
    real a;
    long f;
    long g;
    long h;
    real s;
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
    Accumulator accum;
    
    typedef std::vector<Connector> clist_t;
    typedef std::map<ObjectID, clist_t> map_t;
    
    map_t map;

    for ( Fiber const* fib = fibers.firstID(); fib; fib = fibers.nextID(fib) )
    {
        map.clear();
        // check all connecting Hands and record abscissa, depending on the fiber that is linked
        for ( Hand const* ha = fib->firstHand(); ha; ha = ha->next() )
        {
            Hand const* oh = ha->otherHand();
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
                accum.add(fib->pos(p->a));
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
    accum.subtract_mean();
    accum.print_doc(out);
    accum.print(out, 1);
}


#include "motor_prop.h"

real hand_speed(HandProp const* hp)
{
    if ( hp->activity == "move" )
        return static_cast<MotorProp const*>(hp)->unloaded_speed;
    return 0;
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
    
    const real speedh1 = hand_speed(hp1);
    const real speedh2 = hand_speed(hp2);

    for ( Fiber const* fib = fibers.firstID(); fib; fib = fibers.nextID(fib) )
    {
        map.clear();
        // check all connecting Hands and record abscissa, depending on the fiber that is linked
        for ( Hand * ha = fib->firstHand(); ha; ha = ha->next() )
        {
            Hand * oh = ha->otherHand();
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
    out << LIN << prop->time;
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
void reportSystemSet(std::ostream& out, SET& set, PropertyList const& properties)
{
    for ( Property const* i : properties.find_all(set.title()) )
    {
        size_t points = 0, sup = 0;
        ObjectList objs = set.collect(match_property, i);
        for ( Object * o : objs )
        {
            Mecable * mec = Simul::toMecable(o);
            if ( mec )
            {
                points += mec->nbPoints();
                sup = std::max(sup, mec->nbPoints());
            }
        }
        if ( points > 0 )
        {
            out << LIN << ljust(i->name(), 2);
            out << SEP << objs.size();
            out << SEP << points << SEP << sup;
        }
    }
}

void Simul::reportSystem(std::ostream& out) const
{
    out << COM << ljust("class", 2, 2) << SEP << "count" << SEP << "vertices" << SEP << "largest";
    reportSystemSet(out,  fibers, properties);
    reportSystemSet(out,  solids, properties);
    reportSystemSet(out, spheres, properties);
    reportSystemSet(out,   beads, properties);
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
        if ( obj->tag() == Aster::TAG )
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
void Simul::reportBeadPosition(std::ostream& out, Property const* sel, bool com) const
{
    if ( com )
        out << COM << "class" << SEP << "identity" << SEP << repeatXYZ("pos");
    
    for ( Bead const* obj=beads.first(); obj; obj=obj->next() )
    {
        if ( !sel || sel == obj->prop )
        {
            out << LIN << obj->prop->number();
            out << SEP << obj->identity();
            out << SEP << obj->position();
        }
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
void Simul::reportSolidPosition(std::ostream& out, Property const* sel, bool com) const
{
    if ( com )
    {
        out << COM << "class" << SEP << "identity" << SEP << repeatXYZ("cen");
        out << SEP << repeatXYZ("point0") << SEP << repeatXYZ("point1");
    }
        
    for ( Solid const* obj=solids.first(); obj; obj=obj->next() )
    {
        if ( !sel || sel == obj->prop )
        {
            out << LIN << obj->prop->number();
            out << SEP << obj->identity();
            out << SEP << obj->centroid();
            out << SEP << obj->posP(0);
            if ( obj->nbPoints() > 1 )
                out << SEP << obj->posP(1);

            if ( modulo )
            {
                Vector pos = obj->centroid();
                modulo->fold(pos);
                out << SEP << pos;
            }
        }
    }
}

/**
 Export position of Solids with counts of Hands and attached Hands
 */
void Simul::reportSolidHands(std::ostream& out) const
{
    out << COM << "class" << SEP << "identity" << SEP << repeatXYZ("pos");
    out << SEP << "nb_hand" << SEP << "nb_link";
    
    for ( Solid const* obj = solids.firstID(); obj; obj = solids.nextID(obj) )
    {
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
void Simul::reportSpherePosition(std::ostream& out, Property const* sel, bool com) const
{
    if ( com )
    {
        out << COM << "class" << SEP << "identity";
        out << SEP << repeatXYZ("point0") << SEP << repeatXYZ("point1");
    }
        
    for ( Sphere const* obj=spheres.first(); obj; obj=obj->next() )
    {
        if ( !sel || sel == obj->prop )
        {
            out << LIN << obj->prop->number();
            out << SEP << obj->identity();
            out << SEP << obj->posP(0);
            if ( obj->nbPoints() > 1 )
                out << SEP << obj->posP(1);
        }
    }
}


/**
 Report something about Space (incomplete)
 */
void Simul::reportSpace(std::ostream& out) const
{
    out << COM << "class" << SEP << "identity";
    
    for ( Space const* obj=spaces.first(); obj; obj=obj->next() )
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
    
    // report total substance in each Field
    for ( Field const* obj=fields.first(); obj; obj=obj->next() )
    {
        real vol = obj->cellVolume();
        Field::value_type s, n, x;
        obj->infoValues(s, n, x);
        out << LIN << ljust(obj->prop->name(), 2);
        out << SEP << s;
        out << SEP << s / ( vol * obj->nbCells() );
        out << SEP << n / vol;
        out << SEP << x / vol;
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
void Simul::reportSingleState(std::ostream& out, Property const* sel, bool com) const
{
    if ( com )
    {
        out << COM << "class" << SEP << "identity";
        out << SEP << repeatXYZ("pos") << SEP << repeatXYZ("force");
        out << SEP << "fiber" << SEP << "abscissa" << SEP << "aster";
    }
    
    for ( Single const* obj = singles.firstID(); obj; obj = singles.nextID(obj) )
        if ( !sel || sel == obj->prop )
            write(out, obj, this);
}


/**
 Export details of attached Singles
 */
void Simul::reportSinglePosition(std::ostream& out, Property const* sel, bool com) const
{
    if ( com )
    {
        out << COM << "class" << SEP << "identity" << SEP << repeatXYZ("pos");
        out << SEP << "fiber" << SEP << "abscissa";
    }
        
    for ( Single const* obj = singles.firstID(); obj; obj = singles.nextID(obj) )
    {
        if ( !sel || sel == obj->prop )
        {
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
}

/**
 Export details of attached Singles
 */
void Simul::reportSingleLink(std::ostream& out, Property const* sel, bool com) const
{
    if ( com )
    {
        out << COM << "class" << SEP << "identity";
        out << SEP << repeatXYZ("pos") << SEP << repeatXYZ("force");
        out << SEP << "fiber" << SEP << "abscissa" << SEP << "aster";
    }
        
    for ( Single const* obj = singles.firstID(); obj; obj = singles.nextID(obj) )
    {
        if ( obj->attached()  && ( !sel || sel == obj->prop ))
            write(out, obj, this);
    }
}


/**
 Export number of Single in each state
 */
void Simul::reportSingle(std::ostream& out) const
{
    constexpr size_t SUP = 128;
    
    size_t free[SUP+1] = { 0 }, bound[SUP+1] = { 0 }, based[SUP+1] = { 0 };
    
    for ( Single const* i = singles.firstF(); i ; i = i->next() )
    {
        assert_true(!i->attached());
        ++free[std::min(i->prop->number(), SUP)];
        based[std::min(i->prop->number(), SUP)] += ( i->base() != nullptr );
    }
    
    for ( Single const* i=singles.firstA(); i ; i=i->next() )
    {
        assert_true(i->attached());
        ++bound[std::min(i->prop->number(), SUP)];
        based[std::min(i->prop->number(), SUP)] += ( i->base() != nullptr );
    }
    
    if ( 1 )
    {
        out << COM << ljust("single", 2, 2);
        out << SEP << "total";
        out << SEP << "free";
        out << SEP << "bound";
        out << SEP << "based";
    }
    
    for ( Property const* i : properties.find_all("single") )
    {
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
void Simul::reportSingleForce(std::ostream& out, Property const* sel, bool com) const
{
    constexpr size_t MAX = 8;
    real cnt[MAX+1] = { 0 };
    real avg[MAX+1] = { 0 };
    real sup[MAX+1] = { 0 };
    real len[MAX+1] = { 0 };

    // accumulate counts:
    for ( Single const* i=singles.firstA(); i; i=i->next() )
    {
        if ( i->hasForce() && ( !sel || sel == i->prop ))
        {
            size_t x = std::min(MAX, i->prop->number());
            real f = i->force().norm();
            avg[x] += f;
            cnt[x] += 1;
            sup[x] = std::max(sup[x], f);
            len[x] = std::max(len[x], i->stretch().norm());
        }
    }
    
    if ( com )
        out << COM << ljust("class", 2, 2) << SEP << "avg_force" << SEP << "max_force" << SEP << "max_len";
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
void Simul::reportCoupleState(std::ostream& out, Property const* sel, bool com) const
{
    if ( com )
    {
        out << COM << "class" << SEP << "identity" << SEP << "active" << SEP << repeatXYZ("pos");
        out << SEP << "fiber1" << SEP << "abscissa1" << SEP << "fiber2" << SEP << "abscissa2";
    }
    
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
void Simul::reportCoupleActive(std::ostream& out, Property const* sel, bool com) const
{
    if ( com )
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
void Simul::reportCoupleLink(std::ostream& out, Property const* sel, bool com) const
{
    if ( com )
    {
        out << COM << "class" << SEP << "identity";
        out << SEP << "fiber1" << SEP << "abscissa1";// << SEP << repeatXYZ("pos1");
        out << SEP << "fiber2" << SEP << "abscissa2";// << SEP << repeatXYZ("pos2");
        out << SEP << "force" << SEP << "cos_angle";
    }
        
    for ( Couple const* obj=couples.firstAA(); obj ; obj=obj->next() )
    {
        if ( !sel || sel == obj->prop )
        {
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
}


/**
 Export configuration of bridging couple, as
 P: parallel
 A: antiparallel
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
                                      bool com, Glossary& opt) const
{
    real threshold = 0.010;
    FiberEnd end = PLUS_END;
    
    opt.set(threshold, "threshold");
    opt.set(end, "end", {{"plus_end", PLUS_END}, {"minus_end", MINUS_END}});
    
    size_t T[6] = { 0 };
    for ( Couple const* obj=couples.firstAA(); obj ; obj=obj->next() )
    {
        if ( !sel || sel == obj->prop )
            ++T[obj->configuration(end, threshold)];
    }
    size_t sum = T[0]+T[1]+T[2]+T[3]+T[4]+T[5];
    
    if ( com )
        out << COM << "couples" << SEP << "P" << SEP << "A" << SEP << "X" << SEP << "T" << SEP << "V";
    out << LIN << sum << SEP << T[0] << SEP << T[1] << SEP << T[2] << SEP << T[3] << SEP << T[4];
 }



/**
 Export average properties of Couples forces
 */
void Simul::reportCoupleForce(std::ostream& out, Property const* sel, bool com) const
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
            size_t x = std::min(MAX, i->prop->number());
            real f = i->force().norm();
            avg[x] += f;
            cnt[x] += 1;
            sup[x] = std::max(sup[x], f);
            len[x] = std::max(len[x], i->stretch().norm());
        }
    }
        
    if ( com )
        out << COM << ljust("class", 2) << SEP << "avg_force" << SEP << "max_force" << SEP << "max_len";
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
        out << COM << "force_distribution" << " (`scale` indicates the center of each bin)";
        out << LIN << ljust("scale", 2);
        for ( size_t u = 0; u <= nbin; ++u )
            out << " " << std::setw(5) << delta * ( u + 0.5 );
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
void Simul::reportCouple(std::ostream& out) const
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
    
    if ( 1 )
    {
        out << COM << ljust("couple", 2, 2);
        out << SEP << "total";
        out << SEP << "active";
        out << SEP << "FF";
        out << SEP << "AF";
        out << SEP << "FA";
        out << SEP << "AA";
    }
    
    for ( Property const* i : properties.find_all("couple") )
    {
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
void Simul::reportCoupleAnatomy(std::ostream& out) const
{
    out << COM << "hand_id" << SEP << rjust("hand_name", 2, 1);
    
    for ( Property const* i : properties.find_all("hand") )
    {
        HandProp const* p = static_cast<HandProp const*>(i);
        out << LIN << p->number();
        out << SEP << rjust(p->name(), 2);
    }
    out << '\n';

    out << COM << "class_id" << SEP << rjust("couple_name", 2, 1);
    out << SEP << rjust("hand1", 2) << SEP << rjust("hand2", 2);

    for ( Property const* i : properties.find_all("couple") )
    {
        CoupleProp const* p = static_cast<CoupleProp const*>(i);
        out << LIN << p->number();
        out << SEP << rjust(p->name(), 2);
        out << SEP << rjust(p->hand1_prop->name(), 2);
        out << SEP << rjust(p->hand2_prop->name(), 2);
    }
}

//------------------------------------------------------------------------------
#pragma mark - Clusters


// set flag() to unique value for all fibers
void resetFlags(FiberSet const& set)
{
    for ( Fiber * fib = set.first(); fib; fib=fib->next() )
        fib->matchFlagIdentity();
}


// set flag() for all fibers to `f`
void resetFlags(FiberSet const& set, ObjectFlag f)
{
    for ( Fiber * fib = set.first(); fib; fib=fib->next() )
        fib->flag(f);
}


/**
 Substitute the values of fiber:flag() such that both `a` and `b`
 values are replaced by min(a, b).
 */
void reFlag(FiberSet const& set, ObjectFlag a, ObjectFlag b)
{
    // swap to ensure b < a
    if ( a < b )
        std::swap(a, b);

    // replace a -> b = min(a, b)
    for ( Fiber* fib = set.first(); fib; fib=fib->next() )
    {
        if ( fib->flag() == a )
            fib->flag(b);
    }
}


/**
 equalize flag() when fibers are connected by a Couple:
 */
void Simul::flagClustersCouples() const
{
    for ( Couple const* cx = couples.firstAA(); cx ; cx=cx->next() )
    {
        ObjectFlag f1 = cx->fiber1()->flag();
        ObjectFlag f2 = cx->fiber2()->flag();
        if ( f1 != f2 )
            reFlag(fibers, f1, f2);
    }
}

/**
 equalize flag() when fibers are connected by a Couple of given type:
 */
void Simul::flagClustersCouples(Property const* sel) const
{
    for ( Couple const* cx = couples.firstAA(); cx ; cx=cx->next() )
    {
        if ( cx->prop == sel )
        {
            ObjectFlag f1 = cx->fiber1()->flag();
            ObjectFlag f2 = cx->fiber2()->flag();
            if ( f1 != f2 )
                reFlag(fibers, f1, f2);
        }
    }
}


/**
 equalize flag() when fibers are connected through blobs:
 */
void Simul::flagClustersSolids() const
{
    for ( Solid const* obj=solids.first(); obj; obj=obj->next() )
    {
        SingleList anchored = singles.collectWrists(obj);
        // find the minimun flag value:
        ObjectFlag flg = 0;
        for ( Single const* s : anchored )
        {
            if ( s->attached() )
            {
                ObjectFlag f = s->fiber()->flag();
                if ( flg == 0 || f < flg )
                    flg = f;
            }
        }
        // reflag:
        if ( flg > 0 )
        {
            for ( Single const* s : anchored )
            {
                if ( s->attached() )
                    reFlag(fibers, s->fiber()->flag(), flg);
            }
        }
    }
}


void Simul::flagClusters(bool cop, bool sol, bool mec) const
{
    if ( ! ( cop | sol | mec ) )
        throw InvalidSyntax("missing cluster type [couple|solid|meca]");

    resetFlags(fibers);
    if ( mec ) flagMecaClusters();
    if ( cop ) flagClustersCouples();
    if ( sol ) flagClustersSolids();
    flagLargestCluster(1UL);
}


/// class to store info about a Cluster
struct Cluster
{
    ObjectFlag flg;
    size_t     cnt;
    
    Cluster(ObjectFlag f, size_t n) { flg = f; cnt = n; }
        
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
    typedef std::map<ObjectFlag, size_t> map_t;
    map_t map;

    // extract clusters in 'map' and reset fiber's flag:
    for ( Fiber* fib = fibers.first(); fib; fib = fib->next() )
        ++map[fib->flag()];
    
    size_t size = 0;
    ObjectFlag largest = 0;
    // insert clusters with size information to get them sorted:
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
int Simul::orderClusters(std::ostream& out, size_t threshold, int details) const
{
    typedef std::vector<Fiber*> list_t;
    typedef std::map<ObjectFlag, list_t> map_t;
    // the std::set will keep its elements ordered:
    typedef std::multiset<Cluster> set_t;
    map_t map;
    set_t clusters;

    // extract clusters in 'map' and reset fiber's flag:
    for ( Fiber* fib = fibers.firstID(); fib; fib = fibers.nextID(fib) )
    {
        //std::clog << fib->reference() << " " << fib->flag() << "\n";
        map[fib->flag()].push_back(fib);
        fib->flag(0);
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
    } else if ( details > 0 )
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
            out << LIN << c.flg << SEP << c.cnt;
            
            if ( details > 1 )
            {
                out << " :";
                for ( Fiber * i : list )
                    out << " " << i->identity();
            }
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
    int details = 2;
    bool C = false, S = false, M = false;

    opt.set(details, "details");
    opt.set(C, "couples") || opt.set(C, "couple");
    opt.set(S, "solids") || opt.set(S, "solid");
    opt.set(M, "meca");
    
    flagClusters(C, S, M);
    
    out << COM << "cluster by couples " << C << " solids " << S << " meca " << M;
    orderClusters(out, 2, details);
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
                    Vector pos = fib->posPoint(s, abs);
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
    for ( unsigned a = 0; a <= 360; ++a )
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
                        Vector pos = fib->posPoint(s, abs);
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
    real pol, dev, mn, mx;
    fibers.infoLength(fibers.collect(), nfib, pol, dev, mn, mx);
    pol *= nfib;
    
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

    out << COM << "nb_fiber" << SEP << "polymer" << SEP << "tension" << SEP << "force" << SEP << "length" << SEP << "radius";
    out << LIN << nfib << SEP << std::fixed << std::setprecision(3) << pol << SEP << ten << SEP << force << SEP << len << SEP << rad;
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
    fibers.infoSpindle(ixa, ixp, Vector(1,0,0), 0, 30, 1);
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
    out << COM << "position" << SEP << "left-pointing" << SEP << "right-pointing";
    Vector n(1,0,0);
    real m = 40, dm = 0.5;
    int nr, nl;
    for ( real p = -m ; p <= m ; p += dm )
    {
        fibers.infoPlane(nr, nl, n, -p);
        out << LIN << p;
        out << SEP << nl;
        out << SEP << nr;
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
        out << SEP << obj->posP(0);
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
 Export end-to-end distance of Fiber
 */
void Simul::reportCustom(std::ostream& out) const
{
    for ( Fiber const* fib = fibers.firstID(); fib; fib = fibers.nextID(fib) )
    {
        Vector ee = fib->posEndP() - fib->posEndM();
        out << ee.norm() << " ";
    }
}
