// Cytosim was created by Francois Nedelec. Copyright 2022 Cambridge University
#include "dim.h"
#include "cymdef.h"
#include <fstream>
#include <unistd.h>
#include "messages.h"
#include "filepath.h"
#include "time_date.h"
#include "parser.h"
#include "print_color.h"

// Use second definition to trace execution
#define VLOG(ARG) ((void) 0)
//#define VLOG(ARG) std::clog << ARG;

/**
 \var Simul::currentFormatID
 An integer `currentFormatID` is used to record the format of the trajectory file,
 allowing some backward compatibility as the format of the trajectory file evolves.
 It is stored in the file and accessible upon reading through `Inputter::formatID()`
 
 History of changes in file format:

 - 60: 2/12/2023 NULL_TAG changed from 'v' to '!'
 - 59: 11/01/2023 Chain::birth_time moved to Fiber
 - 58: 26/11/2022 References written on 2 or 4 bytes
 - ??: 11/06/2021 References written always on 4 bytes with tag+identity
 - 56: 23/06/2021 Secondary TAG use capital letters, but ID was not changed
 - 56: 19/01/2021 All fiber dynamic states stored on 16 bytes, really
 - 55: 02/10/2020 Interpolation4 stores coefficients only if mecable!=nullptr
 - 54: 25/05/2020 All fiber dynamic stored on 16 bytes
 - 53: 18/10/2019 Aster, Nucleus and Bundle store their primary Mecable directly
 - 52: 18/10/2019 Space's shape is stored always on 16 characters
 - 51: 03/03/2019 Storing number of Aster links
 - 50: 19/12/2018 Using ASCII's 8th bit for fat references. Fiber's birth time in Filament (now Chain)
 - 49: 12/12/2018 FiberSite writes the Lattice index but not the abscissa
 - 49: 22/11/2018 reference do not include mark, which is writen in object header
 - 48: 04/07/2018 Fiber stores its birth time
 - 47: 13/02/2017 Wrist and Aster built on a local reference frame of Solid
 - 46: 23/10/2015 GrowingFiber writes dynamic states, changed ClassicFiber:write()
 - 45: 18/09/2015 Indices of Properties are numbers starting from one, and not zero
 - 44: 25/07/2015 the dynamic state of fiber ends is stored with a separate TAG
 - 43: 24/04/2014 number of dimension in Space stored in 16 bits
 - 42: 09/11/2013 All fibers store end_state on 32 bits
 -     08/12/2012 FRAME_TAG was changed from "#frame " to "#Cytosim "
 - 41: 17/11/2012 Space stores its shape as a string in objects.cmo
 - 40: 11/09/2012 Aster format simplified
 - 39: 11/07/2012 Object::mark is stored on 32 bits instead of 16 previously
 - 38: 03/05/2012 Fiber stores length instead of segment-length
 - 37: 30/04/2012 Couple::Duo stores its activity state
 - 36: 22/02/2012 All Spaces store their dimensions in objects.cmo
 - 35: 15/09/2011 Some Spaces store their dimensions in objects.cmo
 - 34: 20/12/2010 Moved Fiber::mark to Object::mark
 - 33: 29/04/2010 Added Fiber::mark
 - 32: 15/04/2010 Space became an Object
 - 31: 01/04/2010 Fiber became a Mecable
 - 30: The Tag were reduced to 1 char: saves space & simplifies code
 - 28, 29: 26/05/2009 started cytosim-PI: a revolution!
 - 27: 22/03/2008 new Fiber::write(), called in Tubule::write()
 - 26: 03/11/2007 Hand do not record haEnd flag
 - 24: 14/12/2006 started cytosim 3, lots of changes
 - 23: 10/12/2005 new class Solid
 - 22: modified Sphere
 - 21: modified Sphere
 - 20: 12/07/2004
 - 19: introduced different kinds of Single
*/


//------------------------------------------------------------------------------
#pragma mark - Write Objects

/**
 This writes all objects of the current state to a trajectory file
*/
void Simul::writeObjects(Outputter& out) const
{
    // write a line identifying a new frame:
    out.write("\n#Cytosim  "+std::to_string(getpid())+"  "+TimeDate::date_string());
    
    // record file format:
    out.write("\n#format "+std::to_string(currentFormatID)+" dim "+std::to_string(DIM));
    
    // identify the file as binary, with its endianess:
    if ( out.binary() )
    {
        out.write("\n#binary ");
        out.writeEndianess();
    }
    
    out.write("\n#time "+std::to_string(prop.time)+" sec");
    
    // save first line of text:
    if ( text_.size() )
    {
        size_t n = text_.find('\n');
        out.write("\n#text "+text_.substr(0, n));
    }

    /*
     An object should be written after any other objects that it refers to.
     For example, Aster is written after Fiber, Couple after Fiber...
     This makes it easier to reconstruct the state during input.
     */
    spaces.writeSet(out);
    fields.writeSet(out);
    fibers.writeSet(out);
    solids.writeSet(out);
    beads.writeSet(out);
    spheres.writeSet(out);
    singles.writeSet(out);
    couples.writeSet(out);
    organizers.writeSet(out);
    //events.write(out);

    out.write("\n#section end");
    out.write("\n#end cytosim\n");
}


/**
 This writes the current state to a trajectory file called `name`.
 If this file does not exist, it is created de novo.
 If `append == true` the state is added to the file, otherwise it is cleared.
 If `binary == true` a binary format is used, otherwise a text-format is used.
*/
void Simul::writeObjects(std::string const& name, bool append, int binary) const
{
    if ( prop.clear_system_file && name==prop.system_file )
    {
        VLOG("delete ---> " << name <<"\n");
        std::remove(prop.system_file.c_str());
        prop.clear_system_file = false;
    }
    
    Outputter out(name.c_str(), append, binary);
    
    if ( ! out.good() )
        throw InvalidIO("could not open output file `"+name+"' for writing");
    
    VLOG(name << (append?"":"<--- clear") << " <--- frame @ " << time() << "\n");
    
    try
    {
        out.lock();
        writeObjects(out);
        out.unlock();
    }
    catch( InvalidIO & e )
    {
        out.unlock();
        print_blue(stderr, e.brief());
        std::cerr << ", writing trajectory at t=" << time() << '\n';
    }
}

//------------------------------------------------------------------------------
#pragma mark - Accessory functions to read file

/// InportLock is a helper class used to import a cytosim state from a file
class Simul::ImportLock
{
private:
    
    /// pointer
    Simul * sim;
    
public:
    
    /// mark objects to later be able to tell if they have been updated
    ImportLock(Simul * s)
    : sim(s)
    {
        //Cytosim::log("Simul::ImportLock created with ", sim->nbObjects(), " objects\n");
        sim->couples.freeze();
        sim->singles.freeze();
        sim->fibers.freeze();
        sim->beads.freeze();
        sim->solids.freeze();
        sim->spheres.freeze();
        sim->organizers.freeze();
        sim->fields.freeze();
        sim->spaces.freeze();
        //sim->events.freeze();
    }
    
    /// erase objects which have not been upated
    void prune_all()
    {
        //sim->events.defrost();
        sim->organizers.defrost();
        sim->couples.defrostStore();
        sim->singles.defrostStore();
        sim->beads.defrostMore();
        sim->solids.defrostMore();
        sim->spheres.defrost();
        sim->fibers.defrost();
        sim->spaces.defrost();
        sim->fields.defrost();
        sim = nullptr;
    }
    
    /// move back all objects to normal lists, even if they have not been updated
    void thaw_all()
    {
        /*
         Attention: The order of the thaw() below is important:
         destroying a Fiber will detach any motor attached to it,
         and thus automatically move them to the 'unattached' list,
         as if they had been updated from reading the file.
         Destroying couples and singles before the fibers avoids this problem.
         */
        //sim->events.thaw();
        sim->couples.thaw();
        sim->singles.thaw();
        sim->organizers.thaw();
        sim->beads.thaw();
        sim->solids.thaw();
        sim->spheres.thaw();
        sim->fibers.thaw();
        sim->spaces.thaw();
        sim->fields.thaw();
        sim = nullptr;
    }

    /// reset flags
    ~ImportLock()
    {
        if ( sim )
            thaw_all();
        //Cytosim::log("Simul::ImportLock deleted with ", sim->nbObjects(), " objects\n");
    }
};


static bool isAlpha(int i)
{
    return ( 'a' <= i && i <= 'z' ) || ( 'A' <= i && i <= 'Z' );
}


/** Compatibility function for formats < 50 */
[[ maybe_unused ]]
static ObjectID readObjectID_old(Inputter& in, ObjectTag& tag)
{
    int c;
    do
        c = in.get_char();
    while ( c == ' ' );
    
    if ( c == EOF )
        throw InvalidIO("unexpected end of file");

    tag = c & Object::LOW_BITS;
    // detect fat reference:
    int fat = ( c & Object::HIGH_BIT );
#if BACKWARD_COMPATIBILITY < 50
    // up to format 49, a '$' was added to indicate fat format
    if ( c == '$' )
    {
        c = in.get_char();
        if ( c == EOF )
            throw InvalidIO("unexpected end of file");
        tag = c;
        fat = 1;
    }
#endif
    
    if ( tag == 'v' ) // NULL_TAG before 2.12.2023
        return 0;
    
    ObjectID id = 0;

#if BACKWARD_COMPATIBILITY < 49
    if ( in.binary() && in.formatID() < 49 )
    {
        if ( fat )
        {
            in.readUInt16();         // skip property index
            id = in.readUInt32();
            in.readUInt32();         // skip ObjectMark
        }
        else
        {
            in.readUInt8();          // skip property index
            id = in.readUInt16();
        }
    }
    else
#endif
    if ( in.binary() )
    {
        if ( fat )
            id = in.readUInt32bin();
        else
            id = in.readUInt16bin();
    }
    else
    {
        FILE * file = in.file();
#if BACKWARD_COMPATIBILITY < 49
        // skip property index
        if ( in.formatID() < 49 )
        {
            unsigned u;
            if ( 1 != fscanf(file, "%u", &u) )
                throw InvalidIO("readReference (compatibility) failed");
            if ( in.get_char() != ':' )
                throw InvalidSyntax("missing ':'");
        }
#endif
        if ( 1 != fscanf(file, "%u", &id) )
            throw InvalidIO("readReference failed");
#if BACKWARD_COMPATIBILITY < 49
        if ( in.formatID() < 49 )
        {
            // skip ObjectMark which is not used
            int h = in.get_char();
            if ( h == ':' )
            {
                unsigned long u;
                if ( 1 != fscanf(file, "%lu", &u) )
                throw InvalidIO("readReference (compatibility) failed");
            }
            else
            in.unget_char(h);
        }
#endif
    }
    return id;
}


static ObjectID readObjectID(Inputter& in, ObjectTag& tag)
{
    ObjectID id = 0;

    if ( in.binary() )
    {
#if BACKWARD_COMPATIBILITY < 58
        if ( in.formatID() < 58 ) // 57 before 26.11.2022
        {
            char c = in.get_char();
            tag = c & Object::LOW_BITS;
            if ( c & Object::HIGH_BIT )
                id = in.readUInt32bin();
            else if ( c != 'v' ) // NULL_TAG before 2.12.2023
                id = in.readUInt16bin();
            return id;
        }
#endif
        // binary format 58 (26.11.2022)
        int c = in.get_char();
        tag = c & Object::LOW_BITS;
        //assert_true(isalpha(tag));
        if ( c & Object::HIGH_BIT )
        {
            ObjectID x = in.get_char();
            ObjectID y = in.get_char();
            ObjectID z = in.get_char();
            id = ( x << 16 ) + ( y << 8 ) + z;
        }
        else if ( c != Object::NULL_TAG )
        {
#if BACKWARD_COMPATIBILITY < 60
            if ( in.formatID() < 60 && c == 'v') // NULL_TAG before 2.12.2023
                return 0;
#endif
            ObjectID a = in.get_char();
            ObjectID b = in.get_char();
            id = ( b << 8 ) + a;
        }
    }
    else
    {
        do
            tag = in.get_char();
        while ( tag == ' ' );
#if BACKWARD_COMPATIBILITY < 60
        if ( in.formatID() < 60 && tag == 'v') // NULL_TAG before 2.12.2023
            return 0;
#endif
        if ( tag != Object::NULL_TAG )
            id = in.readUInt();
    }
    return id;
}


/**
 Read a fiber (compatible with format tryout 11.06.2021)
 */
Fiber * Simul::readFiberReference(Inputter& in, ObjectTag& tag, ObjectID& id)
{
#if BACKWARD_COMPATIBILITY < 50
    if ( in.formatID() < 50 )
        id = readObjectID_old(in, tag);
    else
#endif
        id = readObjectID(in, tag);
    
    Fiber* fib = nullptr;
    switch( tag )
    {
        case Object::NULL_TAG:
            break;
#if BACKWARD_COMPATIBILITY < 60
        case 'v': // NULL_TAG before 2.12.2023
            break;
#endif
#if BACKWARD_COMPATIBILITY < 57
        case 'l': // LATTICE_TAG was 'l' before 23/06/2021
            tag = Fiber::LATTICE_TAG;
#endif
        case Fiber::TAG:
        case Fiber::COMPACT_TAG:
        case Fiber::LATTICE_TAG:
            fib = fibers.identifyObject(id);
            if ( !fib )
                std::clog << "unknown fiber ID " << id << " (" << tag << ")\n";
             //throw InvalidIO("unknown fiber ID "+std::to_string(id));
            break;
        default:
            if ( in.eof() )
                std::clog << "unexpected end-of-file\n";
            else
                std::clog << "unknown fiber tag (" << tag << ")\n";
            break;
    }
    return fib;
}


/**
 Read an object reference (compatible with format tryout 11.06.2021)
 */
Object * Simul::readReference(Inputter& in, ObjectTag& tag)
{
    ObjectID id;
#if BACKWARD_COMPATIBILITY < 50
    if ( in.formatID() < 50 )
        id = readObjectID_old(in, tag);
    else
#endif
        id = readObjectID(in, tag);
    
    if ( tag == Object::NULL_TAG || id == 0 )
        return nullptr;

    const ObjectSet * set = findSetT(tag);
    
    if ( !set )
    {
        if ( !isalpha(tag) )
            throw InvalidIO("`"+std::string(1,tag)+"' is not a valid class TAG");
        throw InvalidIO("`"+std::string(1,tag)+"' is not a known class TAG");
    }
    
    Object * res = set->identifyObject(id);
    
    if ( !res )
    {
        //throw InvalidIO("unknown object `"+((char)tag+std::to_string(id))+"' referenced");
        std::clog << "unknown object `"+((char)tag+std::to_string(id))+"' referenced\n";
    }
    
    return res;
}

//------------------------------------------------------------------------------
#pragma mark - Read Objects

/**
 This will update the current state to make it identical to what has been saved
 in the file.
 
 Before reading, all objects are marked with flag().
 Every object found in the file is unflagged as it is updated.
 
 When the read is complete, the objects that are still marked are deleted.
 In this way the new state reflects exactly the system that was stored on file.
 
 @returns
 - 0 = success
 - 1 = EOF
 - 2 and above: ERROR code
 .
 */
int Simul::reloadObjects(Inputter& in, bool prune, ObjectSet* subset)
{
    in.lock();
    ImportLock lock(this);
    try
    {
        VLOG("readObjects starts at `" << in.peek() << "'\n");
        int err = readObjects(in, subset);
        VLOG("readObjects done\n");
        in.unlock();

        primed_ = 0;
        // if no error occurred, process objects that have not been updated
        if ( 0 == err && prune )
            lock.prune_all();
        else
            lock.thaw_all();
        // renew pointers to objects, particularly 'confine_spec'
        prop.complete(*this);
        
        if ( !subset || subset == &fibers )
            fibers.updateFibers();

        assert_false(singles.bad());
        assert_false(couples.bad());
        return err;
    }
    catch(Exception & e)
    {
        VLOG("readObjects failed " << e.brief() << "\n");
        in.unlock();
        throw;
    }
}


int Simul::loadObjects(char const* filename)
{
    Inputter in(DIM, filename, true);

    if ( ! in.good() )
        throw InvalidIO("Could not open specified file for reading");
    if ( in.eof() )
        return 1;

    return reloadObjects(in, 0, nullptr);
}


//------------------------------------------------------------------------------

/**
 This returns:
     - 1 if the start marker of a frame is found
     - 2 if the end marker of a frame is found
     - 0 otherwise
*/
int Simul::readMetaData(Inputter& in, std::string& section, ObjectSet*& objset, ObjectSet* subset)
{
    std::string line = in.get_line();
    VLOG("      #|" << line << "|" << '\n');
    std::istringstream iss(line);
    std::string tok;
    iss >> tok;

    // section heading
    if ( tok == "section" )
    {
        std::string sec;
        iss >> sec >> tok;
        if ( sec != section )
        {
            section = sec;
            VLOG("section <---|" << sec << "|\n");
            if ( sec == "end" )
                return 0;
            objset = findSet(sec);
            if ( !objset )
                throw InvalidIO("unknown section |"+sec+"|");
        }
        if ( subset && objset != subset )
            in.skip_until("#section ");
        else if ( tok == "reheat" )
        {
            bool skip = ( objset != &couples && objset != &singles );
            if ( section == "couple" && prop.skip_free_couple )
                skip = true;
            if ( section == "single" && prop.skip_free_single )
                skip = true;
            if ( skip )
                in.skip_until("#section ");
            else
            {
                PropertyID i = 0;
                size_t cnt[16] = { 0 };
                while ( i < 16 && iss >> cnt[i] )
                    ++i;
                if ( i > 0 )
                {
                    if ( section == "couple" )
                    {
                        couples.reheat(cnt, i);
                        couples.makeCouples(cnt, i);
                    } else if ( section == "single" ) {
                        singles.reheat(cnt, i);
                        singles.makeSingles(cnt, i);
                    }
                }
#if BACKWARD_COMPATIBILITY < 60 // until 2.04.2023
                else if ( in.formatID() < 60 )
                {
                    if ( section == "couple" )
                        couples.reheat();
                    else if ( section == "single" )
                        singles.reheat();
                }
#endif
            }
        }
        else if ( section == "single" && tok == "F" )
        {
            if ( prop.skip_free_single )
                in.skip_until("#section ");
#if BACKWARD_COMPATIBILITY < 58 // until 11.11.2022
            if ( in.formatID() < 58 )
            {
                int mod = 0;
                iss >> mod;
                if ( iss.good() && mod == 1 )
                    singles.reheat();
            }
#endif
            return 0;
        }
        else if ( section == "couple" && tok == "FF" )
        {
            if ( prop.skip_free_couple )
                in.skip_until("#section ");
#if BACKWARD_COMPATIBILITY < 58 // until 11.11.2022
            if ( in.formatID() < 58 )
            {
                int mod = 0;
                iss >> mod;
                if ( iss.good() && mod == 1 )
                    couples.reheat();
            }
#endif
        }
        return 0;
    }
    // optional indication giving the number of objects in the section
    else if ( tok == "record" )
    {
        if ( objset )
        {
            size_t cnt = 0, sup_id = 0;
            if ( iss >> cnt >> sup_id )
                objset->reserve(sup_id);
        }
        return 0;
    }
    // frame start
    else if ( tok == "Cytosim" || tok == "cytosim" || tok == "frame" )
        return 1;
    // binary format signature
    else if ( tok == "binary" )
    {
        in.setEndianess(line.substr(7).c_str());
        return 1;
    }
    // text
    else if ( tok == "text" )
    {
        text_ = line.substr(5);
        return 1;
    }
    // info line "#format 48 dim 2"
    else if ( tok == "format" )
    {
        unsigned f = 0, d = 0;
        iss >> f >> tok >> d;
        if ( f != in.formatID() )
        {
            in.setFormatID(f);
            VLOG("Cytosim is reading format " << f << '\n');
        }
        if ( tok == "dim" )
        {
            if ( in.vectorSize() != d )
                Cytosim::warn("mismatch between file (", d, "D) and executable (", DIM, "D)\n");
            in.vectorSize(d);
        }
        return 1;
    }
    // time data "#time 1.2345"
    else if ( tok == "time" )
    {
        iss >> prop.time;
#if BACKWARD_COMPATIBILITY < 48
        // old format info line "#time 14.000000, dim 2, format 47"
        if ( iss.get() == ',' )
        {
            unsigned i = 0;
            iss >> tok >> i;
            if ( tok == "dim" )
            {
                in.vectorSize(i);
                if ( i != DIM )
                    Cytosim::warn("mismatch between file (", i, "D) and executable (", DIM, "D)\n");
            }
            if ( iss.get() == ',' )
            {
                iss >> tok >> i;
                if ( tok == "format" )
                {
                    if ( i != in.formatID() )
                    {
                        in.setFormatID(i);
                        if ( i < BACKWARD_COMPATIBILITY )
                            std::clog << "Cytosim is attempting to read old format " << i << "\n";
                    }
                }
            }
        }
#endif
        return 1;
    }
    //detect the mark at the end of the frame
    else if ( tok == "end" )
    {
        iss >> tok;
        if ( tok == "cytosim" )
            return 2;
#if BACKWARD_COMPATIBILITY < 50
        if ( tok == "frame" )
            return 2;
#endif
    }
    return 0;
}


/**
 Read file, updating existing objects, and creating new ones for those not 
 already present in the Simul.
 If 'subset!=0' only objects from this class will be imported.
 The Inputter should be locked in a multithreaded application
 
 @returns
 - 0 : success
 - 1 : EOF
 - 2 and above: ERROR code

  */
int Simul::readObjects(Inputter& in, ObjectSet* subset)
{
    ObjectSet * objset = nullptr;
    std::string section;
    size_t nb_objects = 0;
    int has_frame = 0;
    int tag = 0;
    int fat = 0;

    while ( 1 )
    {
        do {
            tag = in.get_char();
            if ( tag == '#' )
            {
                int h = readMetaData(in, section, objset, subset);
                if ( h )
                {
                    if ( h == 2 && has_frame )
                        return 0;
                    if ( h == 1 && nb_objects > 0 )
                    {
                        // found another frame start without finishing current one?
                        return 2;
                    }
                    has_frame = h;
                    fresh_ = 1;
                }
                continue;
            }
            if ( tag == '%' )
            {
                if ( objset )
                    objset->readObjectTypes(in);
                else
                    throw InvalidIO("unset objset with `%`");
                continue;
            }
            if ( tag == '\n' )
                continue;
            if ( tag == EOF )
                return 1;
#if BACKWARD_COMPATIBILITY < 50
            // detect fat header, formatID() < 50
            if ( tag == '$' )
            {
                fat = 1 + in.binary();
                tag = in.get_char();
            }
#endif
            // extract ASCII from character:
            fat = ( tag & Object::HIGH_BIT ) + in.binary();
            tag = ( tag & Object::LOW_BITS );
        } while ( !isAlpha(tag) );
        
        //VLOG("READ '" << (char)tag << "' " << (fat?"fat\n":"\n"));

        try
        {
#if BACKWARD_COMPATIBILITY < 50
            // this is an 'older' code pathway, before 2017?
            if ( !objset )
            {
                ObjectSet * set = findSetT(tag);
                if ( set )
                    set->loadObject(in, tag, fat);
            }
            else
#endif
            {
                objset->loadObject(in, tag, fat);
            }
            ++nb_objects;
        }
        catch( Exception & e )
        {
            print_blue(stderr, e.brief());
            if ( objset )
            {
                in.skip_until("#section ");
                if ( in.eof() )
                {
                    fprintf(stderr, " (eof)\n");
                    return 3;
                }
                std::cerr << e.info() << " (skipping section "+section+")\n";
            }
            else
                std::cerr << e.info() << " (section "+section+")\n";
        }
    }
    return 7;
}


//------------------------------------------------------------------------------
#pragma mark - Write/Read Properties


/**
 The order of the output is important, since properties may depend
 on each other (eg. SingleProp and CoupleProp use HandProp).
 Luckily, there is no circular dependency in Cytosim at the moment.
 
 Thus we simply follow the order in which properties were defined,
 and which is the order in which properties appear in the PropertyList.
 */

void Simul::writeProperties(std::ostream& os, const bool prune) const
{
    //std::clog << "Writing properties" << '\n';
    os << "% Cytosim property file, pid " << getpid() << '\n';

    prop.write(os, prune);
    properties.write(os, prune);
}


/**
At the first call, this will write all properties to file,
and save a copy of what was written to a string `properties_saved`.

The next time this is called, the properties will be compared to the string,
and the file will be rewritten only if there is a difference.
*/
void Simul::writeProperties(bool prune) const
{
    std::ostringstream oss;
    writeProperties(oss, prune);
    if ( oss.str() != properties_saved )
    {
        properties_saved = oss.str();
        // use default property file name
        std::string name = prop.property_file.c_str();
        std::ofstream ofs(name);
        //this should be equivalent to: writeProperties(os, prune);
        ofs << "% " << TimeDate::date_string() << '\n';
        ofs << properties_saved << '\n';
        VLOG(name << " <--- properties @ " << time() << "\n");
    }
}


void Simul::loadProperties(const char file[], int verbose)
{
    std::ifstream is(file, std::ifstream::in);
    if ( !is.good() )
        throw InvalidIO("could not find or read `"+std::string(file)+"'");
    std::streampos ipos(0);
    Parser(this, 1, 1, 0, 0, 0, verbose).evaluate(is, ipos);
}


void Simul::loadProperties()
{
    const char * def = prop.property_file.c_str();
#if BACKWARD_COMPATIBILITY < 57
    const char * old = "properties.cmo";
    if ( ! FilePath::is_file(def) && FilePath::is_file(old) )
        loadProperties(old, 0);
    else
#endif
    loadProperties(def, 0);
}
