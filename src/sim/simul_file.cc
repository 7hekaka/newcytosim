// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "dim.h"
#include "cymdef.h"
#include <fstream>
#include <unistd.h>
#include "messages.h"
#include "parser.h"
#include "print_color.h"


/**
 \var Simul::currentFormatID
 An integer `currentFormatID` is used to record the format of the trajectory file,
 allowing some backward compatibility as the format of the trajectory file evolves.
 It is stored in the file and accessible upon reading through `Inputter::formatID()`
 
 History of changes in file format:

 - ??: 11/06/2021 References written always on 4 bytes with tag+identity
 - 56: 23/06/2021 Secondary TAG use capital letters, but ID was not changed
 - 56: 19/01/2021 All fiber dynamic stored on 16 bytes, really
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
    fprintf(out, "\n\n#Cytosim  %i  %s", getpid(), TicToc::date());
    
    // record file format:
    fprintf(out, "\n#format %i dim %i", currentFormatID, DIM);
    
    // identify the file as binary, with its endianess:
    if ( out.binary() )
    {
        fprintf(out, "\n#binary ");
        out.writeEndianess();
    }
    
    /*
     An object should be written after any other objects that it refers to.
     For example, Aster is written after Fiber, Couple after Fiber...
     This makes it easier to reconstruct the state during input.
     */
    fprintf(out, "\n#time %.6f sec", prop->time);

    spaces.write(out);
    fields.write(out);
    fibers.write(out);
    solids.write(out);
    beads.write(out);
    spheres.write(out);
    singles.write(out);
    couples.write(out);
    organizers.write(out);
    tubules.write(out);
    //events.write(out);
    
    fprintf(out, "\n#section end\n#end cytosim\n");
}


/**
 This writes the current state to a trajectory file called `name`.
 If this file does not exist, it is created de novo.
 If `append == true` the state is added to the file, otherwise it is cleared.
 If `binary == true` a binary format is used, otherwise a text-format is used.
*/
void Simul::writeObjects(std::string const& name, bool append, bool binary) const
{
    Outputter out(name.c_str(), append, binary);
    
    if ( ! out.good() )
        throw InvalidIO("could not open output file `"+name+"' for writing");
    
    try
    {
        out.lock();
        writeObjects(out);
        out.unlock();
    }
    catch( InvalidIO & e )
    {
        print_blue(std::cerr, e.brief());
        std::cerr << ", writing trajectory file" << '\n';
    }
}

//------------------------------------------------------------------------------
#pragma mark - Read from file

/** Compatibility function for formats < 50 */
static ObjectID readOldObjectID(Inputter& in, ObjectTag& tag)
{
    int c;
    do
        c = in.get_char();
    while ( c == ' ' );
    
    if ( c == EOF )
        throw InvalidIO("unexpected end of file");

    tag = c & LOW_BITS;
    // detect fat reference:
    int fat = ( c & HIGH_BIT );
#if BACKWARD_COMPATIBILITY < 50
    // up to format 49, a '$' was added to indicate fat format
    if ( c == '$' )
    {
        tag = in.get_char();
        if ( tag == EOF )
            throw InvalidIO("unexpected end of file");
        fat = 1;
    }
#endif
    
    // Object::TAG is the 'void' reference
    if ( tag == Object::TAG )
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
            in.unget(h);
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
        if ( in.formatID() < 58 )
        {
            char c = in.get_char();
            tag = c & LOW_BITS;
            // Object::TAG is the 'void' reference
            if ( tag == Object::TAG )
                return 0;
            if ( c & HIGH_BIT )
                id = in.readUInt32bin();
            else
                id = in.readUInt16bin();
        }
        else
#endif
        {
            uint32_t u = in.readUInt32bin();
            tag = ( u >> 24 ) & LOW_BITS;
            id = u & 0xFFFFFF;
        }
    }
    else
    {
        do
            tag = in.get_char();
        while ( tag == ' ' );
        if ( tag != Object::TAG )
            id = in.readUInt();
    }
    return id;
}


/**
 Read a fiber (new format 11.06.2021)
 */
Fiber * Simul::readFiberReference(Inputter& in, ObjectTag& tag)
{
    ObjectID id;
#if BACKWARD_COMPATIBILITY < 50
    if ( in.formatID() < 50 )
        id = readOldObjectID(in, tag);
    else
#endif
        id = readObjectID(in, tag);
    
    // Object::TAG is the 'void' reference
    if ( tag == Object::TAG )
        return nullptr;

#if BACKWARD_COMPATIBILITY < 57
    if ( tag == 'l' ) // TAG_LATTICE was 'l' before 23/06/2021
        tag = Fiber::TAG_LATTICE;
#endif
    
    if ( tag != Fiber::TAG && tag != Fiber::TAG_ALT && tag != Fiber::TAG_LATTICE )
        throw InvalidIO("expected reference to a fiber ("+std::string(1,tag)+")");

    return fibers.findID(id);
}


/**
 Read a fiber (new format 11.06.2021)
 */
Object * Simul::readReference(Inputter& in, ObjectTag& tag)
{
    ObjectID id;
#if BACKWARD_COMPATIBILITY < 50
    if ( in.formatID() < 50 )
        id = readOldObjectID(in, tag);
    else
#endif
        id = readObjectID(in, tag);
    
    // Object::TAG is the 'void' reference
    if ( tag == Object::TAG || id == 0 )
        return nullptr;

    const ObjectSet * set = findSetT(tag);
    
    if ( !set )
    {
        if ( !isalpha(tag) )
            throw InvalidIO("`"+std::string(1,tag)+"' is not a valid class TAG");
        throw InvalidIO("`"+std::string(1,tag)+"' is not a known class TAG");
    }
    
    Object * res = set->findID(id);
    
    if ( !res )
        throw InvalidIO("unknown object `"+((char)tag+std::to_string(id))+"' referenced");
    
    return res;
}


/// InputLock is a helper class used to import a cytosim state from a file
class Simul::InputLock
{
private:
    
    /// pointer
    Simul * sim;

    /// state
    bool  frozen;
    
public:
    
    /// flag all objects with FLAG
    InputLock(Simul * s)
    : sim(s)
    {
        //Cytosim::log("Simul::InputLock created with %i objects\n", sim->nbObjects());
        sim->couples.freeze();
        sim->singles.freeze();
        sim->fibers.freeze();
        sim->beads.freeze();
        sim->solids.freeze();
        sim->spheres.freeze();
        sim->tubules.freeze();
        sim->organizers.freeze();
        sim->fields.freeze();
        sim->spaces.freeze();
        //sim->events.freeze();
        frozen = true;
    }
    
    /// erase objects flagged with FLAG
    void prune()
    {
        //sim->events.prune();
        sim->organizers.prune();
        sim->tubules.prune();
        sim->couples.prune();
        sim->singles.prune();
        sim->beads.prune();
        sim->solids.prune();
        sim->spheres.prune();
        sim->fibers.prune();
        sim->spaces.prune();
        sim->fields.prune();
        frozen = false;
    }
    
    void thaw()
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
        sim->tubules.thaw();
        sim->beads.thaw();
        sim->solids.thaw();
        sim->spheres.thaw();
        sim->fibers.thaw();
        sim->spaces.thaw();
        sim->fields.thaw();
        frozen = false;
    }

    /// reset flags
    ~InputLock()
    {
        if ( frozen )
            thaw();
        //Cytosim::log("Simul::InputLock deleted with %i objects\n", sim->nbObjects());
    }
};


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
 .
 */
int Simul::reloadObjects(Inputter& in, bool prune, ObjectSet* subset)
{
    in.lock();
    InputLock lock(this);
    try
    {
        int res = readObjects(in, subset);
        in.unlock();
        if ( 0 == res )
        {
            // if no error occurred, process objects that have not been updated
            if ( prune )
                lock.prune();
            else
                lock.thaw();
            // renew pointers to objects, particularly 'confine_space'
            prop->complete(*this);
        }
        return res;
    }
    catch(Exception & e)
    {
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
 Read file, updating existing objects, and creating new ones for those not 
 already present in the Simul.
 If 'subset!=0' only objects from this class will be imported.
 The Inputter should be locked in a multithreaded application
 
 @returns
 - 0 : success
 - 1 : EOF
 - 2 : the file does not appear to be a valid cytosim archive
 
  */
int Simul::readObjects(Inputter& in, ObjectSet* subset)
{
    ObjectSet * objset = nullptr;
    std::string section, line;
    int has_frame = 0;
    int tag = 0, c = 0;
    int fat = 0;

    while ( in.good() )
    {
        do {
            c = in.get_char();
            if ( c == '#' )
            {
                line = in.get_line();
                break;
            }
            tag = ( c & LOW_BITS );
            fat = ( c & HIGH_BIT );
#if BACKWARD_COMPATIBILITY < 50
            // detect fat header, formatID() < 50
            if ( c == '$' )
            {
                fat = 1;
                tag = in.get_char();
            }
#endif
            if ( c == EOF )
                return 1;
        } while ( !isalpha(tag) );
        
        
        //check for meta-data, contained in lines starting with '#'
        if ( c == '#' )
        {
            //std::clog << "      |#" << line << "|" << '\n';
            std::istringstream iss(line);
            std::string tok;
            iss >> tok;

            // section heading
            if ( tok == "section" )
            {
                iss >> section;
                //std::clog << " section |" << section << "|\n";
                if ( section == "end" )
                    continue;
                else if ( section == "single" )
                {
                    int mod = 0;
                    iss >> tok >> mod;
                    // skip unattached Singles
                    if ( tok == "F" )
                    {
                        if ( prop->skip_free_single > 1 )
                            in.skip_until("#section ");
                        else
                            singles.prune_mode = mod;
                    }
                }
                else if ( section == "couple" )
                {
                    int mod = 0;
                    iss >> tok >> mod;
                    // skip unattached Couples
                    if ( tok == "FF" )
                    {
                        if ( prop->skip_free_couple > 1 )
                            in.skip_until("#section ");
                        else
                            couples.prune_mode = mod;
                    }
                }
                objset = findSet(section);
                if ( !objset )
                    std::clog << " warning: unknown section |" << section << "|\n";
                if ( subset && objset != subset )
                    in.skip_until("#section ");
            }
            // frame start
            else if ( tok == "Cytosim" || tok == "cytosim" || tok == "frame" )
            {
                if ( has_frame )
                    return 2;
                has_frame = 1;
            }
            //binary signature
            else if ( tok == "binary" )
            {
                in.setEndianess(line.substr(7).c_str());
            }
            // info line "#format 48 dim 2"
            else if ( tok == "format" )
            {
                unsigned f = 0, d = 0;
                iss >> f >> tok >> d;
                in.formatID(f);
                //if ( f != currentFormatID )
                //    Cytosim::warn << "Cytosim is reading format ID " << f << '\n';
                if ( tok == "dim" )
                {
                    if ( in.vectorSize() != d )
                        Cytosim::warn << "mismatch between file ("<<d<<"D) and executable ("<<DIM<<"D)\n";
                    in.vectorSize(d);
                }
            }
            // time data "#time 1.2345"
            else if ( tok == "time" )
            {
                iss >> prop->time;
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
                            Cytosim::warn << "mismatch between file ("<<i<<"D) and executable ("<<DIM<<"D)\n";
                    }
                    if ( iss.get() == ',' )
                    {
                        iss >> tok >> i;
                        if ( tok == "format" )
                        {
                            in.formatID(i);
                            //if ( i != currentFormatID )
                            //    std::clog << "Cytosim is reading format ID " << i << "\n";
                        }
                    }
                }
#endif
            }
            //detect the mark at the end of the frame
            else if ( tok == "end" )
            {
                iss >> tok;
                if ( tok == "cytosim" )
                    return 0;
#if BACKWARD_COMPATIBILITY < 50
                if ( tok == "frame" )
                    return 0;
#endif
            }
        }
        else
        {
            //std::clog << "OBJECT |" << (char)tag << "| " << (fat?"fat\n":"\n");
            assert_true( isalpha(tag) );
#if BACKWARD_COMPATIBILITY < 50
            // this is an 'older' code pathway, before 2017?
            if ( !objset )
            {
                ObjectSet * set = findSetT(tag);
                if ( set )
                    set->loadObject(in, tag, fat, true);
                continue;
            }
#endif
            try
            {
                // check that we are using the correct ObjectSet:
                assert_true( objset == findSetT(tag) );
                objset->loadObject(in, tag, fat, true);
            }
            catch( Exception & e )
            {
                print_blue(std::cerr, e.brief());
                std::cerr << e.info() << " ("+section+")\n";
                if ( objset )
                    in.skip_until("#section ");
            }
        }
    }
    return 2;
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
    os << "% " << TicToc::date() << '\n';

    prop->write(os, prune);
    properties.write(os, prune);
}


/**
 At the first call, this will write all properties to file, 
 and save a copy of what was written to a string `properties_saved`.
 
 The next time this is called, the properties will be compared to the string,
 and the file will be rewritten only if there is a difference.
 */
void Simul::writeProperties(char const* name, bool prune) const
{
    std::ostringstream oss;
    writeProperties(oss, prune);
    if ( oss.str() != properties_saved )
    {
        properties_saved = oss.str();

        // use default file name if 'name' is empty or not provided
        if ( !name || *name==0 )
            name = prop->property_file.c_str();
        
        std::ofstream os(name);
        //this should be equivalent to: writeProperties(os, prune);
        os << properties_saved << std::endl;
        //std::clog << "Writing properties at frame " << currentFrame() << '\n';
    }
}


void Simul::loadProperties()
{
    Parser(*this, 1, 1, 0, 0, 0).readConfig(prop->property_file);
}
