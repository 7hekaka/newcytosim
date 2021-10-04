// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University

#include "cymdef.h"
#include "parser.h"
#include "messages.h"
#include "tokenizer.h"
#include "glossary.h"
#include "filepath.h"
#include "stream_func.h"
#include "simul_prop.h"
#include "simul.h"
#include "print_color.h"
#include <fstream>


// Use the second definition to get some reports:
#define VLOG(ARG) ((void) 0)
//#define VLOG(ARG) std::clog << ARG;


//------------------------------------------------------------------------------
/**
 The permission of the parser are set by a number of variables:
 - do_set: can create new Properties
 - do_change: can modify existing Property or Object
 - do_new: can create new Object
 - do_run: can perform simulation steps
 - do_write: can write to disc
 .
 */
Parser::Parser(Simul& sim, bool s, bool c, bool n, bool r, bool w)
: Interface(sim), do_set(s), do_change(c), do_new(n), do_run(r), do_write(w)
{
}

/// check for unused values in Glossary and issue a warning
void check_warnings(Glossary& opt, std::istream& is, std::streampos ipos, size_t cnt = 1)
{
    std::string war;
    if ( opt.has_warning(war, cnt) )
    {
        size_t L;
        Cytosim::log << war << " in " << StreamFunc::get_line(is, ipos, L) << " (line " << L << ")\n";
        // also report to standard error:
        print_yellow(std::cerr, war);
        std::cerr << '\n';
        StreamFunc::print_lines(std::cerr, is, ipos, is.tellg());
    }
}

//------------------------------------------------------------------------------
#pragma mark - Parse

/**
 Create a new Property, which is a set of parameters associated with a class.
 
     set CLASS NAME
     {
       PARAMETER = VALUE
       ...
     }
 
 CLASS should be one of the predefined object class (see @ref ObjectGroup).\n
 NAME should be a string, starting with a letter, followed by alphanumeric characters.
 The underscore character is also allowed, as in `fiber_23`\n
 
 It is also possible to use 'set' to change values of an existing Property:
 
     set NAME
     {
       PARAMETER = VALUE
       ...
     }
 
 This is equivalent to:
 
     set NAME { PARAMETER = VALUE }
 
 or:
 
     set NAME PARAMETER { VALUE }
 
 */

void Parser::parse_set(std::istream& is)
{
    std::streampos ipos = is.tellg();
    std::string cat = Tokenizer::get_symbol(is);
    std::string name, para, blok;
    
#if BACKWARD_COMPATIBILITY <= 50
    {
        /* Read ouput config files anterior to 3.11.2017, which included
         a identity number ('set hand 2 kinesin') */
        size_t inx = 0;
        Tokenizer::get_integer(is, inx);
    }
#endif
    
    Glossary opt;
    Property * pp = nullptr;

    bool spec = ( is.peek() == ':' );
    
    if ( cat == "simul" )
    {
#if BACKWARD_COMPATIBILITY < 50
        // Patch to accept 'set simul:display NAME {}':
        if ( spec )
        {
            is.get(); // skip ':'
            para = Tokenizer::get_symbol(is);
        }
#endif
        name = Tokenizer::get_token(is);
        blok = Tokenizer::get_block(is, '{', true);

        if ( do_change )
        {
#if BACKWARD_COMPATIBILITY < 50
            if ( spec )
                opt.define(para, blok);
            else
#endif
            opt.read(blok);
            execute_change(simul_.prop, opt);
            if ( name != "*" )
                simul_.rename(name);
        }
#if BACKWARD_COMPATIBILITY < 50
        else if ( para == "display" )
        {
            opt.define(para, blok);
            execute_change(simul_.prop, opt);
        }
#endif
    }
    else if ( simul_.isCategory(cat) && !spec )
    {
        /* in this form:
         set CLASS NAME { PARAMETER = VALUE }
         define a new Property
         */
        name = Tokenizer::get_symbol(is);
        blok = Tokenizer::get_block(is, '{', true);
        
        if ( name.empty() )
            throw InvalidParameter("unexpected syntax");

        if ( do_set )
        {
            opt.read(blok);
            pp = execute_set(cat, name, opt);
            
            PropertyID ix;
#if BACKWARD_COMPATIBILITY < 50
            // name changed to `property_number` on 10.12.2017
            if ( opt.set(ix, "identity", "identification") || opt.set(ix, "property_number", "property_index") )
#else
            // name changed to `identification` on 22.06.2021
            // name changed to `identity` on 1.10.2022
            if ( opt.set(ix, "identity", "identification", "property_number") )
#endif
            {
                if ( ix != pp->number() )
                    throw InvalidSyntax("Property identity missmatch");
            }
        }
        else if ( do_change )
        {
            opt.read(blok);
            execute_change(name, opt, false);
        }
    }
    else
    {
        name = cat;
        // in this form, 'set' changes the value of an existing Property
#if BACKWARD_COMPATIBILITY < 50
        if ( spec )
        {
            //set CLASS:PARAMETER NAME { VALUE }
            is.get(); // skip ':'
            para = Tokenizer::get_symbol(is);
            name = Tokenizer::get_token(is);
        }
        else
#endif
        {
            // set NAME PARAMETER { VALUE }
            para = Tokenizer::get_symbol(is);
        }
        
        // set NAME { PARAMETER = VALUE }
        blok = Tokenizer::get_block(is, '{', true);
        
        if ( do_change )
        {
            if ( para.empty() )
                opt.read(blok);
            else
                opt.define(para, blok);
            pp = execute_change(name, opt, do_set);
        }
        else if ( para == "display" )
        {
            opt.define(para, blok);
            execute_change(name, opt, false);
        }
    }

    if ( pp )
        check_warnings(opt, is, ipos);
}

//------------------------------------------------------------------------------
/**
 Change the value of one (or more) parameters for property `NAME`.

     change NAME
     {
       PARAMETER = VALUE
       ...
     }
 
 Short syntax:
 
    change NAME { PARAMETER = VALUE }

The NAME should have been defined previously with the command `set`.
It is also possible to change all properties of a particular class:
 
     change all CLASS
     {
       PARAMETER = VALUE
       ...
     }

Examples:
 
    change system { viscosity = 0.5; }
    change system display { back_color=red; }
    change actin { rigidity = 1; }
    change microtubule { rigidity = 20 }
    change all fiber { confine = inside, 10; }
    change all fiber display { color = white; }

 */

void Parser::parse_change(std::istream& is)
{
    std::streampos ipos = is.tellg();
    bool change_all = false;
    
    std::string name = Tokenizer::get_symbol(is);
    std::string para;
    
    if ( name == "all" )
    {
        change_all = true;
        name = Tokenizer::get_symbol(is);
        if ( !simul_.isCategory(name) )
            throw InvalidSyntax("`"+name+"' is not a known class of object");
    }

#if BACKWARD_COMPATIBILITY < 50
    // Read formats anterior to 3.11.2017
    if ( is.peek() == ':' )
    {
        //change CLASS:PARAMETER NAME { VALUE }
        //change CLASS:PARAMETER * { VALUE }
        is.get();
        para = Tokenizer::get_symbol(is);
        std::string str = Tokenizer::get_token(is);
        if ( str == "*" )
            change_all = true;
        else
            name = str;
    }
    else
#endif
    {
        //change NAME PARAMETER { VALUE }
        para = Tokenizer::get_symbol(is);
        
#if BACKWARD_COMPATIBILITY < 50
        if ( is.peek() == '*' )
        {
            //change CLASS * { VALUE }
            is.get();
            change_all = true;
        }
        else if ( simul_.findProperty(name, para) )
        {
            //change CLASS NAME { VALUE }
            name = para;
            para = "";
        }
#endif
    }

    //change NAME { VALUE }
    std::string blok = Tokenizer::get_block(is, '{', true);
    
    Glossary opt;
    if ( do_change )
    {
        if ( para.empty() )
            opt.read(blok);
        else
            opt.define(para, blok);
        
        if ( change_all )
            execute_change_all(name, opt);
        else
            execute_change(name, opt, do_set);
 
        if ( do_set )
            check_warnings(opt, is, ipos, ~0U);
    }
    else if ( para == "display" )
    {
        opt.define("display", blok);
        if ( change_all )
            execute_change_all(name, opt);
        else
            execute_change(name, opt, false);
    }
}

//------------------------------------------------------------------------------
/**
 The command `new` creates one or more objects with given specifications:
 
     new [MULTIPLICITY] NAME
     {
       position         = POSITION, [SPACE]
       direction        = DIRECTION
       orientation      = ROTATION, [ROTATION]
       mark             = INTEGER
       required         = INTEGER
     }
 
 The NAME should have been defined previously with the command `set`.\n

 The other parameters are:
 
 Parameter        | Type      | Default | Description                          |
 -----------------|-----------|---------|---------------------------------------
 MULTIPLICITY     | INTEGER   |   1     | number of objects.
 `position`       | POSITION  | random  | initial position of the object.
 `orientation`    | ROTATION  | random  | a rotation specified with respect to the object's center of gravity.
 `orientation[1]` | ROTATION  | none    | a rotation specified around the origin.
 `direction`      | DIRECTION | random  | specifies the direction of a fiber.
 `mark`           | INTEGER   |   0     | specifies a mark to be given to all objects created.
 `required`       | INTEGER   |   0     | minimum number of objects that should be created.
 

 Note that `position` only applies to movable objects, and `orientation` will have an effect only on objects that can be rotated. In addition, `position[1]` and `orientation[1]` are relevant only if `(MULTIPLICITY > 1)`, and do not apply to the first object.\n
 
 
 Short syntax:
 
     new [MULTIPLICITY] NAME ( POSITION )
 
 Shorter syntax:
 
     new [MULTIPLICITY] NAME
 
*/

void Parser::parse_new(std::istream& is)
{
    std::streampos ipos = is.tellg();
    Glossary opt;
    size_t cnt = 1;
    Tokenizer::get_integer(is, cnt);
    std::string name = Tokenizer::get_symbol(is);
    
#if BACKWARD_COMPATIBILITY < 50
    // Read formats anterior to 3.11.2017
    if ( simul_.isCategory(name) )
    {
        std::string str = Tokenizer::get_symbol(is);
        if ( !str.empty() )
            name = str;
    }
#endif

    // Syntax sugar: () specify only position
    std::string blok = Tokenizer::get_block(is, '(');
    
    if ( blok.empty() )
    {
        blok = Tokenizer::get_block(is, '{');
        opt.read(blok);
    }
    else {
        opt.define("position", 0, blok);
    }
    
    if ( name.empty() )
        throw InvalidParameter("unexpected syntax");

    if ( do_new & ( cnt > 0 ))
    {
        if ( opt.num_keys() == 0 )
        {
            execute_new(name, cnt);
        }
        else
        {
            // place each object independently from the others:
            execute_new(name, opt, cnt);
            
            if ( opt.has_key("display") )
                throw InvalidParameter("display parameters should be specified within `set'");
            
            check_warnings(opt, is, ipos, ~0U);
        }
    }
}

//------------------------------------------------------------------------------
/**
 Delete objects:

     delete [MULTIPLICITY] NAME
     {
        mark      = INTEGER
        position  = [inside|outside], SPACE
        state     = [0|1], [0|1]
     }
 
 MULTIPLICITY is an integer, or the keyword 'all'.
 
 The parameters (mark, position, state) are all optional.
 All specified conditions must be fulfilled (this is a logical AND).
 The parameter `state` refers to bound/unbound state of Hands for Single and Couple,
 and to dynamic state for Fibers:
 - for Single, `state[0]` refers to the Hand: 0=free, 1=attached.
 - for Couple, `state[0]` refers to the first Hand: 0=free, 1=attached,
           and `state[1]` refers to the second Hand: 0=free, 1=attached.
 - for Fibers, `state[0]` refers to the Dynanic state of the PLUS end,
           and `state[1]` refers to the Dynanic state of the MINUS end.
 .
 
 To delete all objects of specified NAME:
 
     delete all NAME
 
 To delete at most COUNT objects of class NAME:
 
     delete COUNT NAME
 
 To delete all objects with a specified mark:
 
     delete all NAME
     {
        mark = INTEGER
     }
 
 To delete all objects within a Space:
 
     delete NAME
     {
        position = inside, SPACE
     }
 
 The SPACE must be the name of an existing Space.
 Only 'inside' and 'outside' are valid specifications.

 To delete all Couples called NAME that are not bound:
 
     delete all NAME { state1 = 0; state2 = 0; }

 To delete all Couples called NAME that are not crosslinking, use two calls:

     delete all NAME { state1 = 0; }
     delete all NAME { state2 = 0; }
 
 */

void Parser::parse_delete(std::istream& is)
{
    std::streampos ipos = is.tellg();
    size_t cnt = 1;
    bool has_cnt = Tokenizer::get_integer(is, cnt);
    std::string name = Tokenizer::get_symbol(is);
#if BACKWARD_COMPATIBILITY < 50
    // Read formats anterior to 3.11.2017
    if ( simul_.isCategory(name) )
    {
        std::string str = Tokenizer::get_symbol(is);
        if ( !str.empty() )
            name = str;
    }
#endif
    if ( !has_cnt  &&  name == "all" )
    {
        cnt = ~0UL; // this is very large
        name = Tokenizer::get_symbol(is);
    }
    std::string blok = Tokenizer::get_block(is, '{');
    
    if ( name.empty() )
        throw InvalidParameter("unexpected syntax");

    if ( do_new )
    {
        Glossary opt(blok);
        execute_delete(name, opt, cnt);
        check_warnings(opt, is, ipos);
    }
}


/**
 Mark objects:
 
     mark [MULTIPLICITY] NAME
     {
       mark       = INTEGER
       position   = POSITION
     }
 
 or
 
     mark all NAME
     {
         mark       = INTEGER
         position   = POSITION
     }

 NAME can be '*', and the parameter 'position' is optional.
 The syntax is the same as for command `delete`.
 */

void Parser::parse_mark(std::istream& is)
{
    std::streampos ipos = is.tellg();
    size_t cnt = 0;
    bool has_cnt = Tokenizer::get_integer(is, cnt);
    std::string name = Tokenizer::get_symbol(is);
#if BACKWARD_COMPATIBILITY < 50
    // Read formats anterior to 3.11.2017
    if ( simul_.isCategory(name) )
    {
        std::string str = Tokenizer::get_symbol(is);
        if ( !str.empty() )
            name = str;
    }
#endif
    if ( !has_cnt  &&  name == "all" )
        name = Tokenizer::get_symbol(is);
    std::string blok = Tokenizer::get_block(is, '{');
    
    if ( do_new )
    {
        Glossary opt(blok);
        execute_mark(name, opt, cnt);
        check_warnings(opt, is, ipos);
    }
}

//------------------------------------------------------------------------------
/**
 Cut all fibers that intersect a given plane.
 
     cut fiber NAME
     {
        plane = VECTOR, REAL
     }
 
     cut all fiber
     {
        plane = VECTOR, REAL
     }

 The plane is specified by a normal vector `n` (VECTOR) and a scalar `a` (REAL).
 The plane is defined by <em> n.pos + a = 0 </em>
 */

void Parser::parse_cut(std::istream& is)
{    
    std::streampos ipos = is.tellg();
    std::string str = Tokenizer::get_token(is);
 
    if ( str == "all" )
    {
        if ( Tokenizer::get_symbol(is) != "fiber" )
            throw InvalidParameter("only 'fiber' can be cut");
    }
    else if ( str == "fiber" )
    {
        str = Tokenizer::get_symbol(is);
    }
    
    std::string blok = Tokenizer::get_block(is, '{', true);
    
    if ( do_run )
    {
        Glossary opt(blok);
        execute_cut(str, opt);
        check_warnings(opt, is, ipos);
    }
}


//------------------------------------------------------------------------------
/**
 Attach Couple or Single to position on fibers
 
     connect COUPLE_NAME
     {
        
     }
 */
void Parser::parse_connect(std::istream& is)
{
    std::streampos ipos = is.tellg();
    std::string str = Tokenizer::get_token(is);
    std::string blok = Tokenizer::get_block(is, '{');
    
    if ( do_run )
    {
        Glossary opt(blok);
        execute_connect(str, opt);
        check_warnings(opt, is, ipos);
    }
}


//------------------------------------------------------------------------------

/**
 @copydetails Interface::execute_run
 */
void Parser::parse_run(std::istream& is)
{
    std::streampos ipos = is.tellg();
    size_t cnt = 1;
    bool has_cnt = Tokenizer::get_integer(is, cnt);
    std::string name = Tokenizer::get_symbol(is);
    
#if BACKWARD_COMPATIBILITY < 50
    // Read formats anterior to 3.11.2017
    if ( name == "simul" )
    {
        name = Tokenizer::get_symbol(is);
        if ( is.peek() == '*' )
        {
            is.get();
            name = simul_.prop->name();
        }
    }
#endif
    if ( name.empty() )
        throw InvalidSyntax("unexpected syntax (use `run NB_STEPS NAME_OF_SIMUL { }')");
    
    if ( name == "all" )
    {
        if ( Tokenizer::get_symbol(is) != "simul_" )
            throw InvalidSyntax("expected `run all simul_ { }')");
        // There can only be one Simul object:
        name = simul_.prop->name();
    }
    
    if ( name != "*"  &&  name != simul_.prop->name() )
        throw InvalidSyntax("unknown simul name `"+name+"'");

    std::string blok = Tokenizer::get_block(is, '{');
    
    if ( do_run )
    {
        Glossary opt(blok);

        if ( !has_cnt )
        {
            // read `nb_steps' from the option block:
            if ( opt.set(cnt, "nb_steps") )
                opt.clear("nb_steps");
            else
            {
                // instead of `nb_steps', user can specify a duration in seconds:
                real span = 0;
                if ( opt.set(span, "duration", "time") )
                {
                    if ( span <= 0 )
                        throw InvalidParameter("duration must be >= 0'");
                    cnt = (size_t)std::ceil(span/simul_.time_step());
                    opt.clear("duration");
                    opt.clear("time");
                }
            }
        }
        else if ( opt.has_key("nb_steps") )
            throw InvalidSyntax("the number of steps was specified twice");
        
        if ( opt.empty() )
            execute_run(cnt);
        else
            execute_run(cnt, opt, do_write);

        check_warnings(opt, is, ipos);
    }
}

//------------------------------------------------------------------------------
/**
 Read and execute commands from another config file.
 
     read FILE_NAME
     {
       required = BOOL
     }
 
 By default, `required = 1`, and execution will terminate if the file is not found.
 However, if `required=0`, the file will be executed if it is found, but execution
 will continue otherwise.
 
 \todo: able to specify `do_set` and `do_new` for command 'read'
*/

void Parser::parse_read(std::istream& is)
{
    bool required = true;
    std::streampos ipos = is.tellg();
    std::string file = Tokenizer::get_path(is);
    
    if ( file.empty() )
        throw InvalidSyntax("missing/invalid file name after 'read'");
    
    std::string blok = Tokenizer::get_block(is, '{');
    if ( ! blok.empty() )
    {
        Glossary opt(blok);
        opt.set(required, "required");
        check_warnings(opt, is, ipos);
    }
    
    if ( FilePath::is_file(file) )
        readConfig(file);
    else
    {
        if ( required )
            throw InvalidSyntax("could not open file `"+file+"'");
        else
            Cytosim::warn << "could not open file `" << file << "\n";
    }
}

//------------------------------------------------------------------------------
/**
 Import a simulation snapshot from a trajectory file
 
    import WHAT FILE_NAME
    {
        append = BOOL
        frame = INTEGER
    }
 
 The frame to be imported can be specified as an option: `frame=INTEGER`:
 
     import all my_file.cmo { frame = 10 }
 
 By default, this will replace the simulation state by the one loaded from file.
 To add the file objects to the simulation without deleting any of the current 
 object, you should specify `append = 1`:
 
     import all my_file.cmo { append = 1 }
 
 Finally instead of importing all the objects from the file, one can restrict
 the import to a desired class:
 
     import fiber my_file.cmo { append = 1 }
 
 Note that the simulation time will be changed to the one specified in the file,
 but this behavior can be changed by specifying the time:
 
     change system { time = 0 }
 
 ...assuming that the simul is called `system`.
 */

void Parser::parse_import(std::istream& is)
{
    std::streampos ipos = is.tellg();
    std::string what = Tokenizer::get_token(is);
    std::string file = Tokenizer::get_path(is);
    
    if ( what.empty() )
        throw InvalidSyntax("missing class specification (use `import all FILENAME')");

    if ( file.empty() )
        throw InvalidSyntax("missing/invalid file name (use `import all FILENAME')");
    
    std::string blok = Tokenizer::get_block(is, '{');
    
    if ( do_new )
    {
        Glossary opt(blok);
        execute_import(file, what, opt);
        check_warnings(opt, is, ipos);
    }
}


/**
 Export state to file. The general syntax is:
 
     export WHAT FILE_NAME
     {
       append = BOOL
       binary = BOOL
     }
 
 WHAT must be ``objects`` of ``properties``, and by default, both `binary` 
 and `append` are `true`. If `*` is specified instead of a file name,
 the current trajectory file will be used.
 
 Short syntax:
 
     export objects FILE_NAME
 
 Examples:
 
     export all sim_objects.cmo { append=0 }
     export properties properties.txt
 
 Attention: this command is disabled for `play`.
 */

void Parser::parse_export(std::istream& is)
{
    std::streampos ipos = is.tellg();
    std::string what = Tokenizer::get_token(is);
    std::string file = Tokenizer::get_path(is);
    
    if ( what.empty() )
        throw InvalidSyntax("missing class specification (use `export all FILENAME')");
    
    if ( file.empty() )
        throw InvalidSyntax("missing/invalid file name (use `export all FILENAME')");

    std::string blok = Tokenizer::get_block(is, '{');
    
    if ( do_write )
    {
        Glossary opt(blok);
        execute_export(file, what, opt);
        check_warnings(opt, is, ipos);
    }
}


/**
 Export formatted data to file. The general syntax is:
 
     report WHAT FILE_NAME
     {
       append = BOOL
     }
 
 Short syntax:
 
     report WHAT FILE_NAME
 
 WHAT should be a valid argument to `report`:
 @copydetails Simul::report
 
 If `*` is specified instead of a file name, the report is sent to the standard output.
 
 Examples:
 
     report parameters parameters.txt { append=0 }
     report fiber:length fibers_length.txt
 
 Note that writing to a file is normally disabled for `play`.
 */

void Parser::parse_report(std::istream& is)
{
    std::streampos ipos = is.tellg();
    std::string what = Tokenizer::get_polysymbol(is);
    while ( is.peek() == ',' )
        what += char(is.get()) + Tokenizer::get_polysymbol(is);
    std::string file = Tokenizer::get_path(is);
    
    if ( file.empty() )
        throw InvalidSyntax("expected 'report CLASS:REPORT FILE'");
    
    std::string blok = Tokenizer::get_block(is, '{');
    
    if ( do_run && ( do_write || file == "*" ))
    {
        Glossary opt(blok);
        execute_report(file, what, opt);
        check_warnings(opt, is, ipos);
    }
}

/**
Export formatted data to file. The general syntax is:

    write FILE_NAME WHAT [WHAT] ...
    {
      append = BOOL
    }

Short syntax:

    write FILE_NAME WHAT [WHAT] ...

A `*` can be specified instead of a file name, to designate the standard output.
WHAT should be a valid argument to `report`:
@copydetails Simul::report

Examples:

    write * fiber:energy
    write data.txt microtubule:length actin:length { verbose = 0 }

Note that writing to a file is normally disabled for `play`.
*/

void Parser::parse_write(std::istream& is)
{
    std::streampos ipos = is.tellg();
    std::string file = Tokenizer::get_path(is);
    if ( file.empty() )
        throw InvalidSyntax("expected 'write FILE WHAT'");

    std::string what = Tokenizer::get_polysymbol(is);
    std::string more = Tokenizer::get_polysymbol(is);
    while ( !more.empty() )
    {
        what += "," + more;
        more = Tokenizer::get_polysymbol(is);
    }
    
    std::string blok = Tokenizer::get_block(is, '{');
    
    if ( do_run && ( do_write || file == "*" ))
    {
        Glossary opt(blok);
        execute_report(file, what, opt);
        check_warnings(opt, is, ipos);
    }
}

//------------------------------------------------------------------------------

/**
 Call custom function:
 
     call FUNCTION_NAME { OPTIONS }
 
 FUNCTION_NAME should be `equilibrate`, `custom0`, `custom1`, ... `custom9`.

 Note: The Simul::custom() functions need to be modified, to do something!
 */
void Parser::parse_call(std::istream& is)
{
    std::streampos ipos = is.tellg();
    std::string str = Tokenizer::get_symbol(is);
    
    if ( str.empty() )
        throw InvalidSyntax("missing function name after 'call'");
    
    std::string blok = Tokenizer::get_block(is, '{');
    
    if ( do_run )
    {
        Glossary opt(blok);
        execute_call(str, opt);
        check_warnings(opt, is, ipos);
    }
}

//------------------------------------------------------------------------------

/**
 Repeat specified code a number of times.
 
     repeat INTEGER { CODE }
 
 Example:
 
     repeat 100
     {
         run 1000 system
         report microtubule:plus_state states.txt
     }

 */
void Parser::parse_repeat(std::istream& is)
{
    size_t cnt = 1;
    
    if ( ! Tokenizer::get_integer(is, cnt) )
        throw InvalidSyntax("expected positive integer after 'repeat'");

    std::string code = Tokenizer::get_block(is, '{');
    
    for ( size_t c = 0; c < cnt; ++c )
    {
        if ( simul_.prop->verbose )
            Cytosim::log("repeat code %lu/%lu\n", c+1, cnt);
        evaluate(code);
    }
}


/**
 Repeat specified code, with variable substitution
 
     for VAR=INTEGER:INTEGER { CODE }
 
 The two integers are the first and last iterations counters.
 Any occurence of VAR is replaced before the code is executed.
 
 Example:
 
     for CNT=1:10 {
       new 10 filament { length = CNT }
     }
 
 NOTE: This code is a hack, and can be improved vastly!
 */
void Parser::parse_for(std::istream& is)
{
    size_t start = 0, end = 1;
    
    std::string var = Tokenizer::get_symbol(is);
    
    std::string s = Tokenizer::get_token(is);
    if ( s != "=" )
        throw InvalidSyntax("missing '=' in command 'for'");
    
    if ( ! Tokenizer::get_integer(is, start) )
        throw InvalidSyntax("missing number after 'repeat'");

    s = Tokenizer::get_token(is);
    if ( s != ":" )
        throw InvalidSyntax("missing ':' in command 'for'");
    
    if ( ! Tokenizer::get_integer(is, end) )
        throw InvalidSyntax("missing number after 'repeat'");
    
    std::string code = Tokenizer::get_block(is, '{');
    
    for ( size_t c = start; c < end; ++c )
    {
        std::string sub = code;
        // substitute Variable name for this iteration:
        StreamFunc::find_and_replace(sub, var, std::to_string(c));
        // execute code:
        evaluate(sub);
    }
}

//------------------------------------------------------------------------------

/**
 Terminates execution
 
     end
 
 */
void Parser::parse_end(std::istream& is)
{
    if ( do_run )
        throw Exception("terminating program at command 'end'");

    /*
    std::string str = Tokenizer::get_symbol(is);
    
    if ( str == "if" )
    {
        str = Tokenizer::get_token(is);
        ABORT_NOW("the 'if' condition has not been implemented yet");
    }
     */
}

/**
 Save current system's matrix and right-hand-side vector in sub directory
 
     dump DIRECTORY_NAME { mode = {1, 2, 4} }
 
 */
void Parser::parse_dump(std::istream& is)
{
    const std::string str = Tokenizer::get_token(is);
    if ( str.empty() )
        throw InvalidSyntax("missing directory name after 'dump'");
    std::string blok = Tokenizer::get_block(is, '{');

    if ( do_write && do_run )
    {
        Glossary opt(blok);
        int mode = 1;
        opt.set(mode, "mode");

        simul_.sMeca.doNotify = 1;
        simul_.solve_half();
        
        int cwd = FilePath::change_dir(str, true);
        Cytosim::log("Cytosim is dumping a system of size %lu in `%s'...", simul_.sMeca.dimension(), str.c_str());

        if ( mode & 1 ) simul_.sMeca.saveSystem();
        if ( mode & 2 ) simul_.sMeca.dumpSystem();
        if ( mode & 4 ) simul_.sMeca.exportSystem();
        
        FilePath::change_dir(cwd);
        Cytosim::log("done\n");
    }
}

//------------------------------------------------------------------------------
#pragma mark -

/**
 Read and execute the next command to be found in the stream.
 Returns:
 - 0 if successfull
 - 1 if file is exhausted, or has error
 - 2 if 'end' was found.
 Thus parsing should be repeated while the return value is 0.
 */
int Parser::evaluate_one(std::istream& is)
{
    std::string tok = Tokenizer::get_token(is);
    
    if ( tok == "set" )
        parse_set(is);
    else if ( tok == "change" )
        parse_change(is);
    else if ( tok == "new" || tok == "add" )
        parse_new(is);
    else if ( tok == "delete" )
        parse_delete(is);
    else if ( tok == "mark" )
        parse_mark(is);
    else if ( tok == "run" )
        parse_run(is);
    else if ( tok == "read" || tok == "include" )
        parse_read(is);
    else if ( tok == "report" )
        parse_report(is);
    else if ( tok == "write" )
        parse_write(is);
    else if ( tok == "import" )
        parse_import(is);
    else if ( tok == "export" )
        parse_export(is);
    else if ( tok == "call" )
        parse_call(is);
    else if ( tok == "repeat" )
        parse_repeat(is);
    else if ( tok == "for" )
        parse_for(is);
    else if ( tok == "cut" )
        parse_cut(is);
    else if ( tok == "connect" )
        parse_connect(is);
    else if ( tok == "restart" )
    {
        // reset simulation and rewind config file, repeating forever
        simul_.erase_all(1);
        is.clear();
        is.seekg(0);
        return 0;
    }
    else if ( tok == "stop" )
        return 2;
    else if ( tok == ";" )
        return 0;
    else if ( tok == "dump" )
        parse_dump(is);
    else {
        throw InvalidSyntax("unexpected command `"+tok+"'");
    }
    return 0;
}


void Parser::evaluate(std::istream& is)
{
    std::streampos ipos(0);
    try {
        while ( is.good() )
        {
            int c = Tokenizer::skip_space(is, true);
            if ( c == EOF )
                return;
            
            // skip matlab-style comments (%{ })
            if ( c == '%' )
            {
                c = is.get();
                int d = is.get();
                if ( d == '{' )
                    Tokenizer::get_block_text(is, 0, '}');
                else
                    Tokenizer::get_line(is);
                continue;
            }
#if 0
            /*
             skip C-style comments:
             - single-line comment start with '/' and '/'
             - multi-line comments start with '/' and '*'
             */
            if ( c == '/' )
            {
                c = is.get();
                int d = is.get();
                if ( '/' == d )
                    Tokenizer::get_line(is);
                else if ( '*' == d )
                    Tokenizer::get_until(is, "*/");
                else
                    throw InvalidSyntax("unexpected token `/"+std::string(d,1)+"'");
                continue;
            }
#endif
            ipos = is.tellg();
            //StreamFunc::print_lines(std::clog, is, ipos, ipos);
            
            if ( evaluate_one(is) )
                return;
        }
    }
    catch( Exception & e )
    {
        e << "\n" + StreamFunc::get_lines(is, ipos, is.tellg());
        throw;
    }
}


void Parser::evaluate(std::string const& code)
{
    std::istringstream is(code);
    evaluate(is);
}


void Parser::readConfig(std::string const& filename)
{
    std::ifstream is(filename.c_str(), std::ifstream::in);
    if ( !is.is_open() )
        throw InvalidIO("could not open `"+filename+"'");
    if ( !is.good() )
        throw InvalidIO("could not find or read `"+filename+"'");
    VLOG("--Parse `" << filename << "'  set " << do_set << "  change " << do_change);
    VLOG("  new " << do_new << "  run " << do_run << "  write " << do_write << "\n");
    evaluate(is);
}


void Parser::readConfig()
{
    readConfig(simul_.prop->config_file);
}

