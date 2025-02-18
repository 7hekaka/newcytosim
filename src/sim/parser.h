// Cytosim was created by Francois Nedelec.  Copyright 2020 Cambridge University

#ifndef PARSER_H
#define PARSER_H

#include "interface.h"


/// the Parser reads and executes Cytosim config files
/**
 This is also where the syntax of the config file is defined
 */
class Parser : public Interface
{
private:
    
    /// disabled default constructor
    Parser();
    
    /// control switch to enable command 'change' (change a property)
    bool do_change;

    /// control switch to enable command 'set' (creating a property)
    bool do_set;
    
    /// control switch to enable command 'new' and 'delete' (object)
    bool do_new;
    
    /// control switch to enable command 'run' (run simulation)
    bool do_run;
    
    /// control switch to enable command 'write' (write files)
    bool do_write;
    
    /// print warnings
    int do_warn;

protected:
    
    /// counters for the 'restart' command
    size_t restart_;

    /// parse command `set`
    void parse_set(std::istream&);
    
    /// parse command `change`
    void parse_change(std::istream&);
    
    /// parse command `new`
    void parse_new(std::istream&);
    
    /// parse command `delete`
    void parse_delete(std::istream&);
    
    /// parse command `move`
    void parse_move(std::istream&);

    /// parse command `mark`
    void parse_mark(std::istream&);

    /// parse command `cut`
    void parse_cut(std::istream&);
    
    /// parse command `equilibrate`
    void parse_equilibrate(std::istream&);

    /// parse command `run`
    void parse_run(std::istream&);
    
    /// parse command `read`
    void parse_read(std::istream&);
    
    /// parse command `write`
    void parse_write(std::istream&);
    
    /// parse command `report`
    void parse_report(std::istream&);

    /// parse command `import`
    void parse_import(std::istream&);
    
    /// parse command `export`
    void parse_export(std::istream&);
    
    /// parse command `call`
    void parse_call(std::istream&);
    
    /// parse command `repeat`
    void parse_repeat(std::istream&);

    /// parse command `for`
    void parse_for(std::istream&);
    
    /// parse command `end`
    void parse_end(std::istream&);
    
    /// parse command `dump`
    void parse_dump(std::istream&);

    //--------------------------------------------------------------------------

public:
    
    /// construct a Parser with given permissions
    Parser(Simul*, bool Change, bool Set, bool New, bool Run, bool Write, int Warn);
    
    /// Parse next command in stream, advance stream pointer, return 0 if success
    int  evaluate_one(std::istream&);
    
    /// Parse all commands in stream
    void evaluate(std::istream&, std::streampos&);

    /// Parse code in string, and report errors
    void evaluate(std::string const&);
    
    /// Parse the stream. filename is used for error reporting
    void readConfig(std::istream&, std::string const& filename);

    /// Open and parse the config file with the given name
    void readConfig(std::string const& filename);
    
    /// Read the config file, close the file and parse its data
    void readConfigBuffered(std::string const& filename);

    /// Parse the default config file (SimulProp::config)
    void readConfig();
    
    /// execute commands from standard input, return number of lines processed
    size_t read_input(int fildes = 0);

};

#endif

