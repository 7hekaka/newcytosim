// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef INTERFACE_H
#define INTERFACE_H

#include <iostream>
#include "isometry.h"
#include "object.h"

class Glossary;
class Property;
class SimulProp;
class Simul;

/// Cytosim Application Programming Interface
/*
 A reduced set of commands to control and simulate
 a system of objects within cytosim.
 */
class Interface
{
private:

    /// Simul member function pointer
    using SimulFuncPtr = void (Simul::*)();
    
    /// perform `cnt` simulation steps also calling Simul::FUNC at each step
    template <SimulFuncPtr FUNC> void step_simul();
    
    /// create 1 object of type `name`, following options in Glossary
    ObjectList new_object(ObjectSet*, Property const*, Glossary&);
    
    /// change values in given Property as specified in Glossary
    void change_property(Property*, Glossary&);
    
    /// read the specification of position and orientation of an object
    void read_placement(Isometry&, Glossary&);
    
    /// return position and orientation of an object, with verification of 'placement'
    bool find_placement(Isometry&, Glossary&, int placement, size_t& nb_trials);
    
protected:
    
    /// associated Simul object
    Simul * sim_;

public:
    
    /// construct and associates with given Simul
    Interface(Simul*);
    
    //-------------------------------------------------------------------------------
    
    /// `hold()` is called between commands during the execution process
    /**
     This callback provides an opportunity to stop/examine the simulation at regular
     intervals. It does nothing for `sim` and displays the system in `play`.
     */
    virtual void hold() {}
    
    /// erase Simulation
    virtual void erase_simul(bool) const;

    //-------------------------------------------------------------------------------

    /// return Simul object pointer
    Simul* simul() const { return sim_; }
    
    /// change Simul pointer
    void set_simul(Simul* s) { sim_ = s; }

    /// return unmodifiable SimulProp
    SimulProp const& simulProp() const;
    
    /// test if 'name' is a category
    bool isCategory(std::string const& name) const;

    //-------------------------------------------------------------------------------

    /// create a new Property of category `cat` from values specified in Glossary
    Property * execute_set(std::string const& cat, std::string const& name, Glossary&);

    /// change values in Property called `name` as specified in Glossary
    Property * execute_change(std::string const& name, Glossary&, bool strict);
    
    /// change values of all Property of category `cat`
    void execute_change_all(std::string const& cat, Glossary&);
    
    /// change some value in the Simul's property
    void change_simul_property(Glossary& opt);
    
    /// create `cnt` objects of type `name`, following options in Glossary
    ObjectList execute_new(std::string const& cat, std::string const& name, Glossary&, size_t cnt);

    /// create `cnt` objects of type `name`, randomly placed in space (no option)
    ObjectList execute_new(std::string const& name, size_t cnt);
    
    /// delete `cnt` objects of type `name`, following options in Glossary
    void execute_delete(std::string const& name, Glossary&, size_t cnt);
    
    /// move object of type `name`, following options in Glossary
    void execute_move(std::string const& name, Glossary&, size_t cnt);
    
    /// mark `cnt` objects of type `name`, following options in Glossary
    void execute_mark(std::string const& name, Glossary&, size_t cnt);

    /// cut fibers of type `name`, following different options in Glossary
    void execute_cut(std::string const& name, Glossary&, size_t cnt);
    
    /// cut fibers of type `name`, following different options in Glossary
    void execute_connect(std::string const& name, Glossary&);

    /// import objects (or `what`) from a file
    void execute_import(std::string const& filename, std::string const& what, Glossary&);

    /// export objects (or `what`) to a file
    void execute_export(std::string const& filename, std::string const& what, Glossary&);

    /// write information (specified in `what`) to a file
    void execute_report(std::string const& filename, std::string const& what, Glossary&);
    
    /// perform `cnt` simulation steps, following options specified in Glossary
    void execute_run(Glossary&, bool write_permission);
    
    /// perform `cnt` simulation steps
    void execute_run();

    /// execute miscellaneous functions
    void execute_call(std::string& func, Glossary&);

    /// dump system and informations
    void execute_dump(std::string const& path, int mode);
};

#endif

