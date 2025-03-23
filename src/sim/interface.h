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

/// Cytosim's Application Programming Interface
/*
 A reduced set of commands to control and simulate
 a system of objects within Cytosim.
 */
class Interface
{
private:

    /// Simul member function pointer
    using SimulFuncPtr = void (Simul::*)();
    
    /// perform `cnt` simulation steps also calling Simul::FUNC at each step
    template <SimulFuncPtr FUNC> void step_simul();
    
    /// usually create one object, following options in Glossary
    ObjectList new_object(ObjectSet*, Property const*, Glossary&);
    
    /// read specifications of position and orientation for an object
    bool read_placement(Isometry&, Glossary&);
    
    /// set position and orientation of an object, according to 'placement'
    bool find_placement(Isometry&, Glossary&, int placement);
    
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
    
    /// erase all simulation objects, and all properties if 'clear_properties==true'
    virtual void eraseSimul(bool clear_properties) const;

    //-------------------------------------------------------------------------------

    /// return Simul object pointer
    Simul* simul() const { return sim_; }
    
    /// change Simul pointer
    void simul(Simul* s) { sim_ = s; }

    /// return unmodifiable SimulProp
    SimulProp const& simulProp() const;
    
    /// test if 'name' is a category
    bool isCategory(std::string const& name) const;
    
    /// determine Property and find set corresponding to object class `name`
    ObjectSet* findClass(std::string const& name, Property*&);
    
    /// find an object specified as `property:identity`, eg. `microtubule1`
    Object* findObject(std::string const& name, Property*&);

    //-------------------------------------------------------------------------------

    /// create a new Property of category `cat` from values specified in Glossary
    Property * execute_set(std::string const& cat, std::string const& name, Glossary&);

    /// change values in Property called `name` as specified in Glossary
    Property * execute_change(std::string const& name, Glossary&, bool strict);
    
    /// change values of all Property of category `cat`
    void execute_change_all(std::string const& cat, Glossary&);
    
    /// change values in given Property as specified in Glossary
    void change_property(Property*, Glossary&);

    /// create `cnt` objects of type `name`, following options in Glossary
    ObjectList execute_new(std::string const& cat, std::string const& name, Glossary&, size_t cnt);

    /// create `cnt` objects of type `name`, randomly placed in space (no option)
    ObjectList execute_new(std::string const& name, size_t cnt, Space const*, std::string const&);
    
    /// delete `cnt` objects of type `name`, following options in Glossary
    void execute_delete(std::string const& name, Glossary&, size_t cnt);
    
    /// move object of type `name`, following options in Glossary
    size_t execute_move(std::string const& name, Glossary&, size_t cnt);
    
    /// mark `cnt` objects of type `name`, following options in Glossary
    void execute_mark(std::string const& name, Glossary&, size_t cnt);

    /// cut fibers of type `name`, following different options in Glossary
    void execute_cut(std::string const& name, Glossary&, size_t cnt);
    
    /// cut fibers of type `name`, following different options in Glossary
    void execute_equilibrate(std::string const& name, Glossary&);

    /// import objects (or `what`) from a file
    void execute_import(std::string const& filename, std::string const& what, Glossary&);

    /// export objects (or `what`) to a file
    void execute_export(std::string const& filename, std::string const& what, Glossary&);

    /// write information (specified in `what`) to a file
    void execute_report(std::string const& filename, std::string const& what, Glossary&);
    
    /// simulate for `sec` seconds, following options specified in Glossary
    void execute_run(real sec, Glossary&, bool write_permission);
    
    /// perform simulation steps to increment time by `sec` seconds
    void execute_run(real sec);

    /// execute miscellaneous functions
    void execute_call(std::string& func, Glossary&);

    /// dump system and informations
    void execute_dump(std::string const& path, unsigned mode);
};

#endif

