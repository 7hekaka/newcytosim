// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University

#ifndef OBJECT_H
#define OBJECT_H

#include "blinksort.h"
#include "inventoried.h"
#include "movable.h"
#include "random.h"
#include "array.h"

class Simul;
class Modulo;
class Property;
class Inputter;
class Outputter;
class ObjectSet;
class Display;

/// Type for unique class identifier used to store objects in file
typedef unsigned char ObjectTag;

/// Integer type used to mark objects
typedef unsigned int ObjectMark;

/// Interger type used to flag objects
typedef unsigned int ObjectFlag;


/// Parent class for all simulated objects
/**
 This is the interface used for writing / reading from a file.
 
 Three functions identify an Object:
 - tag() [ASCII character] identifies the class of Object.
 - property->number() [integer] identifies its Property.
 - identity() [serial-number] derived from Inventoried identifies unique instantiations.
 .
 These three qualities are concatenated in reference() and writeReference().
 
 Objects are stored in ObjectSet.
 */
class Object : public Movable, public Inventoried
{
    /// to expose next and prev pointers for sorting:
    friend class BlinkSortJob<Object>;

    /// allow container class to access members
    friend class ObjectPool;

protected:
    
    /// the next Object in the list
    Object * next_;
    
    /// the previous Object in the list
    Object * prev_;

private:
    
    /// upstream pointer to container class
    ObjectSet * set_;

    /// integer used for user controlled tasks, recorded to file
    ObjectMark mark_;

    /// integer used for private tasks, not saved to file
    ObjectFlag flag_;

public:
    
    /// the highest bit that is not used by ASCII codes
    static constexpr uint8_t HIGH_BIT = 128;
    
    /// bit mask for those bits which are used by ASCII codes
    static constexpr uint8_t LOW_BITS = 127;

    /// Object::NULL_TAG is the 'void' pointer
    static constexpr ObjectTag NULL_TAG = '!';
    
    /// build a reference string by concatenating (tag, property_number, ObjectID)
    static std::string make_reference(ObjectTag, unsigned, ObjectID);
    
    /// write a reference, but using the provided Tag
    static void writeReference(Outputter&, ObjectTag, ObjectID);

    /// write a reference that identifies the Object uniquely
    static void writeReference(Outputter&, Object const*);
    
    /// write header to object data, using provided tag
    void writeMarker(Outputter&, ObjectTag) const;

public:
    
    /// constructor
    Object() : next_(nullptr), prev_(nullptr), set_(nullptr), mark_(0), flag_(0) { }

    /// copy constructor
    Object(Object const& o) : next_(nullptr), prev_(nullptr), set_(nullptr), mark_(o.mark_), flag_(o.flag_) {}
    
    /// assignment operator
    Object& operator = (const Object& o) { next_=nullptr; prev_=nullptr; set_=nullptr; mark_=o.mark_; flag_=o.flag_; return *this; }
    
    /// destructor
    virtual ~Object();
    
    
    /// a character identifying the class of this object
    virtual ObjectTag tag() const { return NULL_TAG; }
    
    /// Property associated with the Object
    virtual Property const* property() const { return nullptr; };
    
    /// write Object data to file
    virtual void write(Outputter&) const = 0;
    
    /// read Object from file, within the Simul
    virtual void read(Inputter&, Simul&, ObjectTag) = 0;
    
    /// return some characteristics of the object, used for reporting
    virtual void report(std::ostream&) const { }
    
    /// return non-zero value if underlying data is invalid somehow
    virtual int invalid() const { return 0; }

    //--------------------------
    
    /// the next Object in the list, or zero if this is last
    Object * next() const { return next_; }
    
    /// the previous Object in the list, or zero if this is first
    Object * prev() const { return prev_; }

    /// set pointer to next Object
    void next(Object* n) { next_ = n; }
    
    /// set pointer to previous Object
    void prev(Object* n) { prev_ = n; }

    //--------------------------

    /// returns container ObjectSet
    ObjectSet * objset() const { return set_; }
    
    /// returns container Simul
    Simul & simul() const;
    
    /// change container class
    void objset(ObjectSet* s) { set_ = s; }
    
    /// true if Object is registered in a container class
    bool linked() const { return set_ != nullptr; }

    /// concatenation of [ tag(), property()->number(), identity() ] in plain ascii
    std::string reference() const;
    
    //--------------------------

    /// get mark
    ObjectMark mark() const { return mark_; }
    
    /// set mark
    void mark(ObjectMark m) { mark_ = m; }
    
    
    /// retrieve flag value
    ObjectFlag flag() const { return flag_; }
    
    /// set flag (this value is not stored in trajectory files)
    void flag(ObjectFlag f) { flag_ = f; }

};


/// return always 'true'
bool match_all(Object const*, void const*);

/// return 'true' if ( obj->mark() == *mark )
bool match_mark(Object const* obj, void const* mrk);

/// return 'true' if ( obj->property() == val )
bool match_property(Object const* obj, void const* val);


/// a variable list of pointers to Object
typedef Array<Object *, 4> ObjectList;
//typedef std::vector<Object *> ObjectList;


/// output operator
std::ostream& operator << (std::ostream& os, ObjectList const&);


#endif
