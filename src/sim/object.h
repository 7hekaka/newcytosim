// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University

#ifndef OBJECT_H
#define OBJECT_H

#include "inventoried.h"
#include "movable.h"
#include "random.h"
#include "array.h"

class Simul;
class Property;
class Inputter;
class Outputter;
class ObjectSet;
class Display;

/// Type for unique class identifier used to read/write objects from file
typedef char ObjectTag;

/// Type used to mark objects
typedef short ObjectMark;

/// Type used to flag objects
typedef short ObjectFlag;

/// Type used for signature
typedef unsigned ObjectSignature;

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
    /// this class stores the object
    friend class ObjectPool;
    
protected:
    
    /// the next Object in the list
    Object * nextO;
    
    /// the previous Object in the list
    Object * prevO;

private:
    
    /// upstream pointer to container class
    ObjectSet * set_;

    /// integer used for user controlled tasks, recorded to file
    ObjectMark mark_;

    /// integer used for private tasks, not saved to file
    ObjectFlag flag_;
    
    /// a random number associated with this object
    ObjectSignature signature_;
    
public:
    
    /// Object::NULL_TAG = 'v' is the 'void' pointer
    static const ObjectTag NULL_TAG = 'v';
    
    /// build a reference string by concatenating (tag, property_number, ObjectID)
    static std::string reference(ObjectTag, size_t, ObjectID);
    
    /// write a reference, but using the provided Tag
    static void writeReference(Outputter&, ObjectTag, ObjectID);
    
    /// write a reference that does not refer to any Object
    static void writeNullReference(Outputter&);

    /// write a reference that identifies the Object uniquely
    static void writeReference(Outputter&, Object const*);
    
    /// write header to object data, using provided tag
    void writeHeader(Outputter&, ObjectTag) const;

public:
    
    /// constructor
    Object() : nextO(nullptr), prevO(nullptr), set_(nullptr), mark_(0), flag_(0), signature_(RNG.pint32()) { }

    /// copy constructor
    Object(Object const& o) : nextO(nullptr), prevO(nullptr), set_(nullptr), mark_(o.mark_), flag_(o.flag_), signature_(o.signature_) {}
    
    /// assignment operator
    Object& operator =(const Object& o) { nextO=nullptr; prevO=nullptr; set_=nullptr; mark_=o.mark_; flag_=o.flag_; signature_=o.signature_; return *this; }
    
    /// destructor
    virtual ~Object();
    
    
    /// a character identifying the class of this object
    virtual ObjectTag tag() const { return NULL_TAG; }
    
    /// Property associated with the Object
    virtual Property const* property() const = 0;
    
    /// write Object data to file
    virtual void write(Outputter&) const = 0;
    
    /// read Object from file, within the Simul
    virtual void read(Inputter&, Simul&, ObjectTag) = 0;
    
    /// return some characteristics of the object, used for reporting
    virtual void report(std::ostream&) const { }

    //--------------------------
    
    /// the next Object in the list, or zero if this is last
    Object * next() const { return nextO; }
    
    /// the previous Object in the list, or zero if this is first
    Object * prev() const { return prevO; }
    
    /// set next Object
    void next(Object* n) { nextO = n; }
    
    /// set previous Object
    void prev(Object* n) { prevO = n; }

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
    
    /// a random number that makes objects unique
    ObjectSignature signature() const { return signature_; }
    
    /// set signature
    void signature(ObjectSignature s) { signature_ = s; }

};


/// return always 'true'
bool match_all(Object const*, void const*);

/// return 'true' if ( obj->mark() == *mark )
bool match_mark(Object const* obj, void const* mrk);

/// return 'true' if ( obj->property() == val )
bool match_property(Object const* obj, void const* val);


/// a list of pointers to Object
typedef Array<Object *> ObjectList;
//typedef std::vector<Object *> ObjectList;


/// output operator
std::ostream& operator << (std::ostream& os, ObjectList const&);


#endif
