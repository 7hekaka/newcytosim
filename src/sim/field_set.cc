// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "field_set.h"
#include "field_prop.h"
#include "iowrapper.h"
#include "glossary.h"
#include "simul.h"
#include "field.h"


Field * FieldSet::first() const
{
    return static_cast<Field*>(pool_.front());
}


Field * FieldSet::findID(ObjectID n) const
{
    return static_cast<Field*>(inventory_.get(n));
}

//------------------------------------------------------------------------------

void FieldSet::prepare()
{
    for ( Field * f=first(); f; f=f->next() )
    {
        assert_true( f->hasField() );
        f->prepare();
    }
}


void FieldSet::step()
{
    for ( Field * f=first(); f; f=f->next() )
    {
        if ( f->hasField() )
        {
            LOG_ONCE("!!!! Field is active\n");
            f->step(simul_.fibers);
        }
    }
}

//------------------------------------------------------------------------------
#pragma mark -

Property* FieldSet::newProperty(const std::string& cat, const std::string& nom, Glossary&) const
{
    if ( cat == "field" )
        return new FieldProp(nom);
    return nullptr;
}


Object * FieldSet::newObject(const ObjectTag tag, PropertyID pid)
{
    if ( tag == Field::TAG )
    {
        FieldProp * p = simul_.findProperty<FieldProp>("field", pid);
        return new Field(p);
    }
    throw InvalidIO("Warning: unknown Field tag `"+std::to_string(tag)+"'");
    return nullptr;
}


/**
 @ingroup NewObject

 Specify the initial value of the Field:
 
     new field NAME
     {
        value = 0
     }
 
 \todo: read the value of the field from a file, at initialization
 */
void FieldSet::newObjects(ObjectList& res, const std::string& name, Glossary& opt)
{
    Property * p = simul_.properties.find_or_die("field", name);
    FieldProp * fp = static_cast<FieldProp*>(p);
        
    Field * obj = new Field(fp);
        
    // initialize field:
    obj->setField();
        
    // an initial concentration can be specified:
    FieldCell val(0);
    if ( opt.set(val, "value", "initial_value") )
    {
        std::string str;
        if ( opt.set(str, "value", 1) )
        {
            Space const* spc = simul_.findSpace(str);
            if ( !spc )
                spc = obj->prop->field_space_ptr;
            obj->setConcentration(spc, val, 0);
        }
        else
        {
            obj->setConcentration(val);
        }
    }

    res.push_back(obj);
}


void FieldSet::write(Outputter& out) const
{
    if ( size() > 0 )
    {
        out.writeLine("\n#section "+title());
        writeObjects(out, pool_);
    }
}

