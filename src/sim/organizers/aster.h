// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University

#ifndef ASTER_H
#define ASTER_H

#include "object.h"
#include "organizer.h"
#include "aster_prop.h"
#include "solid.h"
#include "fiber.h"

/// A connection between a Fiber and a Solid
//@todo use Interpolation4() to replace coef1[] and coef2[]
class AsterLink
{
    friend class Aster;
    
private:
    
    /// interpolation coefficients on Solid for link 1
    real coef1_[4];
    
    /// interpolation coefficients on Solid for link 2
    real coef2_[4];
    
    /// distance between the two anchoring points
    real len_;
    
    /// rank of link 1
    /**
     0 = no link
     1 = the interpolation corresponds exactly to point 'prime_'
     2 or 3 = link fiber-end with coef1_, fiber-side with coef2_
     */
    index_t rank_;

    /// index of first point on the Solid
    index_t prime_;

#if BACKWARD_COMPATIBILITY < 47
    /// index used for backward compatibility
    index_t alt_;
#endif
    
    /// set coefficients
    void set_coef1(real a, real b, real c) { coef1_[0]=1.0-a-b-c; coef1_[1]=a; coef1_[2]=b; coef1_[3]=c; }
    void set_coef2(real a, real b, real c) { coef2_[0]=1.0-a-b-c; coef2_[1]=a; coef2_[2]=b; coef2_[3]=c; }
    
    /// calculate rank: how many coefficients are not null
    void polish()
    {
        // ensures that sum of coefficients is 1
        coef1_[0] = 1.0 - coef1_[1] - coef1_[2] - coef1_[3];
        coef2_[0] = 1.0 - coef2_[1] - coef2_[2] - coef2_[3];
        
        rank_ = 1;
        for ( index_t i = 1; i <= DIM; ++i )
        {
            if ( abs_real(coef1_[i]) > 0 )
                rank_ = i+1;
        }
    }

public:
    
    /// constructor
    AsterLink()
    {
        reset();
    }

    /// set to zero
    void reset()
    {
        rank_ = 1;
        prime_ = 0;
        len_ = 0;
        set_coef1(0.0, 0.0, 0.0);
        set_coef2(0.0, 0.0, 0.0);
#if BACKWARD_COMPATIBILITY < 47
        alt_ = 0;
#endif
    }

    /// set to interpolate from A to B
    void set(Vector const& A, Vector const& B, unsigned P)
    {
        prime_ = P;
        len_ = ( A - B ).norm();
        
#if ( DIM == 1 )
        set_coef1(A.XX, 0.0, 0.0);
        set_coef2(B.XX, 0.0, 0.0);
#elif ( DIM == 2 )
        set_coef1(A.XX, A.YY, 0.0);
        set_coef2(B.XX, B.YY, 0.0);
#else
        set_coef1(A.XX, A.YY, A.ZZ);
        set_coef2(B.XX, B.YY, B.ZZ);
#endif
        polish();
    }

    /// save coefficient to file (not coef_[0])
    void write(Outputter& out) const
    {
        out.writeUInt16(prime_);
        for ( int d = 1; d < 4; ++d )
            out.writeFloat(coef1_[d]);
        for ( int d = 1; d < 4; ++d )
            out.writeFloat(coef2_[d]);
    }
    
    /// read coefficient from file
    void read(Inputter& in, real rad)
    {
        prime_ = in.readUInt16();
        
        for ( int d = 1; d < 4; ++d )
            coef1_[d] = in.readFloat();
        
        for ( int d = 1; d < 4; ++d )
            coef2_[d] = in.readFloat();
        
        len_ = rad * ( Vector3(coef1_+1) - Vector3(coef2_+1) ).norm();
        polish();
    }
    
#if BACKWARD_COMPATIBILITY < 47
    void readOldFormat(Inputter& in, Solid const* sol)
    {
        reset();
        prime_ = in.readUInt16();
        alt_ = in.readUInt16();
        len_ = ( sol->posPoint(prime_) - sol->posPoint(alt_) ).norm();
        if ( prime_ >= sol->nbPoints() )
            throw InvalidIO("invalid AsterLink index");
        if ( alt_ >= sol->nbPoints() )
            throw InvalidIO("invalid AsterLink index");
    }
#endif
    
    /// help to debugging
    void print(std::ostream& out) const
    {
        const unsigned w = 9;
        out << std::setw(w) << coef1_[0];
        for ( int d = 1; d < 4; ++d )
            out << " " << std::setw(w) << coef1_[d];
        out << "   " << std::setw(w) << coef2_[0];
        for ( int d = 1; d < 4; ++d )
            out << " " << std::setw(w) << coef2_[d];
        out << "\n";
    }

    friend std::ostream& operator << (std::ostream& os, AsterLink const& arg)
    {
        arg.print(os);
        return os;
    }
};


/// A radial configuration of Fiber(s) built around a Solid
/**
 The parameters are defined in AsterProp.
 
 Each Fiber is attached to the Solid:
 - at the end of the Fiber
 - at a secondary point that is tied to the Fiber at some distance from this end.
 .
 This anchors the Fiber to the Solid, both in position and direction.
 The stiffness of the links is defined in AsterProp::stiffness, and can be adjusted independently.
 .
 
 @ingroup OrganizerGroup
 */
class Aster : public Organizer
{
private:
    
    Solid * asSolid;
    
    /// scale of local reference frame
    real asRadius;
    
    /// store the coefficients needed to make the links between Solid and Fiber
    Array<AsterLink> asLinks;
    
    /// Property
    AsterProp const* prop;

    
    /// create and configure the Solid
    unsigned makeSolid(ObjectList&, Simul&, Glossary& opt);

    /// create a Fiber and store additional objects in list
    Fiber * makeFiber(ObjectList&, Simul&, Vector, Vector, FiberProp const*, std::string const&);

    /// define the attachment position of fiber 'inx'
    size_t placeAnchor(Vector, Vector, unsigned origin);

    /// create a new fiber
    void placeFiber(ObjectList&, Simul&, Vector, Vector, unsigned ref, std::string const&, std::string const&);
    
    /// position on Solid of the first link to n-th fiber
    Vector posSolid1(size_t n) const;
    
    /// position on Solid of the second link to n-th fiber
    Vector posSolid2(size_t n) const;

    /// position of Fiber's end involved in first link
    Vector posFiber1(size_t n) const { return fiber(n)->posEnd(prop->pole); }
    
    /// position of Fiber's point involved in second link
    Vector posFiber2(size_t n) const;
    
    void build0(ObjectList&, Glossary& opt, Simul&, unsigned);
    void build1(ObjectList&, Glossary& opt, Simul&, unsigned);
    void build2(ObjectList&, Glossary& opt, Simul&, unsigned);
    void build3(ObjectList&, Glossary& opt, Simul&, unsigned);
    void build4(ObjectList&, Glossary& opt, Simul&, unsigned);
    void build7(ObjectList&, Glossary& opt, Simul&, unsigned);

public:
    
    /// constructor
    Aster(AsterProp const* p) : asSolid(nullptr), asRadius(0), prop(p) {}
    
    /// destructor
    virtual ~Aster();
    
    /// construct all the dependent Objects of the Organizer
    ObjectList build(Glossary&, Simul&);
    
    /// return the scaffolding Solid
    Solid * solid() const { return asSolid; }
    
    /// return the center of the Solid
    Vector position() const { return solid()->posP(0); }
    
    /// return Fiber `n`
    size_t nbFibers() const { return nbOrganized(); }

    /// return Fiber `n`
    Fiber * fiber(size_t n) const { return Fiber::toFiber(organized(n)); }
    
    /// perform one Monte-Carlo step
    void step();
    
    /// add interactions to a Meca
    void setInteractions(Meca&) const;
    
    /// number of links to be displayed using getLink()
    index_t nbLinks() const { return 2 * nbFibers(); }

    /// retrieve link between Solid and end of Fiber number `i`, returning stiffness
    real getLink1(index_t i, Vector&, Vector&) const;
    
    /// retrieve link between Solid and side of Fiber number `i`, returning stiffness
    real getLink2(index_t i, Vector&, Vector&) const;
    
    /// retrieve link of type 1 if `i` is even, of type 2 if `i` is odd
    bool getLink(index_t i, Vector&, Vector&) const;

    //--------------------------------------------------------------------------
    
    /// return unique character identifying the class
    ObjectTag tag() const { return Organizer::ASTER_TAG; }

    /// return associated Property
    Property const* property() const { return prop; }

    /// convert pointer to Aster* if the conversion seems valid; returns 0 otherwise
    static Aster* toAster(Object * obj)
    {
        if ( obj  &&  obj->tag() == Organizer::ASTER_TAG )
            return static_cast<Aster*>(obj);
        return nullptr;
    }
    
    //--------------------------------------------------------------------------

    /// read from IO
    void read(Inputter&, Simul&, ObjectTag);
    
    /// write to IO
    void write(Outputter&) const;

};


#endif

