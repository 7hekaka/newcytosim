// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef ASTER_H
#define ASTER_H

#include "object.h"
#include "organizer.h"
#include "aster_prop.h"
#include "solid.h"
#include "fiber.h"

/// A connection between a Fiber and a Solid
//@todo new Interpolation4() to replace coef1[] and coef2[]
class AsterLink
{
    friend class Aster;
    
private:
    
    /// type of link:
    /**
     0 = no link
     1 = the interpolation corresponds exactly to point 'prime_'
     2 or 3 = link fiber-end with coef1_, fiber-side with coef2_
     */
    size_t   rank_;

    /// index of first point on the Solid
    size_t   prime_;
    
    /// interpolation coefficient for Fiber end
    real     coef1_[4];
    
    /// interpolation coefficient for Fiber side
    real     coef2_[4];
    
    /// distance between the two anchoring points
    real     len_;
    
    /// index used for backward compatibility
    size_t   alt_;
    
public:
    
    /// constructor
    AsterLink()
    {
        reset();
    }
    
    /// set to zero
    void reset()
    {
        rank_ = 0;
        prime_ = 0;
        len_ = 0;
        for ( int i = 0; i < 4; ++i )
        {
            coef1_[i] = 0.0;
            coef2_[i] = 0.0;
        }
        alt_ = 0;
    }
    
    /// calculate rank: how many coefficients are not null
    size_t calcRank() const
    {
        size_t res = 4;
        while ((res > 0) & (abs_real(coef1_[res-1]) < REAL_EPSILON))
            --res;
        return res;
    }

    /// set to interpolate from A to B
    void set(Vector const& A, Vector const& B)
    {
        len_ = ( A - B ).norm();
        
        coef1_[1] = A.XX;
        coef2_[1] = B.XX;
#if ( DIM == 1 )
        coef1_[2] = 0.0;
        coef2_[2] = 0.0;
        coef1_[3] = 0.0;
        coef2_[3] = 0.0;
        coef1_[0] = 1.0 - A.XX;
        coef2_[0] = 1.0 - B.XX;
#elif ( DIM == 2 )
        coef1_[2] = A.YY;
        coef2_[2] = B.YY;
        coef1_[3] = 0.0;
        coef2_[3] = 0.0;
        coef1_[0] = 1.0 - A.XX - A.YY;
        coef2_[0] = 1.0 - B.XX - B.YY;
#elif ( DIM == 3 )
        coef1_[2] = A.YY;
        coef2_[2] = B.YY;
        coef1_[3] = A.ZZ;
        coef2_[3] = B.ZZ;
        coef1_[0] = 1.0 - A.XX - A.YY - A.ZZ;
        coef2_[0] = 1.0 - B.XX - B.YY - B.ZZ;
#endif
        rank_ = calcRank();
    }

    /// save coefficient to file
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
        coef1_[0] = 1.0 - coef1_[1] - coef1_[2] - coef1_[3];
        
        for ( int d = 1; d < 4; ++d )
            coef2_[d] = in.readFloat();
        coef2_[0] = 1.0 - coef2_[1] - coef2_[2] - coef2_[3];
        
        len_ = rad * ( Vector3(coef1_+1) - Vector3(coef2_+1) ).norm();
        rank_ = calcRank();
    }
    
#ifdef BACKWARD_COMPATIBILITY
    void readOldFormat(Inputter& in, Solid const* sol)
    {
        reset();
        prime_ = in.readUInt16();
        alt_ = in.readUInt16();
        coef1_[0] = 1.0;
        len_ = ( sol->posPoint(prime_) - sol->posPoint(alt_) ).norm();
        rank_ = 1;
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
    
    Solid*     asSolid;
    
    /// scale of local reference frame
    real       asRadius;
    
    /// store the coefficients needed to make the links between Solid and Fiber
    Array<AsterLink> asLinks;

    /// create and configure the Solid
    ObjectList makeSolid(Simul&, Glossary& opt, size_t& origin);

    /// create a Fiber for position 'inx'
    ObjectList makeFiber(Simul&, size_t inx, std::string const&, Glossary& opt);

    /// define the attachment position of fiber 'inx'
    void       placeAnchor(Vector const&, Vector const&, size_t origin);

    /// define the anchor points of Fibers
    void       placeAnchors(Glossary& opt, size_t origin, size_t nbf);
    
    /// Property
    AsterProp const* prop;
    
public:
    
    /// constructor
    Aster(AsterProp const* p) : asSolid(nullptr), asRadius(0), prop(p) {}
    
    /// destructor
    virtual      ~Aster();
    
    /// construct all the dependent Objects of the Organizer
    ObjectList    build(Glossary&, Simul&);
    
    /// return the scaffolding Solid
    Solid *       solid() const { return asSolid; }
    
    /// return the center of the Solid
    Vector        position() const { return solid()->posP(0); }
    
    /// return Fiber `n`
    size_t        nbFibers() const { return nbOrganized(); }

    /// return Fiber `n`
    Fiber *       fiber(size_t n) const { return Fiber::toFiber(organized(n)); }
    
    /// perform one Monte-Carlo step
    void          step();
    
    /// add interactions to a Meca
    void          setInteractions(Meca&) const;
    
    /// position of first clamp for Fiber n
    Vector        posLink1(size_t n) const;
    
    /// position of second clamp for Fiber n
    Vector        posLink2(size_t n) const;

    /// position of end on Fiber corresponding to first link
    Vector        posFiber1(size_t n) const { return fiber(n)->posEnd(prop->focus); }
    
    /// position of attachment point on Fiber corresponding to second link
    Vector        posFiber2(size_t n) const;
    
    /// retrieve link between Solid and end of Fiber number `i`, returning stiffness
    real          getLink1(size_t i, Vector&, Vector&) const;
    
    /// retrieve link between Solid and side of Fiber number `i`, returning stiffness
    real          getLink2(size_t i, Vector&, Vector&) const;
    
    /// retrieve link of type 1 if `i` is even, of type 2 if `i` is odd
    bool          getLink(size_t i, Vector&, Vector&) const;

    /// return Solid
    Mecable*      core() const { return asSolid; }
    
    /// return PointDisp of Solid
    PointDisp const* disp() const { if ( asSolid ) return asSolid->prop->disp; return nullptr; }

    //--------------------------------------------------------------------------

    /// a unique character identifying the class
    static const ObjectTag TAG = 'a';
    
    /// return unique character identifying the class
    ObjectTag       tag() const { return TAG; }

    /// return associated Property
    Property const* property() const { return prop; }

    /// convert pointer to Aster* if the conversion seems valid; returns 0 otherwise
    static Aster* toAster(Object * obj)
    {
        if ( obj  &&  obj->tag() == TAG )
            return static_cast<Aster*>(obj);
        return nullptr;
    }
    
    //--------------------------------------------------------------------------

    /// read from IO
    void          read(Inputter&, Simul&, ObjectTag);
    
    /// write to IO
    void          write(Outputter&) const;

};


#endif

