// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef LATTICE_H
#define LATTICE_H

#include <cmath>
#include <iostream>
#include "assert_macro.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "real.h"


/// Array of discrete sites aligned with the abscissa of a Fiber
/**
 A Lattice can have different purposes:
 - to limit the local density of attached Hand,
 - to simulate traffic-jam of motors moving on Fibers,
 - to simulate binding/unbinding of a continuous substance.
 .
 
 The first 2 points are implemented in a specialized Hand: Digit.
 
 The index of a site can be positive or negative, but always corresponds
 to the abscissa taken along the Fiber:
 
     index = (site_t) floor( abscissa / unit );
 
 The Lattice may hold different types of information at each site, eg. a 'number of molecule'.
 The linear density can be derived by dividing by the unit length:
 
     density = value(index) / unit();
 
 */
template <typename CELL>
class Lattice
{
public:
    
    /// type used to index the Lattice
    /** This must be a signed integer type */
    typedef int  site_t;

    /// type stored in each Lattice cell
    typedef CELL cell_t;
    
private:
    
    /// lowest valid index (can be negative or positive)
    site_t     laInf;
    
    /// highest valid index plus one (laSup > laInf)
    site_t     laSup;
    
    /// distance between adjacent sites
    real       laUnit;

    /// Original allocated memory ( laSite = laSite0 - laInf )
    cell_t *   laSite0;
    
    /// Array of sites, of valid range [laInf, laSup[
    cell_t *   laSite;
    
    //------------------------------------------------------------------------------
    #pragma mark - Allocation
    
    /// allocate sites for indices within [inf, sup[, conserving existing indices
    site_t allocateLattice(site_t inf, site_t sup, site_t margin)
    {
        assert_true( inf <= sup );
        
        if ( laSite )
        {
            // return if existing boundaries are sufficient
            if ( laInf <= inf  &&  sup <= laSup )
                return 0;
        }
        
        inf -= margin;
        sup += margin;
        
        if ( laSite )
        {
            // only extend the current coverage:
            inf = std::min(inf, laInf);
            sup = std::max(sup, laSup);
        }
        
        //std::clog << this << " Lattice::allocate [" << inf << ", " << sup << "[\n";
        
        const site_t dim = sup - inf;
        cell_t * laSite0_new = new cell_t[dim];
        cell_t * laSite_new  = laSite0_new - inf;
        
        // reset new array:
        for ( site_t s = 0; s < dim; ++s )
            laSite0_new[s] = 0;
        
        // transfer the current data, from the overlapping zone
        if ( laSite )
        {
            for ( site_t s = laInf; s < laSup; ++s )
                laSite_new[s] = laSite[s];
            delete[] laSite0;
        }
        
        assert_true( inf <= sup );

        laInf = inf;
        laSup = sup;
        laSite0 = laSite0_new;
        laSite  = laSite_new;
        
        //std::clog<<"Lattice "<<this<<" covers ["<<laInf<<","<<laSup<<"[\n";
        return dim;
    }
    
    
    /// free all occupied memory
    void deallocate()
    {
        //printf("Lattice %p realeased\n", this);
        if ( laSite0 )
            delete[] laSite0;
        laSite0 = nullptr;
        laSite  = nullptr;
        laInf   = 0;
        laSup   = 0;
    }
    
    /// forbid access to data via automatic type conversion:
    template<typename T> cell_t   value(const T s) const;
    template<typename T> cell_t&  site(const T s);
    template<typename T> cell_t&  operator[](const T s);
    
    /// disable assignment operator
    Lattice & operator =(const Lattice &);
    
    //------------------------------------------------------------------------------
    #pragma mark - Construction
    
public:
    
    /// Constructor
    Lattice()
    {
        laInf   = 0;
        laSup   = 0;
        laUnit  = 0.0;
        laSite0 = nullptr;
        laSite  = nullptr;
    }
    
    /// Copy constructor
    Lattice(const Lattice & lat)
    {
        std::clog << this << " Lattice::Lattice [" << laInf << ", " << laSup << "[\n";
        
        //copy member variables:
        memcpy(this, &lat, sizeof(Lattice));
        
        // allocate memory:
        laSite0 = new cell_t[laSup-laInf];
        laSite  = laSite0 - laInf;
        
        // copy data:
        for ( site_t s = laInf; s < laSup; ++s )
            laSite[s] = lat.laSite[s];
    }
    
    
    /// Destructor
    ~Lattice()    { deallocate(); }
    
    
    /// set distance betwen adjacent sites (the size of a site)
    void setUnit(real u)
    {
        if ( u < REAL_EPSILON )
            throw InvalidParameter("lattice:unit must be > 0");
        laUnit = u;
    }

    /// change distance betwen adjacent sites
    void changeUnit(real u)
    {
        if ( u < REAL_EPSILON )
            throw InvalidParameter("lattice:unit must be > 0");
        if ( laUnit != u )
        {
            std::clog << "Lattice::unit " << laUnit << " -> " << u << "\n";
            //preserve the same abscissa range, with the new lattice unit
            real s = laUnit * laInf;
            real e = laUnit * laSup + laUnit;
            laUnit = u;
            if ( laSite )
                setRange(s, e);
        }
    }
    
    /// set the range of valid abscissa
    void setRange(real i, real s)
    {
        assert_true( laUnit > REAL_EPSILON );
        //std::clog << this << " Lattice::setRange(" << i << ", " << s << ") " << laUnit << "\n";

        ///\todo use special Lattice value to indicate the ends of the fibers
        /* add here some safety margin */
        allocateLattice(index(i), index(s)+1, 8);
    }
    
    //------------------------------------------------------------------------------
    #pragma mark - Index / Abscissa
    
    /// true if lattice unit size was set
    bool    ready() const { return laUnit > REAL_EPSILON; }

    /// distance between adjacent sites
    real    unit() const { return laUnit; }
    
    /// index of the site containing abscissa `a`
    site_t  index(const real a) const { return (site_t)floor(a/laUnit); }
    
    /// index of the site after the one containing abscissa `a`
    site_t  index_sup(const real a) const { return (site_t)ceil(a/laUnit); }

    /// true if index 'i' is covered by the lattice
    bool    valid(const site_t i) const { return ( laInf <= i  &&  i < laSup ); }
    
    /// true if index 'i' is not covered by the lattice
    bool    invalid(const site_t i) const { return ( i < laInf  || laSup <= i ); }
    
    
    /// the site of index `h` covers the abscissa range `unit * h < s < unit * ( h + 1 )`
    /**
     abscissa(h) returns the abscissa corresponding to the beginning of the site.
     The range covered by site 'h' is [ abscissa(h), abscissa(h+1) ]
     */
    real    abscissa(const real s) const { return s * laUnit; }
    
    /// true if abscissa `a` corresponds to a site that is covered by the lattice
    bool    within(const real a) const { return laInf * laUnit <= a  &&  a <= laSup * laUnit; }

    //------------------------------------------------------------------------------
    #pragma mark - Access

    /// first valid index
    site_t  inf()  const  { return laInf; }
    
    /// last valid index plus one
    site_t  sup()  const  { return laSup; }

    /// pointer to data array
    cell_t* data() const  { return laSite; }

    /// value at index `s`, equivalent to []
    cell_t& data(const site_t& s)     const { assert_true(valid(s)); return laSite[s]; }
    
    /// reference to Site at index s
    cell_t& operator[](const site_t& s)     { assert_true(valid(s)); return laSite[s]; }
    
    /// value at abscissa `a`, with convertion to site index, unlike operator []
    cell_t& cell(real a)             { site_t s=index(a); assert_true(valid(s)); return laSite[s]; }
    
    /// value at abscissa `a`, with convertion to site index, unlike operator []
    cell_t const& cell(real a) const { site_t s=index(a); assert_true(valid(s)); return laSite[s]; }

    /// set all sites to `value`
    void clear(cell_t value = 0)
    {
        for ( site_t s = laInf; s < laSup; ++s )
            laSite[s] = value;
    }
    
    //------------------------------------------------------------------------------
    #pragma mark - Transfer
    
    /// transfer lat->sites within `[s, e[` to *this, and set lat->sites to zero
    void take(Lattice<CELL> & lat, site_t is, site_t ie)
    {
        // check all boundaries
        if ( is < laInf )     is = laInf;
        if ( is < lat.laInf ) is = lat.laInf;
        
        if ( ie > laSup )     ie = laSup;
        if ( ie > lat.laSup ) ie = lat.laSup;
        
        CELL * s = laSite + is;
        CELL * e = laSite + ie;
        CELL * o = lat.laSite + is;

        // transfer values
        while ( s < e )
        {
            *s = *o;
            *o = 0;
            ++s;
            ++o;
        }
    }
    
    
    /// sum all sites in [inf, e[; set sites to zero and return total in `res`
    void takeM(Lattice<CELL> & lat, const site_t e)
    {
        take(lat, laInf, e);
    }
    
    /// sum all sites in [s, sup]; set sites to zero and return total in `res`
    void takeP(Lattice<CELL> & lat, const site_t s)
    {
        take(lat, s, laSup);
    }
    
    //------------------------------------------------------------------------------
    #pragma mark - Collect

    /// sum all sites in `[s, e[`, set them to zero and return total in `res`
    template <typename SUM>
    void collectR(SUM& res, site_t is, site_t ie)
    {
        res = 0;
#if ( 1 )
        // check boundaries
        if ( is < laInf )
            is = laInf;
        
        if ( ie > laSup )
            ie = laSup;
        
        // collect
        for ( ; is < ie; ++is )
        {
            res += laSite[is];
            laSite[is] = 0;
        }
#else
        // check boundaries
        CELL * s = laSite + std::max(is, laInf);
        CELL * e = laSite + std::min(ie, laSup);
        
        // collect
        while ( s < e )
        {
            res += *s;
            *s = 0;
            ++s;
        }
#endif
    }
    
    
    /// sum all sites in [inf, e[; set sites to zero and return total in `res`
    template <typename SUM>
    void collectM(SUM& res, const site_t e)
    {
        collectR(res, laInf, e);
    }
    
    /// sum all sites in ]s, sup]; set sites to zero and return total in `res`
    template <typename SUM>
    void collectP(SUM& res, const site_t s)
    {
        collectR(res, s+1, laSup);
    }
    
    /// sum of values in the entire lattice; set sites to zero
    template <typename SUM>
    void collect(SUM& res) const
    {
        collectR(res, laInf, laSup);
    }
    
    //------------------------------------------------------------------------------
    #pragma mark - Sums

    /// sum of values for all sites in `[s, e[` and return total in `res`
    template <typename SUM>
    void sum(SUM& res, site_t is, site_t ie) const
    {
        // check boundaries
        CELL * s = laSite + std::max(is, laInf);
        CELL * e = laSite + std::min(ie, laSup);
        
        // collect
        res = 0;
        while ( s < e )
        {
            res += *s;
            ++s;
        }
    }
    
    /// sum of values for all sites below `e` (not-included)
    template <typename SUM>
    void sumM(SUM& res, const site_t e) const
    {
        sum(res, laInf, e);
    }

    /// sum of values for all sites above `s` (included)
    template <typename SUM>
    void sumP(SUM& res, const site_t s) const
    {
        sum(res, s, laSup);
    }

    
    /// sum of values in the entire lattice
    template <typename SUM>
    void sum(SUM& res) const
    {
        sum(res, laInf, laSup);
    }

    /// sum of values in the entire lattice using 'real' as accumulator
    real sum() const
    {
        real res;
        sum(res, laInf, laSup);
        return res;
    }

    //------------------------------------------------------------------------------
    #pragma mark - I/O
    
    /// write data within [inf, sup[ to file
    void write_data(Outputter& out, site_t inf, site_t sup) const;
    
    /// write to file
    void write(Outputter& out, site_t inf, site_t sup) const
    {
        if ( sup < inf )
            throw InvalidIO("cannot write incoherent Lattice boundaries");
        if ( laUnit < REAL_EPSILON )
            throw InvalidIO("cannot write incoherent Lattice unit value");

        out.writeInt32(inf);
        out.writeInt32(sup);
        out.writeFloat(laUnit);
        
        if ( laSite )
            write_data(out, inf, sup);
        else
            out.writeUInt32(0);
    }
    
    /// write all data to file
    void write(Outputter& out) const
    {
        write(out, laInf, laSup);
    }

    /// clear all cells and read data from file 
    void read(Inputter& in)
    {
        site_t inf = in.readInt32();
        site_t sup = in.readInt32();
        real uni = in.readFloat();
        in.readUInt16();
        in.readUInt8();
        int nbytes = in.readUInt8();
        
        if ( sup < inf )
            throw InvalidIO("incoherent Lattice boundaries");
        if ( uni < REAL_EPSILON )
            throw InvalidIO("incoherent Lattice unit value");
        
        //only change unit length if the difference is significant
        if ( fabs( laUnit - uni ) > 1e-6 )
            changeUnit(uni);
        
        allocateLattice(inf, sup, 0);
        clear();
        
        if ( nbytes == 1 )
        {
            for ( site_t s = inf; s < sup; ++s )
                laSite[s] = in.readUInt8();
        }
        else if ( nbytes == 2 )
        {
            for ( site_t s = inf; s < sup; ++s )
                laSite[s] = in.readUInt16();
        }
        else if ( nbytes == 4 )
        {
            for ( site_t s = inf; s < sup; ++s )
                laSite[s] = in.readUInt32();
        }
        else
        {
            for ( site_t s = inf; s < sup; ++s )
                laSite[s] = in.readFloat();
        }
    }

    
    /// debug function
    int bad()
    {
        if ( laSite ) {
            if ( laInf > laSup )
                return 1;
        }
        return 0;
    }
    
};

#endif
