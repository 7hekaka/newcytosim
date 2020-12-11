// Cytosim was created by Francois Nedelec.  Copyright 2020 Cambridge University.

#ifndef FIELD_H
#define FIELD_H

#include "dim.h"
#include "real.h"
#include "grid.h"
#include "space.h"
#include "object.h"
#include "iowrapper.h"
#include "messages.h"
#include "exceptions.h"
#include "sparmatsym1.h"
#include "sparmatsym2.h"
#include "field_prop.h"
#include "field_values.h"

class FiberSet;


/// the type of Grid contained in a Field
typedef Grid<FieldScalar, DIM> FieldGrid;


/// value of type VAL defined as a function of position over the simulation Space
/**
 A field represents a value that is varying with position over space.
 It does so by storing a different value for each cell in a regular grid.
 
 Each cell holds the amount of molecules.
 The local concentration can be obtained by dividing by the cell volume:

     concentration = mGrid.cell(position) / mGrid.cellVolume();
 
 Note that the field is build on a Grid with square cells, because diffusion/reaction are then
 easier to implement. The Grid does not necessarily match the edges of the Space exactly,
 but instead extends outside, such as to cover the 'inside' region entirely.
 
 */
class Field : public Object
{
public:

    /// forward the type of value
    typedef FieldGrid::value_type value_type;
    
    /// Grid data object
    FieldGrid  mGrid;
    
    /// property
    FieldProp const* prop;
    
private:
    
    /// disabled default constructor
    Field();
    
    /// duplicate field
    real*    fiTMP;
    
    /// allocated size of fiTMP
    size_t   fiTMPSize;
    
    /// matrix for diffusion
    SparMatSym1 fiDiffusionMatrix;
    
    /// initialize to cover the given Space with squares of size 'step'
    void setGrid(Vector inf, Vector sup, real step, bool tight)
    {
        assert_true( step > REAL_EPSILON );
        // we add a safety border (in micro-meters)
        const real extra = tight ? 0 : 1;
        
        size_t n_cell[3] = { 0, 0, 0 };
        // we use square cells:
        for ( size_t d = 0; d < DIM; ++d )
        {
            n_cell[d] = (size_t)std::ceil( (sup[d]-inf[d]+extra) / step );
            real mid = 0.5 * ( inf[d] + sup[d] );
            inf[d] = mid - 0.5 * step * n_cell[d];
            sup[d] = mid + 0.5 * step * n_cell[d];
        }
        
        mGrid.setDimensions(inf, sup, n_cell);
        
        //verify the cell size:
        for ( size_t d = 0; d < DIM; ++d )
        {
            real dif = abs_real( step - mGrid.cellWidth(d) );
            if ( abs_real(dif) > 1e-3 )
            {
                Cytosim::warn << "Field:step[" << d << "] is not as expected:\n";
                Cytosim::warn << "  field: " << mGrid.cellWidth(d) << "  prop: " << step << "\n";
            }
        }
    }
    
    /// allocate memory for the scalar field (setGrid() must be called before)
    void createCells()
    {
        assert_true( mGrid.hasDimensions() );
        // delete preexisting grid if necessary:
        mGrid.destroy();
        // create the grid using the calculated dimensions:
        mGrid.createCells();
        // set all values to zero (already done in the constructor of FieldValue)
        // mGrid.clear();
        mGrid.printSummary(Cytosim::log, "Field");
    }
    
public:
#pragma mark -
    
    /// constructor
    Field(FieldProp const* p)
    {
        prop      = p;
        fiTMP     = nullptr;
        fiTMPSize = 0;
    }
    
    /// destructor
    ~Field()
    {
        free_real(fiTMP);
    }
    
    /// initialize with squares of size 'step'
    void setField()
    {
        assert_true( prop );
        
        if ( ! mGrid.hasCells() )
        {
            if ( !prop->field_space_ptr )
                throw InvalidParameter("field:space must be defined");
            
            Vector inf, sup;
            prop->field_space_ptr->boundaries(inf, sup);
            
            if ( prop->field_periodic )
            {
                for ( size_t d = 0; d < DIM; ++d )
                    mGrid.setPeriodic(d, true);
                setGrid(inf, sup, prop->step, true);
            }
            else
            {
                setGrid(inf, sup, prop->step, false);
            }
            createCells();
            //std::clog << "Field step "<<prop->step<<" has "<<mGrid.nbCells()<<" cells\n";
            Cytosim::log("Field %lx set with %i cells of size %.3f um\n", this, mGrid.nbCells(), prop->step);
        }
    }
    
    
    /// true if field was set
    size_t hasField() const { return mGrid.hasCells(); }
    
    /// size of cell
    real cellWidth() const { return mGrid.cellWidth(0); }
    
    /// volume of cell
    real cellVolume() const { return mGrid.cellVolume(); }
    
    /// access to data
    value_type& cell(const real w[]) const { return mGrid.cell(w); }
    
    /// access to data
    size_t nbCells() const { return mGrid.nbCells(); }

    /// info
    void infoValues(value_type& s, value_type& n, value_type& x) const { return mGrid.infoValues(s, n, x); }
    
    //------------------------------ simulation --------------------------------
#pragma mark -
    
    /// set all cells to value = volume * conc
    void setConcentration(FieldGrid::value_type conc)
    {
#if 0
        mGrid.setValues( conc * mGrid.cellVolume() );
#else
        FieldGrid::value_type val = conc * mGrid.cellVolume();
        for ( size_t i = 0; i < mGrid.nbCells(); ++i )
            mGrid[i] = val * RNG.exponential();
#endif
    }
    
    
    /// set cells that are inside `spc` to value = volume * conc
    void setConcentration(Space const* spc, FieldGrid::value_type in, FieldGrid::value_type ou)
    {
        real i = in * mGrid.cellVolume();
        real o = ou * mGrid.cellVolume();
        
        for ( size_t c = 0; c < mGrid.nbCells(); ++c )
        {
            Vector w;
            mGrid.setPositionFromIndex(w, c, 0.5);
            if ( spc->inside(w) )
                mGrid.icell(c) = i;
            else
                mGrid.icell(c) = o;
        }
    }
    
    /// initialize Field
    void prepare();

    /// simulation step
    void step(FiberSet&);
    
    /// calculate second derivative of field
    void laplacian(const real*, real*) const;
    
    /// calculate second derivative of field
    void diffuseX(real*, real);
    
    /// set values of field on its edges
    void setEdgesX(real*, real);
    
    /// set values of field on its edges
    void setEdgesY(real*, real);
    
    /// set values of field on its edges
    void setEdgesZ(real*, real);
    
    /// initialize diffusion matrix (only for FieldScalar)
    void prepareDiffusion(real);
    
    /// initialize diffusion matrix (only for FieldScalar)
    void prepareDiffusion(real, unsigned char *);
    
    //------------------------------- object -----------------------------------
#pragma mark -
    
    /// a unique character identifying the class
    static const ObjectTag TAG = 'i';
    
    /// return unique character identifying the class
    ObjectTag       tag() const { return TAG; }
    
    /// return index of 'prop' in corresponding PropertyList
    Property const* property() const { return prop; }
    
    /// a static_cast<> of Object::next()
    Field* next() const { return static_cast<Field*>(nextO); }
    
    /// a static_cast<> of Object::prev()
    Field* prev() const { return static_cast<Field*>(prevO); }
    
    //------------------------------ read/write --------------------------------
#pragma mark -
    
    /// print total, minimum and maximum value
    void   writeInfo(std::ostream& out) const
    {
        real vol = mGrid.cellVolume();
        FieldGrid::value_type sum, mn, mx;
        mGrid.infoValues(sum, mn, mx);
        out << prop->name() << " sum " << sum << " min " << mn/vol << " max " << mx/vol << std::endl;
    }
    
    /// write Field to file using VAL::write()
    /** Some of this should be moved to Grid */
    void   write(Outputter& out) const
    {
        if ( mGrid.hasCells() && prop->save )
        {
            out.writeUInt16(DIM);
            for ( size_t d = 0; d < DIM; ++d )
            {
                out.writeSoftSpace();
                out.writeUInt32(mGrid.breadth(d));
                out.writeFloat(mGrid.inf(d));
                out.writeFloat(mGrid.sup(d));
            }
            out.writeSoftSpace();
            out.writeUInt32(mGrid.nbCells());
            for ( size_t c = 0; c < mGrid.nbCells(); ++c )
                mGrid.icell(c).write(out);
            out.writeSoftNewline();
        }
        
        if ( prop->positive )
        {
            if ( mGrid.hasNegativeValue() )
                throw Exception("Aborting because Field has negative values");
        }
    }
    
    
    /// read Field from file using VAL::read()
    void   readData(Inputter& in, Simul&)
    {
        size_t size[DIM] = { 0 };
        real   minB[DIM] = { 0 }, maxB[DIM] = { 0 };
        
        size_t dim = in.readUInt16();
        if ( dim != DIM )
            throw InvalidIO("cannot read field due to dimensionality mismatch");
        
        for ( size_t d = 0; d < dim; ++d )
        {
            size[d] = in.readUInt32();
            minB[d] = in.readFloat();
            maxB[d] = in.readFloat();
        }
        
        mGrid.setDimensions(minB, maxB, size);
        createCells();
        
        size_t nbc = in.readUInt32();
        if ( nbc != mGrid.nbCells() )
        {
            printf("file: %lu field:%lu\n", nbc, mGrid.nbCells());
            throw InvalidIO("mismatch in Field::size");
        }
        //std::clog << "readData() nb_cells=" << nbc << '\n';
        
        for ( size_t c = 0; c < nbc; ++c )
            mGrid.icell(c).read(in);
    }
    
    /// read Field and checks that the Grid::step has not changed
    void   read(Inputter& in, Simul& sim, ObjectTag)
    {
        readData(in, sim);
        
        if ( prop )
        {
            for ( size_t d = 0; d < DIM; ++d )
            {
                real dif = abs_real( prop->step - mGrid.cellWidth(d) );
                if ( abs_real(dif) > 1e-3 )
                {
                    Cytosim::warn << "Field:step["<<d<<"] has changed:\n";
                    Cytosim::warn << "  file: " << mGrid.cellWidth(d) << " prop: " << prop->step << "\n";
                }
            }
            
            /*
             we should extrapolate the data that were read to a grid with the
             resolution specified by prop->step
             */
        }
    }
    
    /// OpenGL display
    void draw() const;
    
    /// OpenGL display
    void draw(Space const*, Vector3 const& dir, const real pos) const;
    
};


#endif
