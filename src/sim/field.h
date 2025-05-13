// Cytosim was created by Francois Nedelec. Copyright 2020 Cambridge University.

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


/// the type of data contained in a Field
//typedef FieldVector<3> FieldCell;
typedef FieldScalar FieldCell;

/// the type of the Field's underlying Grid
typedef Grid<FieldCell, DIM> FieldGrid;


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
    
    /// Grid data object
    FieldGrid mGrid;
    
    /// property
    FieldProp const* prop;
    
private:
    
    /// disabled default constructor
    Field();
    
    /// duplicate field
    real * fiTMP;
    
    /// allocated size of fiTMP
    size_t fiTMPalc;
    
    /// matrix for diffusion
    SparMatSym1 fiDiffusionMatrix;
    
    /// initialize to cover the given Space with squares of size 'step'
    void setGrid(Vector inf, Vector sup, real step, bool tight)
    {
        assert_true( step > REAL_EPSILON );
        // we add a safety border (in micro-meters)
        const real extra = tight ? 0 : 1;
        
        index_t n_cell[3] = { 0, 0, 0 };
        // we use square voxels:
        for ( int d = 0; d < DIM; ++d )
        {
            n_cell[d] = (index_t)std::ceil( (sup[d]-inf[d]+extra) / step );
            real mid = 0.5 * ( inf[d] + sup[d] );
            inf[d] = mid - 0.5 * step * n_cell[d];
            sup[d] = mid + 0.5 * step * n_cell[d];
        }
        
        mGrid.setDimensions(inf, sup, n_cell);
        
        //verify the cell size:
        for ( int d = 0; d < DIM; ++d )
        {
            real dif = abs_real( step - mGrid.cellWidth(d) );
            if ( abs_real(dif) > 1e-3 )
            {
                Cytosim::warn("Field:step[", d, "] is not as expected:\n",
                              "  field: ", mGrid.cellWidth(d), "  prop: ", step, "\n");
            }
        }
    }
    
    /// allocate memory for the scalar field (setGrid() must be called before)
    void createCells()
    {
        assert_true( mGrid.hasDimensions() );
        // create the grid using the calculated dimensions:
        if ( mGrid.createCells() )
        {
            // set all values to zero (already done in the constructor of FieldValue)
            // mGrid.clear();
            mGrid.printSummary(Cytosim::log, "   Field");
        }
    }
    
public:
#pragma mark -
    
    /// constructor
    Field(FieldProp const* p)
    {
        prop = p;
        fiTMP = nullptr;
        fiTMPalc = 0;
    }
    
    /// destructor
    ~Field()
    {
        free_real(fiTMP);
        prop = nullptr;
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
                for ( int d = 0; d < DIM; ++d )
                    mGrid.setPeriodic(d, true);
                setGrid(inf, sup, prop->step, true);
            }
            else
            {
                setGrid(inf, sup, prop->step, false);
            }
            createCells();
            //std::clog << "Field step "<<prop->step<<" has "<<mGrid.nbCells()<<" cells\n";
            Cytosim::log.print("Field %lx set with %i cells of size %.3f um\n", this, mGrid.nbCells(), prop->step);
        }
    }
    
    
    /// true if field was set
    size_t hasField() const { return mGrid.hasCells(); }
    
    /// size of cell
    real cellWidth() const { return mGrid.cellWidth(0); }
    
    /// volume of cell
    real cellVolume() const { return mGrid.cellVolume(); }
    
    /// access to data
    FieldCell& cell(const real w[]) const { return mGrid.cell(w); }
    
    /// access to data
    index_t nbCells() const { return mGrid.nbCells(); }

    /// info
    void infoValues(FieldCell& s, FieldCell& a, FieldCell& n, FieldCell& x) const { return mGrid.infoValues(s, a, n, x); }
    
    //------------------------------ simulation --------------------------------
#pragma mark -
    
    /// set all cells to value = volume * conc
    void setConcentration(FieldCell conc)
    {
#if 0
        mGrid.setValues( mGrid.cellVolume() * conc );
#else
        FieldCell val = conc * mGrid.cellVolume();
        for ( index_t i = 0; i < mGrid.nbCells(); ++i )
            mGrid[i] = val * RNG.exponential();
#endif
    }
    
    
    /// set cells that are inside `spc` to value = volume * conc
    void setConcentration(Space const* spc, FieldCell val_in, FieldCell val_out)
    {
        FieldCell i = val_in * mGrid.cellVolume();
        FieldCell o = val_out * mGrid.cellVolume();
        
        for ( index_t c = 0; c < mGrid.nbCells(); ++c )
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
    static void laplacian(FieldGrid const&, const real*, real*, bool);
    
    /// implements diffusion of the whole field in direction X
    static void diffuseX(FieldGrid const&, real*, real, real, bool);
    /// implements diffusion of the whole field in direction Y
    static void diffuseY(FieldGrid const&, real*, real, real, bool);
    /// implements diffusion of the whole field in direction Y
    static void diffuseZ(FieldGrid const&, real*, real, real, bool);

    /// set values of field on its edges
    static void setEdgesX(FieldGrid const&, real*, real);
    /// set values of field on its edges
    static void setEdgesY(FieldGrid const&, real*, real);
    /// set values of field on its edges
    static void setEdgesZ(FieldGrid const&, real*, real);
    
    /// initialize diffusion matrix (only for FieldScalar)
    void prepareDiffusion(real);
    
    /// initialize diffusion matrix (only for FieldScalar)
    void prepareDiffusion(real, unsigned char *);
    
    //------------------------------- object -----------------------------------
#pragma mark -
    
    /// a unique character identifying the class
    static const ObjectTag TAG = 'i';
    
    /// return unique character identifying the class
    ObjectTag tag() const { return TAG; }
    
    /// return index of 'prop' in corresponding PropertyList
    Property const* property() const { return prop; }

    //--------------------------------------------------------------------------

    /// a static_cast<> of Object::next()
    Field* next() const { return static_cast<Field*>(next_); }
    
    /// a static_cast<> of Object::prev()
    Field* prev() const { return static_cast<Field*>(prev_); }
    
    //------------------------------ read/write --------------------------------
#pragma mark -
    
    /// write Field to file using VAL::write()
    /** Some of this should be moved to Grid */
    void write(Outputter& o) const
    {
        if ( mGrid.hasCells() && prop->save )
        {
            writeMarker(o, Field::TAG);
            o.writeUInt16(DIM);
            for ( int d = 0; d < DIM; ++d )
            {
                o.writeUInt32(mGrid.breadth(d), ' ');
                o.writeFloat(mGrid.inf(d));
                o.writeFloat(mGrid.sup(d));
            }
            o.writeUInt32(mGrid.nbCells(), ' ');
            for ( index_t c = 0; c < mGrid.nbCells(); ++c )
                mGrid.icell(c).write(o);
            o.writeSoftNewline();
        }
        
        if ( prop->positive )
        {
            if ( mGrid.hasNegativeValue() )
                throw Exception("Aborting because Field has negative values");
        }
    }
    
    
    /// read Field from file using VAL::read()
    void readData(Inputter& in, Simul&)
    {
        index_t size[DIM] = { 0 };
        real minB[DIM] = { 0 };
        real maxB[DIM] = { 0 };
        
        index_t dim = in.readUInt16();
        if ( dim != DIM )
            throw InvalidIO("cannot read field due to dimensionality mismatch");
        
        for ( index_t d = 0; d < dim; ++d )
        {
            size[d] = in.readUInt32();
            minB[d] = in.readFloat();
            maxB[d] = in.readFloat();
        }
        
        mGrid.setDimensions(minB, maxB, size);
        createCells();
        
        index_t nbc = in.readUInt32();
        if ( nbc != mGrid.nbCells() )
        {
            std::cerr << "file: " << nbc << " field: " << mGrid.nbCells() << "\n";
            throw InvalidIO("mismatch in Field::size");
        }
        //std::clog << "readData() num_cells=" << nbc << '\n';
        
        for ( index_t c = 0; c < nbc; ++c )
            mGrid.icell(c).read(in);
    }
    
    /// read Field and checks that the Grid::step has not changed
    void read(Inputter& in, Simul& sim, ObjectTag)
    {
        readData(in, sim);
        
        if ( prop )
        {
            for ( int d = 0; d < DIM; ++d )
            {
                real dif = abs_real( prop->step - mGrid.cellWidth(d) );
                if ( abs_real(dif) > 1e-3 )
                {
                    Cytosim::warn("Field:step[", d, "] has changed:\n",
                                  "  file: ", mGrid.cellWidth(d), " prop: ", prop->step, "\n");
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
