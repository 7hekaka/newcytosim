// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "digit.h"
#include "digit_prop.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "messages.h"
#include "glossary.h"
#include "lattice.h"
#include "simul.h"


Digit::Digit(DigitProp const* p, HandMonitor* h)
: Hand(p,h), prop(p)
{
}


bool Digit::attachmentAllowed(FiberSite & sit) const
{
    if ( Hand::attachmentAllowed(sit) )
    {
        sit.engageLattice();
        if ( ! sit.betweenMP() )
            return false;
        return sit.vacant();
    }
    return false;
}


void Digit::attach(FiberSite const& sit)
{
    assert_true(sit.vacant());
    Hand::attach(sit);
    inc();

    if ( lattice()->unit() != prop->step_size  )
    {
        /*
         This is not necessarily a bug, but we would need to write additional code,
         to be able to handle a step_size that would be a multiple of fiber:lattice_unit.
         */
        throw InvalidParameter("digit:step_size must be defined and equal to fiber:lattice_unit");
    }
}


void Digit::detach()
{
	/*
	 We need to workaround the fact that when reading file, the lattice
	 values can become out of sync with the states of the digits on the fiber.
	 @todo: find a better fix for keeping Digit/Lattice sync
	*/
	if ( haMonitor->parent_flag() != 7 )
	{
		dec();
		assert_true(vacant());
	}
    Hand::detach();
}

//------------------------------------------------------------------------------

/**
 This will attempt to relocate to site `npos`,
 without checking the intermediate sites.

@retval 0 : the move was done
@retval 1 : the move aborted, because the specified destination is occupied
 
 */
int Digit::jumpTo(const site_t npos)
{
    assert_true( npos != site() );
    if ( vacant(npos) )
    {
        hop(npos);
        return 0;
    }
    return 1;
}


int Digit::jumpToEndM()
{
    site_t s = lattice()->index(fiber()->abscissaM()) + 1;
    return jumpTo(s);
}


int Digit::jumpToEndP()
{
    site_t s = lattice()->index(fiber()->abscissaP());
    return jumpTo(s);
}


//------------------------------------------------------------------------------

/**
 Try to jump `n` sites in the PLUS_END direction,
 without checking the intermediate positions.

 @retval 0 : the move was done
 @retval 1 : the move aborted, because the destination is occupied
 @retval 2 : the move aborted, because the destination is outside the Fiber
*/
int Digit::jumpP(const int s)
{
    assert_true( s > 0 );
    assert_true( attached() );
    
    site_t npos = site() + s;

    /*
     If lattice is used, we could have a special value to indicate
     ranges outside of the fibers.
     This would be easier to check than the abscissaM() and abscissaP()
     */
    if ( lattice()->abscissa(npos) > fiber()->abscissaP() )
        return 2;
    
    return jumpTo(npos);
}

/**
 Try to jump to `n` sites in the MINUS_END direction,
 without checking the intermediate positions.

 @retval 0 : the move was done
 @retval 1 : the move aborted, because the destination is occupied
 @retval 2 : the move aborted, because the destination is outside the Fiber
 */
int Digit::jumpM(const int s)
{
    assert_true( s > 0 );
    assert_true( attached() );

    site_t npos = site() - s;
    
    /*
     If lattice is used, we could have a special value to indicate 
     ranges outside of the fibers.
     This would be easier to check than the abscissaM() and abscissaP()
     */
    if ( lattice()->abscissa(npos+1) < fiber()->abscissaM() )
        return 2;
        
    return jumpTo(npos);
}

//------------------------------------------------------------------------------

/**
 Try to move to the adjacent site in the PLUS_END direction.
 
 @retval 0 : the move was done
 @retval 1 : the move aborted, because the destination is occupied
 @retval 2 : the move aborted, because the destination is outside the Fiber
 */
int Digit::stepP()
{
    assert_true( attached() );
    site_t npos = site() + 1;
    
    /*
     If lattice is used, we could have a special value to indicate
     ranges outside of the fibers.
     This would be easier to check than the abscissaM() and abscissaP()
     */
    if ( lattice()->abscissa(npos) >= fiber()->abscissaP() )
        return 2;
    
    return jumpTo(npos);
}

/**
 Try to move to the adjacent site in the MINUS_END direction.
 
 @retval 0 : the move was done
 @retval 1 : the move aborted, because the destination is occupied
 @retval 2 : the move aborted, because the destination is outside the Fiber
 */
int Digit::stepM()
{
    assert_true( attached() );
    site_t npos = site() - 1;
    
    /*
     If lattice is used, we could have a special value to indicate
     ranges outside of the fibers.
     This would be easier to check than the abscissaM() and abscissaP()
     */
    if ( lattice()->abscissa(npos+1) <= fiber()->abscissaM() )
        return 2;
    
    return jumpTo(npos);
}

//------------------------------------------------------------------------------

/**
 Try to move `n` sites in the PLUS_END direction,
 stopping if any intermediate position is already occupied.
 
 @retval 0 : the specified number of steps was done
 @retval 1 : the move stopped at an occupied intermediate site
 @retval 2 : the end of the Fiber was reached
 
 For 1 and 2, the Digit is relocated to the site just before the obstacle.
 */
int Digit::crawlP(const int n)
{
    assert_true( n > 0 );
    assert_true( attached() );
    
    int res = 0, s = 0;
    site_t npos = site();
    
    while ( s < n )
    {
        ++npos;
        
        if ( ! vacant(npos) )
        {
            res = 1;
            break;
        }
        /*
         We could define a special value of Lattice to indicate
         ranges outside of the fibers.
         */
        if ( lattice()->abscissa(npos) > fiber()->abscissaP() )
        {
            res = 2;
            break;
        }
        
        ++s;
    }
    
    if ( s )
    {
        assert_true( vacant(site()+s) );
        hop(site()+s);
    }
    return res;
}


/**
 Try to move `n` sites in the MINUS_END direction,
 stopping if any intermediate position is already occupied.
 
 @retval 0 : the specified number of steps was done
 @retval 1 : the move stopped at an occupied intermediate site
 @retval 2 : the end of the Fiber was reached
 
 For 1 and 2, the Digit is relocated to the site just before the obstacle.
 */
int Digit::crawlM(const int n)
{
    assert_true( n > 0 );
    assert_true( attached() );
    
    int res = 0, s = 0;
    site_t npos = site();
    
    while ( s < n )
    {
        --npos;
        
        if ( ! vacant(npos) )
        {
            res = 1;
            break;
        }
        /*
         We could define a special value of Lattice to indicate
         ranges outside of the fibers.
         */
        if ( lattice()->abscissa(npos+1) < fiber()->abscissaM() )
        {
            res = 2;
            break;
        }
        
        ++s;
    }
    
    if ( s )
    {
        assert_true( vacant(site()-s) );
        hop(site()-s);
    }
    return res;
}

//------------------------------------------------------------------------------


/**
 The Digit normally does not move by itself
 */
void Digit::handleDisassemblyM()
{
    assert_true( attached() );
    
    if ( prop->hold_shrinking_end )
    {
        if ( jumpToEndM() )
            detach();
    }
    else
        detach();
}


/**
 The Digit normally does not move by itself
 */
void Digit::handleDisassemblyP()
{
    assert_true( attached() );
    
    if ( prop->hold_shrinking_end )
    {
        if ( jumpToEndP() )
            detach();
    }
    else
        detach();
}


//------------------------------------------------------------------------------
/**
tests detachment
 */
void Digit::stepUnloaded()
{
    assert_true( attached() );

    testDetachment();
}


/**
(see @ref Stochastic)
 */
void Digit::stepLoaded(Vector const& force, real force_norm)
{
    assert_true( attached() );
    assert_true( nextDetach >= 0 );

    testKramersDetachment(force_norm);
}

