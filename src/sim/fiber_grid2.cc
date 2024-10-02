// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

/** 
 This file implements a 'dummy' grid using STL code, which can be used as a reference.
 For each position, it calculates the geometrical distance to all fiber segments.
 This is algorithmically the slowest method, but it is simple and most likely correct!
 It is useful to get a ground truth and evaluate more advanced methods.
 */

/// a list containing all segments, as a global variable
FiberGrid::SegmentList allSegments;


index_t FiberGrid::setGrid(Vector, Vector, real)
{
    LOG_ONCE("Cytosim is using a crude method to localize fibers!\n");
    return 1;
}


void FiberGrid::createCells()
{
}


index_t FiberGrid::nbCells() const
{
    return 1;
}


size_t FiberGrid::hasGrid() const
{
    return 1;
}


index_t FiberGrid::nbTargets() const
{
    return allSegments.size();
}

//------------------------------------------------------------------------------
#pragma mark - Paint

void FiberGrid::paintGrid(const Fiber * first, const Fiber * last, real)
{
    allSegments.clear();
    // add all segments
    for ( const Fiber * f = first ; f != last ; f=f->next() )
    {
        for ( index_t s = 0; s < f->nbSegments(); ++s )
            allSegments.emplace(f, s);
    }
}


void FiberGrid::tryToAttach(Vector const& place, Hand& ha) const
{
    // randomize the list order
    //std::random_shuffle(allSegments.begin(), allSegments.end());

    // test all segments:
    const uint64_t prob = 0x1p+32 * ha.property()->binding_prob;
    const real sup = square(ha.property()->binding_range);
    for ( FiberSegment const& seg : allSegments )
    {
        if ( RNG.pint32() < prob )
        {
            // Compute the distance between 'place' and segment
            real dis = INFINITY;
            real abs = seg.projectPoint(place, dis);
            
            /*
             Compare to the maximum attachment range of the hand,
             and compare a newly tossed random number with 'prob'
             */
            if ( dis < sup )
            {
                Fiber const* fib = seg.fiber();
                FiberSite sit(fib, seg.abscissa1()+abs);
                
                if ( ha.keyMatch(fib) && ha.attachmentAllowed(sit) )
                {
                    ha.attach(sit);
                    return;
                }
            }
        }
    }
}


FiberGrid::SegmentList FiberGrid::nearbySegments(Vector const& vec) const
{
    return allSegments;
}


FiberGrid::SegmentList FiberGrid::nearbySegments(Vector const& place, const real DD, Fiber const* exclude) const
{
    SegmentList res;
    
    for ( FiberSegment const& seg : allSegments )
    {
        if ( seg.fiber() != exclude )
        {
            // Compute the distance between 'place' and segment
            real dis = INFINITY;
            seg.projectPoint(place, dis);
            
            if ( dis < DD )
                res.push_back(seg);
        }
    }
    
    return res;
}


FiberSegment FiberGrid::closestSegment(Vector const& place) const
{
    FiberSegment res(nullptr, 0);
    real hit = INFINITY;
    
    for ( FiberSegment const& seg : allSegments )
    {
        // Compute the distance between 'place' and segment
        real dis = INFINITY;
        seg.projectPoint(place, dis);
        
        if ( dis < hit )
        {
            hit = dis;
            res = seg;
        }
    }
    
    return res;
}
