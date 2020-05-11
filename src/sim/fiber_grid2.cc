// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

/** 
 This file implements a 'dummy' grid using STL code, which can be used as a reference.
 For each position, it calculates the geometrical distance to all fiber segments.
 This is algorithmically the slowest method, but it is simple and most likely correct!
 It is useful to get a ground truth and evaluate more advanced methods.
 */


typedef std::vector <FiberSegment> SegmentVector;

/// a list containing all segments, as a global variable
SegmentVector allSegments;


unsigned FiberGrid::setGrid(Space const*, real)
{
    LOG_ONCE("Cytosim is using a crude method to localize fibers!\n");
    return 0;
}

void FiberGrid::paintGrid(const Fiber * first, const Fiber * last)
{
    allSegments.clear();
    // add all segments
    for ( const Fiber * f = first ; f != last ; f=f->next() )
    {
        for ( size_t s = 0; s < f->nbSegments(); ++s )
            allSegments.emplace_back(f, s);
    }
}


void FiberGrid::createCells()
{
}

size_t FiberGrid::hasGrid() const
{
    return 1;
}


void FiberGrid::tryToAttach(Vector const& place, Hand& ha) const
{
    // randomize the list order
    std::random_shuffle( allSegments.begin(), allSegments.end() );

    // test all segments:
    for ( FiberSegment const& seg : allSegments )
    {
#if !TRICKY_HAND_ATTACHMENT
        if ( RNG.test(ha.prop->binding_prob) )
#else
        if ( RNG.flip_8th() )
#endif
        {
            // Compute the distance between 'place' and segment
            real dis = INFINITY;
            real abs = seg.projectPoint(place, dis);
            
            /*
             Compare to the maximum attachment range of the hand,
             and compare a newly tossed random number with 'prob'
             */
            if ( dis < ha.prop->binding_range_sqr )
            {
                Fiber * fib = const_cast<Fiber*>(seg.fiber());
                FiberSite pos(fib, seg.abscissa1()+abs);
                
                if ( ha.attachmentAllowed(pos) )
                {
                    ha.attach(pos);
                    return;
                }
            }
        }
    }
}


FiberGrid::SegmentList FiberGrid::nearbySegments(Vector const& place, const real DD, Fiber * exclude) const
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
