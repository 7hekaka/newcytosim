// Cytosim was created by Francois Nedelec. Copyright 2023 Cambridge University.
/*
 Adapted from:
 https://www.geeksforgeeks.org/cpp-program-for-quicksort
 https://www.geeksforgeeks.org/heap-sort
 
 FJN SLCU Cambridge, 4.04.2023
*/

// manage a subtree rooted at node i which is an index in list[]. n is size of heap
template < typename T >
void heapify(T & list, size_t n, size_t i)
{
    size_t largest = i;
    size_t l = 2 * i + 1;
    size_t r = 2 * i + 2;
 
    // If left child is larger than root
    if ( l < n && list.compare(largest, l) )
        largest = l;
 
    // If right child is larger than largest so far
    if ( r < n && list.compare(largest, r) )
        largest = r;
 
    if ( largest != i )
    {
        list.swap(i, largest);
 
        // Recursively heapify the affected sub-tree
        heapify(list, n, largest);
    }
}
 
template < typename T >
void heapsort(T & list, size_t n)
{
    // Build heap (rearrange array)
    for ( size_t i = n / 2; i-- > 0; )
        heapify(list, n, i);
 
    // One by one extract element from heap
    for ( size_t i = n - 1; i > 0; --i )
    {
        // Move current root to end
        list.swap(0, i);
 
        // call max heapify on the reduced heap
        heapify(list, i, 0);
    }
}

//------------------------------------------------------------------------------
// this is not a true quicksort, since there is only using swaps()

template < typename T >
int partition(T & list, size_t start, size_t end)
{
    size_t pivot = start;
    size_t cnt = 0;
    
    for ( size_t i = start + 1; i < end; ++i )
        if ( !list.compare(pivot, i) ) // pivot < list[i]
            ++cnt;
 
    if ( cnt )
    {
        // Giving pivot element its correct position
        pivot += cnt;
        list.swap(pivot, start);
    }
    // Sorting left and right parts of the pivot element
    assert_true( end > 0 );
    size_t i = start, j = end-1;
 
    while ( i < pivot && j > pivot )
    {
        while ( !list.compare(pivot, i) ) // list[i] <= pivot
            ++i;
 
        while ( list.compare(pivot, j) ) // pivot < list[j]
            --j;
 
        if ( i < pivot && j > pivot )
            list.swap(i++, j--);
    }
 
    return pivot;
}

/* argument `end` is the index past the last valid element */
template < typename T >
void quicksort(T & list, size_t start, size_t end)
{
    if ( start < end )
    {
        // partitioning the array
        size_t p = partition(list, start, end);
        
        // Sorting the left part
        quicksort(list, start, p);
        
        // Sorting the right part
        quicksort(list, p+1, end);
    }
}
 
