// Cytosim was created by Francois Nedelec. Copyright 2023 Cambridge University.
/*
 Adapted from:
 https://www.geeksforgeeks.org/introduction-to-smooth-sort/
 FJN SLCU Cambridge, 2.05.2023
*/
 
/// returns Leonardo numbers
static inline int leonardo(int k)
{
    if ( k < 2 )
        return 1;
    return leonardo(k-1) + leonardo(k-2) + 1;
}
 
/// Build the Leonardo heap by merging pairs of adjacent trees
template < typename T >
static void smooth_heapify(T& list, int start, int end)
{
    int i = start;
    int j = 0;
    int k = 0;
 
    while ( k <= end - start )
    {
        if ( k & 0xAAAAAAAA ) {
            j = j + i;
            i = i >> 1;
        }
        else {
            i = i + j;
            j = j >> 1;
        }
 
        ++k;
    }
 
    while ( i > 0 )
    {
        j = j >> 1;
        k = i + j;
        while ( k < end )
        {
            if ( list.compare(k-i, k) ) // list[k] > list[k - i]) {
                break;
            list.swap(k, k-i); //swap(list[k], list[k - i]);
            k = k + i;
        }
        i = j;
    }
}
 
/// Edsger Dijkstra's Smooth Sort algorithm (EWD796a, 16-Aug-1981)
template < typename T >
void smoothsort(T& list, int n)
{
    int p = n - 1;
    int q = p;
    int r = 0;
 
    // Build the Leonardo heap by merging pairs of adjacent trees
    while ( p > 0 )
    {
        if ( (r & 0x03) == 0 )
            smooth_heapify(list, r, q);
 
        if ( leonardo(r) == p ) {
            ++r;
        }
        else {
            r = r - 1;
            q = q - leonardo(r);
            smooth_heapify(list, r, q);
            q = r - 1;
            ++r;
        }
 
        list.swap(0, p);
        --p;
    }
 
    // Convert the Leonardo heap back into an array
    for ( int i = 0; i < n-1; ++i )
    {
        int j = i + 1;
        while ( j > 0 && list.compare(j, j-1) ) //arr[j] < arr[j - 1])
        {
            list.swap(j, j-1);
            --j;
        }
    }
}

