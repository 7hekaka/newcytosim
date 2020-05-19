// Cytosim was created by Francois Nedelec. Copyright Cambridge University 2020
// FJN 18.05.2020

#include "node.h"


Node* Node::split(Node* head)
{
    Node* slow = head;
    Node* fast = head;
    
    while ( fast->nNext && fast->nNext->nNext )
    {
        fast = fast->nNext->nNext;
        slow = slow->nNext;
    }
    
    Node *temp = slow->nNext;
    slow->nNext = nullptr;
    return temp;
}


Node* Node::merge(int (*comp)(const Node*, const Node*), Node *first, Node *second)
{
    // If first linked list is empty
    if ( !first )
        return second;
    
    // If second linked list is empty
    if ( !second )
        return first;
    
    // Pick the smaller value
    if ( comp(second, first) > 0 )
    {
        first->nNext = merge(comp, first->nNext, second);
        first->nNext->nPrev = first;
        first->nPrev = nullptr;
        return first;
    }
    else
    {
        second->nNext = merge(comp, first, second->nNext);
        second->nNext->nPrev = second;
        second->nPrev = nullptr;
        return second;
    }
}


Node* Node::mergesort(int (*comp)(const Node*, const Node*), Node* head)
{
    if ( !head || !head->nNext )
        return head;
    
    Node* tail = split(head);
    
    // Recur for left and right halves
    head = mergesort(comp, head);
    tail = mergesort(comp, tail);
    
    // Merge the two sorted halves
    return merge(comp, head, tail);
}

