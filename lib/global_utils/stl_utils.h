#ifndef _STL_UTILS
#define _STL_UTILS

#include <cassert> //defines assert()
#include <algorithm>
using namespace std;


//When a container reaches it's final size and
//(you do not anticipate adding any elements to it),
//"trimming" the container discards any unnecessary memory
//that was allocated for it (in case of future expansion).
//Restrictions:
//-The copy operator = must be defined.
//-The container must support iterators begin() and end() and
//-The constructor must accept an interval of iterators as an argument.
template<class Container>
inline void TrimContainer(Container& c)
{
  Container trimmed(c.begin(), c.end());
  swap(c, trimmed);
}



//DeleteContents empties the container.  Unlike clear(), this function
//is gauranteed to free up the memory allocated for it's previous contents.
//Restrictions:
//-The copy operator = must be defined.
//-There must be a default constructor
template<class Container>
inline void DeleteContents(Container& c)
{
  Container empty;
  swap(c, empty);
}




// The following function copies the contents of a C array into
// an STL container.
// The container must support the following member functions:
// clear(),
// push_back(T entry)
// begin(), end(),
// A constructor which accepts an interval of iterators as an argument:
// Container(Container::const_iterator start, Container::const_iterator stop)
// and the = operator for the container

template<class Container, class T>
inline void CopyFromCarray(long array_size,
                           T const *array,
                           Container& c)
{
  c.clear();
  for (long i=0; i < array_size; ++i)
    c.push_back(array[i]);
  TrimContainer(c);
}



template<class ContainerIterator, class T>
inline void CopyToCarray(ContainerIterator a,
                         ContainerIterator b,
                         T *array,
                         bool do_realloc = false)
{
  if (do_realloc && array)//(This is dangerous.  I assume the programmer has
    delete [] array;//already to either allocated the array or set it to NULL)
  array = new T [b-a];
  for (ContainerIterator p = a; p != b; ++p)
    array[p-a] = *p;
}




template<class ContainerIterator, class Containee>
//This is a silly function.
//It's just a slightly more conveniant wrapper for the STL find() command.
//This function tests to see if "key" belongs to the set
//specified by the range from "a" to "b"-1.
//Running time is linear in b-a.
bool BelongsTo(Containee & key,
               ContainerIterator a,
               ContainerIterator b)
{
  ContainerIterator p = find(a, b, key);
  assert((p == b) || (*p == key)); //check if found, then *p == key
  return (p != b);
}

template<class ContainerIterator, class Containee>
//Tests to see if key belongs to the set of elements in the range "a" to "b"-1.
//Only works if the elements in the range a...b have been sorted in
//increasing order.
//Running time: O(log(b-a))
bool SortedBelongsTo(Containee const& key,
                     ContainerIterator a,
                     ContainerIterator b)
{
  assert(is_sorted(a, b));
  ContainerIterator p = binary_search(a, b, key);
  assert((p == b) || (*p == key)); //check if found, then *p == key
  return (p != b);
}



//Finds key in the interval [first, last), assuming
//the elements between first and last are in sorted order.
//If no such element exists, it returns "last"
//(Running time: O(log(last-first))
template <class RandomAccessIterator, class T>
RandomAccessIterator 
SortedFind(RandomAccessIterator first,
           RandomAccessIterator last,
           T const& key)
{
  RandomAccessIterator
  p = lower_bound(first,
                  last,
                  key);
  if ((p != last) && (*p != key))
    p = last;
  return p;
}



#endif //#ifndef _STL_UTILS

