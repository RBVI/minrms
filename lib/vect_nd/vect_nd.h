#ifndef _VECT_ND_H
#define _VECT_ND_H


#include <cstring> //contains memcpy()
using namespace std;

#include <simple_numeric_utils.h> //(located in the lib/global_utils directory)



//  Global variables:

// When this .h file is included in other files, either
// 1) the macro G_DIM must be defined beforehand in those files, or
// 2) somewhere in the binary this is being added to, there must be a 
//    global variable named g_dim.

#if (defined G_DIM)
const int g_dim = G_DIM;  //=the number of dimensions
#else
extern int g_dim;         //If "G_DIM" is not defined, then the value of this
                          //variable (#dimensions) can be specified at run time
#endif





namespace vect_nd {



template<class _Real>
inline void
CopyVect(_Real const *aA, _Real *aB)
{
  memcpy( (void*)aB, (void*)aA, sizeof(int)*g_dim);
} // Chain::CopyVect()



template<class _Real>
inline bool
EqualVect(_Real const *aA, _Real const *aB)
{
  int d = 0;
  for(; ((aA[d] == aB[d]) && (d < g_dim)); ++d)
    {}
  return d == g_dim;
}


template<class _Real>
inline _Real
DotProduct(_Real const *A, _Real const *B)
{
  _Real AdotB = 0.0;
  for (int d=0; d < g_dim; ++d)
    AdotB += A[d]*B[d];
  return AdotB;
}



#if ((defined G_DIM) && (G_DIM == 3))
//Takes the cross-product (A x B) and stores the result in dest.
//This function only makes sense in 3-dimensions.
//A somewhat more general version of this function is SubtractFrom2dSubspace().
//(See below.)
template<class _Real>
inline void
CrossProduct(_Real const *A, _Real const *B, _Real *dest)
{
  dest[0] = A[1]*B[2] - A[2]*B[1];
  dest[1] = A[2]*B[0] - A[0]*B[2];
  dest[2] = A[0]*B[1] - A[1]*B[0];
}
#endif




template<class _Real>
inline _Real
LengthSqd(_Real const *v)
{
  return DotProduct(v,v);
}


template<class _Real>
inline _Real
Length(_Real const *v)
{
  return sqrt(LengthSqd(v));
}



template<class _Real>
inline _Real
LengthAbs(_Real const *v)
{
  _Real x_plus_y = 0.0;
  for (int d=0; (d < g_dim); ++d)
    x_plus_y += (v[d]);
  return x_plus_y;
}



template<class _Real>
inline _Real
DistMetricRsqd(_Real const *p1, _Real const *p2)
{
  _Real r_sqd = 0.0;
  for (int d=0; (d < g_dim); ++d)
    r_sqd += SQR(p1[d] - p2[d]);
  return r_sqd;
}


template<class _Real>
inline _Real
DistMetricR(_Real const *p1, _Real const *p2)
{ return sqrt(DistMetricRsqd(p1, p2)); }


template<class _Real>
inline _Real
DistMetricABS(_Real const *p1, _Real const *p2)
{
  _Real x_plus_y = 0;
  for (int d=0; (d < g_dim); ++d)
    x_plus_y += ABS(p1[d] - p2[d]);
  return x_plus_y;
}


template<class _Real>
inline _Real
DistanceSqd(_Real const *p1, _Real const *p2)
{ return DistMetricRsqd(p1, p2);}


template<class _Real>
inline _Real
Distance(_Real const *p1, _Real const *p2)
{ return DistMetricR(p1, p2);}




// Normalize just divides the vector v_source by it's length, and stores the
// new vector in v_dest, unless v_source is the NULL-vector (0,0,0...)
// It returns the length of v_source.
template<class _Real>
inline _Real
Normalize(_Real *v)
{
  _Real length = Length(v);
  if (length != 0.0)
  {
    _Real one_over_length = 1.0 / length;
    for (int d=0; d < g_dim; ++d) 
      v[d] *= one_over_length;
  }
  else
  {
    v[0] = 1.0;
    for (int d=1; d < g_dim; ++d) 
      v[d] = 0.0;
  }
  return length;
}





// ---- ProjectOnto2dSubspace() and ----
// ---- FastProjectOnto2dSubspace() ----
// The following two functions compute the projection
// of vector "V" onto basis vectors "E1" and "E2"
// (not necessarily orthonormal) and stores the result in "dest".
// The version below is fast because many of the dot-products
// have been computed already.
template<class _Real>
inline void FastProjectOnto2dSubspace(_Real const *E1,
                                      _Real const *E2,
                                      _Real *dest,
                                      _Real V_dot_E1,
                                      _Real V_dot_E2,
                                      _Real E1_sqd,
                                      _Real E2_sqd,
                                      _Real E1_dot_E2)
{
  // Formula used:
  //
  //             1      //                     \     /                     \  \
  // dest = ----------- ||(<v,e1>-<v,e2><e1,e2>|e1 + |(<v,e1>-<v,e2><e1,e2>|e2|
  //        1-<e1,e2>^2 \\                     /     \                     /  /
  //
  // where v = V/|V|, e1 = E1/|E1|, e2 = E2/|E2| are unit vectors
  // and <e1,e2> denotes the dot-product.  Slightly more confusing to
  // read, the following formula is less expensive to calculate;
  // it only involves the original vectors, V, E1, and E2 (no square-roots):
  //
  //         //                           \    /                           \  \
  // dest= C ||(<V,E1>|E2|^2-<V,E2><E1,E2>|E1 +|(<V,E2>|E1|^2-<V,E1><E1,E2>|E2|
  //         \\                           /    \                           /  /
  //
  //                       1
  //   and C = -------------------------
  //           |E1|^2*|E2|^2 - <E1,E2>^2
  //
  _Real one_over_outer_coeff =  E1_sqd * E2_sqd - SQR(E1_dot_E2);
  if (one_over_outer_coeff == 0.0)
    return;
  _Real outer_coeff = 1.0 / one_over_outer_coeff;

  _Real E1_coeff = outer_coeff * (V_dot_E1 * E2_sqd  -  V_dot_E2 * E1_dot_E2);
  _Real E2_coeff = outer_coeff * (V_dot_E2 * E1_sqd  -  V_dot_E1 * E1_dot_E2);
  for (int d=0; d < g_dim; ++d)
    dest[d] = E1_coeff * E1[d]  +  E2_coeff * E2[d];
}





// ProjectOnto2dSubspace() is more general than FastProjectOnto2dSubspace().
// But it is slower.  Use this version if you haven't already computed all the
// lengths and the dot-products beforehand.

template<class _Real>
inline void ProjectOnto2dSubspace(_Real const *V,
                                  _Real const *E1,
                                  _Real const *E2,
                                  _Real *dest)
{
  _Real E1_sqd = LengthSqd(E1);
  _Real E2_sqd = LengthSqd(E2);
  _Real V_dot_E1 = DotProduct(V, E1);
  _Real V_dot_E2 = DotProduct(V, E2);
  _Real E1_dot_E2 = DotProduct(E1, E2);

  FastProjectOnto2dSubspace(E1,
                            E2,
                            dest,
                            V_dot_E1,
                            V_dot_E2,
                            E1_sqd,
                            E2_sqd,
                            E1_dot_E2);
} //inline void ProjectOnto2dSubspace()








// ---- SubtractFrom2dSubspace() and ----
// ---- FastSubtractFrom2dSubspace() ----
//
// These two functions compute the difference between the
// vector V, and its projection onto the plane spanned by E1 and E2.
// These function will be used often because they are an alternative
// to the Cross-Product.  Use this function to generate a vector perpindicular
// to two other vectors in dimensions higher than 3.
// 
// The function FastSubtractFrom2dSubspace() requires that
// you have precomputed all of the lengths and dot-products between
// the vectors.

template<class _Real>
inline void FastSubtractFrom2dSubspace(_Real const *V,
                                       _Real const *E1,
                                       _Real const *E2,
                                       _Real *dest,
                                       _Real V_dot_E1,
                                       _Real V_dot_E2,
                                       _Real E1_sqd,
                                       _Real E2_sqd,
                                       _Real E1_dot_E2)
                               
{
  FastProjectOnto2dSubspace(E1,
                            E2,
                            dest,
                            V_dot_E1,
                            V_dot_E2,
                            E1_sqd,
                            E2_sqd,
                            E1_dot_E2);

  for (int d=0; d < g_dim; ++d)
    dest[d] = V[d] - dest[d];
}




// SubtractFrom2dSubspace() is more general than FastSubtractFrom2dSubspace().
// But it is slower.  Use this version if you haven't already computed all the
// lengths and the dot-products beforehand.

template<class _Real>
inline void SubtractFrom2dSubspace(_Real const *V,
                                   _Real const *E1,
                                   _Real const *E2,
                                   _Real *dest)
{
  ProjectOnto2dSubspace(V, E1, E2, dest);
  for (int d=0; d < g_dim; ++d)
    dest[d] = V[d] - dest[d];
}


} // namespace vect_nd;


#endif //#ifndef _VECT_ND_H
