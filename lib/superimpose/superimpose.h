#ifndef _SUPERIMPOSE_H
#define _SUPERIMPOSE_H

// "superimpose.h" 
// This file contains the declaration of a couple functions used for finding
// optimal translation and rotation between two rigid sets of points that
// minimizes the total-square-of-the-distance between each pair.
// The method used is outlined in the article by R. Diamond:
//   "A Note on the Rotational Superposition Problem", (1988)
//   Acta Cryst. A44, pp. 211-216

#include <vect_3d.h>




// (If, for some reason, you want to have run-time debug information
//  printed to the terminal, then set DBG_ROT to 1 instead of 0.)
#define DBG_ROT 0







namespace vect_3d {


class Superimpose
{
  //Private constructor. Do not instantiate or inherit. Namespace only.
  //  Superimpose();
  //static variables inside this namespace:
  // (The next two are used by FindTransformRMSD().  Because
  //  FindTransformRMSD() get's called frequently, I only want to
  //  allocate these guys once.)
  Vect3 *aFixedNew;   //A copy of aFixed    with the centroid subtracted
  Vect3 *aMoveableNew;//A copy of aMoveable with the centroid subtracted
  Real **P; //4x4 pointer array. Passed as argument to tred2()
  Real **eigenvects; //4x4 pointer array.
                              //Stores the eigenvectors returned by tqli()

  // Call the next function before you call FindTransformMinRMS()
  // (The argument "maxNumPoints" should be an upper-bound on the value
  //  of anything passed to FindTransformRMS()'s "numpoints" argument.)
  // (Because FindTransformRMS() get's called frequently, I only wanted
  //  to allocate space for certain variables once.)
  void
  AllocGlobalVariables(int maxNumPoints);
  void
  DeallocGlobalVariables();

public:

  inline Superimpose(int maxNumPoints)
  {
    aFixedNew = aMoveableNew = NULL;
    P = eigenvects = NULL;
    AllocGlobalVariables(maxNumPoints);
  }

  inline ~Superimpose()
  {
    DeallocGlobalVariables();
  }

  // The function Superimpose::FindCentroid() finds the evenly weighted
  // average position of a set of points, and stores the result in "dest"
  // There is no weighting factor ('Wa' in the paper is just 1.0).
  //  void
  //  FindCentroid(short numpoints, ConstVect3 *aPoints, Vect3 dest);

  // The function Superimpose::FindTransformMinRMSD() implements
  // R. Diamond,  'A Note on the Rotational Superposition Problem',
  // Journal of X-ray Crystallography, 1988
  // There is no weighting factor ('Wa' in the paper is just 1.0).
  // It returns the RMSD of the best match, and stores the corresponding
  // rotation and translation matrix (the one that gets applied to the 
  // points in "aMoveable" to produce this best match) in the "dest" argument.)
  Real
  FindTransformMinRMSD(short numpoints,
                       ConstVect3 *aFixed, ConstVect3 *aMoveable,
                       Matrix3x4 dest);

}; //class Superimpose

} //namespace vect_3d

#endif // #ifndef _SUPERIMPOSE_H
