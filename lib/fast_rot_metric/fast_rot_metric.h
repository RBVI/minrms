//
// Copyright (c) 2002 The Regents of the University of California.
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//   1. Redistributions of source code must retain the above copyright
//      notice, this list of conditions, and the following disclaimer.
//   2. Redistributions in binary form must reproduce the above
//      copyright notice, this list of conditions, and the following
//      disclaimer in the documentation and/or other materials provided
//      with the distribution.
//   3. Redistributions must acknowledge that this software was
//      originally developed by the UCSF Computer Graphics Laboratory
//      under support by the NIH National Center for Research Resources,
//      grant P41-RR01081.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE REGENTS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
// OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
// BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
// OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
// EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef _FAST_ROT_METRIC_H
#define _FAST_ROT_METRIC_H

#include <cassert>
using namespace std;

#include <global_utils.h>
#include <vect_3d.h>
#include "principle_axis.h"

namespace minrms {

//SlowRotRMSD() returns the root-mean-squared distance between the positions
//of a set of points, with two different rotations and translations applied
//to it, collectively stored in "O1" and "O2".  The set of points is passed as
//an array called "aPos", and the number of elements in that array is
//passed in "numpoints".

Real
SlowRotSumSqdDist(ConstMatrix3x4 O1, ConstMatrix3x4 O2,
                  ConstVect3 *aPos,
                  int numpoints);



//FastRotSumSqdDist() returns the sum-squared distance between the positions
//of a set of points, with two different rotations and translations applied
//to it, collectively stored in O1 and O2.  Instead of applying the two
//rotations and translations to every point in the set, and calculating the
//rmsd manually, only, the principle axis are needed, and their eigenvalues.
//Not having to loop over every point in the set speeds things up considerably.
//The same answer is obtained, but using a faster method.
//
//Important: The two orientations stored in O1, and O2, are assumed to be
//           acting on a set of points that is already rotated so that
//           its centriod is at the origin, and
//           its 3 principle axis lie along the x, y, and z axis,
//           (with the corresponding eigenvalues are stored in
//            principleAxisEigenvalues[0], [1], [2])
//           If this is not the case, then you must apply additional rotations/
//           translations O1, and O2, to compensate.
//           (Hint: There is a Matrix3x4 returned by FindPrincipleAxis().
//                  If you add a 4th row storing [0 0 0 1] to this matrix,
//                  and multiply them together (with O1 or O2 on the LEFT),
//                  the resulting 3x4 product of the two matrices will have
//                  the desired property.  Pass this matrix to FastRotRMSD())
//
//           (I don't do this automatically, because it's a lot of extra
//            computation, this only needs to be done once per orientation
//            However, this function get's called much more frequently than
//            that.)
//
//Important2: Different kinds of litterature dissagree on their definition of
//           the principle axis.  I diagonalise Qij = sum_over_a {Rai*Raj}.
//           This yeilds eigenvalues that indicate the VARIANCE in the 
//           positions of the points along each axis, which are chosen to
//           minimize this variance.
//           This is as opposed to to finding the axis of minimal rotational
//           energy, which is what physics/mech-E people do.
//           (They turn out to be the same, but the eigenvalues are different.)
//               Using the function FindPrincipleAxis(), gaurantees
//           that the eigenvalues in the correct format and
//           order.  This function is located in
//             ../principle_axis/principle_axis.cc

inline Real
FastRotSumSqdDist(ConstMatrix3x4 O1, ConstMatrix3x4 O2,
                  ConstVect3 principleAxisEigenvalues,
                  int num_points)
{
  assert(O1 && O2 && principleAxisEigenvalues);
  Real sum_sqrd_dist = 0.0f;

  //First, calculate the rotational contribution to the RMSD
  for(int i = 0; i < 3; ++i)
  {
        sum_sqrd_dist += principleAxisEigenvalues[i] * (1 - O1[0][i]*O2[0][i]
                                                          - O1[1][i]*O2[1][i]
                                                          - O1[2][i]*O2[2][i]);
  }
  sum_sqrd_dist *= 2.0f;

  //Now, ad the translational component to the RMSD,
  //(This is stored in the fourth ("3") collumn of O1, and O2)
  sum_sqrd_dist += num_points * (SQR(O2[0][3] - O1[0][3]) +
                                 SQR(O2[1][3] - O1[1][3]) +
                                 SQR(O2[2][3] - O1[2][3]));

  if (sum_sqrd_dist < 0.0f)
  {
    //  cerr << "--err? rmsd orientation metric: sum-squared = "
    //       << sum_sqrd_dist           << "  (setting to 0.0f) --" << endl;
    sum_sqrd_dist = 0.0f;
  }

  //Now sum_sqrd_dist stores the sum of the squared distances between rotated
  //positions (hopefully).  To return RMSD, find the average, and take the
  //sqrt.

  assert(num_points);
  return sum_sqrd_dist;
} //FastRotSumSqdDist()

} //namespace minrms

#endif //#ifndef _FAST_ROT_METRIC_H
