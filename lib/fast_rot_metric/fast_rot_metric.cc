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

#include "fast_rot_metric.h"


//Why not include an RCS initializer string?
//static char rcsid[] = "@(#)$Header: /usr/local/src/minrms/lib/fast_rot_metric/RCS/fast_rot_metric.cc,v 1.3 2002/10/15 21:45:09 conrad Exp $";

namespace minrms
{

//SlowRotSumSqdDist() returns the sum-squared distance between the positions
//of a set of points, with two different rotations and translations applied
//to it, collectively stored in "O1" and "O2".  The set of points is passed as
//an array called "aPos", and the number of elements in that array is
//passed in "numpoints".

float SlowRotSumSqdDist(ConstMatrix3x4 O1, ConstMatrix3x4 O2,
                        ConstVect3 *aPos,
                        int numpoints)
{
  assert(O1 && O2 && aPos);
  float sum_sqrd_dist = 0.0f;

  for (int i=0; i < numpoints; ++i) {
    Vect3 pos1, pos2;
    Mult_Mat3x4_by_Vect3(O1, aPos[i], pos1);
    Mult_Mat3x4_by_Vect3(O2, aPos[i], pos2);
    sum_sqrd_dist += SQR(pos2[0] - pos1[0]) +
                     SQR(pos2[1] - pos1[1]) +
                     SQR(pos2[2] - pos1[2]);
  }
  return sum_sqrd_dist;
} //SlowRotRMSD()

} //namespace minrms
