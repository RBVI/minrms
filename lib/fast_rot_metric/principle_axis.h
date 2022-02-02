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

#ifndef _PRINCIPLE_AXIS_H
#define _PRINCIPLE_AXIS_H

#include <simple_numeric_utils.h>
#include <vect_3d.h>

using namespace vect_3d;

namespace minrms {
//FindPrincipleAxis() takes a list of points in 3 dimensions, and returns
//a 3x4 matrix principleAxisAndCentroid, storing the 3 principle axis of
//those points in the first three collumns,and the centroid in the 4th collumn.
//   All points are weighted equally in the calculation.
//The results are in the following form:
//You can find the principle axis, by examining the rows of T
//The ith collumn of T, store the normalized eigenvector corresponding
//to the ith entry of "eigenvalues".
//The 4th collunn stores the centroid (center-of-mass) of the points.  To
// get the matrix that represents the transformation from original coordinates
//to center-of-mass-principle-axis coordinates, calculate the inverse of the
//"principleAxisAndCentroid" matrix using FindInverseOfRotTransMatrix3x4().
//
// Important:
//           This yeilds eigenvalues that indicate the VARIANCE in the 
//           positions of the points along each axis, and the principle
//           axis are chosen to minimize this variance.
//           This is as opposed to to finding the axis of minimal ROTATIONAL
//           energy (as physics/mech-E people do).  The axis turns out
//           to be the same for both these methods, but the EIGENVALUES
//           are different.  The eigenvalues for the first method yeild
//           degree of variance along each axis (which is what this function
//           returns).
//           The eigenvalues for the second method yeild rotational energy
//           when spinning around each axis (NOT what this function returns).

void
FindPrincipleAxis(ConstVect3 *aR, //<--array of input points
                  int num_points, //<--number of entries in aR
                  Vect3 eigenvalues,
                  Matrix3x4 principleAxisAndCentroid);

} //namespace minrms

#endif //#ifndef _PRINCIPLE_AXIS_H
