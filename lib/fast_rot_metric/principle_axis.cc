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

#include <cmath>
#include <cassert>
using namespace std;

#include <global_utils.h>
#include <eigen.h>
#include "principle_axis.h"

//Why not include an RCS initializer string?
//static char rcsid[] = "@(#)$Header: /usr/local/src/minrms/lib/fast_rot_metric/RCS/principle_axis.cc,v 3.3 2002/10/15 21:45:09 conrad Exp $";

namespace minrms
{

void
FindPrincipleAxis(ConstVect3 *aR,   //<--array of input points
                  int num_points,   //<--number of entries in aR
                  Vect3 eigenvalues,
                  Matrix3x4 principleAxisAndCentroid)
{
  int i,j,a;

  //Allocate a 3x3 Matrix to be passed to eigen-finding function in jacobi.c
  //The format is a pointer array, each of which points to an array.
  //The indexing they use starts at 1 instead of 0, and goes to 3.
  //Hence the call to "new Real [4];" instead of "new Real[3];"
  Real **aaRiRj = new Real * [4];
  Real **eigenvects = new Real * [4];
  if ((! aaRiRj) || (! eigenvects))
    ERR(__FILE__ << __LINE__ <<
        "CalculatePrincipleAxis() Error alloc mem.");
  aaRiRj[0] = NULL; //dosen't hurt to flag enries that shouldn't be used.
  eigenvects[0] = NULL;
  for (i=1; i <= 3; ++i) {
    aaRiRj[i] = new Real [4];
    eigenvects[i] = new Real [4];
    if ((! aaRiRj[i]) || (! eigenvects[i]))
      ERR(__FILE__ << __LINE__ <<
          "CalculatePrincipleAxis() Error alloc mem.");
    aaRiRj[i][0] = -1.0f; //optional: flag memory that shouldn't be accessed
    eigenvects[i][0] = -1.0f;
    //zero fill the useful parts of aaRiRj
    for (j=1; j <= 3; ++j)
      aaRiRj[i][j] = 0.0f;
  }

  //Next, recalculate the coordinates with the center-of-mass (centroid)
  //subtracted
  Vect3 centroid;
  Vect3 *aR_CM = new Vect3 [num_points];
  if (! aR_CM)
    ERR(__FILE__ << __LINE__ <<
        "CalculatePrincipleAxis() Error alloc mem.");

  FindCentroid(num_points, aR, centroid);

  //cout << " num_points = " << num_points << endl;

  for (a=0; a < num_points; ++a)
  {
    aR_CM[a][0] = aR[a][0] - centroid[0];
    aR_CM[a][1] = aR[a][1] - centroid[1];
    aR_CM[a][2] = aR[a][2] - centroid[2];
    //cout << " aR_CM[" << a << "] = ("
    //     << aR_CM[a][0] << "," 
    //     << aR_CM[a][1] << ","
    //     << aR_CM[a][2] << ")" << endl;
  }

  //Now, fill the aaRiRj matrix that will be diagonalized.
  for (a=0; a < num_points; ++a) {
    for (i=0; i<3; ++i) {
      for (j=0; j<3; ++j)
        aaRiRj[i+1][j+1] += aR_CM[a][i] * aR_CM[a][j];
    }
  }

  DEBUG_MSG(DBG_PRINCIPLE_AXIS, " -----------------------\n"
            << "         [" << aaRiRj[1][1] << "  " << aaRiRj[1][2] << "  " << aaRiRj[1][3] << "]\n"
            << "aaRiRj = [" << aaRiRj[2][1] << "  " << aaRiRj[2][2] << "  " << aaRiRj[2][3] << "]\n"
            << "         [" << aaRiRj[3][1] << "  " << aaRiRj[3][2] << "  " << aaRiRj[3][3] << "]");

  assert(principleAxisAndCentroid);

  Real d[4];
  int num_iter_dummy;//Not used.It's required that we pass an int * to jacobi()
  jacobi(aaRiRj, 3, d, eigenvects, &num_iter_dummy);

  DEBUG_MSG(DBG_PRINCIPLE_AXIS, " ----------------------------------\n" 
            << "             ["<<eigenvects[1][1]<<"  "<<eigenvects[1][2]<<"  "<<eigenvects[1][3]<<"]\n"
            << "eigenvects = ["<<eigenvects[2][1]<<"  "<<eigenvects[2][2]<<"  "<<eigenvects[2][3]<<"]\n"
            << "             ["<<eigenvects[3][1]<<"  "<<eigenvects[3][2]<<"  "<<eigenvects[3][3]<<"]"
            );
  //Store the eigenvalues, stored in d, in the "eigenvalues" argument
  eigenvalues[0] = d[1];
  eigenvalues[1] = d[2];
  eigenvalues[2] = d[3];

  DEBUG_MSG(DBG_PRINCIPLE_AXIS, " ----------------------------------\n"
            "               " << eigenvalues[0] << "\n"
            "eigenvalues =  " << eigenvalues[1] << "\n"
            "               " << eigenvalues[2]);

  //Now, store the normalized eigenvector corresponding to the ith eigenvalue
  //in the first three elements of the ith row of the matrix that
  //stores the transformation from original to center-of-mass-
  //principle-axis coordinates.
  for(i=0; i < 3; ++i) {
    for(j=0; j < 3; ++j)
      principleAxisAndCentroid[i][j] = eigenvects[i+1][j+1];
    Normalize3( principleAxisAndCentroid[i] );
  }

  
  //Store the Centroid int the 4th collumn of PrincipleAxisAndCentroid[][]
  principleAxisAndCentroid[0][3] = centroid[0];
  principleAxisAndCentroid[1][3] = centroid[1];
  principleAxisAndCentroid[2][3] = centroid[2];

  //Now, deallocate memory.
  assert(aaRiRj && eigenvects && aR_CM);
  for(i = 1; i <= 3; ++i)
  {
    assert(aaRiRj[i] && eigenvects[i]);
    delete [] aaRiRj[i];
    delete [] eigenvects[i];
  }
  delete [] aaRiRj;
  delete [] eigenvects;
  delete [] aR_CM;
  
} //FindPrincipleAxis()


} //namespace minrms
