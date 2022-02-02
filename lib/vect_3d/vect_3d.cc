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


#include <cassert> //defines assert()
#include <fstream>
using namespace std;

#include <binary_file.h>
#include <simple_numeric_utils.h>
#include "vect_3d.h"



namespace vect_3d {


//global variables:
const Matrix3x4 g_IDENTITY3X4 = { {1.0, 0.0, 0.0, 0.0},
                                  {0.0, 1.0, 0.0, 0.0},
                                  {0.0, 0.0, 1.0, 0.0} };

Real Normalize3(Vect3 x)
{
  Real length = sqrt(SQR(x[0])+SQR(x[1])+SQR(x[2]));
  assert(length != 0.0f);
  Real one_over_length = 1.0/length;
  x[0] *= one_over_length;
  x[1] *= one_over_length;
  x[2] *= one_over_length;
  return length;
} 

Real Normalize4(Vect4 x)
{
  Real length = sqrt(SQR(x[0])+SQR(x[1])+SQR(x[2])+SQR(x[3]));
  assert(length != 0.0f);
  Real one_over_length = 1.0/length;
  x[0] *= one_over_length;
  x[1] *= one_over_length;
  x[2] *= one_over_length;
  x[3] *= one_over_length;
  return length;
} 

void
FindCentroid(short numpoints, ConstVect3 *aPoints, Vect3 dest)
{
  assert( aPoints );
  dest[0] = 0.0f;
  dest[1] = 0.0f;
  dest[2] = 0.0f;
  ConstVect3 *pPoint = aPoints;
  Real one_over_np = 1.0f/(Real)numpoints;
  while(numpoints>0)
  {
    assert( *pPoint );
    dest[0] += (*pPoint)[0];
    dest[1] += (*pPoint)[1];
    dest[2] += (*pPoint)[2];
    ++pPoint;
    --numpoints;
  }
  dest[0] *= one_over_np;
  dest[1] *= one_over_np;
  dest[2] *= one_over_np;
} // FindCentroid()


void
ConvertQuaternionToMatrix3x3(ConstQuaternion input, Matrix3x3 dest)
{
  // This formula was taken from Eq 7 of R. Diamond 1998

  assert(dest && input);
  //first, the diagonal elements:
  dest[0][0] = +SQR(input[0])-SQR(input[1])-SQR(input[2])+SQR(input[3]);
  dest[1][1] = -SQR(input[0])+SQR(input[1])-SQR(input[2])+SQR(input[3]);
  dest[2][2] = -SQR(input[0])-SQR(input[1])+SQR(input[2])+SQR(input[3]);

  //the off-diagonal elements:
  dest[0][1] = 2*(input[0]*input[1] - input[2]*input[3]);
  dest[1][0] = 2*(input[0]*input[1] + input[2]*input[3]);

  dest[1][2] = 2*(input[1]*input[2] - input[0]*input[3]);
  dest[2][1] = 2*(input[1]*input[2] + input[0]*input[3]);

  dest[0][2] = 2*(input[0]*input[2] + input[1]*input[3]);
  dest[2][0] = 2*(input[0]*input[2] - input[1]*input[3]);
} //ConvertQuaternionToMatrix3x3()


void
BinaryWriteMatrix3x4(fstream& binary_file,
                     ConstMatrix3x4 M)
{
  for (int i = 0; i < 3; ++i)
  {
    binary_write(binary_file, M[i][0]);
    binary_write(binary_file, M[i][1]);
    binary_write(binary_file, M[i][2]);
    binary_write(binary_file, M[i][3]);
  }
}



void
BinaryReadMatrix3x4(fstream& binary_file,
                    ConstMatrix3x4 M)
{
  for (int i = 0; i < 3; ++i)
  {
    binary_read(binary_file, M[i][0]);
    binary_read(binary_file, M[i][1]);
    binary_read(binary_file, M[i][2]);
    binary_read(binary_file, M[i][3]);
  }
}
                                         

void
WriteMatrix3x4ToFile(ConstMatrix3x4 M,
                     char const *filename)
{
  assert(filename);
  ofstream mat_file(filename, ios::out);
  if (! mat_file) 
    ERR("Cannot open file \""
        << filename
        << "\" for writing.  Exiting...");
  mat_file.precision(7);
  for (int i = 0; i < 3; ++i)
    mat_file << M[i][0] <<" "<< M[i][1] <<" "<< M[i][2] <<" "<< M[i][3] <<"\n";
}


void
ReadTransFile(char const *filename,
              Matrix3x4 dest)
{
  assert(filename);
  ifstream mat_file(filename, ios::in);
  if (! mat_file) 
    ERR("Cannot open file \""
        << filename
        << "\" for reading.  Exiting...");
  for (int i = 0; i < 3; ++i) {
    mat_file >> dest[i][0];
    mat_file >> dest[i][1];
    mat_file >> dest[i][2];
    mat_file >> dest[i][3];
  }

  if (! mat_file.good())
    ERR("Error reading \"" << filename << "\".  Not a .trans file.");
}


} //namespace vect_3d

