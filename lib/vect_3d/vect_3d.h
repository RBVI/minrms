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


#ifndef _VECT_3D_H
#define _VECT_3D_H

// This file contains a bunch of data types and short (mostly inlined)
// utilities for matrix multiplication, dot/cross-products,
// copying matrices, and rotation-using-quaterions.
//(Note: More complicated routines for reducing and diagonalizing
//       matrices are not declared here, but should be located in
//       ../numrec/include/tred2.h and ../numrec/include/tqli.h)

#include <cassert> //defines assert()
#include <cmath>
#include <fstream>
using namespace std;

#include <global_utils.h>
#include <simple_numeric_types.h>





namespace vect_3d {



//typedefs:

// I use floats.  Doubles are incompatible with the
// code in the eigen library (taken from "Numerical Recipies in C"
// They use floats.)  
// I get tired of trying to remember if I am using floats or doubles.
// I call them "Reals" from now on.

typedef float Real;



//Note: The type "Real" is defined in multiple places
//      It is defined here, and also it may also be defined separately
//      again in the source files of the binaries that use this file.
//      Be sure to make the distinction between the version of "Real"
//      defined there, and "vect_3d::Real" defined here.




typedef Real Vect3[3];
typedef Real const ConstVect3[3];

typedef Real Vect4[4];
typedef Real const ConstVect4[4];

typedef Vect4  Quaternion;
typedef Vect4  const ConstQuaternion;

typedef Real Matrix3x3[3][3];
typedef Real const ConstMatrix3x3[3][3];

typedef Real Matrix3x4[3][4];
typedef Real const ConstMatrix3x4[3][4];

typedef Real Matrix4x4[4][4];
typedef Real const ConstMatrix4x4[4][4];


//global variables
extern const Matrix3x4 g_IDENTITY3X4;


//function prototypes
// (none, for now)


//Takes the cross-product (A x B) and stores the result in dest.
//It returns whether or not A x B == (0,0,0)  (true if A x B == (0,0,0))
inline bool CrossProduct(ConstVect3 A, ConstVect3 B, Vect3 dest)
{
  bool equals_zero = true;
  dest[0] = A[1]*B[2] - A[2]*B[1];
  if (dest[0] != 0.0)
    equals_zero = false;
  dest[1] = A[2]*B[0] - A[0]*B[2];
  if (dest[1] != 0.0)
    equals_zero = false;
  dest[2] = A[0]*B[1] - A[1]*B[0];
  if (dest[2] != 0.0)
    equals_zero = false;
  return equals_zero;
}


// Multiplies the vector X, by the matrix M, and stores the result in "dest".
// (It's okay if X and dest, are the same, -temporary space is used and a
// copy is performed)
inline void Mult_Mat3x3_by_Vect3(ConstMatrix3x3 M, ConstVect3 X, Vect3 dest)
{
  Vect3 temp;

  temp [0] = M[0][0]*X[0] + M[0][1]*X[1] + M[0][2]*X[2];
  temp [1] = M[1][0]*X[0] + M[1][1]*X[1] + M[1][2]*X[2];
  temp [2] = M[2][0]*X[0] + M[2][1]*X[1] + M[2][2]*X[2];
  dest [0] = temp[0];
  dest [1] = temp[1];
  dest [2] = temp[2];
}


// Multiplies the vector "X", by the left-most 3x3 submatrix in "M", and adds
// The remaining collumn of "M" is added to "X" (as a translational offset).
// The result is stored in "dest".
// (It's okay if X and dest, are the same, -temporary space is used and a
// copy is performed)
inline void Mult_Mat3x4_by_Vect3(ConstMatrix3x4 M, ConstVect3 X, Vect3 dest)
{
  Vect3 temp;

  temp [0] = M[0][0]*X[0] + M[0][1]*X[1] + M[0][2]*X[2] + M[0][3];
  temp [1] = M[1][0]*X[0] + M[1][1]*X[1] + M[1][2]*X[2] + M[1][3];
  temp [2] = M[2][0]*X[0] + M[2][1]*X[1] + M[2][2]*X[2] + M[2][3];
  dest [0] = temp[0];
  dest [1] = temp[1];
  dest [2] = temp[2];
}








// The next two functions multiply A by B (from left to right) and
// store the result in "dest".
// It's okay if "dest" is the same address as A or B
inline void
Mult_Mat3x3_by_Mat3x3_equals_Mat3x3(ConstMatrix3x3 A,
                                    ConstMatrix3x3 B,
                                    Matrix3x3 dest)
{
  Real temp[3][3];
  int i,j;
  //calculate the product of A and B
  for (i=0; (i<3); ++i) {
    for (j=0; (j<3); ++j) {
      *((*(temp + i))+j) = (*((*(A+i))+0)) * (*((*(B+0))+j)) +
                       (*((*(A+i))+1)) * (*((*(B+1))+j)) +
                           (*((*(A+i))+2)) * (*((*(B+2))+j));
    }
  }
  //now copy the result into dest
  for (i=0; (i<3); ++i)
    for (j=0; (j<3); ++j)
      /* dest[i][j] = temp[i][j]; */
      *((*(dest + i))+j) =  *((*(temp+i))+j);
}

//The next function effectively adds a fourth row [0 0 0 1] to matrix B,
//to get        __           __
//              |  original   |
//              |contents of B|
//         B' = |-------------| 
//              | 0  0  0  1  |
//              --           --
// Then the product of A B' is calculated and the result is stored in "dest"
// It's okay if "dest" is the same address as A or B
inline void
Mult_Mat3x4_by_Mat3x4_equals_Mat3x4(ConstMatrix3x4 A,
                                    ConstMatrix3x4 B,
                                    Matrix3x4 dest)
{
  Real temp[4][4];
  int i,j;
  //calculate the product of A and B
  for (i=0; (i<3); ++i) {
    for (j=0; (j<3); ++j) {
      *((*(temp + i))+j) = (*((*(A+i))+0)) * (*((*(B+0))+j)) +
                           (*((*(A+i))+1)) * (*((*(B+1))+j)) +
                           (*((*(A+i))+2)) * (*((*(B+2))+j));
      //                 + (*((*(A+i))+3)) * (*((*(B+3))+j));
    }
    //Now, handle the special case where j = 3 (ie, fourth collumn)
    *((*(temp + i))+3) = (*((*(A+i))+0)) * (*((*(B+0))+3)) +
                         (*((*(A+i))+1)) * (*((*(B+1))+3)) +
                         (*((*(A+i))+2)) * (*((*(B+2))+3)) +
                         (*((*(A+i))+3));
  }
  //now copy the result into dest
  for (i=0; (i<3); ++i)
    for (j=0; (j<4); ++j)
      // dest[i][j] = temp[i][j];
      *((*(dest + i))+j) =  *((*(temp+i))+j);
}


//Finds the inverse of a transformation composed solely by a rotation, and
//a translation, formatted as a 3x4 matrix  (That is, it is assumed that the
//rotational part of the transformation is stored in the first three collumns,
//and the translational offset (which is added after the rotation) is stored
//in the last (4th) collumn of RT.
//Note: it's okay if RT is the same address as dest.
inline void
FindInverseOfRotTransMatrix3x4(ConstMatrix3x4 RT, Matrix3x4 dest)
{
  //First the rotational (first 3 collumns) of dest is calculated by
  //taking the transpose of the rotational part of RT.
  //This is because for a 3x3 rotation matrix, R, inverse(R) = transpose(R)
  for(int i = 0; i < 3; ++i) {
    dest[i][0] = RT[0][i];
    dest[i][1] = RT[1][i];
    dest[i][2] = RT[2][i];
  }
  //Now, calculate the translational offset that counteracts the effects of
  //the original translational offset (stored in the 4th collumn of RT).
  //This is done in the following way.
  //Suppose  x'= Rx + t, where R is a rotation matrix, and t is a vector
  //storing the translation.
  //To invert this relationship, suptract "t" from both sides, and 
  //multiply by the left on both sides by the transpose of R. (same as inverse)
  //         x = transpose(R) (x' - t)
  //So, the translational offset is  (transpose(R)*(-t)), not -t.

  Vect3 temp;
  temp[0] = -dest[0][0]*RT[0][3] - dest[0][1]*RT[1][3] - dest[0][2]*RT[2][3];
  temp[1] = -dest[1][0]*RT[0][3] - dest[1][1]*RT[1][3] - dest[1][2]*RT[2][3];
  temp[2] = -dest[2][0]*RT[0][3] - dest[2][1]*RT[1][3] - dest[2][2]*RT[2][3];

  //now copy "temp" into the fourth collumn of "dest"
  dest[0][3] = temp[0];
  dest[1][3] = temp[1];
  dest[2][3] = temp[2];
} //FindInverseOfRotTransMatrix3x4()



// returns the dot-product of vectors A . B
inline Real DotProduct3(Vect3 A, Vect3 B)
{ return ((A)[0]*(B)[0] + (A)[1]*(B)[1] + (A)[2]*(B)[2]); }

inline void SubtractVect3(ConstVect3 const A, ConstVect3 const B, Vect3 dest)
{
  dest[0] = (A)[0]-(B)[0];
  dest[1] = (A)[1]-(B)[1];
  dest[2] = (A)[2]-(B)[2];
}

inline void AddVect3(ConstVect3 A, ConstVect3 B, Vect3 dest)
{
  dest[0] = (A)[0]+(B)[0];
  dest[1] = (A)[1]+(B)[1];
  dest[2] = (A)[2]+(B)[2];
}

inline void ScalarMultVect3(Real a, ConstVect3 const v, Vect3 dest)
{
  dest[0] = a * v[0];
  dest[1] = a * v[1];
  dest[2] = a * v[2];
}


inline void ResetToIdentity3x3(Matrix3x3 dest)
{
  for (int i = 0; i<3; ++i)
  {
    for (int j = 0; j<3; ++j)
      dest[i][j] = ((i==j) ? 1.0f : 0.0f);
  }
}


inline void ResetToIdentity3x4(Matrix3x4 dest)
{
  for (int i = 0; i<3; ++i)
  {
    for (int j = 0; j<4; ++j)
      dest[i][j] = ((i==j) ? 1.0f : 0.0f);
  }
}


inline void ResetToIdentity4x4(Matrix4x4 dest)
{
  for (int i = 0; i<4; ++i)
  {
    for (int j = 0; j<4; ++j)
      dest[i][j] = ((i==j) ? 1.0f : 0.0f);
  }
}

inline void ResetToZero3x3(Matrix3x3 dest)
{
  assert(dest);
  for (int i = 0; i<3; ++i)
  {
    for (int j = 0; j<3; ++j)
      dest[i][j] = 0.0f;
  }
}

inline void ResetToZero3x4(Matrix3x4 dest)
{
  assert(dest);
  for (int i = 0; i<3; ++i)
  {
    for (int j = 0; j<4; ++j)
      dest[i][j] = 0.0f;
  }
}


inline void ResetToZero4x4(Matrix4x4 dest)
{
  assert(dest);
  for (int i = 0; i<4; ++i)
  {
    for (int j = 0; j<4; ++j)
      dest[i][j] = 0.0f;
  }
}

inline Real FindTrace3x3(ConstMatrix3x3 M)
{
  assert(M);
  //  Real t=0.0f;
  //  for (int i = 0; i<3; ++i)
  //    t += M[i][i];
  return M[0][0] + M[1][1] + M[2][2];
}

inline Real FindTrace4x4(ConstMatrix4x4 M)
{
  assert(M);
  //  Real t=0.0f;
  //  for (int i = 0; i<4; ++i)
  //    t += M[i][i];
  return M[0][0] + M[1][1] + M[2][2] + M[3][3];
}

inline void CopyMatrix3x3(ConstMatrix3x3 source, Matrix3x3 dest)
{
  for(int i=0; i<3; ++i)
  {
    dest[i][0] = source[i][0];
    dest[i][1] = source[i][1];
    dest[i][2] = source[i][2];
  }
}

inline void CopyMatrix3x4(ConstMatrix3x4 source, Matrix3x4 dest)
{
  for(int i=0; i<3; ++i)
  {
    dest[i][0] = source[i][0];
    dest[i][1] = source[i][1];
    dest[i][2] = source[i][2];
    dest[i][3] = source[i][3];
  }
}

inline void CopyMatrix4x4(ConstMatrix4x4 source, Matrix4x4 dest)
{
  for(int i=0; i<4; ++i)
  {
    dest[i][0] = source[i][0];
    dest[i][1] = source[i][1];
    dest[i][2] = source[i][2];
    dest[i][3] = source[i][3];
  }
}

inline bool EqualMatrices3x4(ConstMatrix3x4 A, ConstMatrix3x4 B)
{
  //I'm sure there are faster implementations, but it doesn't matter _too_
  //much for it to be fast in the situation I will be using it in.
  for (int i=0; i < 3; ++i) 
    for (int j=0; j < 4; ++j)
      if (A[i][j] != B[i][j]) return false;
  return true;
}

Real Normalize3(Vect3 x);
Real Normalize4(Vect4 x);


//ConvertQuaternionToMatrix3x3 akes a normalized quaternion "input" and
//outputs a 3x3 matrix "dest", corresponding to the rotation produced
//by that quaternion.
//  The results are meainingless unless "input" is normalized
//ie:
//1.0f=input[0]*input[0]+input[1]*input[1]+input[2]*input[2]+input[3]*input[3]

void
ConvertQuaternionToMatrix3x3(ConstQuaternion input, Matrix3x3 dest);


void
BinaryWriteMatrix3x4(fstream& binary_file,
                     ConstMatrix3x4 M);

void
BinaryReadMatrix3x4(fstream& binary_file,
                    ConstMatrix3x4 M);

void
WriteMatrix3x4ToFile(ConstMatrix3x4 M,
                     char const *filename);

void
ReadTransFile(char const *filename,
              Matrix3x4 dest);

void
FindCentroid(short numpoints, ConstVect3 *aPoints, Vect3 dest);

} //namespace vect_3d

#endif //#ifndef _VECT_3D_H
