#include <cassert>
#include <cmath>
using namespace std;

#include <global_utils.h>
#include <simple_numeric_utils.h>
#include <vect_3d.h>
#include <eigen.h>
#include "superimpose.h"



using namespace vect_3d;



Real
Superimpose::FindTransformMinRMSD(short numpoints,
                                ConstVect3 *aFixed, ConstVect3 *aMoveable,
                                Matrix3x4 dest)
{

  assert(aFixed && aMoveable);

  //assert(numpoints>=3);//Need at least 3 points to do calculation accuratley
                       //(This dosen't insure everything's okay, but we might
                       // as well check for it anyway.)

  //The paper describing the method (see "superimpose.h")
  //says that the optimal translation between the two sets of
  //points has their centroids overlap.  So, for each sequence, subtract off
  //its centroid from each point, and store the result in the global arrays
  //aFixedNew and aMoveableNew.  (I made these arrays global instead
  //of local to avoid calling new() and delete())
  Vect3 centFixed, centMoveable;
  FindCentroid(numpoints, aFixed,    centFixed);
  FindCentroid(numpoints, aMoveable, centMoveable);
  DEBUG_MSG(DBG_ROT, "Centroids: fixed=["
            << centFixed[0] << ","
            << centFixed[1] << ","
            << centFixed[2] << "]"
            << ",  moveable=["
            << centMoveable[0] << ","
            << centMoveable[1] << ","
            << centMoveable[2] << "]");

  int a; //used as an index that loops over 0..numpoints-1
  assert(aFixedNew && aMoveableNew);
  for (a=0; a < numpoints; ++a)
  {
    SubtractVect3(aFixed[a],    centFixed,     aFixedNew[a]);
    SubtractVect3(aMoveable[a], centMoveable,  aMoveableNew[a]);
  }
  //Okay, cool.  Now I will use aFixedNew, and aMoveableNew in all 
  //future calculations.

  int i, j; //indices into a 3x3 matrix. loop over 0..2;
  Real M[3][3];
  ResetToZero3x3(M); //fill with zeros
  
  for(a=0; a<numpoints; ++a)
  {
    for(i=0; i<3; ++i)    {
      for(j=0; j<3; ++j)  {
        M[i][j] += aMoveableNew[a][i] * aFixedNew[a][j];
      }
    }
  }

  DEBUG_MSG(DBG_ROT, " -----------------------\n"
            << "    [" << M[0][0] << "," << M[0][1] << "," << M[0][2] << "]\n"
            << "M = [" << M[1][0] << "," << M[1][1] << "," << M[1][2] << "]\n"
            << "    [" << M[2][0] << "," << M[2][1] << "," << M[2][2] << "]");

  Real traceM = FindTrace3x3(M);
    
  Real Q[3][3];
  for(i=0; i<3; ++i)  {
    for(j=0; j<3; ++j)  {
      Q[i][j] = M[i][j] + M[j][i] - ((i == j) ? (2 * traceM) : 0.0f);
    }
  }

  DEBUG_MSG(DBG_ROT, " -----------------------\n"
            << "    [" << Q[0][0] << "," << Q[0][1] << "," << Q[0][2] << "]\n"
            << "Q = [" << Q[1][0] << "," << Q[1][1] << "," << Q[1][2] << "]\n"
            << "    [" << Q[2][0] << "," << Q[2][1] << "," << Q[2][2] << "]");

  Real V[3];
  V[0] = M[1][2] - M[2][1];
  V[1] = M[2][0] - M[0][2];
  V[2] = M[0][1] - M[1][0];

  DEBUG_MSG(DBG_ROT, "V = [" << V[0] << "," << V[1] << "," << V[2] << "]");

  //  Commenting out next line. P is a static member.
  // Real **P = new Real* [4];
  // (P is declared in find_rotation.h, and should be previously allocated)

  //Now copy Q into the 3x3 upper-right corner of P (P is a 4x4 pointer array)
  //(Actually, P is a 5x5 pointer array.
  // The 0 indices aren't used. Indexing starts at 1.
  // I'm not sure why the authors of Numerical Recipies in C did it this way.)
  for(i=0; i<3; ++i)  {
    for(j=0; j<3; ++j)  {
      P[i+1][j+1] = Q[i][j];
    }
  }

  // Now fill in the remaining lower-right row and collumn of P
  P[1][4] = P[4][1] = V[0];
  P[2][4] = P[4][2] = V[1];
  P[3][4] = P[4][3] = V[2];
  P[4][4] = 0.0f;

  DEBUG_MSG(DBG_ROT, " ----------------------------------\n" 
        << "    ["<<P[1][1]<<","<<P[1][2]<<","<<P[1][3]<<","<<P[1][4]<<"]\n"
        << "P = ["<<P[2][1]<<","<<P[2][2]<<","<<P[2][3]<<","<<P[2][4]<<"]\n"
        << "    ["<<P[3][1]<<","<<P[3][2]<<","<<P[3][3]<<","<<P[3][4]<<"]\n"
        << "    ["<<P[4][1]<<","<<P[4][2]<<","<<P[4][3]<<","<<P[4][4]<<"]");

  assert(dest);
  assert(dest[0] && dest[1] && dest[2]);



#if 0
  // *******************************************************************
  //To calculate the maximum eigenvalue/eigenvector, I use the functions
  //from sections 11.2 and 11.3 of "Numerical Recipies in C".
  //These functions calculate all the eigenvalues/eigenvectors, which is 
  //a lot more information than what we need, and there is a bit more
  //additional overhead doing it this way because it's optimized for general
  //nxn matrices (where n is large).  I'll change it if profiling reveals 
  //that there is a noticeable slowdown here.

  //The next two variables are required by "tred2()" and "tqli()"
  Real d[5];
  Real e[5];
  //Reduce P to "tri-diagonal" form (Num. Rec. in C, section 11.2)
  tred2(P, 4, d, e);

  //We're required to set eigenvects to the identity matrix before calling
  //tqli()
  for (i = 1; i<=4; ++i)
  {
    for (j = 1; j<=4; ++j)
      eigenvects[i][j] = ((i==j) ? 1.0f : 0.0f);
  }
  /* no, okay try this instead: */
  //We're required to set "eigenvects" to the matrix output by tred2
  //before calling tqli()
  for (i = 1; i<=4; ++i)
  {
    for (j = 1; j<=4; ++j) {
      if (i==j)
        eigenvects[i][j] = d[i];
      else if (i==j+1)
        eigenvects[i][j] = e[i];
      else if (j==i+1)
        eigenvects[i][j] = e[j];
      else 
        eigenvects[i][j] = 0.0f;
    }
  }

  //Find the eignevalues/vectors of the tri-diagonal matrix specified
  //by d[] and e[]. Store the eigenvalues in d, the eigenvectors in eigenvects.
  tqli(d, e, 4, eigenvects);
#endif //woops, confusing.  looks like I was wrong.  Instead use:




  Real d[5];
  int num_iter_dummy;//Not used.It's required that we pass an int * to jacobi()
  jacobi(P, 4, d, eigenvects, &num_iter_dummy);

  DEBUG_MSG(DBG_ROT, " ----------------------------------\n" 
        << "             ["<<eigenvects[1][1]<<","<<eigenvects[1][2]<<","<<eigenvects[1][3]<<","<<eigenvects[1][4]<<"]\n"
        << "eigenvects = ["<<eigenvects[2][1]<<","<<eigenvects[2][2]<<","<<eigenvects[2][3]<<","<<eigenvects[2][4]<<"]\n"
        << "             ["<<eigenvects[3][1]<<","<<eigenvects[3][2]<<","<<eigenvects[3][3]<<","<<eigenvects[3][4]<<"]\n"
        << "             ["<<eigenvects[4][1]<<","<<eigenvects[4][2]<<","<<eigenvects[4][3]<<","<<eigenvects[4][4]<<"]");

  //The eigenvalues returned by tqli/jacobi don't seemed to be in sorted order,
  //so I'll have to search for the largest one.
  int whichIsLargest = 1;
  Real largestEigenvalue = d[1];
  for (i=2; i<=4; ++i)
  {
    if (d[i] > largestEigenvalue)
    {
      largestEigenvalue = d[i];
      whichIsLargest = i;
    }
  }
  //Now, find the eigenvector corresponding to the largest eigenvalue.
  Quaternion  p; //(p corresponds to the optimal rotation)
  p[0] = eigenvects[1][ whichIsLargest ];
  p[1] = eigenvects[2][ whichIsLargest ];
  p[2] = eigenvects[3][ whichIsLargest ];
  p[3] = eigenvects[4][ whichIsLargest ];
  Normalize4(p);

  DEBUG_MSG(DBG_ROT, " -----------------------\n"
            << "largest eigenvalue = " << largestEigenvalue << "\n"
            << "corresponding eigenvector (quaternion of best angle) = ["
            << p[0] << ","
            << p[1] << ","
            << p[2] << ","
            << p[3] << "]");

  //**************finished calculating eigenstuff*************************

  Matrix3x3 R;

  //Now find the rotational part of the transformation matrix.
  //Store the resulting rotation matrix in "R"
  ConvertQuaternionToMatrix3x3(p, R);

  DEBUG_MSG(DBG_ROT, "         -----\n"
            << " associated rotation matrix is:\n"
            << "    [" << R[0][0] << "," << R[0][1] << "," << R[0][2] << "]\n"
            << "R = [" << R[1][0] << "," << R[1][1] << "," << R[1][2] << "]\n"
            << "    [" << R[2][0] << "," << R[2][1] << "," << R[2][2] << "]");

  //Now, find the translational part of the transform matrix
  //This is the vector
  //
  //  c1[] - R[][]*c2[]
  //       (where R[][] is the rotation matrix {which we just calculated
  //                                            and stored in "R"}
  //         c1=centroid of first (fixed) sequence (called "centFixed")
  //         c2=centroid of first (moveable) sequence (called "centMoveable")
  //
  // You derive this result by considering that the total transformation is:
  // You want to rotate X around its centroid, and center the result around
  // the centroid of the other sequence.  This looks like:
  //
  // Xnew[] = R[][]*(X-c2) + c1 (simplification yeilds R[][]*X + c1-R[][]*c2)
  //
  //  ..Where X=position of a vector expressed in the second ("moveable")
  // sequences coordinate system,

  Vect3 translation;
  Mult_Mat3x3_by_Vect3(R, centMoveable, translation);
  SubtractVect3(centFixed, translation, translation);

  //Now, copy the Rotation and translation into "dest"
  for(i=0; i<3; ++i) {
    dest[i][0] = R[i][0];
    dest[i][1] = R[i][1];
    dest[i][2] = R[i][2];
    dest[i][3] = translation[i];
  }
  
  //Now, calculate the smallest rms_distance (as derived in Eq(23) of Diamond)
  Real E0 = 0.0f;
  for(a=0; a<numpoints; ++a)  {
    for(i=0; i<3; ++i)    {
      E0 += SQR(aFixedNew[a][i] - aMoveableNew[a][i]);
    }
  }
  //The quantity that Diamond calculates in Eq(23) is the "minimum of E"
  //where "E" is the sum of the squares of the distances.
  Real min_sum_squared_distance = E0 - 2*largestEigenvalue;


  if (min_sum_squared_distance < 0.0f)
  {
    //When minimizing identical structures, min_sum_squared_distance should 
    //be 0.0, but sometimes I get a very small negative value from roundoff.
    //So reset it to zero (but print out an error message to the user)
    //  cerr << "--err?:" <<__FILE__<<": sum-squared = "
    //  <<min_sum_squared_distance << "  (setting to 0.0f) --" << endl;
    min_sum_squared_distance = 0.0f;
  }

  DEBUG_MSG(DBG_ROT, " -----------------------\n"
            << "min-sum-squared-distance = "
            << min_sum_squared_distance << ",   min_rms = "
            << sqrt( min_sum_squared_distance / numpoints ));
  
  //Now return the square-root-of-the-average
  return sqrt( min_sum_squared_distance / numpoints );

} // FindTransformMinRMS()



void
Superimpose::AllocGlobalVariables(int maxNumPoints)
{
  //I would like to include the following lines, but this makes purify think
  //I am freeing unallocate memory (I can't see the error)
  //Please include later or there will be a memory leak.
  //  if (aFixedNew)
  //      delete [] aFixedNew;
  assert(aFixedNew == NULL);
  aFixedNew = new Vect3 [maxNumPoints];
  //I would like to include the following lines, but this makes purify think
  //I am freeing unallocate memory (I can't see the error)
  //Please include later or there will be a memory leak.
  //  if (aMoveableNew)
  //    delete [] aMoveableNew;
  assert(aMoveableNew == NULL);
  aMoveableNew = new Vect3 [maxNumPoints];
  if ((!aFixedNew) || (!aMoveableNew))
    ERR(__FILE__<<":"<<__LINE__<<": can't alloc mem.");
  if (!P || !eigenvects)
  {
    P = new Real * [5];
    eigenvects = new Real * [5];
    if ((!P) || (!eigenvects))
      ERR(__FILE__<<":"<<__LINE__<<": can't alloc mem.");
    for(int i=0; i<=4; ++i)
    {
      P[i] = new Real [5];
      eigenvects[i] = new Real [5];
      if ((!P[i]) || (!eigenvects[i]))
        ERR(__FILE__<<":"<<__LINE__<<": can't alloc mem.");
    }
  }
}


// Call the next function after you are through using FindTransformMinRMS()
void
Superimpose::DeallocGlobalVariables()
{
  if (aFixedNew) delete [] aFixedNew;
  if (aMoveableNew) delete [] aMoveableNew;
  if (P)
  {
    for(int i=0; i<=4; ++i) {
      if (P[i]) delete [] P[i];
    }
    delete [] P;
  }
  if (eigenvects)
  {
    for(int i=0; i<=4; ++i) {
      if (eigenvects[i]) delete [] eigenvects[i];
    }
    delete [] eigenvects;
  }
}

