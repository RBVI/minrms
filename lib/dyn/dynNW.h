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

#ifndef _DYNNW_H
#define _DYNNW_H

#define NW_PRECOMPUTE_CIJ


#include <cassert>
using namespace std;

#include <pair_alignment.h>


namespace minrms {


class DynNW
{
  DynNW(DynNW const& ); //no implicit copying
  DynNW& operator = (DynNW const& ); //no explicit copying

public:

  DynNW(int numRes1, int numRes2, int set_min_num_matches=0)
  {
    int aNumResidues[2];
    aNumResidues[0] = numRes1;
    aNumResidues[1] = numRes2;
    Alloc(aNumResidues, set_min_num_matches);
  }

  ~DynNW()
  { Dealloc(); }

  //Call this next function once before invoking Solve().
  //This is where you specify the gap-penalty you want to use.

  void InitTable(Real gap_penalty//Gap penalty you want to use. This parameter
                                  //is simply the maximum permissable increase
                                  //in total sum-squared-distance when
                                  //increasing the number of matches made
                                  //by one, the quantity then divided by 2.
                                  //  In most cases, this is the same
                                  //as the square of the maximum permissable
                                  //distance between any pair of residues
                                  //(divided by 2).
                                  //  It has units of distance-squared.
                 );

  //Solve() generates a structural alignment between molecules m1 and m2, (and
  // these molecules should match the dimensions specified earlier in Alloc())
  //Solve() can be called multiple times without having to call
  //Alloc(), InitTable(), and Dealloc() each time.
  void
  Solve(Biopolymer const& m1, //input the molecules you want to align
        Biopolymer const& m2,       
        PairwiseAlignment& dest   //Data structure where the resulting
                                  //alignment/solution will be stored.

        #ifdef NW_PRECOMPUTE_CIJ
        //Set the "lookup_cij" parameter to true only if you have allready
        //invoked ComputePairwiseDistances() (see below).
        //This saves it from having to recompute inter-residue distances
        //over again.
        ,bool lookup_cij = false
        #endif //#ifdef NW_PRECOMPUTE_CIJ
        );




  #ifdef NW_PRECOMPUTE_CIJ
  //ComputePairwiseDistances() must be called once before
  // Solve(lookup_cij=true) is invoked.  ComputePairwiseMatchCosts()
  // precomputes the distances between all pairs
  // of residues from opposite proteins, and stores them in
  // a table for later use.
  // (If "pMinCost" & "pMaxCost" are non-NULL, the minimum and maximum distance
  //  btwn pairs are returned to the caller by storing them in these addresses)
  void ComputePairwiseMatchCosts(Biopolymer const& m1,
                                 Biopolymer const& m2,
                                 Real *pMinCost = NULL,
                                 Real *pMaxCost = NULL);
  #endif //#ifdef NW_PRECOMPUTE_CIJ




private:
  Real  *aAij;
  enum PrevCellSpecifier {CASE_I1_J1, //previous cell is i-1 j-1
                          CASE_I1_J,  //previous cell is i-1 j
                          CASE_I_J1,  //previous cell is i   j-1
                          PREV_PTR_UNINITIALIZED}; //previous cell is 0,0

  PrevCellSpecifier *aPrevPointers; //(backtrace to find the alignment)

  #ifdef NW_PRECOMPUTE_CIJ
  Real  *aCij;
  #endif

  int     max_num_matches;
  int    *aNumRes;
  int    *aTableSize;
  int    *aWindowSize;
  int     n_min;

  void Alloc(int const *aSetNumResidues,
             int set_n_min);
  void Dealloc();


  inline void GetPrevCell(int i, int j, int *pPrevI, int *pPrevJ)
  {
    switch(aPrevPointers[i*aTableSize[1]+j])
    {
    case CASE_I1_J1:
      *pPrevI = i-1;
      *pPrevJ = j-1;
      break;
    case CASE_I1_J:
      *pPrevI = i-1;
      *pPrevJ = j;
      break;
    case CASE_I_J1:
      *pPrevI = i;
      *pPrevJ = j-1;
      break;
    case PREV_PTR_UNINITIALIZED:
      *pPrevI = 0;
      *pPrevJ = 0;
      break;
    default:
      assert(0);
      break;
    }
  }

  inline bool IsMatch(int i, int j)
  { return (aPrevPointers[i*aTableSize[1]+j] == CASE_I1_J1); }


  //ExportSolution takes a reference to a
  //PairwiseAlignment variable, "dest", where
  //it will store the solution to the most recent
  //alignment calculated by Solve().
  void ExportSolution(PairwiseAlignment& dest);


}; //class DynNW


} //namespace minrms

#endif // #ifndef _DYNNW_H
