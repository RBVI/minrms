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

#include <cassert>
using namespace std;

#include <global_utils.h>
#include <simple_numeric_utils.h>
#include "distance_metric.h"
#include "dynNW.h"



//Hopefully, the following line includes an RCS-version-string in the binary.
static char rcs_id[] = "@(#)$Header: /usr/local/src/minrms/lib/dyn/RCS/dynNW.cc,v 3.4 2005/02/24 00:48:21 conrad Exp $";


using namespace minrms;




void
DynNW::Dealloc()
{
  delete [] aWindowSize;
  delete [] aNumRes;
  delete [] aTableSize;
  delete [] aAij;
  delete [] aPrevPointers;
  #ifdef NW_PRECOMPUTE_CIJ
  delete [] aCij;
  #endif //#ifdef NW_PRECOMPUTE_CIJ
}





void
DynNW::Alloc(int const *aSetNumResidues,
             int set_n_min)
{
  int const num_sequences = 2;

  assert(aSetNumResidues);
  DEBUG_MSG(DBG_DYNNW, "number of sequences = " << num_sequences);

  // First, calculate the range over which n can vary.
  // (this is 0 to the maximum number of matches possible)
  max_num_matches = MIN( aSetNumResidues[0],
                         aSetNumResidues[1] );
  DEBUG_MSG(DBG_DYNNW, "max_num_matches = " << max_num_matches );
  if ((set_n_min < 1) || (set_n_min > max_num_matches))
    ERR("Minimum and/or maximum number of matches allowed does not "
        "make sense.  The range of allowed matches must be >=1, and <="
        << max_num_matches
        << "\n"
        "Presently min_n = " << set_n_min);

  n_min = set_n_min;

  aWindowSize  = new int [num_sequences];
  aNumRes = new int [num_sequences];
  aTableSize = new int [num_sequences];
  if ((!aNumRes) || (! aTableSize) || (! aWindowSize))
  {
    ERR("Error  \"" << __FILE__ << "\":" << __FILE__ << "\n"
        "Failed to allocate memory.\n"
        "Please report this error to the developer.");
  }
  for (short seq=0; seq< num_sequences; ++seq)
  {
    aNumRes[seq] = aSetNumResidues[seq];
    DEBUG_MSG(DBG_DYNNW, "num residues in seq " << seq+1 << " = "
              << aNumRes[seq]);
    aTableSize[seq] = aNumRes[seq] + 1;
  }

  aWindowSize[0] = aNumRes[0] - n_min + 1;
  aWindowSize[1] = aNumRes[1] - n_min + 1;
  DEBUG_MSG(DBG_DYNNW,
            "minimum number of matches: " << n_min << "\n"
            "window size for sequence #1: " << aWindowSize[0] << "\n"
            "window size for sequence #2: " << aWindowSize[1]);
  aAij = new Real [ aTableSize[0] * aTableSize[1] ];
  #ifdef NW_PRECOMPUTE_CIJ
  aCij = new Real [ aTableSize[0] * aTableSize[1] ];
  #endif //#ifdef NW_PRECOMPUTE_CIJ
  aPrevPointers = new PrevCellSpecifier [ aTableSize[0] * aTableSize[1] ];

  if ((! aAij)
      || (! aPrevPointers)
  #ifdef NW_PRECOMPUTE_CIJ
      || (! aCij)
  #endif //#ifdef NW_PRECOMPUTE_CIJ
      )
  {
    ERR("Error allocating memory for Needleman & Wunsch cost table.\n");
  }

  DEBUG_MSG(DBG_DYNNW,
            "Alocated Needleman & Wunsch Table.\n"
            "  &(A[0][0]) = " << aAij << "\n"
            "  &(A[" << aTableSize[0]-1 << "][" << aTableSize[1]-1 
            << "]) = " << aAij + aTableSize[0]*aTableSize[1] -1);

} //DynNW::Alloc()






//... All the function InitTable() does is fill in the table with (i+j)*G
//values (where i,j are the indices into the table and G=gap_penalty).
//Initializing the table saves us from having to complicate
//the dynamic programming code by checking for special
//"base cases" during the search.  I thought this way would be easier.
void
DynNW::InitTable(Real gap_penalty)
{
  DEBUG_MSG(DBG_DYNNW,
            "Initializing the special edge-section of the table.");

  for(int i = 0; i <= aNumRes[0]; ++i)
  {
    for(int j = 0; j <= aNumRes[1]; ++j)
    {
      //Initialize whole table with an impossible "warning"
      //value. The good cells will be overwritten with good
      //values.  Upon reading, we allways make sure we
      //are not reading uninitialized memory.
      aAij[ i*aTableSize[1]+j ] = -1.0f;
    } 
  }

  //to only initialize the "border cells" which will actually
  //be referenced.  These cells are indicated by a small square
  //  _
  // |_|
  //
  //in the figure below:
  //
  //         j ^
  //           :             a = windowSize[0] = m - Nmin
  //           :                  <--a->
  //           :................_._._._
  //           :              _|_|     | /|\
  //           :            _|_|       |  | 
  //           :          _|_|         |  b = windowSize[1] = n - Nmin
  //           :        _|_|          x|  |   (where "m" and "n" are the # of
  //           :      _|_|          x  |  |    residues in each molecule, and
  //           :    _|_|          x   _| \|/   Nmin = minumum # of matches made
  //           :  _|_|          x   _|_|       in the "windowing" optimization.
  //           :_|<-----b-----x   _|_| :       In this case m=11, n=14, Nmin=8)
  //           |_|          x a _|_|   :
  //           |_|        x   ||_|     :
  //         . |_|      x   _|v|       :
  //         . |_|    x   _|_|         :
  //         . |_|  x   _|_|           :
  //        j=1|_|x____|_|             :
  //        j=0|_|_|_|_|...............:......\ i
  //            i i                           /
  //            " " ...
  //            0 1
  //
  //The reason the "border cells" lie in this funny shape is because of
  //the "windowing" optimization.  (If windowing is not used, then
  //border cells just lie on the lower edge or left hand edge of the table,
  //that is all the cells with i=0, or j=0.)

  //First, initialize the value of A00
  aAij[0] = 0.0f;

  //Then we fill in the values of  A i,(i-a)
  //  (the lower diagonal line)
  for(int i = 1; i <= aNumRes[0]; ++i)
  {
    int j = i - aWindowSize[0];
    if (j < 0)
      j = 0;

    aAij[ i*aTableSize[1]+j ] = (i+j) * gap_penalty;
    aPrevPointers[i*aTableSize[1]+j] = PREV_PTR_UNINITIALIZED;
    DEBUG_MSG(DBG_DYNNW,
              "&A[" << i << "]["<< j << "] = "
              << &(aAij[ i*aTableSize[1]+j ])
              << " cell initialized.");
  } //for(int i = 1; i <= aNumRes[0]; ++i)

  //Then we fill in the values of  A (j-b),j
  //  (the upper diagonal line in the figure)
  for(int j = 1; j <= aNumRes[1]; ++j)
  {
    int i = j - aWindowSize[1];
    if (i < 0)
      i = 0;

    aAij[ i*aTableSize[1]+j ] = (i+j) * gap_penalty;
    aPrevPointers[i*aTableSize[1]+j] = PREV_PTR_UNINITIALIZED;

    DEBUG_MSG(DBG_DYNNW,
              "&A[" << i << "]["<< j << "] = "
              << &(aAij[ i*aTableSize[1]+j ])
              << " cell initialized.");
  } //for(int j = 1; j <= aNumRes[1]; ++j)

} //DynNW::InitTable()






void
DynNW::Solve(Biopolymer const& m1,
             Biopolymer const& m2,
             PairwiseAlignment& dest
             #ifdef NW_PRECOMPUTE_CIJ
             , bool lookup_cij
             #endif
             )
{
  assert(aAij);
  //first, find the gap penalty (this should have been set by a previous
  //call to InitTable())
  Real G = aAij[1];

  if (G == 0.0f)
    ERR("Error occured during the initial psuedo-Needleman & Wunsch\n"
        "  RMSD minimization stage.\n"
        "       Using a gap penalty of 0.0 will produce meaningless\n"
        "       alignments containing 0 matches.  Aborting...\n");
  if (G < 0.0f)
    ERR("Error occured during the initial psuedo-Needleman & Wunsch\n"
        "  RMSD minimization stage.\n"
        "         The \"gap penalty\" is supposed to signify the largest\n"
        "       acceptable squared-distance between two residues x 2.\n"
        "          The user requested"
        " a negative gap_penalty (" << G << ").\n"
        "       However the squared physical distance between\n"
        "       any two points should be positive.\n"
        "       Also note:\n"
        "          Unlike traditional sequence alignment, our object\n"
        "       is to minimize a quantity (RMSD), not maximize it,\n"
        "       so our gap penalty should be a positive value.\n"
        "       Please try again.\n"
        "       Aborting...");

  Real Cij; //distance between residues i and j
#ifdef NW_PRECOMPUTE_CIJ
  Real *pCij; //pointer into the array of distances
#endif //#ifdef NW_PRECOMPUTE_CIJ
  Real case_ij_match; //lowest cost if residues i and j match
  Real case_i_no_match;//lowest cost if residue i does not match
  Real case_j_no_match;//lowest cost if residue j does not match

  Real *pAij;   //pointer into aAij at current position (i,j)
  Real *pAi1_j; //pointer into aAij at offset position (i-1,j)
  Real *pAi_j1; //pointer into aAij at offset position (i,j-1)
  Real *pAi1_j1;//pointer into aAij at offset position (i-1,j-1)
  PrevCellSpecifier *pPrev;//A pointer into the "aPrevPointers" array
                           //at the current position(i,j)

  for(int i = 1; i <= aNumRes[0]; ++i)
  {
    //j varies from [j_lower_bound, j_upper_bound)
    int j, j_lower_bound, j_upper_bound;

    if (i <= aWindowSize[0])
      j_lower_bound = 1;
    else
      j_lower_bound = i - aWindowSize[0] + 1;

    if (i + aWindowSize[1] <= aNumRes[1])
      j_upper_bound = i + aWindowSize[1]-1;
    else
      j_upper_bound = aNumRes[1];

    pAij = aAij + i*aTableSize[1] + j_lower_bound;
    #ifdef NW_PRECOMPUTE_CIJ
    pCij = aCij + i*aTableSize[1] + j_lower_bound;
    #endif //#ifdef NW_PRECOMPUTE_CIJ
    pAi1_j = aAij + (i-1)*aTableSize[1] + j_lower_bound;
    pAi_j1 = aAij + i*aTableSize[1] + (j_lower_bound-1);
    pAi1_j1 = aAij + (i-1)*aTableSize[1] + (j_lower_bound-1);
    pPrev =  aPrevPointers + i*aTableSize[1] + j_lower_bound;

    // ------------ INNER LOOP -------------
    for(j = j_lower_bound; j <= j_upper_bound; ++j)
    {

      #ifdef NW_PRECOMPUTE_CIJ
      if (lookup_cij)
      {
        //First, I make sure I have not screwed up and am not reading
        //what would normally be unitialized memory in the table
        //("-1.0f" is an impossible value, and is used to signal an error.
        // With the DEBUG flag turned off, the memory marked -1.0 would
        // normally be uninitialized.)
        assert(*pCij != -1.0f);
        Cij = *pCij; //do a table lookup instead of recomputing the distance.
      }
      else
      #endif //#ifdef NW_PRECOMPUTE_CIJ
        Cij = DistanceMetric::
          InlinedSquaredDistanceBetweenIandJ(m1.begin() + i-1,
                                             m2.begin() + j-1);

      //Now, make sure I have not screwed up and am not reading
      //what would normally be unitialized memory in the table
      //("-1.0f" is an impossible value, and is used to signal an error.
      // With the DEBUG flag turned off, the memory marked -1.0 would
      // normally be uninitialized.)
      //So the objective of the next three asserts is to verify that
      //the parts of the table I will be reading from contain values
      //computed from previous iterations.
      assert(*pAi1_j  != -1.0f);
      assert(*pAi_j1  != -1.0f);
      assert(*pAi1_j1 != -1.0f);

      case_ij_match   = *pAi1_j1 + Cij;
      case_i_no_match = *pAi1_j + G;
      case_j_no_match = *pAi_j1 + G;

      *pAij = case_i_no_match;
      *pPrev = CASE_I1_J;
      if (case_j_no_match < *pAij) {
        *pAij = case_j_no_match;
        *pPrev = CASE_I_J1;
      }                
      if (case_ij_match < *pAij) {
        *pAij = case_ij_match;
        *pPrev = CASE_I1_J1;
      }


      #if 0
      DEBUG_MSG(DBG_DYNNW,
                "&Aij = " << pAij
                << "; A["        << i << "]["<< j << "] = " << *pAij
                << "(choices: "
                << case_ij_match << ","
                << case_i_no_match << ","
                << case_j_no_match << ")\n"
                "previous cell at ["
                << i-(*pPrev != CASE_I_J1) << ","
                << j-(*pPrev != CASE_I1_J) << "]");
      #endif //#if 0


      ++pAij;

      #ifdef NW_PRECOMPUTE_CIJ
      ++pCij;
      #endif //#ifdef NW_PRECOMPUTE_CIJ

      ++pAi1_j;
      ++pAi_j1;
      ++pAi1_j1;
      ++pPrev;
    } //for(j = j_lower_bound+1; j < j_upper_bound;  ++j)
  } //for(i = 0; i < aTableSize[0]; ++i)

  ExportSolution(dest);

} // DynNW::Solve()




#ifdef NW_PRECOMPUTE_CIJ
void DynNW::ComputePairwiseMatchCosts(Biopolymer const& m1,
                                      Biopolymer const& m2,
                                      Real *pMinCost,
                                      Real *pMaxCost)
{
  for(int i = 0; i <= aNumRes[0]; ++i)
  {
    for(int j = 0; j <= aNumRes[1]; ++j)
    {
      //Initialize whole table with an impossible "warning"
      //value. The good cells will be overwritten with good
      //values.  Upon reading, we allways make sure we
      //are not reading uninitialized memory.
      aCij[ i*aTableSize[1]+j ] = -1.0f;
    } 
  }
  Real min_cost, max_cost;
  bool first_iter = true;
  for(int i = 1; i <= aNumRes[0]; ++i)
  {
    //j varies from [j_lower_bound, j_upper_bound)
    int j, j_lower_bound, j_upper_bound;

    if (i <= aWindowSize[0])
      j_lower_bound = 1;
    else
      j_lower_bound = i - aWindowSize[0] + 1;

    if (i + aWindowSize[1] <= aNumRes[1])
      j_upper_bound = i + aWindowSize[1]-1;
    else
      j_upper_bound = aNumRes[1];

    Real *pCij //pointer into the array of squared distances
      = aCij + i*aTableSize[1] + j_lower_bound;

    // ------------ INNER LOOP -------------
    for(j = j_lower_bound; j <= j_upper_bound; ++j)
    {
      Real Cij = DistanceMetric::
        InlinedSquaredDistanceBetweenIandJ(m1.begin() + i-1,
                                           m2.begin() + j-1);
      if (first_iter)
      {
        min_cost = max_cost = Cij;
        first_iter = false;
      }
      else if (Cij < min_cost)
        min_cost = Cij;
      else if (Cij > max_cost)
      {
        max_cost = Cij;
        if ( max_cost > DistanceMetric::UPPER_BOUND )
          ERR("Error: \"" << __FILE__ << "\":" << __LINE__ << "\n"
              "       Somehow the pairwise-cost (usually distance, in angstroms,\n"
              "       squared) of matching residues "
              << (m1.begin() + i-1)->id <<" and "<< (m2.begin() + j-1)->id <<"\n"
              "       from the two proteins has exceeded "<< DistanceMetric::UPPER_BOUND << "\n"
              "       It seems very unlikely that they could be that far away.\n"
              "       Aborting\n");
      }
      //Finally, store the squared-distance between these two residues
      //in the table.
      DEBUG_MSG(DBG_LOOKUP_TABLE_NW,
                "Distance lookup table &aCij[" << i << "]["<< j << "] = "
                << pCij << ", and aCij[" << i << "]["<< j << "] = " << Cij);
      *pCij = Cij;
      ++pCij;
    } // for(j = j_lower_bound; j <= j_upper_bound; ++j)
  } // for(int i = 1; i <= aNumRes[0]; ++i)
  if (pMinCost)
    *pMinCost = min_cost;
  if (pMaxCost)
    *pMaxCost = max_cost;
} // DynNW::ComputePairwiseMatchCosts()
#endif //#ifdef NW_PRECOMPUTE_CIJ



void
DynNW::ExportSolution(PairwiseAlignment& dest)
{
  //Scan through it once, just to see how many matches are in the alignment.
  int num_matches = 0;
  int cur_i = aNumRes[0];
  int cur_j = aNumRes[1];
  //trace backwards through the alignment and count the matches
  while ((cur_i > 0) && (cur_j > 0))
  {
    assert((cur_i>0) && (cur_j>0));
    assert((cur_i <= aNumRes[0]) && (cur_j <= aNumRes[1]));
    if (IsMatch(cur_i, cur_j))
      ++num_matches;

    //cerr << "ExportSoln: GetPrevCell("<<cur_i<<","<<cur_j<<") = (";
    GetPrevCell(cur_i, cur_j, &cur_i, &cur_j);
    //cerr << cur_i<<","<<cur_j<<")" << endl;
  }

  dest.resize(num_matches);

  //Now store the matches in the "dest" argument
  cur_i = aNumRes[0];
  cur_j = aNumRes[1];
  int cur_n = num_matches;
  //trace backwards through the alignment
  while (cur_n > 0)
  {
    assert((cur_i > 0) && (cur_j > 0));
    assert((cur_i <= aNumRes[0]) &&
           (cur_j <= aNumRes[1]) &&
           (cur_n <= num_matches));
    if (IsMatch(cur_i, cur_j))
    {
      assert((cur_i>0) && (cur_j>0));

      DEBUG_MSG(DBG_PAIR_ALIGN_DATA_STRUCT,
                cur_n << "th match at i = " <<
                cur_i << ", j = " <<
                cur_j);
      dest[cur_n-1][0] = cur_i-1;
      dest[cur_n-1][1] = cur_j-1;
      --cur_n;
    }

    GetPrevCell(cur_i, cur_j, &cur_i, &cur_j);

    assert((cur_i <= aNumRes[0]) && (cur_j <= aNumRes[1]));
  } // while (cur_n > 0)
  assert(cur_n == 0);
} // DynNW::ExportSolution()

