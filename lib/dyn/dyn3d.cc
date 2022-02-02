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

#include <cstdio>  //needed for "sprintf()"
#include <fstream> //needed for ofstream in Dyn3d::DisplayInMidas()
using namespace std;


#include <global_utils.h>
#include <simple_numeric_utils.h>
#include <biopolymer.h>
#include <pair_alignment.h> //needed because
                            //DisplayUsingMidas() accesses
                            //PairwiseAlignment::aColorNames
#include "distance_metric.h"  //needed by Dyn3d::FillLookupTable()
#include "dyn3d.h"


//Hopefully, the following line stores an RCS-version-string in the binary.
static char rcs_id[] = "@(#)$Header: /usr/local/src/minrms/lib/dyn/RCS/dyn3d.cc,v 3.7 2002/11/18 17:23:02 conrad Exp $";


using namespace minrms;


const Dyn3d::NLayerTable::PrevCellSpecifier
  Dyn3d::NLayerTable::CASE_N1_I1_J1 = 0;
const Dyn3d::NLayerTable::PrevCellSpecifier
  Dyn3d::NLayerTable::CASE_N_I1_J   = 1;
const Dyn3d::NLayerTable::PrevCellSpecifier
  Dyn3d::NLayerTable::CASE_N_I_J1   = 2;



int Dyn3d::NLayerTable::GetPrevCell(int n, int i, int j, 
                                    int which_molecule) const
{
  //assert()s below make sure n,i,j refers to a vallid cell in the table.
  assert(n>=0); assert(n <= max_num_matches);
  assert(i>=n); assert(i <= numResA);
  assert(j>=n); assert(j <= numResB);

  PrevCellSpecifier c = aaaTable[n][i][j];

  //The assert()s below make sure the prev-cell pointer is not
  //not filled with garbage, but actually points to a meaningful cell.
  assert(! ((c == CASE_N_I1_J) && (i == n)));
  assert(! ((c == CASE_N_I_J1) && (j == n)));

  switch(which_molecule)
  {
  case 0:
    return ((c == CASE_N_I_J1) ? i : (i-1));
  case 1:
    return ((c == CASE_N_I1_J) ? j : (j-1));
  default:
    assert(0);
    return -1;//need to put something here to avoid getting a compiler warning.
  }
} // Dyn3d::NLayerTable::GetPrevCell()





void
Dyn3d::Solve()
{
  int n; //for looping over number of matches

  // The indexing is weird, because zero is included, ...since it is
  // possible to have an alignment with zero matches,
  // and this needs to be a base-case for the recursion.
  // So the jth elemeht of the table is table[j], not table[j-1].

  //First, initialize pathelogical cases for when n = 0
  //loop over all n=0 and initialize
  InitTable_IncrN();

  for (n=2; n<=n_max; ++n)
  {
        // The next function does the dirty work...
        // It calculates the optimal alignment for n matches,
        // (it runs in time O(n^2) and can take a fairly long time
        // (Eg: 7-sec on a 233MHz Dec, for a sequence w/ 4500 residues)
        Solve_IncrN(n);

        cout
          << "n="
          << n
          << ",  min_cost="
          << " "
          << ((pApply2Results)
                  ? 
                  (*pApply2Results)(aOutcomes[n-1].cost, n)
                  :
                  aOutcomes[n-1].cost)
          //  << ",  wordst_pair = (i=" << worst_i << ",j=" << worst_j << ",Cij="
          //  << worst_cost
          //  << ")"
          << endl;
        
  } // for (n = 1 ... loop over n
} // Dyn3d::Solve()



void
Dyn3d::InitTable_IncrN()
{
  int i; //for looping over residues in the  first sequence
  int j; //for looping over residues in the second sequence

  assert((windowSizeA <= aNumResidues[0]) &&
         (windowSizeB <= aNumResidues[1]));
  //initialize the special base cases for n=1 (no previous matches possible)
  for (i=1; i<=windowSizeA; ++i)
  {
    for (j=1; j<=windowSizeB; ++j)
        {
          prev.SetPrevCell(1,i,j, NLayerTable::CASE_N1_I1_J1);
          Real best_cost =
                Cij[ (i-1)*(aNumResidues[1]) + j-1 ];

          if ((i>1) && (costs.Get(1,i-1,j) <= best_cost))
          {
                best_cost = costs.Get(1,i-1,j);
                prev.SetPrevCell(1,i,j, NLayerTable::CASE_N_I1_J);
          }
          if ((j>1) && (costs.Get(1,i,j-1) <= best_cost))
          {
                best_cost = costs.Get(1,i,j-1);
                prev.SetPrevCell(1,i,j, NLayerTable::CASE_N_I_J1);
          }
          costs.Set(1,i,j,best_cost);

          DEBUG_MSG(DBG_DYN3D, "initializing(1) cell[1]["
                                << i << "]["
                                << j << "]  prev=["
                                << prev.GetPrevCell(1,i,j,0) << "]["
                                << prev.GetPrevCell(1,i,j,1) << "] "
                                << "cost=" << best_cost
                                );

    } //for (j=1; j<=aNumResidues[1]; ++j)
  } //for (i=1; i<=aNumResidues[0]; ++i)
  assert( (windowSizeA <= aNumResidues[0]) &&
          (windowSizeB <= aNumResidues[1]) );
  aOutcomes[0].cost = costs.Get(1,
                                windowSizeA,
                                windowSizeB);
} //Dyn3d::InitTable_IncrN()



void
Dyn3d::Solve_IncrN(int n)
{
  int i; //for looping over residues in the  first sequence
  int j; //for looping over residues in the second sequence
  int i_max = MIN(aNumResidues[0], (n + windowSizeA - 1)); //loop over i,j= n
  int j_max = MIN(aNumResidues[1], (n + windowSizeB - 1)); //to i_max,j_max

  // The indexing is weird, because zero is included, ...since it is
  // possible to have an alignment with zero matches,
  // and this needs to be a base-case for the recursion.
  // So the jth elemeht of the table is table[j], not table[j-1].

  //initialize special cases for when i = n and j>=n
  for (j=n; j<=j_max; ++j)
  {
        prev.SetPrevCell(n,n,j, NLayerTable::CASE_N1_I1_J1);
        Real best_cost = costs.Get(n-1, n-1, j-1) +
          Cij[ (n-1)*aNumResidues[1] + (j-1) ];

        //Note: The next case is only considered if (j != n)
        if ((j>n) && (costs.Get(n,n,j-1) <= best_cost))
        {
          best_cost = costs.Get(n,n,j-1);
          prev.SetPrevCell(n,n,j, NLayerTable::CASE_N_I_J1);
        }
        costs.Set(n,n,j,best_cost);

        DEBUG_MSG(DBG_DYN3D, "initializing(2) cell["
                          << n << "]["
                          << n << "]["
                          << j << "] prev=["
                          << prev.GetPrevCell(n,n,j,0) << "]["
                          << prev.GetPrevCell(n,n,j,1) << "] "
                          "cost=" << best_cost
                          );

  } //for (j=n; j<=aNumResidues[1]; ++j) {


  //initialize special case for when i > n and j=n
  for (i=n+1; i<=i_max; ++i)
  {
    prev.SetPrevCell(n,i,n, NLayerTable::CASE_N1_I1_J1);
    Real best_cost = costs.Get(n-1, i-1, n-1) +
      Cij[ (i-1)*aNumResidues[1] + (n-1) ];
    if ((i>n) && (costs.Get(n,i-1,n) <= best_cost))
    {
      best_cost = costs.Get(n,i-1,n);
      prev.SetPrevCell(n,i,n, NLayerTable::CASE_N_I1_J);
    }
    costs.Set(n,i,n,best_cost);

    DEBUG_MSG(DBG_DYN3D, "initializing(3) cell["
                          << n << "]["
                          << i << "]["
                          << n << "]  prev=["
                          << prev.GetPrevCell(n,i,n,0) << "]["
                          << prev.GetPrevCell(n,i,n,1) << "] "
                          "cost=" << best_cost
                          );

    //variables needed for the inner-loop:
    Real lowest_cost;

    //I will explicitly use pointer arithmatic.
    //The next variables point to the current entry(s) in the dynamic
    //programming arrays, which are in use for the current value of j.
    //The naming convention is illustrated by:
    //     "*pCosts_n_i_j1" is the lowest cost of matching n pairs of
    //of residues when considering only the first i residues from
    //sequence A, with the first j-minus-one residues from sequence B.
    //(Note: Consequently, I start out at the beginning of the loop
    //       with j = n+1  (and so (j-1) = n).
    Real *pCosts_n_i_j = costs.GetAddressOfNIJ(n, i, n+1);

    Real *pCosts_n1_i1_j1 = costs.GetAddressOfNIJ(n-1, i-1, n);
    Real *pCosts_n_i1_j = costs.GetAddressOfNIJ(n, i-1, n+1);
    Real *pCosts_n_i_j1 = costs.GetAddressOfNIJ(n, i, n);

    //pPrev_n_i_j stores an indicator (equivalent to a pointer)
    //that specifies which cell was the precursor to this cell.
    // { Note: the choices are:
    //   cell# (n-1,i-1,j-1), (n,i-1,j), or (n,i,j-1) }
    //   provided that this cell is (n,i,j) }
    //In other words, it indicates which of the three values was smallest
    //  A(n-1),(i-1),(j-1) + C(i-1),(j-1)
    //  An,(i-1),j
    //  An,i,(j-1)
    //...and in so doing leaves a trail of pointers that the
    //   user will use in order to reconstruct the alignment
    //   which generated the lowest value of An,i,j.
    //   (This is, of course, the ultimate goal.)

    NLayerTable::PrevCellSpecifier
      *pPrev_n_i_j = prev.GetAddressOfNIJ(n, i, n+1);
    Real *pCij_i_j = &(Cij[ (i-1)*aNumResidues[1] + ((n-1)+1) ]);


    //Inner-loop:
    for (j=n+1; j<=j_max; ++j)
    {
      // In case 1, the j-th residue from sequence B matches with
      // the i-th residue from sequence A.
      // (This is the only time we have to consider anything new.)
      // Calculate the incremental cost for matching these two residues,
      // and add it to the cost for matching (n-1) matches up to residues
      // (i-1) and (j-1).
      // (In a lot of case, the contribution to the total cost of
      //  matching these two residues is independent of the current
      //  alignment-so-far and has been calculated in advance
      //  and stored in a lookup-table.
      //   If this is the case, then the third argument to the cost-
      //  function will be completely ignored, the operation will be just
      //  a table lookup and hopefully will be very fast indeed.)
      lowest_cost =   *pCosts_n1_i1_j1 + *pCij_i_j;
      *pPrev_n_i_j = NLayerTable::CASE_N1_I1_J1;

      // In case 2, the residue from sequence B at position j matches
      // one of the residues from sequence A at position k, where k < i;
      //    If this is the case, then since residue i from sequence A
      // is not used, we can throw it away and consider the minimal
      // solution up to i-1 (n and j are the same).
      //    But things are easy than you might expect, because,
      // of the order in which we have filled in the table.  
      // Because we scan over all possible values of j before we increment
      // i, we have allready calculated the cost for this allignment.
      // It is in located in *pCosts_n_i1_j.
      if (*pCosts_n_i1_j <= lowest_cost) {
        *pPrev_n_i_j = NLayerTable::CASE_N_I1_J;
        lowest_cost = *pCosts_n_i1_j;
      }

      // In case 3, the residue from sequence B at position j does not
      // match with any residues from sequence A.
      // We can throw it away, and use whatever solution we allready
      // had for j-1 (keeping n and i the same).
      // If newly computed cost is lower, then update
      if (*pCosts_n_i_j1 <= lowest_cost) {
        *pPrev_n_i_j = NLayerTable::CASE_N_I_J1;
        lowest_cost = *pCosts_n_i_j1;
      }

      // now that we know what case this is, make the new choice
      // accordingly, and update the table at this(n,i,j) point
      // to reflect that choice.
      *pCosts_n_i_j = lowest_cost;
      ++pCosts_n_i_j;

      ++pCosts_n1_i1_j1;
      ++pCosts_n_i1_j;
      ++pCosts_n_i_j1;

      #ifdef DEBUG
      //Print some debugging information to the user (if enabled)
      //Update some variables used for debugging
      Real solutions[NLayerTable::CASE_N_I_J1 + 1];
      solutions[NLayerTable::CASE_N1_I1_J1]=*pCosts_n1_i1_j1 + *pCij_i_j;
      solutions[NLayerTable::CASE_N_I1_J] = *pCosts_n_i1_j;
      solutions[NLayerTable::CASE_N_I_J1] = *pCosts_n_i_j1;
      int prev_n, prev_i, prev_j;
      switch (*pPrev_n_i_j)
      {
      case NLayerTable::CASE_N1_I1_J1:
        DEBUG_MSG(DBG_DYN3D, "["
                             << n << "]["
                             << i << "]["
                             << j<< "] current entry stores a match");
        prev_n = n-1;
        prev_i = i-1;
        prev_j = j-1;
        break;
      case NLayerTable::CASE_N_I1_J:
        prev_n = n;
        prev_i = i-1;
        prev_j = j;
        break;
      case NLayerTable::CASE_N_I_J1:
        prev_n = n;
        prev_i = i;
        prev_j = j-1;
        break;
      default:
        assert(0);
        break;
      }
      DEBUG_MSG(DBG_DYN3D, "[" << n << "][" << i << "][" << j << "]  prev=["
                           << prev_n << "]["
                           << prev_i << "]["
                           << prev_j << "] "
                           "(cost="
                           << solutions[NLayerTable::CASE_N1_I1_J1] << ","
                           << solutions[NLayerTable::CASE_N_I1_J] << ","
                           << solutions[NLayerTable::CASE_N_I_J1] << ")");
      #endif // #ifdef DEBUG

      ++pPrev_n_i_j;
      ++pCij_i_j;

    } // for (j = 1 ... loop over j  (Inner-loop)
  } // for (i = 1 ... loop over i

  aOutcomes[n-1].cost = costs.Get(n, i_max, j_max);

}  //Dyn3d::Solve_IncrN()



void
Dyn3d::ExportSolution(int n,
                      PairwiseAlignment& dest)
{
  assert((n_min <= n) && (n <= max_num_matches));
  dest.resize(n);

  int cur_n = n;
  int cur_i = GetNumRes(0);
  int cur_j = GetNumRes(1);

  //trace backwards through the alignment
  while (cur_n > 0)
  {
    int new_i = GetPrevCell(cur_n,cur_i,cur_j, 0);
    int new_j = GetPrevCell(cur_n,cur_i,cur_j, 1);
    if (IsMatch(cur_n, cur_i, cur_j))
    {
      assert((cur_i>0) && (cur_j>0));
      //#define DEBUG_MSG( type, msg ) cerr <<" MSG: " << msg << endl
      DEBUG_MSG(DBG_PAIR_ALIGN_DATA_STRUCT,
                "match at cell [" << cur_n << "][" << cur_i << "][" << cur_j
                << "]" );
                
      dest[cur_n-1][0] = cur_i-1;
      dest[cur_n-1][1] = cur_j-1;
          
      --cur_n;
    }
    else {
      DEBUG_MSG(DBG_PAIR_ALIGN_DATA_STRUCT,
                "No match at cell [" << cur_n << "][" << cur_i << "]["
                << cur_j << "]" );
    }
    cur_i = new_i;
    cur_j = new_j;
    DEBUG_MSG(DBG_PAIR_ALIGN_DATA_STRUCT,
              "...previous case of this cell is ["
              << cur_n << "][" << cur_i << "][" << cur_j << "]" );
  } // (cur_n >= 0)
} // Dyn3d::ExportSolution()


Dyn3d::~Dyn3d()
{
  Dealloc();
} // ~Dyn3d()



// Alloc() allocates storage space for the the table, 
// and performs some initialization.
void
Dyn3d::Alloc(int const *aSetNumResidues,
             int set_n_min,
             int set_n_max)
{
  int const num_sequences = 2;

  assert(aSetNumResidues);
  //  num_sequences = set_num_sequences;
  DEBUG_MSG(DBG_ALLOC, "number of sequences = " << num_sequences);

  // First, calculate the range that n can vary over.
  // (this is 0 to the maximum number of matches possible)
  max_num_matches = MIN( aSetNumResidues[0],
                         aSetNumResidues[1] );
  DEBUG_MSG(DBG_ALLOC, "max_num_matches = " << max_num_matches );
  if ((set_n_min < 1) || (set_n_min > max_num_matches) ||
      (set_n_max < 1) || (set_n_max > max_num_matches) ||
      (set_n_max < set_n_min))
        ERR("Minimum and/or maximum number of matches allowed does not "
            "make sense.  The range of allowed matches must be >=1, and <="
            << max_num_matches
            << "\n"
            "Presently min_n = " << set_n_min 
            <<", max_n = " << set_n_max);

  n_min = set_n_min;
  n_max = set_n_max;

  windowSizeA = aSetNumResidues[0] - n_min + 1;
  windowSizeB = aSetNumResidues[1] - n_min + 1;

  aNumResidues = new int [num_sequences];
  if (! aNumResidues) {
    ERR("dyn.cc:Dyn3d::Alloc()  failed to allocate memory.");
  }
  for (int seq=0; seq<num_sequences; ++seq)
  {
    aNumResidues[seq] = aSetNumResidues[seq];
    DEBUG_MSG(DBG_ALLOC, "num residues in seq " << seq+1 << " = "
              << aNumResidues[seq]);
  }


  prev.Alloc( aSetNumResidues[0], aSetNumResidues[1],
              set_n_min, set_n_max );
  DEBUG_MSG(DBG_ALLOC,"got past dyn_table.prev.Alloc()");
  costs.Alloc( aSetNumResidues[0], aSetNumResidues[1] );
  DEBUG_MSG(DBG_ALLOC,"got past dyn_table.costs.Alloc()");

  DEBUG_MSG(DBG_ALLOC,"Allocating Lookup table, Cij.\n");
  Cij = new Real [ aNumResidues[0]*aNumResidues[1] ];
  DEBUG_MSG(DBG_LOOKUP_TABLE_DYN3D,
            "  &(C[0][0]) = " << Cij << "\n"
            "  &(C[" << aNumResidues[0]-1 << "][" << aNumResidues[1]-1 
            << "]) = " << Cij + aNumResidues[0]*aNumResidues[1] -1);

  aOutcomes = new Outcome [ max_num_matches ];

  #ifdef CREATE_OVERALL_PAIR_HISTOGRAM
  // The next array is needed if you wan't to do "pair-histograms"
  // of the final results.
  pair_count = new long* [ aNumResidues[0] ];
  //pair_count = new Real* [ aNumResidues[0] ];
  for (int i=0; i < aNumResidues[0]; ++i)
    pair_count[i] = new long [ aNumResidues[1] ];
    //pair_count[i] = new Real [ aNumResidues[1] ];
  ResetPairHistogram();
  #endif //#ifdef CREATE_OVERALL_PAIR_HISTOGRAM
} // Dyn3d::Alloc()



void
Dyn3d::Alloc(Biopolymer const& m1,
             Biopolymer const& m2,
             int set_n_min,
             int set_n_max)
{
  int aTempNumRes[2];
  aTempNumRes[0] = m1.size();
  aTempNumRes[1] = m2.size();
  Alloc(aTempNumRes,
        set_n_min, set_n_max);
}


void
Dyn3d::Dealloc()
{
  prev.Dealloc();
  costs.Dealloc();
  delete [] Cij;

  #ifdef CREATE_OVERALL_PAIR_HISTOGRAM
  if (pair_count && aNumResidues) {
    for (int i=0; i < aNumResidues[0]; ++i)
      delete [] pair_count[i];
  }
  #endif //#ifdef CREATE_OVERALL_PAIR_HISTOGRAM

  if (aNumResidues)
    delete [] aNumResidues;

  if (aOutcomes)
    delete [] aOutcomes;
} // Dyn3d::Dealloc()




void
Dyn3d::NLayerTable::Alloc(int setNumResA, int setNumResB,
                                                          int set_n_min, int set_n_max)
{
  numResA = setNumResA;
  numResB = setNumResB;
  //numResB = setNumResB;
  max_num_matches = MIN( setNumResA, setNumResB );
  if ((set_n_min < 1) || (set_n_min > max_num_matches) ||
      (set_n_max < 1) || (set_n_max > max_num_matches) ||
      (set_n_max < set_n_min))
    ERR("Minimum and/or maximum number of matches allowed does not "
        "make sense.  The range of allowed matches must be > 0, and < "
        << max_num_matches
        << "\n"
        "Presently min_n = " << set_n_min
        <<", max_n = " << set_n_max);
  n_min = set_n_min;
  n_max = set_n_max;

  if (setNumResA == max_num_matches) {
    assert( setNumResB >= setNumResA );
    windowSizeA = max_num_matches - n_min + 1;
    windowSizeB = windowSizeA + setNumResB-setNumResA;
  }
  else if (setNumResB == max_num_matches) {
    assert( setNumResA >= setNumResB );
    windowSizeB = max_num_matches - n_min + 1;
    windowSizeA = windowSizeB + setNumResA-setNumResB;
  }
  else {
    assert(0);
  }


  int n,i;

  aaaTable = new PrevCellSpecifier** [max_num_matches+1];
  DEBUG_MSG(DBG_ALLOC_DETAIL, "allocating table[][][] = " << aaaTable);
  CHECK_ALLOC(aaaTable);
  for (n=1; n<=n_max; ++n)
  {
    //Since the nth match cannot involve residues from either sequence 
    //that occur before the nth position in that sequence, we will never
    //set i or j to a value less than n.  So why allocate space for cells
    //whose i and j coordinates are less than n?
    //Avoiding this allocation saves a Lot of space (by a factor of 3).
    //
    //   Additionally, if we impose a maximum "window-size", we do not have
    //to allocate space for cells whose i > n+window_sizeA.
    //(The use of a moving window is justified if we are sure that we can
    //demand that at least n_min matches are made (whenever n_min > 1).
    //This should be explained in the Jewett and Huang paper somewhere.)
    long layer_sizeA = MIN(windowSizeA, (setNumResA - n + 1));
    aaaTable[n] = new PrevCellSpecifier* [ layer_sizeA ];
    DEBUG_MSG(DBG_ALLOC_DETAIL, "allocating table[" << n << "][][] = "<<aaaTable[n]);
    CHECK_ALLOC(aaaTable[n]);
    aaaTable[n] -= n;
    for (i=n; i<(n+layer_sizeA); ++i)
    {
      long layer_sizeB = MIN(windowSizeB, (setNumResB - n + 1));
      aaaTable[n][i] = new PrevCellSpecifier [ layer_sizeB ];
      DEBUG_MSG(DBG_ALLOC_DETAIL,
                "allocating table["
                << n << "]["
                << i << "][] = " << aaaTable[n][i]);
      CHECK_ALLOC(aaaTable[n][i]);
      aaaTable[n][i] -= n;
    }
  }
} // NLayerTable::Alloc()



void
Dyn3d::NLayerTable::Dealloc()
{
  //Have to change this part if num_sequences != 2.  I'll worry about
  //this later. (probably using some kind of recursive thing.
  // but still don't know how to get around the type-checking.)
  if (aaaTable)
  {
    for (int n=1; n<=n_max; ++n)
    {
      for (int i=n; i < MIN(numResA+1,(n+windowSizeA)) ; ++i)
      {
        aaaTable[n][i] += n;
        DEBUG_MSG( DBG_ALLOC_DETAIL, 
                   "deleting " << aaaTable[n][i]
                   << " = table["
                   << n << "][" << i << "][]");
        delete [] aaaTable[n][i];
      }
      aaaTable[n] += n;
      DEBUG_MSG( DBG_ALLOC_DETAIL, 
                 "deleting " << aaaTable[n]
                 << " = table["
                 << n << "][][]");
      delete [] aaaTable[n];
    }
    DEBUG_MSG( DBG_ALLOC_DETAIL, 
               "deleting " << aaaTable
               << " = table[][][]");
    delete [] aaaTable;
  }
  aaaTable = NULL;
} // Dyn3d::NLayerTable::Dealloc()



void
Dyn3d::FillLookupTable(Biopolymer const& m1,
                       Biopolymer const& m2)
{
  for(int i = 0;
      i < aNumResidues[0];
      ++i)
  {
    //j varies from [j_lower_bound, j_upper_bound)
    int j, j_lower_bound, j_upper_bound;

    if (i < windowSizeA)
      j_lower_bound = 0;
    else
      j_lower_bound = i - windowSizeA + 1;

    if (i + windowSizeB < aNumResidues[1])
      j_upper_bound = i + windowSizeB;
    else
      j_upper_bound = aNumResidues[1];

    Real *pCost = Cij + i*aNumResidues[1] + j_lower_bound;

    for(j = j_lower_bound;
        j < j_upper_bound;
        ++j)
    {
      *pCost = (*pCostFunc)(m1.begin() + i, m2.begin() + j);

      DEBUG_MSG(DBG_LOOKUP_TABLE_DYN3D,
                "&(Cij) = " << pCost
                << ";   Cij["
                << i + 1 << "]["
                << j + 1 << "] = " << *pCost );
      ++pCost;
    }
  }
} //Dyn3d::FillLookupTable()



#ifdef CREATE_OVERALL_PAIR_HISTOGRAM

void
Dyn3d::ResetPairHistogram()
{
  for (int i=0; i < aNumResidues[0]; ++i)
    for (int j=0; j < aNumResidues[1]; ++j)
      pair_count[i][j] = 0;
}



void
Dyn3d::AccumPairHistogram(int n_start,
                          int n_end,
                          Biopolymer const& m1, //(m1 and m2 are ignored)
                          Biopolymer const& m2)
                          //bool minimize_rmsd_beforehand=false)
{
  PairwiseAlignment nth_alignment;
  for(int n=n_start; n<=n_end; ++n)
  {
    //If there is a valid alilgnment with n matches, include
    //it in the histogram.
    if (aOutcomes[n-n_min].cost < DistanceMetric::UPPER_BOUND)
    {
      ExportSolution(n, nth_alignment);
      nth_alignment.AccumPairHistogram(pair_count);
    }
  }
} // Dyn3d::AccumPairHistogram()




#if 0
// *******   Alternate version of AccumPairHistogram()  ******
//    The following version of AccumPairHistogram() produces an array with
//one entry for every pair of residues (one residue in one chain,
//the other in the opposite chain).
//This array keeps track of the best Gerstein&Levitt score of any alignment
//that matches this pair of residues.
//Residues that belong to alignments with good Gerstein&Levitt scores
//will stand out more noticably in the array.
//So, we are not really taking a "histogram" at all, but this code
// is nearly identical code I wrote with this same function name
// and I decided to keep the name.
void
Dyn3d::AccumPairHistogram(int n_start,
                          int n_end,
                          Biopolymer const& m1,
                          Biopolymer const& m2)
                          // bool minimize_rmsd_beforehand)
{
  PairwiseAlignment nth_alignment;
  if (n_start < 30)
        n_start = 30;

  //loop over the list of alignments (each with different numbers of matches)
  //For each pair of matched residues in the alignment, update
  //that pairs current score to reflect if the new score is lower.
  //(Recall lower means the structures correspond better.)
  for(int n=n_start; n<=n_end; ++n)
  {
    ExportSolution(n, nth_alignment);
    
    //Calculate the LevittGerstein structural probability score
    Real Sstr, mu_str, delta_str, Z; //some dummy parameters to pacify syntax
    Real Pstr = LevittGerstein98::Sstr_Probability(nth_alignment,
                                                   m1,
                                                   m2,
                                                   &Sstr,
                                                   &mu_str,
                                                   &delta_str,
                                                   &Z);
    for(int k=0; k < n; ++k)
    {
      assert(nth_alignment.NumMatches() == n);
      int i = nth_alignment[k][0];
      int j = nth_alignment[k][1];

      #if DBG_PAIR_HISTOGRAM
      if ((Pstr < pair_count[i][j])
          && (Pstr < 0.001))
      {
        if ((i == 45) && (j == 15))
        {
          char comment_str[512];
          sprintf(comment_str,
                  "This alignment matched i=45 j=15 with Pstr=%e",
                  Pstr);
          nth_alignment.ExportMSF("align_weird.msf",
                                  m1,
                                  m2,
                                  "Mol1",
                                  "Mol2",
                                  comment_str);

          sprintf(comment_str,
                  "i=45 j=15 Pstr=%e",
                  Pstr);

          nth_alignment.ExportGFX(m1,
                                  m2,
                                  "align_weird.gfx",
                                  comment_str,
                                  true,
                                  true,
                                  0.25,
                                  NULL,
                                  1,
                                  45,
                                  15);

          //cerr <<
          //  "New align_weird.msf/gfx file generated.\n"
          //  "  Press <Enter> to continue..." << flush;
          //char s[32];
          //cin >> s;
          //cerr << "continuing..."<<endl;

        }
      } //if ((Pstr < pair_count[i][j]) && (Pstr < 0.001))
      #endif // #if DBG_PAIR_HISTOGRAM
      pair_count[i][j] = MIN(Pstr, pair_count[i][j]);

      ++pair_count[i][j];

    } //for(int k=0; k < n; ++k)
  } //for(int n=n_start; n<=n_end; ++n)
} // Dyn3d::AccumPairHistogram()
#endif //#if 0



void
Dyn3d::PrintPairHistogram(char *filename) const
{
  FILE *histogram_file;
  //  fstream histogram_file("pair_hist_overall", ios::out);
  histogram_file = fopen("pair_hist_overall.vrml", "w");

  if (! histogram_file)
    ERR("Can't open file \"pair_hist_overall.vrml\" for writing.");

  fprintf(histogram_file, "#VRML V2.0 utf8\n\n");
  fprintf(histogram_file, "ElevationGrid { \n");
  fprintf(histogram_file, "  xDimension %d\n", aNumResidues[0]);
  fprintf(histogram_file, "  zDimension %d\n", aNumResidues[1]);
  fprintf(histogram_file, "  xSpacing 1.0\n");
  fprintf(histogram_file, "  zSpacing 1.0\n");
  fprintf(histogram_file, "  height [");

  //  histogram_file << setw(4);
  for (int i=0; i < aNumResidues[0]; ++i)
  {
    for (int j=0; j < aNumResidues[1]; ++j)
    {
      //histogram_file << pair_count[j][i];
      //fprintf(histogram_file, "%ld", pair_count[i][j]);
      Real max_val = 15.0;
      Real min_val = 0;
      Real x = (Real)pair_count[i][j];
      Real y = ((x == 0.0) ? max_val : -log10(x));
      y = MIN(y, max_val);
      y = MAX(y, min_val);
      fprintf(histogram_file, "%f", y);
      //histogram_file << "\n";
      if (((j+1) != aNumResidues[1])
          || ((i+1) != aNumResidues[0]))
        fprintf(histogram_file, ",");
          fprintf(histogram_file, "\n");
    }
    //histogram_file << "\n";
    //        fprintf(histogram_file, "#row %d\n", i+1);
    fprintf(histogram_file, "\n");
  }

  fprintf(histogram_file, "  ]\n}\n");
  fclose(histogram_file);



  //fstream histogram_file("pair_hist_overall", ios::out);
  histogram_file = fopen("pair_hist_overall.ppm", "w");

  if (! histogram_file)
    ERR("Can't open file \"pair_hist_overall.ppm\" for writing.");

  int const NUM_COLOR_LEVELS = 256;

  fprintf(histogram_file, "P3\n");
  //  fprintf(histogram_file, "# pair_hist_overall.ppm\n");
  fprintf(histogram_file, "%d %d\n", aNumResidues[0], aNumResidues[1]);
  fprintf(histogram_file, "%d\n", NUM_COLOR_LEVELS-1);

  //  histogram_file << setw(4);
  for (int j=0; j < aNumResidues[1]; ++j)
  {
    for (int i=0; i < aNumResidues[0]; ++i)
    {
      //histogram_file << pair_count[j][i];
      Real max_val = 15.0;
      Real min_val = 0;
      Real x = (Real)pair_count[i][aNumResidues[1]-1-j];
      Real y = ((x == 0.0) ? max_val : -log10(x));
      y = MIN(y, max_val);
      y = MAX(y, min_val);
      int value = floor((y-min_val)*(NUM_COLOR_LEVELS-1) / (max_val-min_val));
      //value = (value * 2)/3;
      //  if (value)
      //  value += (NUM_COLOR_LEVELS-1)/3;

      //if (value) {
      if (value != 1.0) {
        fprintf(histogram_file, "%d %d %d",
                value,
                value*2/3,
                (NUM_COLOR_LEVELS-1-value)*1 / 2
                );
      }
      else 
        fprintf(histogram_file, "0 0 0");
          
      //histogram_file << "\n";
      //          if (((j+1) != aNumResidues[1])
      //                  || ((i+1) != aNumResidues[0]))
      //                fprintf(histogram_file, ",");
      fprintf(histogram_file, "\n");
    }
    //histogram_file << "\n";
    //        fprintf(histogram_file, "#row %d\n", j+1);
    //        fprintf(histogram_file, "\n", j+1);
  }
  fclose(histogram_file);

  //  fstream histogram_file("pair_hist_overall", ios::out);
  histogram_file = fopen(filename, "w");
  if (! histogram_file)
    ERR("Can't open file \"" << filename << "\" for writing.");
  
  fprintf(histogram_file,
          "num rows equals length of sequence A:%d\n",
          aNumResidues[0]);
  fprintf(histogram_file,
          "num columns equals length of sequence B:%d\n",
          aNumResidues[1]);

  //  histogram_file << setw(4);
  for (int i=0; i < aNumResidues[0]; ++i)
  {
    for (int j=0; j < aNumResidues[1]; ++j)
    {
      //histogram_file << pair_count[j][i];
      fprintf(histogram_file, "%ld", pair_count[i][j]);
      //fprintf(histogram_file, "%e", pair_count[i][j]);
      //histogram_file << "\n";
      if ((j+1) != aNumResidues[1])
        fprintf(histogram_file, " ");
    }
    //histogram_file << "\n";
    //fprintf(histogram_file, "#row %d\n", i+1);
    fprintf(histogram_file, "\n");
  }

  fclose(histogram_file);
} // Dyn3d::PrintPairHistogram()

#endif //#ifdef CREATE_OVERALL_PAIR_HISTOGRAM

