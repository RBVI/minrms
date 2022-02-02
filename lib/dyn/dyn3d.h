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

#ifndef _DYN3D_H
#define _DYN3D_H

#include <iostream>
#include <cassert>
using namespace std;

#include <biopolymer.h>
#include <pair_alignment.h>
#include "distance_metric.h"


//rcs_id: "@(#)$Header: /usr/local/src/minrms/lib/dyn/RCS/dyn3d.h,v 3.7 2005/02/24 00:48:21 conrad Exp $"


namespace minrms {


#ifdef COST_FUNC_DEPENDS_ON_ALIGNMENT
#undef COST_FUNC_DEPENDS_ON_ALIGNMENT
//    You may try to use this method to minimize things like
//intramolecular similarity RMS (comparisons of difference in distance
//between analagous residues _within_ each molecule),
//instead of intermolecular RMSD  (difference in position
//_between_ the two molecules)
//    It seems deceptively reasonable to allow the marginal cost of
//matching an additional pair of residues to an alignment (Cij),
//to depend on the alignment made so far.  This destroys the
//assumption necessary to justify dynamic programming in the first place:
//that "Optimal solutions must have optimal sub-solutions".
//You will get really weird answers if you do this.
//This is a bad idea.
//Leaving this comment in here as a cautionary tale. -AJ03.31.98
#endif //#ifdef COST_FUNC_DEPENDS_ON_ALIGNMENT




class Dyn3d
{


  class TwoLayerTable
  {

    TwoLayerTable(const TwoLayerTable &); //dissable, no implicit copying allowed
    TwoLayerTable &operator =(const TwoLayerTable &); //no implicit copying

    long size_of_n_layer;
    long size_of_ni_row;
    //T *aT;
    Real *aT;

  public:

    inline TwoLayerTable() {
      size_of_n_layer = size_of_ni_row = 0;
    }
    inline TwoLayerTable(int numResA, int numResB) {
      Alloc(numResA, numResB);
    }
    inline ~TwoLayerTable() { Dealloc(); }

    //memory managment:
    inline void Alloc(int numResA, int numResB)
    {
      size_of_ni_row  =  numResB;
      size_of_n_layer = (numResA) * size_of_ni_row;
      //aT = new T [ size_of_n_layer * 2 ];
      aT = new Real [ size_of_n_layer * 2 ];
      DEBUG_MSG( DBG_ALLOC, 
                 "alocating " << aT
                 << " = 2LayerTable[]");
      CHECK_ALLOC(aT);
    }
        
    inline void Dealloc()
    {
      size_of_n_layer = size_of_ni_row = 0;
      if (aT!=NULL)
      {
        DEBUG_MSG( DBG_ALLOC, 
                   "deleting " << aT
                   << " = 2LayerTable[]");
        delete [] aT;
        aT = NULL;
      }
    }

    //        inline T Get(int n, int i, int j) const
    inline Real Get(int n, int i, int j) const
    {
      assert((1<=n) && (1<=i) && (1<=j) && (j <= size_of_ni_row));
      return aT[
                (n&1) * size_of_n_layer +
                (i-1) * size_of_ni_row  +
                (j-1)
                ];
    }
  
    //        inline void Set(int n, int i, int j, T set_value)
    inline void Set(int n, int i, int j, Real set_value)
    {
      assert((1<=n) && (1<=i) && (1<=j) && (j <= size_of_ni_row));
      aT[
         (n&1) * size_of_n_layer +
         (i-1) * size_of_ni_row  +
         (j-1)
         ] = set_value;
    }

    //As an alternative, for faster access, you can find the address of the
    //first element of the row by calling the function below, and in the
    //inner-loop do strictly additive pointer arithmetic.
    //Here, you can read or modify the (private) af data, so be careful.
    //    GetAddressOfNIthRow()
    //returns the address of the first element of the ith row of the
    //(n&1 even/odd)th layer of the cost-table.
    //(Note1:The indexing for n and i begins at 1, not 0, although space
    //       for the zeroith element is allocated...
    // Note2:This address is located at an offset of 1 from the
    //       actual beginning of space allocated for that row,
    //       no thanks to the confusing indexing method I'm using.)
    //        inline T *GetAddressOfNIJ(int n, int i, int j) const
    inline Real *GetAddressOfNIJ(int n, int i, int j) const
    {
      return &( aT[
                   (n&1) * size_of_n_layer + 
                   (i-1) * size_of_ni_row +
                   (j-1)
                   ]
                );
    }
  }; //class TwoLayerTable



// (Note: DEC compiler doesn't like private nested classes)
public:

  class NLayerTable
  {
    NLayerTable(const NLayerTable &); // dissable, no implicit copying allowed
    NLayerTable &operator =(const NLayerTable &); //no implicit copying
    int numResA, numResB, max_num_matches;
    //window-specific parameters:
    int n_min, n_max;
    int windowSizeA, windowSizeB;

  public:

    typedef unsigned char PrevCellSpecifier;//Pointer to the previous recursion
    static const PrevCellSpecifier CASE_N1_I1_J1;//previous cell is n-1,i-1,j-1
    static const PrevCellSpecifier CASE_N_I1_J;  //previous cell is   n,i-1,j
    static const PrevCellSpecifier CASE_N_I_J1;  //previous cell is   n, i, j-1

  private:
    PrevCellSpecifier ***aaaTable; //Backtrace through table 2find the alignment

  public:

    inline NLayerTable() { aaaTable = NULL; }
    inline NLayerTable(int numResA, int numResB,
                                           int n_min, int n_max)
    { Alloc(numResA, numResB, n_min, n_max); }

    inline ~NLayerTable() { Dealloc(); }

    //memory managment:
    void Alloc(int numResA, int numResB,
                           int n_min, int n_max);

    void Dealloc();
        
    //data access/retrieval:
    int GetPrevCell(int n, int i, int j, int which_protein) const;

    inline bool IsMatch(int n, int i, int j) const
    {
      return (aaaTable[n][i][j] == CASE_N1_I1_J1);
    }

    inline void SetPrevCell(int n, int i, int j, PrevCellSpecifier which_case)
    { aaaTable[n][i][j] = which_case; }

    //As an alternative, for faster access, you can find the address of the
    //first element of the row by calling the function below, and in the
    //inner-loop do strictly additive pointer arithmetic.
    //Here, you can read or modify the (private) af data, so be careful.
    //    GetAddressOfNIthRow()
    //returns the address of the first element of the ith row of the
    //(n&1 even/odd)th layer of the cost-table.
    //(Note1:The indexing for n and i begins at 1, not 0, although space
    //       for the zeroith element is allocated...
    // Note2:This address is located at an offset of 1 from the
    //       actual beginning of space allocated for that row,
    //       no thanks to the joyous indexing method I'm using.)
    inline PrevCellSpecifier *GetAddressOfNIJ(int n, int i, int j) const
    { return &(aaaTable[n][i][j]); }
  }; //class NLayerTable
private:




  Dyn3d(const Dyn3d &); // dissable, no implicit copying allowed
  Dyn3d &operator =(const Dyn3d &); //no implicit copying




  // Alloc() allocates storage space for the the table, 
  // and performs some initialization.
  void Alloc(int const *aSetNumResidues,
             int set_n_min,
             int set_n_max);

  // Or, if you prefer, you can call Alloc() by specifying the Loaded
  // sequence data. Whichever you prefer.
  void Alloc(Biopolymer const& m1,
             Biopolymer const& m2,
             int set_n_min,
             int set_n_max);

  // The next function deallocates the memory allocated by Alloc().
  void Dealloc();




public:

  //void SetMinN(int set_n_min);
  //inline void SetMaxN(int set_n_max)  {
  //if ((set_n_max < 1) || (set_n_max > max_num_matches))
  //ERR("You cannot set the minimum number of matches\n"
  //"made to less than 1, or greater than the length\n"
  //"of the shorter sequence.");
  //n_max = set_n_max;
  //}
  //inline int GetMinN() const { return n_min; }
  //inline int GetMaxN() const { return n_max; }

  inline static void SetCutoffWorstPair(Real set_cutoff_worst_pair)
  { cutoff_worst_pair = set_cutoff_worst_pair; }
  inline static Real GetCutoffWorstPair()
  { return cutoff_worst_pair; }
  static const Real NO_CUTOFF_WORST_PAIR;

  //No default constructor
  //Dyn3d(); 

  //The next constructor allocates memory for comparing two sequences,
  //one, with num_res_A residues, and another with num_res_B residues.
  inline Dyn3d(int num_res_A, int num_res_B,
               int n_min, int n_max)
  {
    int num_residues[2];
    num_residues[0] = num_res_A;
    num_residues[1] = num_res_B;
    Alloc(num_residues, n_min, n_max);
  }

  ~Dyn3d();

  //  inline PairwiseMatch
  //  *GetPtr2Entry(int n, int i, int j) const {
  //        assert(n>=0); assert(n <= max_num_matches);
  //        assert(i>=n); assert(i <= aNumResidues[0]);
  //        assert(j>=n); assert(j <= aNumResidues[1]);
  //        return prev.GetAddressOfNIJ(n,i,j);
  //  }

  inline
  int GetPrevCell(int n, int i, int j, int which_sequence) const {
    return prev.GetPrevCell(n,i,j, which_sequence);
  }

  inline
  bool IsMatch(int n, int i, int j) const {
    return prev.IsMatch(n,i,j);
  }

  // The following function specifies the individual contribution
  // to the cost from matching residue i, from sequence A, with
  // residue j from sequence B.  It should be supplied by the user.
  //   An additional third argument is passed to to this function.
  // You can pass information about the allignment-so-far with this argument.
  // It is only necessary if the incremental cost
  // for the match between residues i and j depends on present state
  // of the alignment before this match is made.  (Such would be the case,
  // if you were using intra-distance matrices as described by Holm&Sanders,
  // because the nth additional match, causes n additional terms to be added
  // to the overall summation.)
  //   If you don't need this (as would be if we were just using RMS or
  // average-distance between sequences whose rotation is fixed), then
  // you can just pass NULL.
  Real (*pCostFunc)(Biopolymer::const_iterator i,
                    Biopolymer::const_iterator j);
  void FillLookupTable(Biopolymer const& m1,
                       Biopolymer const& m2);

  //      You can do final manipulations on the total-sum-of-the-costs,
  // you get after alignment by setting the pApply2Results() function.
  //      I wanted to add this because, if you want to find the root-mean-
  // squared-distance between matching residues in an allignment, you need
  // something like this to compute the average, and then the square-root.
  //      In this example, you'd choose the cost-function to be:
  //       "pCostFunc = &SquareDistanceBetween(i,j)",
  // and the pApply2Results function to be:
  //       "pApply2Results = &DivideByNthenSqrt(tot_cost, n)"
  // The composition of these two functions returns the Root-Mean-Squared.
  //      Lastly, if you do not set this function, the default setting
  // The default setting (NULL) just returns the identity (raw cost).

  Real (*pApply2Results)(Real total_cost, int n);

  // Solve() finds the (set of) alignment(s) that minimize the cost function
  // defined above.
  void Solve();

  // ..If you prefer to, you can call SolveIncrN(n) over all n
  // from 1 to max_num_matches.
  // This way you can make function calls in between each step.
  // is done between each step.
  void Solve_IncrN(int n);
  // However, you must first call InitTable_IncrN()
  // before you call Solve_IncrN().  (Not required when calling Solve())
  void InitTable_IncrN();

  // Output:
  //  The next function returns one of the alignments calculated
  // by Solve() or Solve_IncrN() back to the caller, by saving it
  // in the form of a "PairwiseAlignment" data structure, "dest".
  // The integer "n" indicates the number of equivalenced-residues
  // in the solution being sought (of course indexing starts at 1, not 0).
  // (Note: The "dest" argument should be pre-allocated so that
  //        it is large enough to store "n" matches before calling
  //        this function.)
  void ExportSolution(int n,
                      PairwiseAlignment& dest);

  // GetMaxNumMatches returns the maximum possible number pairs that can
  // be matched, based on the number of residues in each sequence, and
  // the number of sequences.  (I will have to generalize what I mean by this
  // if NUM_SEQUENCES != 2.  Worry about that later.)
  inline int
  GetMaxNumMatches() const { return max_num_matches; }

  inline int GetNumRes(int seq) const {
    assert(seq < 2); 
    return aNumResidues[seq];
  }

  //GetCost() is called after Solve() or Solve_IncrN() is called.
  //GetCost() retrieves the "cost" (or "score", whichever term you prefer)
  //of the alignment containing n matches, after it has been calculated.
  //This value happens to be the minimal sum-squared-
  //distance of aligning n residues from either protein, at this orientation.
  //This function is undefined unless "n" lies in the range [n_min,n_max].
  inline Real
  GetCost(int n) const {
    assert((n_min <= n) && (n <= n_max));
    return aOutcomes[n-1].cost;
  }

  //MaxDistanceExceeded(), checks to see if an alignment generated
  //by Solve_IncrN() violated the maximum distance constraint
  //set previously (using DistanceMetric::SetMaxDistance()).
  //Note:
  //It is valid for, n = 1 ... n_max, not just n= n_min...n_max.
  //If (n < n_min), then it predicts whether or not any alignment with
  //at least n_min matches will violate the distance constraint,
  //even before those alignments are calculated.
  //   It is useful for terminating a calculation really early if the
  //partial alignments generated so far are already looking really poor.
  bool MaxDistanceExceeded(int n)
  { return (aOutcomes[n-1].cost >= DistanceMetric::UPPER_BOUND); }

  #ifdef CREATE_OVERALL_PAIR_HISTOGRAM

  //The histogram-table (pair_count) is initialized at creation to zero.
  //Unless reset by manually calling ClearPairHistogram(), successive calls
  //to this function will continue to accumulate onto the previous counts.
  void AccumPairHistogram(int n_start,
                          int n_end,
                          Biopolymer const& m1,
                          Biopolymer const& m2);
                          //bool minimize_rmsd_beforehand=false);
  inline void AccumPairHistogram(Biopolymer const& m1,
                                 Biopolymer const& m2)
                               //bool minimize_rmsd_beforehand=false);
  { AccumPairHistogram(1, max_num_matches, m1, m2); }

  void ResetPairHistogram();


  inline long * const *GetPairHistogram() const { return pair_count; }
  //inline Real * const *GetPairHistogram() const { return pair_count; }

  //This next function assumes you have previously called AccumPairHistogram().
  //The next function prints out a big two-dimensional array of
  //whitespace-delimited (long) integers.  Each one being the number of times
  //in the set of alignments (that were accepted in the final results)
  //that a particular residue (i) from the first sequence (A),
  //was matched with a particular residue (j) from the second sequence (B)
  //  The first two lines of this file are:
  //"length_sequence_A:"(which is the number of rows in the array)
  //"length_sequence_B:" (which is the number of collumns).
  //The name of the file is specified by the "filename" argument.
  void PrintPairHistogram(char *filename) const;
  #endif //#ifdef CREATE_OVERALL_PAIR_HISTOGRAM

private:

  //The next member stores a history of the decisions made at every
  //step of the algorithm. (whether a match was made, if so, who gets
  //matched to who)  This array is used afterwards, to figure out 
  //which residues get matched to who in the 'optimal' alignment
  //it calculates.
  //  (Incidentally, the costs associated with the alignment up to this
  //point are stored seperately in a (much smaller) array. see below.)
  NLayerTable        prev;
  
  Real *Cij;

  //     The next member stores the minimum total-cost-so-far
  // (The sum of the Cij's) required to match n of the first
  // i residues of sequence A, with n of the first
  // j residues f sequence B.
  //    Note: Unlike "prev",  "most_recent" and "costs" only store the last
  // two layers of this (the last two layers being the data relevant to
  // the most recently calculated value of n, and n-1).  This is because 
  // this data is never referenced later than two iterations of 
  // the outer loop (over which n varies) from the time it is originally
  // set.  So there's no need to keep track of it for longer than that.
  // This saves space by a factor of 2 or 3.
  //  TwoLayerTable<Real> costs;
  TwoLayerTable costs;


  //    Once the final results are calculated, they are stored in
  // a single-dimensional "aOutcomes[]" array (indexed by "n", the
  // number-of-matches-made) which stores the best cost for matching
  // n residues from sequence A with n residues from sequence B.
  struct Outcome
  {
    Real    cost;
    //(I had more data members here originally, but I removed them
    // to get rid of "#ifdef COST_FUNC_DEPENDS_ON_ALIGNMENT" AJ08.11.98)
  };
  Outcome *aOutcomes;

  //Stores the number of residues in each sequence
  int *aNumResidues;

  //Stores the maximum number of matches possible, given the number of
  //residues in each sequence, and the number of sequences.
  int max_num_matches;

  #ifdef CREATE_OVERALL_PAIR_HISTOGRAM
  //The next array is only used to store "pair-histogram" data.
  //(Stores a count of how many times residue i from seqA was matched
  // with residue j from seqB over all the 'best' alignments with 
  // various "n"s (numbers of matches made) for this run, but with
  // all other factors (like physical orientation) held constant.)
  long **pair_count;
  //Real **pair_count;
  #endif //#ifdef CREATE_OVERALL_PAIR_HISTOGRAM

  //specifies whether or not to apply the additional constraint of a 
  //worst-pair cutoff.
  static Real cutoff_worst_pair;

  //window-specific parameters:
  int n_min, n_max;
  int windowSizeA, windowSizeB;

}; //class Dyn3d

} //namespace minrms

#endif // #ifndef _DYN3D_H
