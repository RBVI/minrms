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

#include <simple_numeric_utils.h>
#include <apply_trans.h>
#include <distance_metric.h>
#include "orientation_filter.h"

using namespace minrms;


OrientationFilter::
OrientationFilter(Real set_max_phys_distance,
                  int  set_min_num_matches,
                  int  numRes1,
                  int  numRes2):
  dyn_nw(numRes1, numRes2, set_min_num_matches),
  solution(MIN(numRes1, numRes2)),
  max_phys_distance(set_max_phys_distance),
  min_num_matches(set_min_num_matches)

{
  assert(min_num_matches <= MIN(numRes1, numRes2));

  cout <<
      "Needleman & Wunsch orientation filter selected.  Parameters:\n"
      "   A minimum of " << set_min_num_matches <<
      " matches required with\n"
      "   no match exceeding a distance of: "
         << set_max_phys_distance
         << " angstroms" << endl;


  //It's too hard to explain in detail why I use the
  //following gap penalty.  I chose it to be a
  //value I can be sure is larger than the
  //difference in sum-squared-distance between
  //two successive minimal RMSD alignments
  //with N, and N+1 matched pairs of residues.
  //This means, alignments generated later on
  //will never have gaps in them if that's possible.
  //We try to use as small a value as possible
  //because we don't want any number of gaps to
  //equal the penalty incured by violating the distance
  //constraint (which is DistanceMetric::max_allowed_dist_metric,
  //see "distance_metric.h")
  gap_penalty = 0.5 * SQR(set_max_phys_distance)
                    * MIN(numRes1, numRes2)
                    * Biopolymer::Residue::NumBackboneAtoms();
  assert(gap_penalty >= 0.0f);

  dyn_nw.InitTable(gap_penalty);

}//OrientationFilter::OrientationFilter()
  

  
bool OrientationFilter::
GoodOrientation(Biopolymer const& m1_rotated,
                Biopolymer const& m2_rotated)
{
  //We need to set our environment up so that orientation filtering
  //will work correctly.
  //So we temporarily set the "maximum distance cutoff" to the distance
  //the user specified for use with orientation filtering.
  Real old_max_dist = DistanceMetric::GetMaxPhysicalDistance();
  DistanceMetric::SetMaxPhysicalDistance(max_phys_distance);




  #if (! defined NW_PRECOMPUTE_CIJ)
  // Find the alignment with the maximum number of matched pairs possible,
  // given the maximum-allowed-distance-constraint.
  dyn_nw.Solve(m1_rotated, m2_rotated, solution);

  #else

  // First, figure out distances between all pairs of residues in advance...
  dyn_nw.ComputePairwiseMatchCosts(m1_rotated, m2_rotated);
  // Once the distances have been computed, we can invoke dyn_nw.Solve()
  // to find the alignment with the maximum number of matched pairs possible,
  // given the maximum-allowed-distance-constraint.
  dyn_nw.Solve(m1_rotated, m2_rotated, solution, true);

  #endif // #if (! defined NW_PRECOMPUTE_CIJ)





  DEBUG_MSG(DBG_NW_ORIENTATION_FILTER,
            "N.W. solution contains " << solution.NumMatches()
            << " matches.");

  //and reset the original distance limits.
  DistanceMetric::SetMaxPhysicalDistance(old_max_dist);

  return (solution.NumMatches() >= min_num_matches);
}//OrientationFilter::GoodOrientation()

