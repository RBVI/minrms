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

#include <map> //not <multimap>, at least on the dec-alpha's
#include <set>
#include <sstream>
#include <vector>
#include <cstdio>  //needed for "sprintf()"
#include <cassert>
using namespace std;

#include <unistd.h>       //needed for "unlink()"
#include <sys/time.h>     //required for "struct timezone" and "gettimeofday()"
#include <sys/resource.h> //required for "struct rlimit" and "getrlimit()"
#include <global_utils.h>
#include <short_string.h>
#include <simple_numeric_utils.h> //required for reference to FindMedian()
#include <vect_3d.h>      //required for referrence to Vect3 and Matrix3x4
#include <fast_rot_metric.h>
#include <binary_file.h>
#include <biopolymer.h>
#include <mol2sequence.h>
#include <apply_trans.h>
#include <superimpose.h>
#include <or_gen.h>
#include <orientation_filter.h>
#include <dyn.h>
#include "search_or.h"
#include "search_or_dyn3d.h"

using namespace minrms;
//RCS-version-string:
//static char rcs_id[] = "@(#)$Header: /usr/local/src/minrms/bin/minrms/RCS/search_or_dyn3d.cc,v 3.9 2005/02/24 00:49:03 conrad Exp $";



//some constants
const char *SearchOrDyn3d::BASE_FILE_NAME = "align";
const Real SearchOrDyn3d::Settings::dyn3d_NO_MAX_RMSD = -1.0;


#ifdef COMPILER_BUG2

SearchOrDyn3d::Settings::Settings()
 :SearchOr::Settings(),
  n_max(NO_MAX_NUM_MATCHES),
  dyn3d_max_rmsd(dyn3d_NO_MAX_RMSD),
  dyn3d_max_dist(DistanceMetric::NO_MAX_PHYS_DISTANCE),
  orientation_filter_nw(false),
  orientation_filter_nw_max_dist(8.0),
  orientation_filter_nw_min_num_matches(0.33),
  display_stats_interval(DISSABLE_DISPLAY_STATS_INTERVAL)
{
  //cerr << "SearchOrDyn3d::Settings, this = " << this << endl;
  //cerr << "default constructor called." << endl;
}


  //copy constructor:
  //(The alpha compiler produces buggy
  // binaries unless I define one explicitely)

SearchOrDyn3d::Settings::Settings(Settings const& s)
 :SearchOr::Settings(s),
  n_max(s.n_max),
  dyn3d_max_rmsd(s.dyn3d_max_rmsd),
  dyn3d_max_dist(s.dyn3d_max_dist),
  orientation_filter_nw(s.orientation_filter_nw),
  orientation_filter_nw_max_dist(s.orientation_filter_nw_max_dist),
  orientation_filter_nw_min_num_matches(s.orientation_filter_nw_min_num_matches),
  display_stats_interval(s.display_stats_interval)
{
  //cerr << "SearchOrDyn3d::Settings, this = " << this << endl;
  //cerr << "s.aMol = " << s.aMol << ";  aMol = " << aMol << endl;
  //cerr << "s.aMol_f = " << s.aMol_f << ";  aMol_f = " << aMol_f << endl;
}

#endif //#ifdef COMPILER_BUG2



SearchOrDyn3d::SearchOrDyn3d(SearchOrDyn3d::Settings const& s):
  settings(s),
  dyn3d(s.aMol_f[0].size(), s.aMol_f[1].size(), s.n_min, s.n_max)
{
  assert(settings.aMol);
  assert(settings.aMol_f);

  //Assign the cost-function criteria used by the dynamic programming algorithm
  //(The following line sets RMSD-between-corresponding-atoms,
  // as the cost-function used to select the alignments.)
  dyn3d.pCostFunc = &DistanceMetric::SquaredDistanceBetweenIandJ;

  //Now, set the maximum distance between residues that are allowed
  //to be matched.  Pairs of residues further than this distance
  //from eachother will not be matched.
  DistanceMetric::SetMaxPhysicalDistance(settings.dyn3d_max_dist);

  //The next line is not really necessary:
  //dyn3d.pApply2Results = &DistanceMetric::DivideByNthenSqrt; commenting out

  aMol_c = new Biopolymer[2];
  for (int m = 0; m < 2; ++m)
  {
    aMol_c[m] = settings.aMol_f[m];
  }
  numRes1 = settings.aMol_f[0].size();
  numRes2 = settings.aMol_f[1].size();
  max_num_matches = MIN(numRes1, numRes2);
  n_min = settings.n_min; //shorthand
  n_max = settings.n_max; //shorthand
  range_of_n = n_max - n_min + 1;
  alt_method_requires_alt_solution_table =
    ((settings.alt_method == Settings::ALT_METHOD_3D)
     ||
     (settings.alt_method == Settings::ALT_METHOD_3D_ORIG)
     ||
     ((settings.alt_method == Settings::ALT_METHOD_PAIRS)
      &&
      (settings.alt_rmsd_tolerance
       !=
       Settings::ALTERNATES_NOT_LIMITED_BY_RMSD)
      )
     );

  #ifdef DEBUG
  int n1_orig = settings.aMol[0].size();
  int n2_orig = settings.aMol[1].size();
  #endif // #ifdef DEBUG
  int max_num_matches_unfiltered = MIN(settings.aMol[0].size(),
                                       settings.aMol[1].size());
  assert(max_num_matches_unfiltered >= max_num_matches);
  alignment_temp.resize(max_num_matches_unfiltered);

  assert(settings.num_solutions >= 1);

  apAltOrientationFiles = NULL;
  aaSumSqdDists = new Real* [settings.num_solutions];
  aaTransforms  = new Matrix3x4* [settings.num_solutions];
  aaAlignments  = new PairwiseAlignment* [settings.num_solutions];
  CHECK_ALLOC(aaSumSqdDists && aaTransforms && aaAlignments);

  for (int a = 0; a < settings.num_solutions; ++a)
  {
    aaSumSqdDists[a] = new Real [range_of_n];
    aaTransforms[a] = new Matrix3x4 [range_of_n];
    aaAlignments[a] = new PairwiseAlignment [range_of_n];
    CHECK_ALLOC(aaSumSqdDists[a] && aaTransforms[a] && aaAlignments[a]);

    for (int n=n_min; n<=n_max; ++n)
    {
      aaAlignments[a][n-n_min].resize(n);//The nth alignment has (n+1) matches in it.
      aaSumSqdDists[a][n-n_min] = DistanceMetric::UPPER_BOUND; //Flag: Impossibly high value.
      CopyMatrix3x4(g_IDENTITY3X4,
                    aaTransforms[a][n-n_min]);
    }
  }

  aRotateUsingTheseCoords1 =new Vect3[ numRes1 *
                                     Biopolymer::Residue::NumBackboneAtoms()];
  aRotateUsingTheseCoords2 =new Vect3[ numRes2 *
                                     Biopolymer::Residue::NumBackboneAtoms()];
  CHECK_ALLOC(aRotateUsingTheseCoords1 && aRotateUsingTheseCoords2);

  int num_points = 0;
  for(Biopolymer::const_iterator ps = settings.aMol_f[0].begin();
      ps < settings.aMol_f[0].end();
      ++ps)
  {
    for(int q=0; q < Biopolymer::Residue::NumBackboneAtoms(); ++q)
    {
      Biopolymer::Residue::const_iterator pa;
      pa = ps->GetBackboneAtom(q);
      assert(pa != ps->end());
      aRotateUsingTheseCoords1[num_points][0] = (*pa).second.xyz[0];
      aRotateUsingTheseCoords1[num_points][1] = (*pa).second.xyz[1];
      aRotateUsingTheseCoords1[num_points][2] = (*pa).second.xyz[2];
      ++num_points;
    }
  }

  num_points = 0;
  for(Biopolymer::const_iterator ps = settings.aMol_f[1].begin();
      ps < settings.aMol_f[1].end();
      ++ps)
  {
    for(int q=0; q < Biopolymer::Residue::NumBackboneAtoms(); ++q)
    {
      Biopolymer::Residue::const_iterator pa;
      pa = ps->GetBackboneAtom(q);
      assert(pa != ps->end());
      aRotateUsingTheseCoords2[num_points][0] = (*pa).second.xyz[0];
      aRotateUsingTheseCoords2[num_points][1] = (*pa).second.xyz[1];
      aRotateUsingTheseCoords2[num_points][2] = (*pa).second.xyz[2];
      ++num_points;
    }
  }

#ifdef CREATE_PROFILE_HISTOGRAM
  //This next array was only created for the purpose of diagnosing
  //the effectiveness of an optimization we wrote.
  //It is not something the user would generally be interested in.
  aNumLayersHistogram = new long[max_num_matches];
  for(int n=0; n < max_num_matches; ++n)
    aNumLayersHistogram[n] = 0;
#endif //#ifdef CREATE_PROFILE_HISTOGRAM
} //SearchOrDyn3d::SearchOrDyn3d()



SearchOrDyn3d::~SearchOrDyn3d()
{
  CleanUpAltSolutionTable();

  if (aMol_c) delete [] aMol_c;

  if (aaSumSqdDists)
  {
    for (int a=0; a < settings.num_solutions; ++a)
    {
      assert(aaSumSqdDists[a]);
      delete [] aaSumSqdDists[a];
    }
    delete [] aaSumSqdDists;
  }

  if (aaAlignments)
  {
    for (int a=0; a < settings.num_solutions; ++a)
    {
      assert(aaAlignments[a]);
      delete [] aaAlignments[a];
    }
    delete [] aaAlignments;
  }

  if (aaTransforms)
  {
    for (int a=0; a < settings.num_solutions; ++a)
    {
      assert(aaTransforms[a]);
      delete [] aaTransforms[a];
    }
    delete [] aaTransforms;
  }

  if (aRotateUsingTheseCoords1) delete [] aRotateUsingTheseCoords1;
  if (aRotateUsingTheseCoords2) delete [] aRotateUsingTheseCoords2;
} //SearchOrDyn3d::~SearchOrDyn3d()






void SearchOrDyn3d::Search(OrientationGenerator const& og,
                           long start,
                           long stop)
{
  // (Necessary for finding alternate solutions)
  if ((settings.num_solutions > 1) && (alt_method_requires_alt_solution_table))
    InitAltSolutionTable();


  //***    This next function does all the dirty-work.
  //*** FirstPass() finds the best alignments (the ones with
  //*** the minimal RMSD) over the space of orientations og.
  //*** It does not calculate any alternate alignments,
  //*** however, it does store promissing looking alignments
  //*** in a table for consideration later.

  FirstPass(og, start, stop);


  // (Necessary for finding alternate solutions)  What a mess!
  if (settings.num_solutions > 1)
  {
    switch (settings.alt_method)
    {
    case Settings::ALT_METHOD_3D:
      CalcAltSolutionsUsing_ALT_METHOD_3D(og, start, stop);
      break;
    case Settings::ALT_METHOD_3D_ORIG:
      CalcAltSolutionsUsing_ALT_METHOD_3D_ORIG(og, start, stop);
      break;
    case Settings::ALT_METHOD_PAIRS:
      CalcAltSolutionsUsing_ALT_METHOD_PAIRS(og, start, stop);
      break;
    default:
      assert(0);
      break;
    }
  }

  // Clean Up (associated with finding alternate solutions)
  if ((settings.num_solutions > 1) && (alt_method_requires_alt_solution_table))
    CleanUpAltSolutionTable();
} // SearchOrDyn3d::Search()




void SearchOrDyn3d::
CalcAltSolutionsUsing_ALT_METHOD_3D(OrientationGenerator const& og,
                                    long start,
                                    long stop)
{
  CalcAltSolutionsUsingOrientation(og);
}



void SearchOrDyn3d::
CalcAltSolutionsUsing_ALT_METHOD_3D_ORIG(OrientationGenerator const& og,
                                         long start,
                                         long stop)
{
  CalcAltSolutionsUsingOrientation(og);
}



namespace CalcAltSolutionsNamespace
{
  int max_num_common_pairs;
  bool NotTooManyPairsInCommon(PairwiseAlignment const& a1,
                               PairwiseAlignment const& a2)
  {
    int e = NumMatchesInBothAssumingMonotonicity(a1, a2);
    int E = max_num_common_pairs;
    //cout << "&max_pairs = "
    //<< &max_num_common_pairs << ";  max_pairs = "
    //<< max_num_common_pairs <<";  e = " << e << ";  E = " << E << endl;
    return (e <= E);
  }
} //namespace CalcAltSolutionsNamespace



void SearchOrDyn3d::
CalcAltSolutionsUsing_ALT_METHOD_PAIRS(OrientationGenerator const& og,
                                       long start,
                                       long stop)
{
  //FindBest() has just
  //completed it's calculation
  //As a courtesy, print out the alignments at this stage to the user.
  //(The user may have to wait for a long time for the final
  // results.)

  //WriteMSF(0); <-- This may change a lot, just print the list of RMSD vs. N
  WriteChimeraInfo(0);

  //Initialize a static global variable needed by the function we will pass to
  //CalcAltSolutionsMultiPass(), the function which finds the alternate solutions.
  CalcAltSolutionsNamespace::max_num_common_pairs = settings.alt_max_pairs_common;

  set<long> use_these_orientations;

  if (settings.alt_rmsd_tolerance != Settings::ALTERNATES_NOT_LIMITED_BY_RMSD)
  {
    //If possible, exploit the alt-solutions-table to reduce the computation
    //time by knowing which orientations don't yeild any low RMSD alignments.
    assert(alt_method_requires_alt_solution_table);
    GetOrientationsOfLowRmsdFromAltTable(og, use_these_orientations);
  }
  else
  {
    //Otherwise, just re-use the entire list of orientations.
    for (long i = start; i != stop; ++i)
      use_these_orientations.insert(i);
  }
  for (int r = 1; r < settings.num_solutions; ++r)
  {
    cout << "----- Searching for the " << r+1 
         << "th-best alignments -----" << endl;

    MultiPass(r,
              og,
              use_these_orientations,
              &CalcAltSolutionsNamespace::NotTooManyPairsInCommon);

    //As a courtesy, print out results while you calculate.
    //WriteMSF(r); <-- This may change a lot, just print the list of RMSD vs. N
    WriteChimeraInfo(r);
  }
} //SearchOrDyn3d::CalcAltSolutionsUsing_ALT_METHOD_PAIRS();






void
SearchOrDyn3d::CalcAltSolutionsUsingOrientation(OrientationGenerator const& og)
{
  //The following variable stores when each orientation produced an optimal,
  //or second-best, or third-best (etc...) solution
  multimap<long, OrientationOwner> orientation_owners;

  //Now mine the contents of the database of alternative orientations
  //and their corresponding RMSD that we have stored in the files
  //stored in the "apAltOrientationFiles[]" array.
  //  We must first through these files to find the orientations which
  //orientations at which the second, third, fourth, etc.. best alignments,
  //were created.
  //(This is done separately for each different number of matches made, n.)
  ReadAltSolutionFiles(og,
                       orientation_owners);

  //Now, once the 'winning' orientations have been found, re-compute
  //the corresponding winning alignments made at those orientations.
  //Note: We have allready calculated these alignments before, but
  //      we did not save them because to save every aligmnet attempted
  //      would have taken up way to much space.  Saving just the
  //      orientations (in the form of an integer index) that produced
  //      these alignments saves a _lot_ of space.
  //  This should fill in the contents of the aaAlignments[][], aaSumSqdDists[][],
  //  and aaTransforms[][] arrays.

  CalcAltSolutionsFromAltSolutionTable(og, orientation_owners);

}//SearchOrDyn3d::CalcAltSolutionsUsingOrientation()




bool
SearchOrDyn3d::MaxDistanceExceeded(int n)
{
  return dyn3d.MaxDistanceExceeded(n);
}



bool
SearchOrDyn3d::MaxRMSDExceeded(Real sum_sqd_dist, int n)
{
  return ((settings.dyn3d_max_rmsd != settings.dyn3d_NO_MAX_RMSD)
          &&
          (sum_sqd_dist > n
                          *
                          Biopolymer::Residue::NumBackboneAtoms()
                          *
                          SQR(settings.dyn3d_max_rmsd)));
}



void SearchOrDyn3d::
KeepTrackOfAltSolutions(int n,
                        long orig_orientation_id,
                        Real sum_sqd_dist,
                        Matrix3x4 orientation_reduced,
                        PairwiseAlignment const& alignment)
{
  assert(settings.num_solutions > 1);

  Real lowest_sum_sqd_dist_so_far = aaSumSqdDists[0][n-n_min];

  // Then if rmsd is within the limit set by the "alt_rmsd_tolerance" ratio,
  // then some statistics about this solution will be saved to disk.
  // Later on, they will be loaded again and this solution will be
  // re-examined to see if it is a potential alternate solution
  // (for example, a 2nd-best or 3rd-best solution).
  if (alt_method_requires_alt_solution_table
      &&
      ((sum_sqd_dist
        <=
        lowest_sum_sqd_dist_so_far * SQR(settings.alt_rmsd_tolerance))
       ||
       (settings.alt_rmsd_tolerance ==
        Settings::ALTERNATES_NOT_LIMITED_BY_RMSD)))
  {
    assert(apAltOrientationFiles);
    assert(apAltOrientationFiles[n-n_min]);
    assert(aNumAltOrientations);

    //What kind of information we keep track of (if any)
    //to describe this alignment, depends on settings.alt_method.
    //The various alt_methods have different trade-offs over
    //the amount of information that must be stored, versus:
    //the ease with which the alternate solutions can be interpreted.
    //(For example, ALT_METHOD_3D requires more temporary file space than
    // ALT_METHOD_3D_ORIG, but is easier to undertand.)
    DEBUG_MSG(DBG_ALT_FILES_WRITE,
              " found potential alternate. Writing to file.\n"
              "  o_id = " << orig_orientation_id <<
              ",  n = " << n <<
              ", rmsd = " <<
              sqrt(sum_sqd_dist
                   /
                   (n * Biopolymer::Residue::NumBackboneAtoms())) <<
              ", rmsd_best = " << 
              sqrt(lowest_sum_sqd_dist_so_far
                   /
                   (n * Biopolymer::Residue::NumBackboneAtoms())));
    aNumAltOrientations[n-n_min]++;
    Real rmsd = sqrt(sum_sqd_dist
                         /
                     (n * Biopolymer::Residue::NumBackboneAtoms()));
    binary_write(*(apAltOrientationFiles[n-n_min]), orig_orientation_id);
    binary_write(*(apAltOrientationFiles[n-n_min]), rmsd);
    if (settings.alt_method == Settings::ALT_METHOD_3D)
    {
      assert(settings.refine_while_searching);
      //If alt.method == Settings::ALT_METHOD_3D,
      //then save the orientation matrix that minimizes
      //RMSD of that alignment.  This matrix will be used later
      //in order to compare the "3d" similarity of the different solutions.
      BinaryWriteMatrix3x4(*(apAltOrientationFiles[n-n_min]),
                           orientation_reduced);
    }
  } //if (sum_sqd_dist <= lowest_sum_sqd_dist * SQR(settings.alt_rmsd_tolerance))
} //SearchOrDyn3d::KeepTrackOfAltSolutions()







void
SearchOrDyn3d::ExtractDynTableSolution(int               n,
                                       PairwiseAlignment &dest_alignment,
                                       Real&             sum_sqd_dist_orig,
                                       Real&             sum_sqd_dist_reduced,
                                       Matrix3x4         orientation_reduced,
                                       Superimpose       *pSuperimpose)
{
  assert((settings.n_min <= n) && (n <= settings.n_max));
  assert(pSuperimpose || (! settings.refine_while_searching));
  assert(orientation_reduced);

  // *** A solution is composed of an alignment, and orientation,
  // *** and the RMSD of that alignment at that orientation.

  // *** First, extract the alignment from the dyn table.
  dest_alignment.resize(n); //make sure there's a enough space
  dyn3d.ExportSolution(n, dest_alignment);

  // *** Now, extract the other information associated with
  // *** this solution (the RMSD, and the orientation).
  sum_sqd_dist_orig = dyn3d.GetCost(n);


  // *** If refine_while_searching == true, also
  // *** calculate the RMSD after optimal superposition.
  if (settings.refine_while_searching)
  {
    assert(pSuperimpose);
    assert(aRotateUsingTheseCoords1 && aRotateUsingTheseCoords2);
    Real rmsd_reduced = dest_alignment.CalcMinRMSD(settings.aMol_f[0],
                                                   settings.aMol_f[1],
                                                   orientation_reduced,
                                                   pSuperimpose,
                                                   aRotateUsingTheseCoords1,
                                                   aRotateUsingTheseCoords2);
    //(translate from RMSD to sum-squared-distance)
    sum_sqd_dist_reduced =
      SQR(rmsd_reduced) * n * Biopolymer::Residue::NumBackboneAtoms();
  } // if (refine_while_searching)
} // SearchOrDyn3d::ExtractDynTableSolution()





  
void SearchOrDyn3d::RecordSolution(int r,
                                   int n,
                                   Real sum_sqd_dist,
                                   ConstMatrix3x4 orientation,
                                   PairwiseAlignment const& alignment)
{
  assert(alignment.NumMatches() == n);
  aaAlignments[r][n-n_min] = alignment;
  aaSumSqdDists[r][n-n_min] = sum_sqd_dist;
  CopyMatrix3x4(orientation, aaTransforms[r][n-n_min]);

#if 0
  // The following code was commented out because I dont
  // keep track of "pivots" anymore.
  // (the pairs of blocks that got used to generate the orientation used)
  // It's still nice to keep this code around.
  if (pAssignPivot) //If the caller bothered to specify "pAssignPivot", then
    aaPivots[r][n-n_min] = *pAssignPivot; //keep track of which pivot was used
                                      //to generate the orientation that
                                      //resulted in this alignment.
  else { //otherwise store (an impossible) default value
    aaPivots[n-n_min].i = -num_consec_matches;
    aaPivots[r][n-n_min].j = -num_consec_matches;
  }
#endif
} //void SearchOrDyn3d::RecordSolution()



void SearchOrDyn3d::
ReadAltSolutionFiles(OrientationGenerator const& og,
                     multimap<long, OrientationOwner>& orientation_owners)
{
  cout <<
    "Choosing the best alternate alignments from database of \".tmp\" files.\n"
    //"  Note: Each alternate alignment must have been made at \n"
    //"        orientation which placed the 2nd molecule in a\n"
    //"        position no closer than "
    //   << alt_o_min_delta_orientation_RMSD <<
    //"angstroms (RMSD)\n"
    //"        to any of its other positions in alternate alignments...\n"
       << flush;

  for (int n = n_min; n <= n_max; ++n)
  {
    
    if (aaSumSqdDists[0][n-n_min] != DistanceMetric::UPPER_BOUND)//(make sure that 
    {     //there exists at least one valid alignment at this value of n first)
      ReadAltSolutionFile(n,
                          *(apAltOrientationFiles[n-n_min]),
                          aNumAltOrientations[n-n_min],
                          og,
                          orientation_owners);
    }
  }
} //SearchOrDyn3d::ReadAltSolutionFiles()





//UGH!
// Because of the two different modes for generating alternate alignments
// (ALT_METHOD_3D, and ALT_METHOD_3D_ORIG), this function is a mess.

void SearchOrDyn3d::
ReadAltSolutionFile(int n,
                    fstream& alt_file,
                    int num_alt_orientations,
                    OrientationGenerator const& og,
                    multimap<long, OrientationOwner>&
                    orientation_owners)
{
  assert(alt_file);
  assert(num_alt_orientations >= 1);
  alt_file.seekg(0); //might as well rewind (unnecessary?)

  //vInstances is a vector that stores information for
  //all the solutions that have ever been calculated containing
  //n matched pairs of residues!
  vector<OrientationUsageInstance> vInstances;
  vInstances.reserve(num_alt_orientations);

  //*** Load all the solutions from the file into the vInstances[]
  //*** vector so the solutions can be sorted according their RMSD.
  //  Once we have read all the orientations into the vInstances[] array,
  //we will sort the array by the RMSD of the alignments (with n matches)
  //created at this orientation so that the orientations which yield
  //the best alignments are at the beginning.  Alignments that are
  //"too close together" will be discarded from this array.
  //Then, later on, we will call CalcAltSolutionsFromAltSolutionTable()
  //to calculate all the alignments corresponding to orientations
  //stored in the vInstances[i] array, where i is an integer >= 1.
  //We will skip the first element vInstances[0], because this was
  //suposedly calculated in FirstPass.
  //   ------  Weird Quirk: ------
  //  This means we need to make sure that the orientation
  //that yields the "best" alignment with this number of
  //matches (stored in aaTransforms[0][n-n_min])
  //is actually located at the beginning of this array.
  //The problem is after sorting the the array, using STL's sort function,
  //this might not be the case.  If there are multiple alignments with
  //this same lowest RMSD, their relative ordering in the array is random.
  //To insure that aaTransforms[0][n-n_min] is stored in the first element
  //of the vInstances[] array, I use the "best_orientation" and
  //"best_orientation_found" variables.

  bool best_orientation_found = false;
  OrientationUsageInstance best_orientation;

  DEBUG_MSG(DBG_ALT_ORIENTATION,
            "Loading Alt file for n=" << n << " into vector of size "
            << num_alt_orientations);

  while (alt_file)
  {
    OrientationUsageInstance o;
    binary_read(alt_file, o.orig_orientation_id);
    binary_read(alt_file, o.rmsd);
    if (settings.alt_method == Settings::ALT_METHOD_3D)
    {
      assert(settings.refine_while_searching);
      BinaryReadMatrix3x4(alt_file, o.orientation_reduced);
    }

    if (alt_file) //make sure last read was successful
    {
      bool just_found_best;
      switch (settings.alt_method)
      {
      case Settings::ALT_METHOD_3D:
        just_found_best = EqualMatrices3x4(aaTransforms[0][n-n_min],
                                           o.orientation_reduced);
        break;
      case Settings::ALT_METHOD_3D_ORIG:
        just_found_best = EqualMatrices3x4(aaTransforms[0][n-n_min],
                                           og[o.orig_orientation_id].transform);
        break;
      default:
        assert(0);
        break;
      }
      just_found_best = (just_found_best && (! best_orientation_found));
      if (just_found_best)
      {
        best_orientation = o;
        best_orientation_found = true;
      }
      else
        vInstances.push_back(o);
    } //if (alt_file)
  } //while (alt_file)

  DEBUG_MSG(DBG_ALT_ORIENTATION, "n=" << n << ": "
            "Size of vector is " << vInstances.size() << "\n"
            "                  Sorting vector by RMSD\n");

  //Now, sort the list of solutions by RMSD.
  //(Solutions with lowest RMSD will be first, highest last.)
  sort(vInstances.begin(), vInstances.end());
  vInstances.insert(vInstances.begin(), best_orientation);

  assert((settings.alt_method != Settings::ALT_METHOD_3D)
         ||
         (EqualMatrices3x4(vInstances[0].orientation_reduced,
                           aaTransforms[0][n-n_min])));
  assert((settings.alt_method != Settings::ALT_METHOD_3D_ORIG)
         ||
         (EqualMatrices3x4(og[vInstances[0].orig_orientation_id].transform,
                           aaTransforms[0][n-n_min])));

  assert(vInstances.size() <= num_alt_orientations);

  long *aAltOrigOrientationIDs = NULL;
  Matrix3x4 *aAltReducedOrientations = NULL;
  switch (settings.alt_method)
  {
  case Settings::ALT_METHOD_3D:
    assert(settings.refine_while_searching);
    aAltReducedOrientations = new Matrix3x4 [settings.num_solutions];
    CHECK_ALLOC(aAltReducedOrientations);
    //The first entry in the vInstances array has the solution with lowest RMSD
    //The first element of the array below should refer to this solution.
    CopyMatrix3x4(vInstances[0].orientation_reduced,
                  aAltReducedOrientations[0]);
    assert(EqualMatrices3x4(aAltReducedOrientations[0],
                            aaTransforms[0][n-n_min]));
    break;
  case Settings::ALT_METHOD_3D_ORIG:
    assert(! settings.refine_while_searching);
    aAltOrigOrientationIDs = new long [settings.num_solutions];
    CHECK_ALLOC(aAltOrigOrientationIDs);
    //The first entry in the vInstances array has the solution with lowest RMSD
    //The first element of the array below should refer to this solution.
    aAltOrigOrientationIDs[0] = vInstances[0].orig_orientation_id;
    assert(EqualMatrices3x4(og[aAltOrigOrientationIDs[0]].transform,
                            aaTransforms[0][n-n_min]));
    break;
  default:
    assert(0);
  }
  Real best_rmsd = vInstances[0].rmsd;

  int alt_count = 1; //Indicates whether we are presently looking for the
                     //first best(0), second best(1), third best(2), etc..
                     //alignment.
                     //
                     //We set it to 1, not 0, because we should have
                     //allready calculated the first-best alignments allready,
                     //so, we don't want to recalculate them!
                     //(This is also why we begin with current_orientation=1
                     // instead of 0.)
  
  for(vector<OrientationUsageInstance>::iterator pOUI = vInstances.begin()+1;
      ((pOUI != vInstances.end())
       && (alt_count < settings.num_solutions)
       && ((pOUI->rmsd <= best_rmsd * settings.alt_rmsd_tolerance)
           ||
           (settings.alt_rmsd_tolerance == Settings::ALTERNATES_NOT_LIMITED_BY_RMSD)
          )
      );
      ++pOUI)
  {
    //Now, check to see if all the previous alternate-best solutions are
    //sufficiently different than this one to be considered.
    //If this candidate solution is not different enough than all of them,
    //it will be discarded.
    DEBUG_MSG(DBG_ALT_ORIENTATION,
              "n=" << n <<
              ", o_id=" << pOUI->orig_orientation_id <<
              ":\n"
              "The following alt solutions are sufficiently different from candidate:");

    int r;
    for(r = 0;
        ((r < alt_count) &&
         (SufficientlyDifferentOrientationFromRthAlternate(
                     *pOUI,
                     aAltOrigOrientationIDs,
                     aAltReducedOrientations,
                     r,
                     og)));
        ++r)
    {
      #ifdef DEBUG
      #if DBG_ALT_ORIENTATION
      if (r != 0)
        cerr << ", ";
      cerr << "(n=" << n << ", r=" << r+1 << ")"
           << flush;
      #endif //#if DBG_ALT_ORIENTATION
      #endif //#ifdef DEBUG
    }
    #ifdef DEBUG
    #if DBG_ALT_ORIENTATION
    cerr << endl;
    #endif //#if DBG_ALT_ORIENTATION
    #endif //#ifdef DEBUG
    
    //If they are, then this is orientation itself, produces
    //an alternate-best-solution, so append it to the list.
    if (r == alt_count)
    {
      OrientationOwner entry;
      entry.n = n;
      entry.r = alt_count;
      orientation_owners.insert(make_pair(pOUI->orig_orientation_id, entry));
      DEBUG_MSG(DBG_ALT_ORIENTATION,
                "successful: "
                << r+1
                << "th alternate solution for "
                "n=" << n << " is candidate#" << pOUI->orig_orientation_id);

      switch (settings.alt_method)
      {
      case Settings::ALT_METHOD_3D:
        CopyMatrix3x4(pOUI->orientation_reduced,
                      aAltReducedOrientations[alt_count]);
        break;
      case Settings::ALT_METHOD_3D_ORIG:
        aAltOrigOrientationIDs[alt_count] = pOUI->orig_orientation_id;
        break;
      default:
        assert(0);
      }

      ++alt_count;
    } //if (r == alt_count)
    #ifdef DEBUG
    #if DBG_ALT_ORIENTATION
    else
    {
      cerr << pOUI->orig_orientation_id <<
        " is too similar to: ";
      switch (settings.alt_method)
      {
      case Settings::ALT_METHOD_3D_ORIG:
        cerr << "RMSD("
             << pOUI->orig_orientation_id << ","
             << aAltOrigOrientationIDs[r] << ") = " <<
          sqrt(og.SumSqdDistBetweenOrientations(pOUI->orig_orientation_id,
                                                aAltOrigOrientationIDs[r])
               /
               (aMol_c[1].size() * Biopolymer::Residue::NumBackboneAtoms()))
           << endl;
        break;
      case Settings::ALT_METHOD_3D:
        cerr << "RMSD(n=" << n << ",r=" << r+1 << ")= " <<
          sqrt(og.SumSqdDistBetweenOrientations(pOUI->orientation_reduced,
                                                aAltReducedOrientations[r])
               /
               (aMol_c[1].size() * Biopolymer::Residue::NumBackboneAtoms()))
           << endl;
        
      } //switch (settings.alt_method)
    } //else clause for "if (r == alt_count)"
    #endif //#if DBG_ALT_ORIENTATION
    #endif //#ifdef DEBUG

  } //loop over orientations (or "instances") read from alt_file

  //Clean up:
  if (aAltOrigOrientationIDs) delete [] aAltOrigOrientationIDs;
  if (aAltReducedOrientations) delete [] aAltReducedOrientations;
} //SearchOrDyn3d::ReadAltSolutionFile()




//The following function is only ever called by ReadAltSolutionFile().
//See ReadAltSolutionFile() for usage context.
bool SearchOrDyn3d::
SufficientlyDifferentOrientationFromRthAlternate(OrientationUsageInstance const& oui,
                                                 long      *aAltOrigOrientationIDs,
                                                 Matrix3x4 *aAltReducedOrientations,
                                                 int  r,
                                                 OrientationGenerator const& og)
{
  Real min_delta_orientation_sum_sqd_dist = 
    SQR(settings.alt_min_3d_difference)
    * aMol_c[1].size()
    * Biopolymer::Residue::NumBackboneAtoms();

  Real sumSqdDistBetweenSolutions = 0.0;
  switch (settings.alt_method)
  {
  case Settings::ALT_METHOD_3D:
    assert(aAltReducedOrientations);
    sumSqdDistBetweenSolutions =
      og.SumSqdDistBetweenOrientations(oui.orientation_reduced,
                                       aAltReducedOrientations[r]);
    break;
  case Settings::ALT_METHOD_3D_ORIG:
    assert(aAltOrigOrientationIDs);
    sumSqdDistBetweenSolutions =
      og.SumSqdDistBetweenOrientations(oui.orig_orientation_id,
                                       aAltOrigOrientationIDs[r]);
    break;
  default:
    assert(0);
  }
  return (sumSqdDistBetweenSolutions >= min_delta_orientation_sum_sqd_dist);
}





void
SearchOrDyn3d::GetOrientationsOfLowRmsdFromAltTable(OrientationGenerator const& og,
                                                    set<long>& desirable_orientations)
{
  assert(settings.num_solutions > 1);
  assert(apAltOrientationFiles);
  assert(aNumAltOrientations);

  cout <<
    "Calculating alternate alignments by searching over all orientations again\n"
    "  Filtering orientations based on RMSDs of alignments calculated\n"
    "  Originally the search consisted of " << og.size() <<
    " orientations." << endl;

  //Read in all the orientations in all the orientation temp-files.
  //Check to see if the RMSD of any of the alignments at that orientation
  //fall within alt_rmsd_tolerance of any of the the lowest-RMSD-alignments
  //with the same number of matches.
  //If so, then that orientation is appended to the list of orientations.
  //Otherwise it is discarded.

  //first loop over all orientation files
  for(int n = n_min; n <= n_max; ++n)
  {
    //(The next if statement may be unnecessary.
    //It just insures that valid 1st-best solutions
    //with this number of matches (n) have been already been calculated that
    //satisfy all the constraints.  If not, then theres no sense examining
    //the other orientations that might yeild alternate suboptimal alignments
    //with the same number of matches.  They will surely not satisfy the
    //constraints either.)
    if (aaSumSqdDists[0][n-n_min] != DistanceMetric::UPPER_BOUND)
    {
      assert(apAltOrientationFiles[n-n_min]);
      apAltOrientationFiles[n-n_min]->seekg(0); //might as well rewind (unnecessary?)
      //then loop over all orientations in each file
      while (*apAltOrientationFiles[n-n_min])
      {
        Real best_rmsd = sqrt(aaSumSqdDists[0][n-n_min]/
                               (n * Biopolymer::Residue::NumBackboneAtoms()));

        OrientationUsageInstance o;
        binary_read(*apAltOrientationFiles[n-n_min], o.orig_orientation_id);
        binary_read(*apAltOrientationFiles[n-n_min], o.rmsd);

        Real tolerable_rmsd = settings.alt_rmsd_tolerance * best_rmsd;

        if (settings.alt_method == Settings::ALT_METHOD_3D)
        {  
          BinaryReadMatrix3x4(*apAltOrientationFiles[n-n_min],
                              o.orientation_reduced);
        }

        if ((tolerable_rmsd >= o.rmsd)
            ||
            (settings.alt_rmsd_tolerance ==
             Settings::ALTERNATES_NOT_LIMITED_BY_RMSD))
        {
          DEBUG_MSG(DBG_ALT_FILES_READ,
                    "Desireable orientation found:"
                    " retaining orientationID #" <<o.orig_orientation_id<< "\n"
                    "  ...produced an alignment with n=" << n <<
                    " matches and RMSD=" << o.rmsd <<
                    "; where best_RMSD=" << best_rmsd);
          desirable_orientations.insert(o.orig_orientation_id);
        }
      } //while (*apAltOrientationFiles[n-n_min])
    } //if (aaSumSqdDists[0][n-n_min] != DistanceMetric::UPPER_BOUND)
  } //for(int n = n_min; n <= n_max; ++n)

  cout <<
    "  ...After the RMSD-filter, the number of orientations\n"
    "     we have to consider is: "
       << desirable_orientations.size() << endl;

} //GetOrientationsOfLowRmsdFromAltTable()




bool SearchOrDyn3d::
AltSolutionsSufficientlyDissimilar(int n,
                                   int r,
                                   int r_prev)
{
  if ((settings.alt_method == Settings::ALT_METHOD_3D) ||
      (settings.alt_method == Settings::ALT_METHOD_3D_ORIG))
  {
    int num_points = aMol_c[1].size()
      *
      Biopolymer::Residue::NumBackboneAtoms();

    Real min_delta_orientation_sum_sqd_dist = 
      SQR(settings.alt_min_3d_difference) * num_points;

    //Just in case, somebody's been messing with the contents of
    //the aRotateUsingTheseCoords2[] array, I fill it again
    int n_count = 0;
    for(Biopolymer::const_iterator ps = settings.aMol_f[1].begin();
        ps < settings.aMol_f[1].end();
        ++ps)
    {
      for(int q=0; q < Biopolymer::Residue::NumBackboneAtoms(); ++q)
      {
        Biopolymer::Residue::const_iterator pa;
        pa = ps->GetBackboneAtom(q);
        assert(pa != ps->end());
        aRotateUsingTheseCoords2[n_count][0] = (*pa).second.xyz[0];
        aRotateUsingTheseCoords2[n_count][1] = (*pa).second.xyz[1];
        aRotateUsingTheseCoords2[n_count][2] = (*pa).second.xyz[2];
        ++n_count;
      }
    }

    Real sum_sqd_distance =
      SlowRotSumSqdDist(aaTransforms[r_prev][n-n_min],
                        aaTransforms[r][n-n_min],
                        aRotateUsingTheseCoords2,
                        num_points);

    return (sum_sqd_distance >= min_delta_orientation_sum_sqd_dist);
  }
  else 
  {
    assert(settings.alt_method == Settings::ALT_METHOD_PAIRS);
    int num_identical_matches =
        NumMatchesInBothAssumingMonotonicity(aaAlignments[r_prev][n-n_min],
                                             aaAlignments[r][n-n_min]);
    return num_identical_matches <= settings.alt_max_pairs_common;
  }
} //SearchOrDyn3d::AltSolutionsSufficientlyDissimilar()



void
SearchOrDyn3d::CheckAlternatesForConsistency()
{
  cout << "Checking the list of alternate solutions for consistency." << endl;

  int num_points = aMol_c[1].size() *
    Biopolymer::Residue::NumBackboneAtoms();

  Real min_delta_orientation_sum_sqd_dist = 
    SQR(settings.alt_min_3d_difference) * num_points;
  
  for(int n = n_min;
      n <= n_max;
      ++n)
  {
    SortAltList(n);

    //Now loop through list of alternate alignments (for this number of
    //matches made, n) looking for alternate alignments
    //"that are too close together" and remove them.
    for(int r = 1;
        ((r < settings.num_solutions) &&
         (aaSumSqdDists[r][n-n_min] < DistanceMetric::UPPER_BOUND));
        ++r)
    {
      bool delete_this_solution = false; //if true, solution is deleted later



      //Solutions can be discarded if they are too "similar"
      //to any of the other solutions.  I check for this here.
      //If an alignment is "too similar" to any of the other
      //alignments with lower rmsd, and the same number of matches,
      //I delete it from the list.

      int r_prev;
      for(r_prev = 0;
          (r_prev < r)
            &&
           AltSolutionsSufficientlyDissimilar(n, r, r_prev);
          ++r_prev)
      {}

      if (r_prev < r)
      {
        delete_this_solution = true;
        cout <<
          " ...The " << r+1 << "th-best solution with " << n << " matches was too\n"
          "    similar to other, lower-RMSD solutions.  It will be eliminated.\n"
             << flush;

        //print out some debug information:
        
        #ifdef DEBUG
        switch(settings.alt_method)
        {
        case Settings::ALT_METHOD_3D:
        case Settings::ALT_METHOD_3D_ORIG:
          DEBUG_MSG(DBG_ALTERNATES," CheckAlternatesForConsistency(): \n"
                    "sum_sqd_dist(aaTransforms["<< r_prev <<"]["<<n-n_min<<"],\n"
                    "             aaTransforms["<< r <<"]["<<n-n_min<<"]) = "
                    << SlowRotSumSqdDist(aaTransforms[r_prev][n-n_min],
                                         aaTransforms[r][n-n_min],
                                         aRotateUsingTheseCoords2,
                                         num_points) <<
                    " < " << min_delta_orientation_sum_sqd_dist << "\n"
                    "This solution will be replaced by it's successors.");
          

          break;
        case Settings::ALT_METHOD_PAIRS:
          DEBUG_MSG(DBG_ALTERNATES," CheckAlternatesForConsistency(): \n"
                    "NumIdenticalMatches(aaAlignments["<< r_prev <<"]["<<n-n_min<<"],\n"
                    "                    aaAlignments["<< r <<"]["<<n-n_min<<"]) = "
                    << NumMatchesInBothAssumingMonotonicity(
                                             aaAlignments[r_prev][n-n_min],
                                             aaAlignments[r][n-n_min])
                    << " > "
                    << settings.alt_max_pairs_common << "\n"
                    "This solution will be replaced by it's successors.");
          break;
        default:
          assert(0);
          break;
        } //switch(settings.alt_method)
        #endif //#ifdef DEBUG
      } // if (r_prev < r)

      //If an upper bound was placed on the RMSD of the
      //solutions, make sure the rmsd of these solutions is below
      //that rmsd.  Otherwise this solution should be discarded.
      if (settings.alt_rmsd_tolerance !=
          Settings::ALTERNATES_NOT_LIMITED_BY_RMSD)
      {
        //Check to see if the RMSD exceeds the
        Real rmsd = sqrt(aaSumSqdDists[r][n-n_min]
                          /
                         (n*Biopolymer::Residue::NumBackboneAtoms()));

        Real best_rmsd = sqrt(aaSumSqdDists[0][n-n_min]
                              /
                              (n*Biopolymer::Residue::NumBackboneAtoms()));
        Real tolerable_rmsd = best_rmsd * settings.alt_rmsd_tolerance;

        if (rmsd > tolerable_rmsd)
        {
          delete_this_solution = true;

          DEBUG_MSG(DBG_ALTERNATES," CheckAlternatesForConsistency(): \n"
                    "      "
                    << r+1 <<
                    "th best solution with "
                    << n <<
                    " matches deleted.\n"
                    "      RMSD = "
                    << rmsd <<
                    "  >  rmsd_ratio * best_rmsd = "
                    << settings.alt_rmsd_tolerance << " * " << best_rmsd);
        } // if (rmsd > tolerable_rmsd)
      } // if (alt_rmsd_tolerance != ALTERNATES_NOT_LIMITED_BY_RMSD)


      if (delete_this_solution)
      {
        //...then erase this solution by copying over it with the solutions
        //that come later.
        for(int r_copy = r;
            r_copy+1 < settings.num_solutions;
            ++r_copy)
        {
          assert((aaSumSqdDists[r_copy][n-n_min] != DistanceMetric::UPPER_BOUND)
                 ||
                 (aaSumSqdDists[r_copy+1][n-n_min] == DistanceMetric::UPPER_BOUND));
          DEBUG_MSG(DBG_ALTERNATES,
                    "  replacing aaAlignments["<< r_copy <<"]["<<n-n_min<<"]"
                    " with aaAlignments["<< r_copy+1 <<"]["<<n-n_min<<"])\n");
          aaSumSqdDists[r_copy][n-n_min] = aaSumSqdDists[r_copy+1][n-n_min];
          aaAlignments[r_copy][n-n_min] = aaAlignments[r_copy+1][n-n_min];
          CopyMatrix3x4(aaTransforms[r_copy+1][n-n_min],aaTransforms[r_copy][n-n_min]);
        }
        //void the solution at the end.
        aaSumSqdDists[settings.num_solutions-1][n-n_min] = DistanceMetric::UPPER_BOUND;
        CopyMatrix3x4(g_IDENTITY3X4, aaTransforms[settings.num_solutions-1][n-n_min]);
        --r; //I do this so I won't skip over the next orientation we have
             //copied into the current slot.
      } //if (delete_this_solution)
    } //for(int r = 1; (r < settings.num_solutions) &&...
  } //for(int n = n_min; n <= n_max; ++n)
} // SearchOrDyn3d::CheckAlternatesForConsistency()




void
SearchOrDyn3d::CalcAltSolutionsFromAltSolutionTable(
                                    OrientationGenerator const& og,
                                    multimap<long, OrientationOwner>&
                                           orientation_owners)
{
  int orientation_counter = 0; //These two variables are used for printing
  int alignment_counter   = 0; //out debug information to the user.

  cout <<
    "Recalculating the candidates for alternate alignments\n"
    "from the data in the \".tmp\" files."
    << endl;

  Superimpose superimpose(max_num_matches
                          * Biopolymer::Residue::NumBackboneAtoms());

  //Now scan the list of orientations, and for each orientation,
  //search for the n,r pairs (ie. n = number of matches, r = rank of solution)
  //that each orientation generates.
  //For each n,r-pair save an alignment in the appropriate place in
  //the aaSumSqdDists[r][n-n_min], aaAlignments[r][n-n_min], and aaTransforms[r][n-n_min] arrays.
  multimap<long,OrientationOwner>::iterator pMMOO
    = orientation_owners.begin();
  while(pMMOO != orientation_owners.end())
  {
    long o_id = pMMOO->first;
    DEBUG_MSG(DBG_ALT_ORIENTATION,
              " Calculating alternate-alignments from orientation id# "
              << o_id);

    Matrix3x4 const& orientation_orig
      = og[o_id].transform;

    //Rotate the second molecule to the desired orientation.
    ApplyTransformBackbone(orientation_orig,
                           settings.aMol_f[1],
                           aMol_c[1]);

    //The relative positions are different so we have to recalculate the
    //squared-distances between the residues in the two sequences.
    dyn3d.FillLookupTable(aMol_c[0], aMol_c[1]);

    // Initialize the first layer of the dynamic-programming-table
    // (required before calling Solve_IncrN())
    dyn3d.InitTable_IncrN();
    //    DEBUG_MSG(DBG_SEARCH_OR_DYN3D, "Iteration #" << orientation_counter
    //              << " got passed Dyn3d->InitTable_IncrN()");

    #ifdef DEBUG
    pair<multimap<long,OrientationOwner>::iterator,
         multimap<long,OrientationOwner>::iterator>
      solutions_for_this_orientation =
      orientation_owners.equal_range(o_id);
    assert(pMMOO == solutions_for_this_orientation.first);
    #endif //#ifdef DEBUG

    //Figure out how many alignments of interest (or "owners",
    //as I've taken to calling them), were generated at this
    //orientation.  Sort this subset, by number-of-matches-made, n.
    //To do this, we must copy this to another container that allows sorting.

    vector<OrientationOwner> vSolutionsForThisOrientation;
    vSolutionsForThisOrientation.reserve(max_num_matches);
    //copy into the vector.
    for(;
        ((pMMOO != orientation_owners.end()) && (pMMOO->first == o_id));
        ++pMMOO)
    {
      #ifdef DEBUG
      OrientationOwner oo = pMMOO->second;
      #endif //#ifdef DEBUG
      vSolutionsForThisOrientation.push_back(pMMOO->second);
    }
    assert(pMMOO == solutions_for_this_orientation.second);

    DEBUG_MSG(DBG_ALT_ORIENTATION,
              " At orientation " << o_id << " there are "
              << vSolutionsForThisOrientation.size() << " solutions:");

    //sort the vector (by number-of-matches-made, n)
    sort(vSolutionsForThisOrientation.begin(),
         vSolutionsForThisOrientation.end());

    #ifdef DEBUG
    #if DBG_ALT_ORIENTATION
    //For debugging purposes, print out stats for all the solutions
    //(alternate or otherwise) at this orientation.
    //"n" indicates the number of matches in the solution
    //"r" indicates the rank of the solution (1st-best, 2nd-best, 3rd, etc...)
    for(vector<OrientationOwner>::iterator p= vSolutionsForThisOrientation.begin();
        p != vSolutionsForThisOrientation.end();
        ++p)
    {
      if (p != vSolutionsForThisOrientation.begin())
        cerr << ", ";
      cerr << "(n=" << p->n << ", r=" << p->r+1 << ")" << flush;
    }
    cerr << endl;
    #endif //#if DBG_ALT_ORIENTATION
    #endif //#ifdef DEBUG

    //   Now generate the alignments for this orientation.
    //Successive alignments will have increasing numbers of matches made
    //(ie. increasing values of n).  If we calculate an alignment whose
    //number of matches, n, equals the "n" field for one of the entries of,
    //the vector "vSolutionsForThisOrientation", then we save this solution
    //in the aaAlignments[][] array. (Note: the rmsd, and associated
    //orientation also get saved in the aaSumSqdDists[][], and aaTransforms[][]
    //arrays.)  So as we are increasing n, we step through the
    //"vSolutionsForThisOrientation" vector, which we can do because it
    //is sorted by n.

    int n;
    vector<OrientationOwner>::iterator pVOO
      = vSolutionsForThisOrientation.begin();
    for(n = 2;
        ((n <= n_max)
         &&
         (pVOO != vSolutionsForThisOrientation.end()));
        ++n)
    {
      dyn3d.Solve_IncrN(n);

      #ifdef DEBUG
      if (n >= n_min)
        DEBUG_MSG(DBG_ALT_ORIENTATION,
                  "Orientation=" << o_id <<
                  " recalc alt solutions: n = " << n << ", rmsd = "
                  << sqrt(dyn3d.GetCost(n)
                          /
                          (n * Biopolymer::Residue::NumBackboneAtoms())));
      
      #endif //#ifdef DEBUG

      //If the alignment we just calculated has one of the values of n
      //(number of matches made) that we're looking for, save it
      //in the aaAlignments[][] array.

      if (pVOO->n == n)
      {
        assert((n_min <= n) && (n <= n_max));
        int r = pVOO->r;
        //print a debug message, and do some debug-checking
        DEBUG_MSG(DBG_ALT_ORIENTATION,
                  r+1 << "th-best alignment with "<<
                  n << " matches found at orientation#"
                  << o_id+1);

        //The assert() below will check to see if the solution for r,n has been
        //assigned yet.  During the life of this function, each solution
        //should only be calculated once.  This will check to make sure it is.
        assert(aaSumSqdDists[r][n-n_min] == DistanceMetric::UPPER_BOUND);

        alignment_temp.resize(n);
        dyn3d.ExportSolution(n, alignment_temp);

        Real sum_sqd_dist;
        Real rmsd;
        Matrix3x4 orientation_reduced;
        if (settings.refine_while_searching)
        {
          rmsd = alignment_temp.CalcMinRMSD(settings.aMol_f[0],
                                            settings.aMol_f[1],
                                            orientation_reduced,
                                            &superimpose,
                                            aRotateUsingTheseCoords1,
                                            aRotateUsingTheseCoords2);

          sum_sqd_dist = SQR(rmsd)*n *Biopolymer::Residue::NumBackboneAtoms();
                             //(translate from RMSD to sum-squared-distance)
          DEBUG_MSG(DBG_REFINE_OR_DYN3D,
                    " ..transform found.  New cost = "<< sum_sqd_dist);
        } // if (refine_while_searching)
        else {
          sum_sqd_dist = dyn3d.GetCost(n);
        } //else clause for if (refine_while_searching)

        //Now, (finally) write to the array of solutions.
        aaSumSqdDists[r][n-n_min] = sum_sqd_dist;
        aaAlignments[r][n-n_min] = alignment_temp;
        if (settings.refine_while_searching)
          CopyMatrix3x4(orientation_reduced, aaTransforms[r][n-n_min]);
        else
          CopyMatrix3x4(orientation_orig, aaTransforms[r][n-n_min]);

        //now, find the next owner (ie. (n,r) pair) of this orientation.
        ++pVOO;
        //for debugging purposes, count the number of
        //alternate alignments computed so far.
        ++alignment_counter;

      }//if (pVOO->n == n)..ie. if this is one of the alignments we're seeking.
    } //for(n=2; (n<=n_max) && (! terminate_early); ++n)
    ++orientation_counter;
  } //while(pMMOO != orientation_owners.end())
  cout << "During re-calculation to find the alternate alignments,\n"
       << orientation_counter << " orientations were tried, and\n"
       << alignment_counter << " alternate alignments were found\n" << flush;

} //SearchOrDyn3d::CalcAltSolutionsFromOrientationFileData()





void
SearchOrDyn3d::InitAltSolutionTable()
{
  assert(settings.num_solutions > 0);

  if (settings.num_solutions == 1)
  {
    apAltOrientationFiles = NULL;
    aNumAltOrientations = NULL;
    return;
  }

  //Now, deal with setting up the files necessary to handle
  //alternative-best alignments.
  rlimit num_allowed_files;
  getrlimit(RLIMIT_NOFILE, &num_allowed_files);
  DEBUG_MSG(DBG_ALT_FILES_WRITE,
            "Number of simultaneous open file descriptors allowed\n"
            "soft_limit:   " << num_allowed_files.rlim_cur <<
            "  hard_limit: " << num_allowed_files.rlim_max);

  if (range_of_n > num_allowed_files.rlim_max)
    ERR("ERROR: In order to use the \"--alt-orientation\" option,\n"
              "       (or the \"--alt-residues\" option with the third argument specified),\n"
                  "       the number of different \"1st-choice\" alignments that will\n"
                  "       generated (usually this is just the number of residues\n"
                  "       be in the shorter seqence.  In this case, it is currently: "
                  << range_of_n << "),\n"
                  
                  "       cannot exceed the number of simultaneous open\n"
                  "       file-descriptors you are allowed in your operating\n"
                  "       system.  Your operating system presently allows you to open\n"
                  "       "<< num_allowed_files.rlim_cur <<" files simultaneously.\n"
                  "       This can be extended to "
                  << num_allowed_files.rlim_max << ", by a superuser, but it still\n"
                  "       falls short of the " << range_of_n <<
                  " descriptors you need.\n"
                  "           You will either have to\n" 
                  "       1) stop using the \"--alt-orientation\" option.\n"
                  "       2) reduce the size of your proteins, perhaps by removing\n"
                  "          some non-essential residues\n"
                  "       3) reduce the number of aligmnents generated by using the\n"
                  "          \"--cutoff-min-n\" and \"cutoff-max-n\" command-line\n"
                  "          options\n"
                  "       4) try this program on a different operating system\n"
                  "\n"
                  "       Exiting...\n");
  else if (range_of_n > num_allowed_files.rlim_cur)
    ERR("ERROR: In order to use the \"--alt-orientation\" option,\n"
              "       (or the \"--alt-residues\" option with the third argument specified),\n"
                  "       the number of different \"1st-choice\" alignments that will\n"
                  "       be generated (usually this is just the number of residues\n"
                  "       in the shorter seqence.  In this case, it is currently: "
                  << range_of_n << ")\n"
                  
                  "       cannot exceed the number of simultaneous open\n"
                  "       file-descriptors you are allowed in your operating\n"
                  "       system.  Your operating system grants you access to\n"
                  "       "<< num_allowed_files.rlim_cur <<" simultaneous files\n"
                  "       which falls short of the " << range_of_n <<
                  " descriptors you need.\n"
                  "             However, a superuser can extend this to "
                  << num_allowed_files.rlim_max << ".\n"
                  "           You will either have to\n" 
                  "       1) use the \"unlimit\" command to increase the number of files you are\n"
                  "          allowed to open at one time.\n"
                  "       2) stop using the \"--alt-orientation\" option.\n"
                  "       3) reduce the size of your proteins, perhaps by removing\n"
                  "          some non-essential residues\n"
                  "       4) reduce the number of aligmnents generated by using the\n"
                  "          \"--cutoff-min-n\" and \"cutoff-max-n\" command-line\n"
                  "          options\n"
                  "       5) try this program on a different operating system\n"
                  "\n"
                  "       Exiting...\n");

  apAltOrientationFiles = new fstream* [range_of_n];
  aNumAltOrientations = new long [range_of_n];
  CHECK_ALLOC(apAltOrientationFiles && aNumAltOrientations);

  for (int n = n_min; n <= n_max; ++n)
  {
    stringstream alt_filename(BASE_FILE_NAME);
    alt_filename << n << ".tmp";
    apAltOrientationFiles[n-n_min] = new
      fstream(alt_filename.str().c_str(), ios::in | ios::out | ios::trunc);
    CHECK_ALLOC(apAltOrientationFiles[n-n_min]);
    aNumAltOrientations[n-n_min] = 0;
  }
} //SearchOrDyn3d::InitAltSolutionTable()



void
SearchOrDyn3d::CleanUpAltSolutionTable()
{
  if (apAltOrientationFiles)
  {
    assert(settings.num_solutions > 1);
    //First, close the temporary files.
    for(int n = n_min; n <= n_max; ++n)
    {
      if (apAltOrientationFiles[n-n_min])
        delete apAltOrientationFiles[n-n_min];
    }
    delete [] apAltOrientationFiles;
    apAltOrientationFiles = NULL;

    //now delete (unlink) the temporary files from the filesystem
    for(int n = n_min; n <= n_max; ++n)
    {
      stringstream alt_filename(BASE_FILE_NAME);
      alt_filename << n << ".tmp";
      if (unlink(alt_filename.str().c_str()))
        cerr << "Warning: There was an error trying to delete temporary file \""
             << alt_filename.str() << "\"" << endl;
    }
  }
  //  else 
  //    assert(settings.num_solutions == 1);
} //SearchOrDyn3d::CleanUpAltSolutionTable()






void
SearchOrDyn3d::SortAltList(int n)
{
  //Save the the alignment that was spawned by the lowest-RMSD solution
  //so far with n matches, to see if it get's disturbed by the sorting.
  // (This is interesting but unnecessary. I may want to take this out later.)
  alignment_temp.resize(n);
  alignment_temp = aaAlignments[0][n-n_min];

  //Why not prealocate some things we will use a lot later.
  vector<SolutionInfo> vSolutions;
  vSolutions.reserve(settings.num_solutions);
  SolutionInfo s;
  s.alignment.reserve(n);

  //In order to sort all the solutions according to RMSD,
  //copy the solutions arrays into a temporary C++-vector
  //(Using an STL vector to do the sorting saves me writing a lot of code.)
  for(int r=0;
      (r < settings.num_solutions) // let r range from 0...(settings.num_solutions-1)
        &&
      (aaSumSqdDists[r][n-n_min] != DistanceMetric::UPPER_BOUND);//terminate early if
                                                  //this solution has not been
                                                  //initialized.  It means
                                                  //all solutions with higher r
                                                  //will not be either.
      ++r)
  {
    s.sum_sqd_dist = aaSumSqdDists[r][n-n_min];
    s.alignment    = aaAlignments[r][n-n_min];
    CopyMatrix3x4(aaTransforms[r][n-n_min], s.transform);
    vSolutions.push_back(s);
  }
  //sort this vector by increasing rmsd
  sort(vSolutions.begin(), vSolutions.end());
  //copy them back into the solutions array.
  for(int r=0; r < vSolutions.size(); ++r)
  {
    aaSumSqdDists[r][n-n_min]      = vSolutions[r].sum_sqd_dist;
    aaAlignments[r][n-n_min] = vSolutions[r].alignment;
    CopyMatrix3x4(vSolutions[r].transform, aaTransforms[r][n-n_min]);
  }


  //(Again This is interesting but unnecessary.I may want to take the next lines out)
  if (alignment_temp != aaAlignments[0][n-n_min])
    cout <<
      " ...Interesting, after iterative refinement, the RMSD of formerly the\n"
      "    best alignment _exceeded_ the RMSD of one of the alternate alignments\n"
      "    with " << n << " matches.  These alignments will be reordered.\n"
         << flush;
} //SearchOrDyn3d::SortAltList(int n)






void SearchOrDyn3d::
FirstPass(OrientationGenerator const& og,
          long start,
          long stop)
{
  assert((start >= 0) && (start < og.size()));
  assert((stop > 0) && (stop <= og.size()));
  assert(stop >= start);

  //initialize some variables:

  //For the purpose of predicting running time:
  struct timezone dummy_tz;//we only need this to pacify the syntax of gettimeofday
  struct timeval start_time, current_time, prev_time;
  gettimeofday(&start_time, &dummy_tz);
  prev_time = start_time;

  //Generate the class used for generating superpositions.
  //(This is only necessary if intermediate refinement was selected.)
  Superimpose *pSuperimpose = NULL;
  if (settings.refine_while_searching)
    pSuperimpose = new Superimpose(max_num_matches *
                                   Biopolymer::Residue::NumBackboneAtoms());

  //Generate the class used for filtering out bad orientations
  OrientationFilter *pFilter = NULL;
  if (settings.orientation_filter_nw)
    pFilter= new OrientationFilter(settings.orientation_filter_nw_max_dist,
                                   settings.orientation_filter_nw_min_num_matches,
                                   numRes1,
                                   numRes2);

  //Now, loop over all orientations:
  for(long o_id = start;
      o_id != stop;
      ++o_id)
  {
    DEBUG_MSG(DBG_SEARCH_OR_DYN3D, "Generating alignments for orientation ID#"
              << o_id);

    Matrix3x4 const& orientation_orig   //The original orientation
      = og[o_id].transform;             //retrieved from the
                                        //OrientationGenerator og.
                                        //The alignments will be
                                        //generated at this
                                        //orientation.

    //Rotate the second molecule to the desired orientation.
    ApplyTransformBackbone(orientation_orig,
                           settings.aMol_f[1],
                           aMol_c[1]);

    //If the orientation passes the filter, then use it
    if ((! settings.orientation_filter_nw)
        ||
        pFilter->GoodOrientation(aMol_c[0], aMol_c[1]))
    {

      //The relative positions are different so we have to recalculate the
      //squared-distances between the residues in the two sequences.
      dyn3d.FillLookupTable(aMol_c[0], aMol_c[1]);

      // Initialize the first layer of the dynamic-programming-table
      // (required before calling Solve_IncrN())
      dyn3d.InitTable_IncrN();
    

      // Now fill in the rest of the table, updating the best alignments so
      // far for each different value of n.
      int count_num_updates=0;
      int n;
      for(n = 2; n <= n_max; ++n)
      {

        // calling dyn3d.Solve_IncrN(n) creates an alignment with n matches
        dyn3d.Solve_IncrN(n);

        if (MaxDistanceExceeded(n))
          break; //Terminate early if the distance constraint was violated
                 //...that is, if any residues matched in the alignment exceeded
                 //the maximum distance as set by DistanceMetric::SetMaxDistance().

        //We only consider storing the alignments that
        //contain at least "n_min" residues.  Otherwise,
        //we discard the solution.
        if (n >= n_min)
        {
          // ***    We keep track of three pieces of information for
          // *** every solution we calculate during the search:
          // *** 1) the candidate alignment
          // *** 2) the orientation between the two structures,
          // *** 3) the RMSD of the alignment (at that orientation)
          // *** The third piece of information #3) is used to judge
          // *** the quality of this alignment.  Alignments with high RMSD
          // *** are discarded in favor of those with low RMSD.  We want to
          // *** use an accurate method of computing RMSD from this alignment.
          // ***   But which RMSD do we use?
          // *** -The RMSD between matched residues at the
          // ***  orientation in which the alignment was calculated, or
          // *** -the RMSD between matched residues at an orientation
          // ***  that has been optimized to minimize the RMSD of the
          // ***  alignment we just calculated.
          // *** The second RMSD is likely to be closer to the RMSD
          // *** you would get if we were to pick this solution and use
          // *** it as a starting point to begin iterative refinement.
          // *** (by calling the Refine() function).
          // *** Thus it is a better indicator of the quality of what the final
          // *** solution will be at the end.
          // ***    On the down side, it takes a little extra computation
          // *** time to compute the second kind of RMSD.
          // *** (It is hardly noticeable for large proteins, though)
          // *** The decision as to which kind of RMSD to use is optional.
          // *** The variable "settings.refine_while_searching"
          // *** determines this.  If it is true, then we use the second
          // *** the RMSD after optimal superposition.  If it is false,
          // *** we use the first kind of RMSD, the RMSD at the original
          // *** orientation generated by matching fragments together.
          // ***    No matter which kind of RMSD we use, the orientation
          // *** we associate with this solution should be consistent
          // *** with the RMSD.  In the code that follows I save
          // *** the two versions of the rmsd and the orientation
          // *** associated with this alignment in different variables:
          // *** The first set of variables is:
          Real sum_sqd_dist_orig;  // The sum_sqd_dist of the
                                   // alignments (with the structures placed
                                   // in the orientation produced by the
                                   // OrientationGenerator, og).

          // *** And the second set of variables is:
          // *** (Note: These variables are not calculated unless
          // ***        settings.refine_while_searching == true.)
          Matrix3x4 orientation_reduced;//the optimal orientation of minimal RMSD for
                                        //that alignment.

          Real sum_sqd_dist_reduced;   //The sum-squared-
                                       //distance for that alignment
                                       //with the two structures in this
                                       //new orientation.  (This reflects
                                       //the "quality" of the solution
                                       //more accurately than
                                       //sum_sqd_dist_orig.)

          ExtractDynTableSolution(n,
                                  alignment_temp,
                                  sum_sqd_dist_orig,
                                  sum_sqd_dist_reduced,
                                  orientation_reduced,
                                  pSuperimpose);

          // *** Finally, for ease of notation,
          // *** we denote the sum_sqd_dist and orientation.
          // *** to refer to whichever quantities will be actually
          // *** be associated with this solution.
          Real sum_sqd_dist;
          Matrix3x4 orientation;

          if (settings.refine_while_searching)
          {
            sum_sqd_dist = sum_sqd_dist_reduced;
            CopyMatrix3x4(orientation_reduced, orientation);
          }
          else
          {
            sum_sqd_dist = sum_sqd_dist_orig;
            CopyMatrix3x4(orientation_orig, orientation);
          }

          if (MaxRMSDExceeded(sum_sqd_dist, n))
            break; //Terminate early if an RMSD constraint was violated,
                   //...that is, if sqrt(sum_sqd_dist/n) > settings.dyn3d_max_rmsd

          // *** Now we have a solution
          DEBUG_MSG(DBG_SEARCH_OR_DYN3D,
                    " alignment found at orientation#"<< o_id <<", " << 0+1 <<
                    "th-best, n=" << n <<
                    ", sum_sqd_dist=" << sum_sqd_dist <<
                    ", rmsd="
                    << sqrt(sum_sqd_dist
                            /
                            (n * Biopolymer::Residue::NumBackboneAtoms())));

          // *** If it's the best solution, update the best-so-far table
          if (sum_sqd_dist < aaSumSqdDists[0][n-n_min])
          {
            DEBUG_MSG(DBG_SEARCH_OR_DYN3D,
                      "Updating " << 0+1 << "th-best solution with "<< n << " matches.");
            RecordSolution(0,
                           n,
                           sum_sqd_dist,
                           orientation,
                           alignment_temp);
            ++count_num_updates;
          } //if (sum_sqd_dist < aaSumSqdDists[0][n-n_min])

          // *** Whether it is the best solution or not, we will need
          // *** to hold on to this solution until later,
          // *** when searching for alternate solutions.
          // *** At that time we will need to compare it
          // *** against potential alternate solutions.
          if (settings.num_solutions > 1)
          {
            KeepTrackOfAltSolutions(n,
                                    o_id,
                                    sum_sqd_dist,
                                    orientation_reduced,
                                    alignment_temp);
          }
        } // if (n >= n_min)
      } // for(n = n_min; (n <= n_max); ++n)


#ifdef CREATE_PROFILE_HISTOGRAM
      // **** Optional stuff: ****
      // Why not record the number of layers of table calculated?
      // Explanation:
      //   Because maximum-distancer and maximum-rmsd cutoffs, terminate the
      //   calculation early, I thought it would be interesting to find out how
      //   much of the dynamic programming table gets calculated, on average,
      //   and what is the distribution?  Do these optimizations
      //   save us much time?
      int last_layer_calculated = n-1;
      ++(aNumLayersHistogram[ last_layer_calculated-1 ]);
#endif //#ifdef CREATE_PROFILE_HISTOGRAM

#ifdef CREATE_OVERALL_PAIR_HISTOGRAM
      // **** More optional stuff: ****
      //Now, accumulate a histogram of which residues got matched
      //from the alignments for this orientation.
      //(Actually, I may take this out later.  These pair-histograms
      // are definitely not turning out to be very meaningful at least so far.)

      dyn3d.AccumPairHistogram(n_min,
                               n-1,  //<-Stop accumulating at the last
                               //  iteration we sucessfully computed
                               //  (which occured at n-1, not n)
                               aMol_c[0],
                               aMol_c[1]);
#endif //#ifdef CREATE_OVERALL_PAIR_HISTOGRAM

      //If this orientation generated alignments that were better
      //than the best alignments so far, this is a good time
      //to print out some stats/information to the user:
      if (count_num_updates > 0)
      {
        gettimeofday(&current_time, &dummy_tz);
        cout << count_num_updates << " updates for orientation "
             << o_id + 1;
        if (current_time.tv_sec - start_time.tv_sec != 0) 
        {
          float fraction_completed= ((float)(o_id - start + 1)/(stop - start));
          float fraction_remaining= ((float)(stop - o_id) / (stop - start));
          float time_remaining = (current_time.tv_sec - start_time.tv_sec) *
                                 (fraction_remaining / fraction_completed);

          if (settings.num_solutions > 1)
          cout << "   time left(1st iter): ";
            else
          cout << "   time left: ";

          cout << time_remaining / 60 << " min.";
        }
        cout << endl;
      }
    } // if (pFilter->GoodOrientation(aMol_c[0], aMol_c[1]))


    //Now check for ellapsed time.
    //If enough of it has ellapsed, print out a status message to the user.
    gettimeofday(&current_time, &dummy_tz);
    if ((settings.display_stats_interval != Settings::DISSABLE_DISPLAY_STATS_INTERVAL)
        &&
        (current_time.tv_sec - prev_time.tv_sec >= settings.display_stats_interval))
    {
      float fraction_completed= ((float)(o_id - start + 1) / (stop - start));
      float fraction_remaining= ((float)(stop - o_id) / (stop - start));
      float time_remaining = (current_time.tv_sec - start_time.tv_sec) *
                             (fraction_remaining / fraction_completed);
      cout << floor(fraction_completed*100 + 0.5) << "% complete"
        "   time left: " << time_remaining / 60 << " min."
           << endl;
      prev_time = current_time;
    }
  } // loop over all orientations "for (o_id = start; o_id != end; ++o_id)"


  //If no solutions were found, print out an error message. (aren't I nice?)
  {
    int n;
    for (n = n_min;
         ((n <= n_max) &&   //Condition below is true only if no solution was
          (aaSumSqdDists[0][n-n_min] == DistanceMetric::UPPER_BOUND)); //found.
         ++n)
      {}
    if (n > n_max)
      ERR("Error: No solutions were found satisfying your criteria.\n"
          "       The restrictions you placed on accepting solutions were\n"
          "       too strict.  Try using a lower minimum-number-of-matches,\n"
          "       a larger maximum distance cutoff, for example.");
  }

  //clean up
  if (pSuperimpose) delete pSuperimpose;
  if (pFilter)      delete pFilter;
} //SearchOrDyn3d::FirstPass()





void
SearchOrDyn3d::
MultiPass(OrientationGenerator const& og,
          set<long> const& use_these_orientations,
          bool (*pDissimilar)(PairwiseAlignment const &,
                              PairwiseAlignment const &))
{
  for (int r = 1; r < settings.num_solutions; ++r)
  {
    cout << "----- Searching for the " <<
      r+1 << "th-best alignments -----" << endl;
    MultiPass(r, og, use_these_orientations, pDissimilar);
  }
}



void
SearchOrDyn3d::
MultiPass(int r,
          OrientationGenerator const& og,
          set<long> const& use_these_orientations,
          bool (*pDissimilar)(PairwiseAlignment const &,
                              PairwiseAlignment const &))
{
  assert(pDissimilar); //Make sure a similarity function has been specified.

  //initialize some variables:

  //For the purpose of predicting running time:
  struct timezone dummy_tz;//we only need this to pacify the syntax of gettimeofday
  struct timeval start_time, current_time, prev_time;
  gettimeofday(&start_time, &dummy_tz);
  prev_time = start_time;

  //Generate the class used for generating superpositions.
  //(This is only necessary if intermediate refinement was selected.)
  Superimpose *pSuperimpose = NULL;
  if (settings.refine_while_searching)
    pSuperimpose = new Superimpose(max_num_matches *
                                   Biopolymer::Residue::NumBackboneAtoms());

  //Generate the class used for filtering out bad orientations
  OrientationFilter *pFilter = NULL;
  if (settings.orientation_filter_nw)
    pFilter= new OrientationFilter(settings.orientation_filter_nw_max_dist,
                                   settings.orientation_filter_nw_min_num_matches,
                                   numRes1,
                                   numRes2);

  //Now, loop over all orientations:
  //Loop through the set of orientations
  long orientation_counter = 1;
  for(set<long>::const_iterator pO_id = use_these_orientations.begin();
      pO_id != use_these_orientations.end();
      ++pO_id)
  {
    long o_id = *pO_id;
    
    DEBUG_MSG(DBG_SEARCH_OR_DYN3D, "Generating alignments for orientation ID#"
              << o_id);

    Matrix3x4 const& orientation_orig   //The original orientation
      = og[o_id].transform;             //retrieved from the
                                        //OrientationGenerator og.
                                        //The alignments will be
                                        //generated at this
                                        //orientation.

    //Rotate the second molecule to the desired orientation.
    ApplyTransformBackbone(orientation_orig,
                           settings.aMol_f[1],
                           aMol_c[1]);

    //If the orientation passes the filter, then use it
    if ((! settings.orientation_filter_nw)
        ||
        pFilter->GoodOrientation(aMol_c[0], aMol_c[1]))
    {

      //The relative positions are different so we have to recalculate the
      //squared-distances between the residues in the two sequences.
      dyn3d.FillLookupTable(aMol_c[0], aMol_c[1]);

      // Initialize the first layer of the dynamic-programming-table
      // (required before calling Solve_IncrN())
      dyn3d.InitTable_IncrN();
    

      // Now fill in the rest of the table, updating the best alignments so
      // far for each different value of n.
      int count_num_updates=0;
      int n;
      for(n = 2; n <= n_max; ++n)
      {

        // calling dyn3d.Solve_IncrN(n) creates an alignment with n matches
        dyn3d.Solve_IncrN(n);

        if (MaxDistanceExceeded(n))
          break; //Terminate early if the distance constraint was violated
                 //...that is, if any residues matched in the alignment exceeded
                 //the maximum distance as set by DistanceMetric::SetMaxDistance().
        if (n >= n_min)
        {

          Real sum_sqd_dist_orig;  //The sum_sqd_dist of the
                                   // alignments (with the structures placed
                                   // in the orientation produced by the
                                   // OrientationGenerator, og).

          // *** And the second set of variables is:
          // *** (Note: These variables are not calculated unless
          // ***        settings.refine_while_searching == true.)
          Matrix3x4 orientation_reduced;//the optimal orientation of minimal RMSD for
                                        //that alignment.

          Real sum_sqd_dist_reduced;               //the new sum-squared-
                                                   //distance for that alignment
                                                   //with the two structures in this
                                                   //new orientation.  (This reflects
                                                   //the "quality" of the solution
                                                   //more accurately than 
                                                   //sum_sqd_dist_orig.)

          ExtractDynTableSolution(n,
                                  alignment_temp,
                                  sum_sqd_dist_orig,
                                  sum_sqd_dist_reduced,
                                  orientation_reduced,
                                  pSuperimpose);

          // *** Finally, for ease of notation,
          // *** we denote the sum_sqd_dist, and orientation.
          // *** to refer to whichever quantities will be actually
          // *** be associated with this solution.
          Real sum_sqd_dist;
          Matrix3x4 orientation;
          if (settings.refine_while_searching)
          {
            sum_sqd_dist = sum_sqd_dist_reduced;
            CopyMatrix3x4(orientation_reduced, orientation);
          }
          else
          {
            sum_sqd_dist = sum_sqd_dist_orig;
            CopyMatrix3x4(orientation_orig, orientation);
          }



          if (MaxRMSDExceeded(sum_sqd_dist, n))
            break; //Terminate early if an RMSD constraint was violated,
                   //...that is, if sqrt(sum_sqd_dist/n) > settings.dyn3d_max_rmsd

          if (sum_sqd_dist < aaSumSqdDists[r][n-n_min])
          {
            DEBUG_MSG(DBG_ALT_MULTI_PASS,
                      "Candidate alignment is sufficiently different than:");
            int r_alt;
            for (r_alt = 0;
                 (r_alt < r)
                   &&
                   (*pDissimilar)(alignment_temp, aaAlignments[r_alt][n-n_min]);
                 ++r_alt)
            {
              #ifdef DEBUG
              #if DBG_ALT_MULTI_PASS
              cerr << "(r=" << r_alt+1 << ",n=" << n << ") ";
              #endif //#if DBG_ALT_MULTI_PASS
              #endif //#ifdef DEBUG
            }
            #ifdef DEBUG
            #if DBG_ALT_MULTI_PASS
            cerr << endl;
            #endif //#if DBG_ALT_MULTI_PASS
            #endif //#ifdef DEBUG

            if (r_alt == r) //i.e. if all the previous alignments, from
                            //r_alt = {0 .. r-1} are sufficiently different
                            //to this one, then replace the rth-best alignment
                            //with this one.
            {
              RecordSolution(r,
                             n,
                             sum_sqd_dist,
                             orientation,
                             alignment_temp);
              DEBUG_MSG(DBG_ALT_MULTI_PASS,
                        "Updating " << r+1 << "th-best solution with "<< n << " matches.");
              ++count_num_updates;
            }
          } //if (sum_sqd_dist < aaSumSqdDists[r][n-n_min])
        } // if (n >= n_min)
      } // for(n = n_min; (n <= n_max); ++n)


#ifdef CREATE_PROFILE_HISTOGRAM
      // **** Optional stuff: ****
      // Why not record the number of layers of table calculated?
      // Explanation:
      //   Because maximum-distancer and maximum-rmsd cutoffs, terminate the
      //   calculation early, I thought it would be interesting to find out how
      //   much of the dynamic programming table gets calculated, on average,
      //   and what is the distribution?  Do these optimizations
      //   save us much time?
      int last_layer_calculated = n-1;
      ++(aNumLayersHistogram[ last_layer_calculated-1 ]);

      // **** More optional stuff: ****
      //Now, accumulate a histogram of which residues got matched
      //from the alignments for this orientation.
      //(Actually, I may take this out later.  These pair-histograms
      // are definitely not turning out to be very meaningful at least so far.)
      dyn3d.AccumPairHistogram(n_min,
                               n-1,  //<-Stop accumulating at the last
                               //  iteration we sucessfully computed
                               //  (which occured at n-1, not n)
                               aMol_c[0],
                               aMol_c[1]);
#endif //#ifdef CREATE_PROFILE_HISTOGRAM

      //If this orientation generated alignments that were better
      //than the best alignments so far, this is a good time
      //to print out some stats/information to the user:
      if (count_num_updates > 0)
      {
        gettimeofday(&current_time, &dummy_tz);
        cout << count_num_updates << " updates for orientation "
             << o_id + 1;
        if (current_time.tv_sec - start_time.tv_sec != 0) 
        {
          float fraction_completed= ((float)(orientation_counter + 1) /
                                             use_these_orientations.size());
          float fraction_remaining= ((float)(use_these_orientations.size() -
                                             orientation_counter)
                                              /
                                              use_these_orientations.size());
          float time_remaining = (current_time.tv_sec - start_time.tv_sec) *
                                 (fraction_remaining / fraction_completed);
          cout << "   time left(" << r+1 << "th iter): "
               << time_remaining / 60 << " min."
               << endl;
        }
        else
          cout << endl;
      }
    } // if (pFilter->GoodOrientation(aMol_c[0], aMol_c[1]))


    //Now check for ellapsed time.
    //If enough of it has ellapsed, print out a status message to the user.
    gettimeofday(&current_time, &dummy_tz);
    if ((settings.display_stats_interval != Settings::DISSABLE_DISPLAY_STATS_INTERVAL)
        &&
        (current_time.tv_sec - prev_time.tv_sec >= settings.display_stats_interval))
    {
      float fraction_completed= ((float)(orientation_counter + 1) /
                                         use_these_orientations.size());
      float fraction_remaining= ((float)(use_these_orientations.size() -
                                         orientation_counter)
                                         /
                                         use_these_orientations.size());
      float time_remaining = (current_time.tv_sec - start_time.tv_sec) *
                             (fraction_remaining / fraction_completed);
      cout << floor(fraction_completed*100 + 0.5) << "% complete"
        "   time left(" << r+1 <<"th iter): " << time_remaining / 60 << " min."
           << endl;
      prev_time = current_time;
    }

    orientation_counter++;
  } //loop orientations "for(pO_id != != use_these_orientations.end(); ++pO_id)"


  //If no solutions were found, print out an warning message. (aren't I nice?)
  {
    int n;
    for (n = n_min;
         ((n <= n_max) &&   //Condition below is true only if no solution was
          (aaSumSqdDists[r][n-n_min] == DistanceMetric::UPPER_BOUND)); //found.
         ++n)
      {}
    if (n > n_max)
      cout <<
        "-------------------------------------------------------------\n"
        " WARNING:  No "
         << r+1 << "th-best solutions were found\n"
        "           satisfying your criteria.  The restrictions you\n"
        "           placed on selecting alternate solutions were too\n"
        "           strict to find this many alternate solutions.\n"
        "-------------------------------------------------------------\n";
  }


  //clean up
  if (pSuperimpose) delete pSuperimpose;
  if (pFilter)      delete pFilter;
} //SearchOrDyn3d::CalcAltSolutionsMultiPass(r)










void
SearchOrDyn3d::Refine()
{
  // At this point, we have a set of optimal alignments based on a guess
  // of the angle orientation which was obtained by minimizing the RMS
  // distance between num_consec_matches arbitrarily chosen consecutive
  // residues from the two sequences...
  //    Now, if the "settings.refine_max_iters" option is set, we can use these
  // positions of the residues that got matched
  // to make better guesses as to how to orient the two sequences.
  // If for a given number-of-matches, n, a best alignment is found, 
  // generate a new rotation/translation that optimizes the RMS for
  // that alignment.  Now, find the best alignments for this new
  // orientation.  In so doing, perhaps a new alignment is discovered
  // for containing n-matches, which is the best-so-far (for that n).
  //    Repeat the process iteravely, until either there are no new
  // updates, or too-many iterations have occured and we want to bail-out.
  //    The algorithm I described above would potentially seem to spawn
  // many children with each iteration so you have to be a little
  // careful to avoid having the recursion grow exponentially.
  // (see below)
  //              Details:
  //    We hold the orientation fixed generate a list of alignments.
  // We keep track of the ones that are better (for any value of n)
  // than the best alignments so far (for that n).
  // 1) Starting with a new alignment with with a number of matched
  //   pairs of residues, n, somewhere in the middle. (see "Ugly Details")
  //   new orientation to optimize that alignment, and calculate
  //   and even newer set of best alignments using this new orientation.
  // 2) If any of these are better than the best alignments so-far (for
  //   any value of n), update the list of best-alignments-so-far
  //   (possibly overwriting alignments that were generated in the
  //    last recursion and have not been iterated over yet) and return to
  //    step 1)
  //   Otherwise go to step 3)
  // 3) Exit iteration.
  //
  //    * Ugly Details:
  //       Instead of starting with the alignment with the highest number
  //   of matched-pairs,n, or lowest number of matched pairs of residues,
  //   I chose an alignment with n as the median over the set of all the
  //   alignments, that have unoptimized alignments
  //   (for which aNeed2Refine[n-n_min] == true).  I do this because orientations
  //   that are optimal for matching very few residues (low n) tend to be
  //   very poor orientations for matching the sequences as a whole.
  //   Similarly, this is also true about orientations that try
  //   to match almost all the residues (high n).  Choosing a value of
  //   n that is in the middle, produces orientations that
  //   tend not to be as spurious, and more often are similar
  //   to the orientations used for matching many of the other alignments
  //       So, hopefully using the middle values of n first, will
  //   yeild the overall best orientation guesses as earlier as possible,
  //   which I hope will speed up the convergence of this method.

  cout << "Entering refinement stage." << endl;

  assert(settings.refine_method != Settings::NO_REFINE);
  int num_solutions_refined;
  switch (settings.refine_method)
  {
  case Settings::REFINE_BEST:
    num_solutions_refined = 1; //ie, only refine the first best,
                            //do not refine alternates (ie 2nd-best, 3rd-best)
    break;
  case Settings::REFINE_ALL:
    num_solutions_refined = settings.num_solutions; //ie, refine all the solutions,
                                                    //the best, 2nd best, etc...
    break;
  default:
    assert(0);
    break;
  }
  
  int num_new_solutions = 0;
  int num_iterations;

  Superimpose superimpose(max_num_matches
                      * Biopolymer::Residue::NumBackboneAtoms());
  bool  *aNeed2Refine = new bool [max_num_matches];
  CHECK_ALLOC(aNeed2Refine);

  PairwiseAlignment oldNthAlignment(max_num_matches);

  for(int r = 0; r < num_solutions_refined; ++r)
  {
    int n; //n is an index for looping over number-of-matches-made

    for(n = n_min; n <= n_max; ++n)
    { //At this point, ALL valid alignments are unrefined and
      //need to be refined, so set aNeed2Refine[n-n_min] to true.
      //An alignment is considered "valid" if it does not violate any
      //maximum distance limits set by the user.
      //This can be tested by checking to see if "Sum-Squared-Distance" is
      //less than DistanceMetric::UPPER_BOUND ( <--- a value which only
      //exists to signal that a distance constraint has been violated).
      aNeed2Refine[n-n_min] =
        (aaSumSqdDists[r][n-n_min] < DistanceMetric::UPPER_BOUND);
    }

    if (settings.num_solutions != 1)
      cout << "Currently refining the " << r+1 << "th-best alignments.\n"
           << flush;

    for(num_iterations = 1;
        (num_iterations <= settings.refine_max_iters)
          || (settings.refine_max_iters == Settings::ITERATE_UNTIL_CONVERGENCE);
        ++num_iterations)
    {
      //--- Which of our alignments do we refine next? ----
      // The one whose number of matches, n, is the median of the numbers
      // of matches in all alignments that still need to be refined.
      // (See "Ugly Details" discussion above.)
      n = FindMedian(aNeed2Refine, n_min, n_max); //See <simple_numeric_utils.h>

      // Figure out if it's time to exit the loop.
      // -Exit if there are no more unrefined alignments (num_in_set==0).
      if (n == (n_min-1))
      {
        DEBUG_MSG(DBG_REFINE_OR_DYN3D,
                  "No more unrefined "
                  << r+1 << "th-best alignments left.\n"
                  "                    Exiting Refinement for "
                  << r+1 << "th-best...");
        break; //exit the loop.  we're done.
      }

      assert((aaSumSqdDists[r][n-n_min] < DistanceMetric::UPPER_BOUND)
             &&
             (aNeed2Refine[n-n_min] == true));

      //At this point "n" should store the # matches in the alignment we want to refine!
      DEBUG_MSG(DBG_REFINE_OR_DYN3D,
                "Refine() going into "
                << num_iterations << "th orientation, optimized for "
                << r+1 << "th-best, n=" << n);

      //Otherwise, find the Minimum RMSD of this alignment, as well as
      //the orientation between the molecules (this will be stored in
      //"orientation_of_refinement_for_n") that minimizes it.

      Matrix3x4 orientation_of_refinement_for_n;
      aaAlignments[r][n-n_min].CalcMinRMSD(settings.aMol_f[0],
                                           settings.aMol_f[1],
                                           orientation_of_refinement_for_n,
                                           &superimpose,
                                           aRotateUsingTheseCoords1,
                                           aRotateUsingTheseCoords2);

      //Okay, now that we have the optimal transformation for this
      //alignment, apply it to the residues in the second sequence.
      ApplyTransformBackbone(orientation_of_refinement_for_n,
                             settings.aMol_f[1],
                             aMol_c[1]);

      //The molecule has been rotated, so the distances between atoms
      //have changed.  We have to recalculate the inter-distances
      dyn3d.FillLookupTable(aMol_c[0], aMol_c[1]);
          
      //Now, actually calculate the best alignments for this orientation
      // for all numbers-of-matched-pairs, keeping track of any alignments
      // which are record-setters (to be iterated over later).

      //backup the current best alignment for n matches
      oldNthAlignment.resize(n);
      oldNthAlignment = aaAlignments[r][n-n_min];

      // Initialize the first layer of the dynamic-programming-table
      //  (required before calling Solve_IncrN())
      dyn3d.InitTable_IncrN();

      // Now fill in the rest of the table.
      for(int m = 2;
          ((m <= n_max)
           &&
         //If settings.refine_method==Settings::REFINE_ALL,
         //then we only want to improve the alignment
         //with n matches.  We are not interested in calculating
         //any other alignments in particular, those with m>n.
         //So, after we get to the alignment with m=n, stop.
           (! ((r > 0) &&
               (settings.refine_method == Settings::REFINE_ALL) &&
               (m > n))));
          ++m)
      {
        // calling dyn3d.Solve_IncrN(n) creates an alignment with n matches
        dyn3d.Solve_IncrN(m);

        if (MaxDistanceExceeded(m))
          break; //Terminate early if the distance constraint was violated
                 //...that is, if any residues matched in the alignment exceeded
                 //the maximum distance as set by DistanceMetric::SetMaxDistance().
        //We only consider storing the alignments that
        //contain at least "n_min" residues.  Otherwise,
        //we discard the solution.
        if (m >= n_min)
        {

          Real sum_sqd_dist_orig;  //The sum_sqd_dist of the
                                   // alignments (with the structures placed
                                   // in the orientation produced by the
                                   // OrientationGenerator, og).

          // *** (Note: The next three variables are not calculated unless
          // ***        settings.refine_while_searching == true.)
          Matrix3x4 orientation_reduced;//the optimal orientation of minimal
                                        //RMSD for that alignment.

          Real sum_sqd_dist_reduced;         //the new sum-squared-
                                             //distance for that alignment
                                             //with the two structures in this
                                             //new orientation.  (This reflects
                                             //the "quality" of the solution
                                             //more accurately than
                                             //sum_sqd_dist_orig.)


          ExtractDynTableSolution(m,
                                  alignment_temp,
                                  sum_sqd_dist_orig,
                                  sum_sqd_dist_reduced,
                                  orientation_reduced,
                                  &superimpose);

          // *** Finally, for ease of notation,
          // *** we denote the sum_sqd_dist and orientation.
          // *** to refer to whichever quantities will be actually
          // *** be associated with this solution.
          Real sum_sqd_dist;
          Matrix3x4 orientation;
          if (settings.refine_while_searching)
          {
            sum_sqd_dist = sum_sqd_dist_reduced;
            CopyMatrix3x4(orientation_reduced, orientation);
          }
          else
          {
            sum_sqd_dist = sum_sqd_dist_orig;
            CopyMatrix3x4(orientation_of_refinement_for_n, orientation);
          }


          if (MaxRMSDExceeded(sum_sqd_dist, m))
            break; //Terminate early if an RMSD constraint was violated,
                   //...that is, if sqrt(sum_sqd_dist/n) > settings.dyn3d_max_rmsd

          DEBUG_MSG(DBG_REFINE_OR_DYN3D,
                    " alignment generated: refining for " << r+1 <<
                    "th-best alignments; m=" << m <<
                    ", sum_sqd_dist=" << sum_sqd_dist <<
                    ", rmsd="
                    << sqrt(sum_sqd_dist
                            /
                           (m * Biopolymer::Residue::NumBackboneAtoms())));

          if (sum_sqd_dist < aaSumSqdDists[r][m-n_min])
          {
            //     In the (r == 0) case we were looking for solutions that
            //are global minima in RMSD, so anytime an alignment with any
            //number of matches (say "m"), has a lower RMSD than the
            //best-one-so far (also with m matches), discard the old solution
            //in favor of the new.  I call this "cross-n-referencing":
            //updating alignments which are not the solutions we originally
            //set out to optimize during the course of this iteration.
            //(Recall that at this point, the orientation between the two
            // molecules has been optimized for the best alignment so-far
            // with with "n" matches, which is a different solution to the
            // best one containing "m" matches.)
            //     I make this distinction because, when refining 2nd and
            //3rd-best alignments (ie. when r != 0), we cannot do this.
            //In that case, we have to be more-careful.  (See below.)

            if (r == 0)
            {
              DEBUG_MSG(DBG_REFINE_OR_DYN3D,
                        "Updating " << r+1 <<
                        "th-best solution with " << m << " matches.");

              //Update the best-so-far-list
              RecordSolution(r,
                             m,
                             sum_sqd_dist,
                             orientation,
                             alignment_temp);
              if (m != n) aNeed2Refine[m-n_min] = true;
              ++num_new_solutions;

              //     The current orientation was chosen to minimize the RMSD
              //between the pairs of residues matched in the best alignment
              //so-far with n matches.
              //     There are two ways you can reach convergence.
              //Either the alignment doesn't change:
              //(oldNthAlignment == alignment_temp)
              //or it changes, but it doesn't get any better:
              //(sum_sqd_dist >= aaSumSqdDists[r][m-n_min]).
              //Either condition means we have reached convergence.
              //This is confusing, and I check for both because if I don't,
              //in some cases floating point roundoff error can cause it
              //to go into an infinite loop (by oscillating back and forth 
            } //if (r == 0)
            else
            {
              //If we are looking for alternate solutions (r != 0), instead of
              //the best solutions, then only enable refinement if
              //"settings.refine_method" is set to "Settings::REFINE_ALL".
              assert(settings.refine_method == Settings::REFINE_ALL);

              //   If (r != 0), we are interested in refining
              //the 2nd-best or 3rd-best solutions, (or rth-best).
              //We are not seeking the overalllowest-rmsd solution.  Instead,
              //we are looking for local, not global minima in RMSD.
              //    Recall that in the (r == 0) case we were looking for
              //solutions that were global minima in RMSD, so anytime an
              //alignment with any number of matches (say "m"), had a lower
              //RMSD than the best-one-so far with m matches, we would
              //discard the old solution in favor of the new.  We don't
              //want to do that anymore.
              //    Also recall that the orientation we are using was
              //chosen minimize the RMSD between paired residues in the
              //alignment with "n" matches.  However, the alignment we just
              //calculated at this orientation contains "m" matches.  These
              //two alignments are not directly comparable unless m==n.
              //Note that the rth-best alignment with "n" matches
              //does not necessarilly fall within the same
              //local-minima in RMSD that the rth-best alignment
              //with "m" matches belongs to.  And it is _very_ likely
              //that, for some value of "n", at least one of the
              //rth-best-solutions with n matches orients the molecules
              //in a similar way to the way they were oriented in the
              //_best_ alignments, (at some different number of matches).
              //Overwriting the rth-best alignments with "m" matches with
              //the alignment made at this orientation with m matches
              //in this case, would replace the rth-best solutions with _the_
              //best solutions.  This is not what we want.
              //When seeking alternate alignments, (ie. r != 0),
              //we are looking for local, not global minima.
              //     So in this case (when r != 0), we only update if m == n.
              //This will isolate the refinement process separately
              //for each solution.
              if (m == n)
              {
                DEBUG_MSG(DBG_REFINE_OR_DYN3D,
                      "Updating " << r+1 <<
                      "th-best solution with " << m << " matches.");

                //Update the best-so-far-list
                RecordSolution(r,
                       m,
                       sum_sqd_dist,
                       orientation,
                       alignment_temp);

                ++num_new_solutions;
              } //if (m == n)
            } //else clause for "if (r == 0)"
          } //if (sum_sqd_dist < aaSumSqdDists[r][n-n_min])
        } //if (m >= n_min)
      } //for(m=2; (m <= n_max); ++m)

      if (oldNthAlignment == aaAlignments[r][n-n_min])
      { 
        DEBUG_MSG(DBG_REFINE_OR_DYN3D,
              " ..Reached convergence for " << r+1
              << "th-best, n = " << n);
        aNeed2Refine[n-n_min] = false;
      }
      else
      {
        aNeed2Refine[n-n_min] = true;
      } //else clause for "if (oldNthAlignment == alignment)"

    } // for (num_iters=1; num_iters<=refine_max_iters; num_iters++)

    #ifdef DEBUG
    //make sure every alignment has been refined
    for(n=n_min;
        (n <= max_num_matches) && (aaSumSqdDists[r][n-n_min] < DistanceMetric::UPPER_BOUND);
        ++n)
      assert(! aNeed2Refine[n-n_min]);
    #endif
  } // for(int r = 0; r < num_solutions_refined; ++r)

  //print out statistics
  cout << "(During refinement, " << num_iterations
       << " new orientations were tested,\n"
       " and " << num_new_solutions
       << " instances of better alignments were found.)" << endl;

  if (settings.num_solutions > 1)
  {
    //During refinement, the solutions may have been modified.
    //Make sure that after refinement, none of the alternate solutions
    //are "too similar" to eachother.
    CheckAlternatesForConsistency();
  }

  delete [] aNeed2Refine;

} //SearchOrDyn3d::Refine()



void SearchOrDyn3d::
WritePlotFile(string plot_filename_suffix) const
{
  for(int r=0; r < settings.num_solutions; ++r)
    WritePlotFile(r, plot_filename_suffix);
}


void SearchOrDyn3d::
WritePlotFile(int r, string plot_filename_suffix) const
{
  assert(numRes2 == aMol_c[1].size());
  assert((n_min <= max_num_matches) && (n_min >= 1) &&
         (n_max <= max_num_matches) && (n_max >= 1) &&
         (n_min <= n_max));
  assert(aaTransforms[r]);
  assert(aaAlignments[r]);
  assert(aaSumSqdDists[r]);

  stringstream plot_filename;
  if (settings.num_solutions == 1)
    plot_filename << SearchOrDyn3d::BASE_FILE_NAME
                  << "_" << plot_filename_suffix;
  else
    plot_filename << SearchOrDyn3d::BASE_FILE_NAME
                  << r+1 << "_" << plot_filename_suffix;

  ofstream plot_file(plot_filename.str().c_str(), ios::out);
  if (! plot_file)
    ERR_INTERNAL("WritePlotFile() cannot create file \"" << plot_filename.str()
                    << "\" for writing.");

  //Now label the columns:
  plot_file << "# Matches,RMSD,Longest Distance"
    //",Orientation Increment" Commented out on 9/8/1999.we don't print this now
#ifdef REPORT_ALL_RESULTS
            << ",S_str,mu_str,delta_str,Z,"
#endif //#ifdef REPORT_ALL_RESULTS
            << ",-log(probability)\n";


  for (int n = n_min;
       n <= n_max;
       ++n)
  {
    if (aaSumSqdDists[r][n-n_min] < DistanceMetric::UPPER_BOUND)
    {
      //save how many matches were made in the current alignment
      plot_file << n;
      
      //save the RMS_D of the current alignment
      plot_file // << " n = "<< n << ",   min_cost = "
        << " " << sqrt(aaSumSqdDists[r][n-n_min] 
                /
                (n * Biopolymer::Residue::NumBackboneAtoms()));
      
      //Now, print out the match of this alignment which has the largest 
      //distance between it's residues.  First, print out this distance
      //(in the third collumn),
      //The next two collumns are not being written right now. 
      //Please the ignore next 3 lines of comments...
      //                     ...and then print the residue number from
      //sequence A (in the fourth collumn) and the residue number from
      //sequence B (in the fifth collumn) that were involved in this match.
      ApplyTransformBackbone(aaTransforms[r][n-n_min],
                             settings.aMol_f[1],
                             aMol_c[1]);



      //Find the pair of residues whose quick-atoms are furthest apart.
      Biopolymer::const_iterator worst_i, worst_j;

      plot_file <<" "<< aaAlignments[r][n-n_min].FindFurthestResPair(worst_i,
                                                                     worst_j,
                                                                     aMol_c[0],
                                                                     aMol_c[1]);

      //Now print Levitt & Gerstein's Probability statistics
      //based on S_str
      Real S_str, mu_str, delta_str, Z;
      Real probability = LevittGerstein98::Sstr_Probability(aaAlignments[r][n-n_min],
                                                            aMol_c[0],
                                                            aMol_c[1],
                                                            &S_str,
                                                            &mu_str,
                                                            &delta_str,
                                                            &Z);
      plot_file
#ifdef REPORT_ALL_RESULTS
        << " " << S_str 
        << " " << mu_str
        << " " << delta_str
        << " " << Z
#endif //#ifdef REPORT_ALL_RESULTS
        << " " << ((probability == 0.0) ? 0.0 : -log10(probability))
        << "\n";
    } //if (aaSumSqdDists[r][n-n_min] < DistanceMetric::UPPER_BOUND
  } //for (n = n_min; n <=  n_max; ++n)
} // SearchOrDyn3d::WritePlotFile()



void
SearchOrDyn3d::WriteMidasGraphics(int  r,
                                  bool show_residue_markers,
                                  bool connect_the_dots) const
{
  // Now create some alignment files (called BASE_FILE_NAME++<r>+<n>+".gfx", 
  // where <n> is an integer specifying how many matches were made,
  // and <r> specifies the "rank" of the solutions (ie 1st-best, 2nd-best, etc)
  // Each file corresponds to a different number of matches in the alignment
  // (different value of n) and all it stores is a set of red line segments,
  // one for each match in the alignment.
  for(int n=n_min; n<=n_max; ++n)
  {
    ShortString filename;
    ShortString caption;
    if (aaSumSqdDists[r][n-n_min] < DistanceMetric::UPPER_BOUND)
    {
      //Choose a filename and a caption
      if (settings.num_solutions == 1)
      {
        sprintf(filename, "%s%d.gfx", BASE_FILE_NAME, n);
        sprintf(caption,
                "n = %d, RMS_D = %f",
                n,
                DistanceMetric::DivideByNthenSqrt(aaSumSqdDists[r][n-n_min],
                                n*Biopolymer::Residue::NumBackboneAtoms()));
      }
      else
      {
        sprintf(filename, "%s%d_%d.gfx", BASE_FILE_NAME, r+1, n);
        sprintf(caption,
                "rank = %d, n = %d, RMS_D = %f",
                r+1,
                n,
                DistanceMetric::DivideByNthenSqrt(aaSumSqdDists[r][n-n_min],
                                n*Biopolymer::Residue::NumBackboneAtoms())
                );
      }

      //Now that we have the optimal transformation for this
      //alignment, apply it to the residues in the second sequence.
      ApplyTransformBackbone(aaTransforms[r][n-n_min],
                             settings.aMol_f[1],
                             aMol_c[1]);
        
      aaAlignments[r][n-n_min].ExportGFX(aMol_c[0],
                                         aMol_c[1],
                                         filename,
                                         caption,
                                         false,                //show_residue_markers
                                         connect_the_dots,
                                         PairwiseAlignment::MIDAS_MARKER_SIZE);
    } // if (aaSumSqdDists[0][n-n_min] < DistanceMetric::UPPER_BOUND)
  } // for(int n=n_min; n<=n_max; ++n)
} // SearchOrDyn3d::WriteMidasGraphics()




void SearchOrDyn3d::WriteMidasGraphics(bool show_residue_markers,
                                       bool connect_the_dots) const
{
  for(int r = 0; r < settings.num_solutions; ++r)
  {
    WriteMidasGraphics(r,
                       show_residue_markers,
                       connect_the_dots);
  }
} //WriteMidasGraphics()



void
SearchOrDyn3d::WriteMSF(bool show_number_lines,
                        bool compress_using_lower_case) const
{
  for(int r = 0; r < settings.num_solutions; ++r)
  {
    WriteMSF(r, show_number_lines, compress_using_lower_case);
  }
}


void
SearchOrDyn3d::WriteMSF(int  r,
                        bool show_number_lines,
                        bool compress_using_lower_case) const
{
  assert(settings.vPdbFileNames.size() == 2);
  assert(settings.vMsfOutputLabels.size() == 2);
  PairwiseAlignment alignment(max_num_matches);

  bool unrecognized_residue_names =
    (ContainsUnknownResidues(settings.aMol[0]) ||
     ContainsUnknownResidues(settings.aMol[1]));

  for(int n = n_min; n <= n_max; ++n)
  {
    assert(n == aaAlignments[r][n-n_min].NumMatches());
    if (aaSumSqdDists[r][n-n_min] < DistanceMetric::UPPER_BOUND)
    {
      stringstream msf_filename;
      Real rmsd = sqrt(aaSumSqdDists[r][n-n_min]
                       /
                       (n * Biopolymer::Residue::NumBackboneAtoms()));
      if (settings.num_solutions == 1)
        msf_filename << BASE_FILE_NAME << n << ".msf";
      else
        msf_filename << BASE_FILE_NAME << r+1 << "_" << n << ".msf";

      stringstream comments(ios::out);
      if (r == 0)
      {
        comments
          << "Chimera minimal RMSD structural alignment with "
          << n << " equivalences.\n"
            "RMSD = " << rmsd << "\n"
                  "-----\n"
            "Transform Matrix to apply to structure: "
          << settings.vPdbFileNames[1] << "\n"
            << aaTransforms[r][n-n_min][0][0] << " " << aaTransforms[r][n-n_min][0][1] << " " << aaTransforms[r][n-n_min][0][2] << " " << aaTransforms[r][n-n_min][0][3] << "\n"
            << aaTransforms[r][n-n_min][1][0] << " " << aaTransforms[r][n-n_min][1][1] << " " << aaTransforms[r][n-n_min][1][2] << " " << aaTransforms[r][n-n_min][1][3] << "\n"
            << aaTransforms[r][n-n_min][2][0] << " " << aaTransforms[r][n-n_min][2][1] << " " << aaTransforms[r][n-n_min][2][2] << " " << aaTransforms[r][n-n_min][2][3] << "\n";
        }
        else
        {
          comments
            <<
            "Chimera " << r+1<< "th-lowest RMSD structural alignment with "
            << n << " equivalences\n"
            "RMSD = " << rmsd << ".\n"
            "(subject to the constraint that all alternate\n"
            " versions of the alignments must be sufficiently \"dissimilar\"\n"
            " from eachother according to preferences set by the user.)\n"
            "-----\n"
            "Transform Matrix to apply to structure: "
            << settings.vPdbFileNames[1] << "\n"
            << aaTransforms[r][n-n_min][0][0] << " " << aaTransforms[r][n-n_min][0][1] << " " << aaTransforms[r][n-n_min][0][2] << " " << aaTransforms[r][n-n_min][0][3] << "\n"
            << aaTransforms[r][n-n_min][1][0] << " " << aaTransforms[r][n-n_min][1][1] << " " << aaTransforms[r][n-n_min][1][2] << " " << aaTransforms[r][n-n_min][1][3] << "\n"
            << aaTransforms[r][n-n_min][2][0] << " " << aaTransforms[r][n-n_min][2][1] << " " << aaTransforms[r][n-n_min][2][2] << " " << aaTransforms[r][n-n_min][2][3] << "\n";

        }

        alignment.resize(n);
        TranslateAlignmentBetweenMolecules(aaAlignments[r][n-n_min],
                                           alignment,
                                           settings.aMol_f[0],
                                           settings.aMol_f[1],
                                           settings.aMol[0],
                                           settings.aMol[1]);

        alignment.ExportMSF(msf_filename.str(),
                            settings.aMol[0],
                            settings.aMol[1],
                            settings.vMsfOutputLabels[0],
                            settings.vMsfOutputLabels[1],
                            comments.str(),
                            show_number_lines,
                            compress_using_lower_case
                            && (! unrecognized_residue_names));
    } //if (aaSumSqdDists[r][n-n_min] < DistanceMetric::UPPER_BOUND)
  } // for(int n=n_min; n<=n_max; ++n)

  if (unrecognized_residue_names)
  {
    cerr <<"Warning: Residues are present with unrecognized names.\n"
         <<"         These residues will be named 'X'" << endl;
    if (compress_using_lower_case)
      cerr <<"         This will dissable the ability to compress the\n"
           <<"         msf-files using upper/lower-case to\n"
           <<"         distinguished matched/unmatched residues." << endl;
  }
} // SearchOrDyn3d::WriteMSF()



void
SearchOrDyn3d::WriteChimeraInfo() const
{
  for (int r = 0; r < settings.num_solutions; ++r)
  {
    WriteChimeraInfo(r);
  }

  #ifdef CREATE_PAIR_HISTOGRAM
  PrintPairHistogram("pair_hist");
  #endif //#ifdef CREATE_PAIR_HISTOGRAM
  #ifdef CREATE_PROFILE_HISTOGRAM
  PrintNumLayersHistogram("num_layers_hist.txt");
  #endif //#ifdef CREATE_PROFILE_HISTOGRAM
}

void
SearchOrDyn3d::WriteChimeraInfo(int r // the rank (ie, 1st, 2nd, 3rd solutions)
                                ) const
{
  char const *info_filename_suffix = "chimera.info";
  char const *plot_filename_suffix = "chimera.plot";

  WritePlotFile(r, plot_filename_suffix);

  stringstream info_filename;
  //stringstream plot_filename;
  if (settings.num_solutions == 1)
    info_filename << SearchOrDyn3d::BASE_FILE_NAME << "_" << info_filename_suffix;
  else
    info_filename << SearchOrDyn3d::BASE_FILE_NAME << r+1 << "_" << info_filename_suffix;
  ofstream info_file(info_filename.str().c_str(), ios::out);
  if (! info_file)
    ERR_INTERNAL("Cannot open file \""
                 << info_filename.str()
                 << "\" for writing.  Exiting...");
  info_file << "fixed pdb:" << settings.vPdbFileNames[0] << "\n";
  info_file << "moveable pdb:" << settings.vPdbFileNames[1] << "\n";
  if (settings.num_solutions == 1)
  {
    info_file << "rmsd:" << SearchOrDyn3d::BASE_FILE_NAME << "_"
               << plot_filename_suffix << "\n";
    info_file << "base:" << SearchOrDyn3d::BASE_FILE_NAME << "\n";
  }
  else
  {
    info_file << "rmsd:" << SearchOrDyn3d::BASE_FILE_NAME << r+1 << "_"
               << plot_filename_suffix << "\n";
    info_file << "base:" << SearchOrDyn3d::BASE_FILE_NAME << r+1 << "_\n";
  }

  info_file << "number of atoms used:"
             << Biopolymer::Residue::NumBackboneAtoms() << "\n";
  for(int q=0; q < Biopolymer::Residue::NumBackboneAtoms(); ++q)
    info_file << Biopolymer::Residue::LookupBackboneAtomSymbol(q) << "\n";
} //SearchOrDyn3d::PrintChimeraInfo()



void
SearchOrDyn3d::WriteTransform(int r) const
{
  assert(aaTransforms);
  assert(aaTransforms[r]);

  assert((n_min <= max_num_matches) && (n_min >= 1) &&
         (n_max <= max_num_matches) && (n_max >= 1) &&
         (n_min <= n_max));

  stringstream trans_filename;

  for(int n=n_min; n <= n_max; ++n)
  {
    if (aaSumSqdDists[r][n-n_min] < DistanceMetric::UPPER_BOUND)
    {
      assert(aaTransforms[r][n-n_min]);
      if (settings.num_solutions == 1)
        trans_filename << BASE_FILE_NAME << n << ".trans";
      else
        trans_filename << BASE_FILE_NAME << r+1 << "_" << n << ".trans";
      WriteMatrix3x4ToFile(aaTransforms[r][n-n_min],
                           trans_filename.str().c_str());
    }
  }// for(int n=n_min; n<=n_max; ++n)
}//SearchOrDyn3d::WriteTransform(int r)

void
SearchOrDyn3d::WriteTransform() const
{
  for(int r = 0; r < settings.num_solutions; ++r)
  {
    WriteTransform(r);
  }// for(int r = 0; r < settings.num_solutions; ++r)
}// SearchOrDyn3d::WriteTransform()





#ifdef CREATE_PAIR_HISTOGRAM

void
SearchOrDyn3d::PrintPairHistogram(string hist_filename_base) const
{
  if (settings.num_solutions == 1)
  {
    stringstream hist_filename;
    hist_filename << hist_filename_base << ".txt";
    PrintPairHistogram(0, 0, hist_filename.str());
  }
  else
  {
    //Create individual histogram files for each alternate solution
    for (int r = 0; r < settings.num_solutions; ++r)
    {
      stringstream hist_filename;
      hist_filename << hist_filename_base << r+1 << ".txt";
      PrintPairHistogram(r, r, hist_filename.str());
    }

    //Now create a histogram file that accumulates pairs of residues
    //matched in all of these alternate solutions.
    stringstream hist_filename;
    hist_filename << hist_filename_base << ".txt";
    PrintPairHistogram(0,
                       settings.num_solutions-1,
                       hist_filename.str());
  } // if (settings.num_solutions == 1)


  #ifdef CREATE_OVERALL_PAIR_HISTOGRAM
  //If the dynamic-programming table has been allocated then presumably
  //it has been used, and the pair_histogram table has been accumulating.
  assert(dyn3d.GetMaxNumMatches() == max_num_matches);

  //If so, print out the histogram of these matched-pairs to a separate file.
  //  Unlike the one we just printed, this table accumulates matched-pairs
  //over all the alignments we calculated, not just the ones that got
  //accepted.
  //(Actually this is not quite true: FindBestAlignmentsOverAllAngles()
  // accumulates with each iteration, Refine() dosen't.)
  //Anyhow, print out this histogram.

  dyn3d.PrintPairHistogram("pair_hist_overall.txt");
  #endif //#ifdef CREATE_OVERALL_PAIR_HISTOGRAM
} //SearchOrDyn3d::PrintPairHistogram()





void
SearchOrDyn3d::PrintPairHistogram(int start_rank,
                                  int end_rank,
                                  string hist_filename) const
{
  long **aaHistogram;
  aaHistogram = new long* [ numRes1 ];
  CHECK_ALLOC(aaHistogram);
  int i,j;
  for(i=0; i<numRes1; ++i)
  {
    aaHistogram[i] = new long [ numRes2 ];
    CHECK_ALLOC(aaHistogram[i]);
    for(j=0; j<numRes2; ++j)
      aaHistogram[i][j] = 0;
  }
  
  ofstream histophile(hist_filename.c_str(), ios::out);
  if (! histophile)
    ERR_INTERNAL("Cannot open file \""
                  << hist_filename
                  << "\" for writing.  Exiting...");
  
  histophile << "num rows equals length of sequence A:"
             << numRes1 << "\n";
  histophile << "num columns equals length of sequence B:"
             << numRes2 << "\n";

  assert((n_min >= 1) && (n_max <= max_num_matches) 
         && (n_min <= n_max));

  for (int rank = start_rank; rank <= end_rank; ++rank)
  {
    for(int n=n_min; (n <= n_max); ++n)
    {
      //If there is a valid rank'th-best alignment with n matches,
      //the consider this alignment when making the histogram.
      if (aaSumSqdDists[rank][n-n_min] < DistanceMetric::UPPER_BOUND)
        aaAlignments[rank][n-n_min].AccumPairHistogram(aaHistogram);
    }
  }

  for(i=0; i < numRes1; ++i)
  {
    assert(aaHistogram[i]);
    for(j=0; j < numRes2; ++j)
      histophile << aaHistogram[i][j] << " ";
    histophile << "\n";
  }

  for(i=0; i<numRes1; ++i)
    delete [] aaHistogram[i];
  delete [] aaHistogram;
} //SearchOrDyn3d::PrintPairHistogram()

#endif //#ifdef CREATE_PAIR_HISTOGRAM






#ifdef CREATE_PROFILE_HISTOGRAM

void
SearchOrDyn3d::PrintNumLayersHistogram(const char *filename) const
{
  assert(filename);
  assert(max_num_matches >= 1);
  int windowsize1 = numRes1 - n_min + 1;
  int windowsize2 = numRes2 - n_min + 1;
  assert((windowsize1 >= 1) && (windowsize2 >= 1));


  ofstream hist_file(filename, ios::out);
  if (! hist_file)
    ERR_INTERNAL("PrintNumLayersHistogram() cannot open file: \""
                 << filename << "\" for writing.");

  long total_num_cells = 0;
  for(int n=1; n <= max_num_matches; ++n)
  {
    long layer_size1 = MIN(windowsize1, (numRes1 - n + 1));
    long layer_size2 = MIN(windowsize2, (numRes2 - n + 1));
    long layer_size = layer_size1 * layer_size2;
    total_num_cells += layer_size;
    hist_file << n << " "
              << aNumLayersHistogram[n-1] << " "
              << total_num_cells
              << "\n";
  }
} //SearchOrDyn3d::PrintNumLayersHistogram()

#endif //#ifdef CREATE_PROFILE_HISTOGRAM


