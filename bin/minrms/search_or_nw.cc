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

#include <fstream>
#include <sstream>
using namespace std;

#include <sys/time.h> //required for "struct timezone" and "gettimeofday()"
#include <global_utils.h>
#include <simple_numeric_utils.h>
#include <vect_3d.h>  //required for referrence to Vect3
#include <short_string.h>
#include <apply_trans.h>
#include <biopolymer.h>
#include <superimpose.h>
#include <dynNW.h>
#include <distance_metric.h>
#include "search_or_nw.h"


//The following line includes an RCS-version-string in the binary (I think).
//static char rcs_id[] = "@(#)$Header: /usr/local/src/minrms/bin/minrms/RCS/search_or_nw.cc,v 3.4 2003/03/31 21:28:54 conrad Exp $";


using namespace minrms;



#ifndef NW_PRECOMPUTE_CIJ
Then generate a syntax error because if this flag is;
not defined, then it probably was also not defined;
at the time "dynNW.cc" and "dynNW.h" were compiled;
If so, then the linker will give us error saying that;
DynNW::ComputePairwiseMatchCosts();
(among other functions), is not defined;
(Just make sure that this flag was defined when compiling "dynNW.cc");
#endif //#ifndef NW_PRECOMPUTE_CIJ



#ifdef COMPILER_BUG2
SearchOrNW::Settings::
Settings()
 :SearchOr::Settings(),
  nw_which_criteria(SINGLE_ANSWER_MANUAL_GAP),
  nw_fixed_rmsd(0.0),
  nw_fixed_n(n_min),
  nw_max_dist(0.0),
  nw_gap(0.0),
  nw_max_iters(nw_DEFAULT_MAX_ITERS)
{}

  //copy constructor:
  //(The alpha compiler produces buggy
  // binaries unless I define one explicitely)

SearchOrNW::Settings::
Settings(Settings const& s)
 :SearchOr::Settings(s),
  nw_which_criteria(s.nw_which_criteria),
  nw_fixed_rmsd(s.nw_fixed_rmsd),
  nw_fixed_n(s.nw_fixed_n),
  nw_max_dist(s.nw_max_dist),
  nw_gap(s.nw_gap),
  nw_max_iters(s.nw_max_iters)
{}
#endif //#ifdef COMPILER_BUG2



//some constants
const char *SearchOrNW::BASE_FILE_NAME = "alignNW";



SearchOrNW::SearchOrNW(SearchOrNW::Settings const& s):
  settings(s),
  dynNW(s.aMol_f[0].size(), s.aMol_f[1].size(), s.n_min)
{
  assert(settings.aMol);
  assert(settings.aMol_f);

  aMol_c = new Biopolymer[2];
  for (int m = 0; m < 2; ++m)
  {
    aMol_c[m] = settings.aMol_f[m];
  }
  //use "numRes1" and "numRes2" for shorthand.
  numRes1 = settings.aMol_f[0].size();
  numRes2 = settings.aMol_f[1].size();
  max_num_matches = MIN(numRes1, numRes2);
  best_alignment.reserve(max_num_matches);

  aRotateUsingTheseCoords1 = new Vect3[aMol_c[0].size() *
                                       Biopolymer::Residue::NumBackboneAtoms()];
  aRotateUsingTheseCoords2 = new Vect3[aMol_c[1].size() *
                                       Biopolymer::Residue::NumBackboneAtoms()];
  CHECK_ALLOC(aRotateUsingTheseCoords1 && aRotateUsingTheseCoords2);
} // SearchRotSpace::SearchOrNW()




void SearchOrNW::
Search(OrientationGenerator const& og,
       long start,
       long stop)
{
  assert((start >= 0) && (start < og.size()));
  assert((stop > 0) && (stop <= og.size()));
  assert(stop >= start);

  //"lowest_rmsd" stores the lowest RMSD of any of the subset
  //of alignments made so far containing the maximum number of
  //matches, (given the constraints) of any of the alignments
  //made so far.
  const Real FIRST_ALIGNMENT = -1.0f;
  Real lowest_rmsd = FIRST_ALIGNMENT;
  //(give it an initial (impossible) default value (-1.0f).)
  
  //This following variables are only necessary
  //if the "refine_while_searching" argument is set to "true"
  Superimpose superimpose(max_num_matches
                          * Biopolymer::Residue::NumBackboneAtoms());


  //Now, loop over all orientations
  int num_iters_since_last_update = 0;

  for(long o_id = start;
      o_id != stop;
      ++o_id)
  {

    Matrix3x4 const& orientation_orig
      = og[o_id].transform;

    PairwiseAlignment alignment; //Stores the alignment 
                                 //with the most matches,
                                 //and then, of secondary
                                 //importance, the lowest RMSD
                                 //at this orientation,
                                 //(hopefully).
    alignment.reserve(max_num_matches);

    Matrix3x4 orientation; //the orientation after possible optimization.

    Real rmsd =
      FindBestAlignmentForSingleOrientation(alignment,
                                            orientation_orig,
                                            orientation,
                                            &superimpose,
                                            aRotateUsingTheseCoords1,
                                            aRotateUsingTheseCoords2);

    DEBUG_MSG(DBG_SEARCH_OR_NW,
              "At current orientation, n = "
              << alignment.NumMatches() <<
              ", and RMSD = " << rmsd);

    // ***** Now that we've calculated the new "best" alignment at this
    // ***** orientation, see if it is better than the "best" alignment
    // ***** at ALL orientations.
    // ***** (Best according to number-of-matches-made, given
    // *****  the constraints, and the RMSD).
    if ((rmsd != ORIENTATION_TRIVIALLY_REJECTED)
        &&
        (
         (lowest_rmsd == FIRST_ALIGNMENT)//if this is the first iteration
         ||                              //so far in the loop over orientations,
         (alignment.NumMatches() > best_alignment.NumMatches())
         ||                              //or this alignment has more matches
         ((alignment.NumMatches() == best_alignment.NumMatches())
          &&                             //or this alignment has the same number
          rmsd < lowest_rmsd)            //of matches and a lower RMSD.
                                         //THEN, the new alignment is the "best"
         )
        )
    {
      best_alignment = alignment;
      CopyMatrix3x4(orientation, best_orientation);
      lowest_rmsd = rmsd;

      cout << " A new best alignment was found, orientation #"
           << o_id << "\n"
           << " with n = " << alignment.NumMatches()
           << " and RMSD = " << rmsd
           << ";   orientations since last update: "
           << num_iters_since_last_update << endl;

      num_iters_since_last_update = 1;
    }
    else 
      ++num_iters_since_last_update;

  } // loop over all orientations "for (o_id = start; o_id != end; ++o_id)"

} // SearchOrNW::Search()





const Real SearchOrNW::ORIENTATION_TRIVIALLY_REJECTED = -1.0f;

Real
SearchOrNW::FindBestAlignmentForSingleOrientation(
                                                  PairwiseAlignment& alignment,
                                                  ConstMatrix3x4 orientation_orig,
                                                  //the following four args are only necessary
                                                  //if settings.refine_while_searching == true
                                                  Matrix3x4 orientation_final,
                                                  Superimpose *pSuperimpose,
                                                  Vect3 *aRotateUsingTheseCoords1,
                                                  Vect3 *aRotateUsingTheseCoords2)
{
  assert(! (settings.refine_while_searching && !pSuperimpose));
  assert(! (settings.refine_while_searching && !aRotateUsingTheseCoords1));
  assert(! (settings.refine_while_searching && !aRotateUsingTheseCoords2));

  Real old_cutoff_dist = DistanceMetric::GetMaxPhysicalDistance();
  if (settings.nw_which_criteria == Settings::SINGLE_ANSWER_MAX_DIST)
  {
    DistanceMetric::SetMaxPhysicalDistance(settings.nw_max_dist);
  }

  Real min_distance_sqd;
  Real max_distance_sqd;
  ApplyTransformBackbone(orientation_orig,
                         settings.aMol_f[1],
                         aMol_c[1]);


  //To speed up the process of performing multiple iterations,
  //of dynamic programming, precompute all the
  //pairwise-distances (squared) between
  //all the residues from either sequence.
  dynNW.ComputePairwiseMatchCosts(aMol_c[0],
                                  aMol_c[1],
                                  &min_distance_sqd,
                                  &max_distance_sqd);

  //Paranoid checking for dangerous condition:
  //Suppose RMSD is the criteria being used (ie. SINGLE_ANSWER_FIXED_RMSD)
  //If the orientation doesn't even bring any single pair of
  //residues to within the desired RMSD, then we should not
  //try to binary search the interval 
  //Instead, skip to the next orientation.
  //(This will hardly ever happen, but if it did happen,
  // it would cause this code to go into an infinite loop.)
  if ((settings.nw_which_criteria == Settings::SINGLE_ANSWER_FIXED_RMSD)
      &&
      (sqrt(min_distance_sqd / Biopolymer::Residue::NumBackboneAtoms())
       >
       settings.nw_fixed_rmsd))
  {
    DEBUG_MSG(DBG_SEARCH_OR_NW,
              "At current orientation all residues of either\n"
              "                 molecule are farther apart than: "
              << settings.nw_fixed_rmsd << "\n"
              "                 This orientation will be ignored.");
    return ORIENTATION_TRIVIALLY_REJECTED;
  }


  Real G;
  Real G_min = min_distance_sqd / 2.0;
  Real G_max = max_distance_sqd / 2.0;
  assert(G_max >= G_min);

  //Now choose the value of the gap penalty "G":
  //(depends on "which_criteria")
  switch(settings.nw_which_criteria)
  {
  case Settings::SINGLE_ANSWER_FIXED_RMSD:
  case Settings::SINGLE_ANSWER_FIXED_N:
        G = 0.5 * (G_min + G_max); //(Initially, start half way
                                   // between the min/max bounaries)
        break;
  case Settings::SINGLE_ANSWER_MAX_DIST:
        //G = DistanceMetric::UPPER_BOUND/2.01;
        G = 0.5
                * SQR(settings.nw_max_dist)
            * max_num_matches
                * Biopolymer::Residue::NumBackboneAtoms();
        break;
  case Settings::SINGLE_ANSWER_MANUAL_GAP:
        G = settings.nw_gap;
        break;
  default:
        assert(0);
        break;
  } // switch(settings.nw_which_criteria)

  // ******   Now iterate over different GAP PENALTIES, G,  *****
  // ******      until the desired conditions are met.      *****
  // First, we need to set some things up
  int iter = 0;
  Real rmsd;  //the rmsd of the most recently calculated alignment.

  //The following variables are required for SINGLE_ANSWER_FIXED_RMSD mode
  //(to help it figure out when to stop iterating).
  int  lowest_n_with_rmsd_too_high = max_num_matches + 1;
  int highest_n_with_rmsd_too_low = 0;

  bool exit_loop_early = false;
  while ((iter < settings.nw_max_iters) && (! exit_loop_early))
  {
    DEBUG_MSG(DBG_SEARCH_OR_NW_G_ITER, " iter=" << iter
              << "; G = " << G << "; [G_min,G_max] = ["
              << G_min << "," << G_max << "]");

    dynNW.InitTable(G); //Set the gap parameter to "G"
    dynNW.Solve(aMol_c[0],
                aMol_c[1],
                alignment
                #ifdef NW_PRECOMPUTE_CIJ
                ,true
                #endif
                );
    int n = alignment.NumMatches();

    rmsd = RMSD(alignment,
                aMol_c[0],
                aMol_c[1]);

    DEBUG_MSG(DBG_SEARCH_OR_NW_G_ITER, "      " << iter
              << "; n = " << alignment.NumMatches()
              << "; RMSD = " << rmsd);

    //Now decide whether we want to keep iterating, and if
    //so, pick a new value of "G":
    switch(settings.nw_which_criteria)
    {
    case Settings::SINGLE_ANSWER_FIXED_RMSD:
      //-If our alignment has an RMSD higher than the target RMSD, then in the
      // next iteration, binary search the interval from [G, G_max].
      //-If it is below the target RMSD, search the interval from [G_min, G].
      //-Special case: However, if it is only _slightly_ lower
      // than the target RMSD, that is _within_one_match_ of exceeding
      // the target RMSD, then stop iterating (since we know we cannot
      // get any closer to the target RMSD without exceeding it).
      if (rmsd <= settings.nw_fixed_rmsd)
      {
        // If the rmsd is as close to the target rmsd as possible
        // (that is, if we add one more match,we will exceed the target rmsd)
        // then we're done.  exit the loop early.
        if (n == lowest_n_with_rmsd_too_high - 1)
        {
          exit_loop_early = true;
          DEBUG_MSG(DBG_SEARCH_OR_NW_G_ITER, " RMSD just below (within one match) of\n"
                    "                           exceeding the target RMSD. exiting loop...");
        }
        else //otherwise, binary search the interval from [G,G_max]
          G_min = G;
        //now update "highest_n_with_rmsd_too_low"
        if (n > highest_n_with_rmsd_too_low)
          highest_n_with_rmsd_too_low = n;
      }
      else // then rmsd > settings.nw_fixed_rmsd
      {
        if ((n == highest_n_with_rmsd_too_low + 1) //If we are only one match
            ||                                     //above the target rmsd,
            (iter+1 == settings.nw_max_iters))     //or if we have run out
                                                   //of iterations, then
                                                   //just set G=G_min for
                                                   //the final iteration.
        {
          DEBUG_MSG(DBG_SEARCH_OR_NW_G_ITER, "--Ran out of iterations, or\n"
                    "                            --RMSD within one match of the target RMSD\n"
                    "                            --Employing an old Gap penalty\n"
                    "                            --and repeating this iteration.");
          G_max = G_min;
          //This effectively sets G = G_min in the next iteration,
          //since, later, we will set G=0.5*(G_min+G_max).Note that
          //G_min should equal the highest gap penalty used so far
          //whose resulting alignment's RMSD was below the target.
          //RMSD.  At this lower gap penalty, the alignment should
          //contain one fewer matches than this alignment, which
          //is as many matches as possible without exceeding
          //the target RMSD.

          --iter; //During the next iteration, we want to repeat the
                  //gap penalty we tried during a previous iteration (G_min).
                  //However, by then "iter" may have exceeded "sa_max_iters",
                  //so we decrement iter by one, to keep from leaving the loop
                  //early.
        }
        else //otherwise, binary-search the interval from [G_min,G]
          G_max = G;
        //now update "lowest_n_with_rmsd_too_high"
        if (n < lowest_n_with_rmsd_too_high)
          lowest_n_with_rmsd_too_high = n;
      }
      G = (G_max + G_min)/2.0;
      break;

    case Settings::SINGLE_ANSWER_FIXED_N:
      if (n < settings.nw_fixed_n)
        G_min = G; // then binary-search the interval [G_min,G]
      else if (n > settings.nw_fixed_n)
      {
        if (iter+1 == settings.nw_max_iters) //If we have run out
                                             //of iterations, then
                                             //just set G=G_min for
                                             //the final iteration.
                                             //This will insure the number
                                             //of matches in our alignment
                                             //does not exceed our target n!
        {
          DEBUG_MSG(DBG_SEARCH_OR_NW_G_ITER, "--Ran out of iterations, or\n"
                    "                            --Employing an old Gap penalty\n"
                    "                            --and repeating this iteration.");
          G_max = G_min;
          //This effectively sets G = G_min in the next iteration,
          //since, later, we will set G=0.5*(G_min+G_max).Note that
          //G_min should equal the highest gap penalty used so far
          //whose resulting alignment's RMSD was below the target.
          //RMSD.  At this lower gap penalty, the alignment should
          //contain one fewer matches than this alignment, which
          //is as many matches as possible without exceeding
          //the target RMSD.
          --iter; //During the next iteration, we want to repeat the
                  //gap penalty we tried during a previous iteration (G_min)
                  //However, by then "iter" may have exceeded "sa_max_iters"
                  //so we decrement iter by one, to keep from leaving 
                  //the loop early.
        }//if (iter+1 == sa_max_iters)
        else //otherwise, binary-search the interval from [G_min,G]
          G_max = G;
      }//else if (n > settings.nw_fixed_n)
      else {
        //if our alignment contains the desired number of matches, exit early
        exit_loop_early = true;
        DEBUG_MSG(DBG_SEARCH_OR_NW_G_ITER, " n = n_goal.  exiting loop...");
      }
      G = (G_max + G_min)/2.0;
      break;
    case Settings::SINGLE_ANSWER_MAX_DIST:
      exit_loop_early = true; //only one iteration in this mode.
      break;
    case Settings::SINGLE_ANSWER_MANUAL_GAP:
      exit_loop_early = true; //only one iteration in this mode.
      break;
    default:
      assert(0);
      break;
    } // switch(settings.nw_which_criteria)
    ++iter;
  } // while ((iter < settings.nw_max_iters) && (! exit_loop_early))

  if (settings.refine_while_searching) //This is slower, but more accurate.
  {
    //Then, pick a different orientation which
    //minimizes the rmsd using Diamond88's method.
    //Replace the old "rmsd" and "orientation" with the new.
    rmsd = alignment.CalcMinRMSD(settings.aMol_f[0],
                                 settings.aMol_f[1],
                                 orientation_final,
                                 pSuperimpose,
                                 aRotateUsingTheseCoords1,
                                 aRotateUsingTheseCoords2);
    DEBUG_MSG(DBG_REFINE_OR_NW," ..transform found.  New rmsd = "<< rmsd);
  }
  else
    CopyMatrix3x4(orientation_orig, orientation_final);
  
  DistanceMetric::SetMaxPhysicalDistance(old_cutoff_dist);
  return rmsd;
} //SearchOrNW::FindBestAlignmentForSingleOrientation()







void
SearchOrNW::Refine()
{
  //Dissable "settings.refine_while_searching" for the duration of
  //this function.  It shouldn't matter if we do this, but it's easier
  //to understand what's going this way, and it's a little faster too.
  //Since we are in the process of refinment, every new solution
  //that get's calculated will note be refined regardless of the sate
  //of settings.refine_whiles_searching
  bool old_refine_while_searching = settings.refine_while_searching;
  settings.refine_while_searching = false;


  //First, reset the DistanceMetric's maximum distance cutoff
  //according to settings.nw_which_criteria.
  //This time, we are refining until we reach convergence, so
  //there's no point to choose really accurate values of RMSD,
  //since they will be calculated eventually anyway.
  //(This argument does not hold in the case of SearchOrNW::Refine().)

  Real old_cutoff_dist = DistanceMetric::GetMaxPhysicalDistance();
  if (settings.nw_which_criteria == Settings::SINGLE_ANSWER_MAX_DIST)
    DistanceMetric::SetMaxPhysicalDistance(settings.nw_max_dist);

  //Then, create some variables we will need later.
  Superimpose superimpose(max_num_matches
                          * Biopolymer::Residue::NumBackboneAtoms());

  PairwiseAlignment improved_alignment;
  improved_alignment.reserve(max_num_matches);

  Matrix3x4 improved_orientation;
  Real      improved_rmsd;

  const Real FIRST_ITERATION = -1.0f;
  Real best_rmsd = FIRST_ITERATION;
  bool exit_loop_early = false;
  for(int   iter = 0;
      (((iter < settings.refine_max_iters)
        ||
        (settings.refine_max_iters == Settings::ITERATE_UNTIL_CONVERGENCE))
       &&
       (! exit_loop_early));
          ++iter)
  {
    DEBUG_MSG(DBG_REFINE_OR_NW, "Refine() going into "
              << iter+1 << "th orientation.");

    Real rmsd = best_alignment.CalcMinRMSD(settings.aMol_f[0],
                                           settings.aMol_f[1],
                                           improved_orientation,
                                           &superimpose,
                                           aRotateUsingTheseCoords1,
                                           aRotateUsingTheseCoords2);
    if (best_rmsd == FIRST_ITERATION)
      best_rmsd = rmsd;
    DEBUG_MSG(DBG_REFINE_OR_NW," ..transform found.  New rmsd = "<< rmsd);


    improved_rmsd = FindBestAlignmentForSingleOrientation(improved_alignment,
                                                          improved_orientation,
                                                          improved_orientation,
                                                          &superimpose,
                                                          aRotateUsingTheseCoords1,
                                                          aRotateUsingTheseCoords2);
    //TAKE THIS OUT! best_rmsd = improved_rmsd;
    DEBUG_MSG(DBG_REFINE_OR_NW,
              "  ..At new orientation, n = "
              << improved_alignment.NumMatches()
              << ",  and RMSD = " << improved_rmsd);

    // ***** Now that we've calculated the new "best" alignment at this
    // ***** orientation, see if it is better than the "best" alignment
    // ***** at ALL orientations.
    // ***** (Best according to number-of-matches-made, given
    // *****  the constraints, and the RMSD).
    if ((improved_alignment.NumMatches() < best_alignment.NumMatches())
        ||                              //or this alignment has more matches
        ((improved_alignment.NumMatches() == best_alignment.NumMatches())
         &&
         improved_rmsd >= best_rmsd))
    {
      exit_loop_early = true;
      DEBUG_MSG(DBG_REFINE_OR_NW,
                " ..No improvements made between current and previous iterations.\n"
                "   Exiting loop...");
    }
    else
    {
      best_rmsd = improved_rmsd;
      best_alignment = improved_alignment;
      CopyMatrix3x4(improved_orientation, best_orientation);
    }
  } //while ((iter < refine_max_iters) && (! exit_loop_early))


  //Now, restore some old values:
  settings.refine_while_searching = old_refine_while_searching;
  DistanceMetric::SetMaxPhysicalDistance(old_cutoff_dist);
} //SearchOrNW::Refine()





void
SearchOrNW::WriteMSF(bool show_number_lines,
                     bool compress_using_lower_case) const
{
  stringstream msf_filename;
  msf_filename << SearchOrNW::BASE_FILE_NAME
               << best_alignment.NumMatches()
               << ".msf";

  if (best_alignment.NumMatches() == 0)
  {
    cerr <<
      "--- Warning: no alignment was found satisfying the ---\n"
      "---          criteria specified by the user.       ---\n"
      "---          Consequently, no \"" << msf_filename.str() << "\"\n"
      "---          file will be created."
         << endl;
    return;
  }
  ApplyTransformBackbone(best_orientation,
                         settings.aMol_f[1],
                         aMol_c[1]);
  Real rmsd = RMSD(best_alignment,
                   aMol_c[0],
                   aMol_c[1]);
  stringstream comments;
  comments
    <<
    "Chimera pseudo-minimal RMSD structural alignment with "
    << best_alignment.NumMatches() <<
    " equivalences.\n"
    "(This alignment was generated by a variant of minrms that employs\n"
    " the algorithm by Needleman & Wunsch to generate all alignments.\n"
    " Details can be found in appendices B and C of the paper that\n"
    " describes minrms.)\n"
    "RMSD = " << rmsd << "\n"
    "-----\n"
    "Transform Matrix to apply to structure: "
    << settings.vPdbFileNames[1] << "\n"
    << best_orientation[0][0] << " " << best_orientation[0][1] << " " << best_orientation[0][2] << " " << best_orientation[0][3] << "\n"
    << best_orientation[1][0] << " " << best_orientation[1][1] << " " << best_orientation[1][2] << " " << best_orientation[1][3] << "\n"
    << best_orientation[2][0] << " " << best_orientation[2][1] << " " << best_orientation[2][2] << " " << best_orientation[2][3] << "\n";


  PairwiseAlignment alignment_temp(best_alignment.NumMatches());
  TranslateAlignmentBetweenMolecules(best_alignment,
                                     alignment_temp,
                                     settings.aMol_f[0],
                                     settings.aMol_f[1],
                                     settings.aMol[0],
                                     settings.aMol[1]);

  if (compress_using_lower_case)
    cerr << "Warning: MSF-files cannot be compressed by using upper/lower case\n"
      "       (at least in single-answer mode).\n"
      "       This feature is not supported yet.\n"
      "       (Note: this warning message could be printed merely because of\n"
      "              an internal programming blooper, not the user's fault."
      "              This message is mostly for debugging purposes.\n"
      "              Please ignore it.)\n"
         << endl;

  alignment_temp.ExportMSF(msf_filename.str(),
                           settings.aMol[0],
                           settings.aMol[1],
                           settings.vMsfOutputLabels[0],
                           settings.vMsfOutputLabels[1],
                           comments.str(),
                           show_number_lines,
                           compress_using_lower_case);

} //SearchOrNW::WriteMSF()



void
SearchOrNW::WriteMidasGraphics(bool show_residue_markers,
                               bool connect_the_dots) const
{
  stringstream gfx_filename;
  gfx_filename << SearchOrNW::BASE_FILE_NAME
               << best_alignment.NumMatches()
               << ".gfx";

  if (best_alignment.NumMatches() == 0)
  {
    cerr <<
      "------------------------------------------------------\n"
      "--- Warning: no alignment was found satisfying the ---\n"
      "---          criteria specified by the user.       ---\n"
      "---          Consequently, no \"" << gfx_filename.str() << "\"\n"
      "---          file will be created."
         << endl;
    return;
  }
  ApplyTransformBackbone(best_orientation,
                         settings.aMol_f[1],
                         aMol_c[1]);

  stringstream caption;
  caption << "n = " << best_alignment.NumMatches()
          << ", RMSD = " << RMSD(best_alignment, aMol_c[0], aMol_c[1]);
                  
  
  best_alignment.ExportGFX(aMol_c[0],
                           aMol_c[1],
                           gfx_filename.str(),
                           caption.str(),
                           show_residue_markers,
                           connect_the_dots,
                           PairwiseAlignment::MIDAS_MARKER_SIZE);

} //SearchOrNW::DisplayUsingMidas()





void
SearchOrNW::WriteTransform() const
{
  stringstream trans_filename;
  trans_filename << SearchOrNW::BASE_FILE_NAME
                 << best_alignment.NumMatches()
                 << ".trans";

  if (best_alignment.NumMatches() == 0)
  {
    cerr <<
      "--- Warning: no alignment was found satisfying the ---\n"
      "---          criteria specified by the user.       ---\n"
      "---          Consequently, no \"" << trans_filename.str() << "\"\n"
      "---          file will be created."
         << endl;
    return;
  }
  WriteMatrix3x4ToFile(best_orientation, trans_filename.str().c_str());
} //SearchOrNW::WriteTransform()




void
SearchOrNW::WriteChimeraInfo() const
{
  WritePlotFile();

  string plot_filename("align_chimera.plot");
  string info_filename("align_chimera.info");

  ofstream info_file(info_filename.c_str(), ios::out);
  if (! info_file)
    ERR_INTERNAL("Cannot open file \""
                 << info_filename
                 << "\" for writing.  Exiting...");
  info_file << "fixed pdb:" << settings.vPdbFileNames[0] << "\n";
  info_file << "moveable pdb:" << settings.vPdbFileNames[1] << "\n";
  info_file << "rmsd:" << plot_filename << "\n";
  info_file << "base:" << BASE_FILE_NAME << "\n";

  info_file << "number of atoms used:"
            << Biopolymer::Residue::NumBackboneAtoms() << "\n";
  for(int q=0; q < Biopolymer::Residue::NumBackboneAtoms(); ++q)
    info_file << Biopolymer::Residue::LookupBackboneAtomSymbol(q) << "\n";

} //SearchOrNW::WriteChimeraInfo()



void
SearchOrNW::WritePlotFile() const
{
  string plot_filename("align_chimera.plot");

  ofstream plot_file(plot_filename.c_str(), ios::out);
  if (! plot_file)
    ERR_INTERNAL("WritePlotFile() cannot create file \"" << plot_filename
                                        << "\" for writing.");
 //Now label the columns:
  plot_file << "# Matches,RMSD,Longest Distance"
    //",Orientation Increment" Commented out on 9/8/1999.we don't print this now
            << ",-log(probability)\n";

  ApplyTransformBackbone(best_orientation,
                         settings.aMol_f[1],
                         aMol_c[1]);

  Biopolymer::const_iterator worst_i, worst_j;//needed to pacify syntax.

  Real probability = LevittGerstein98::Sstr_Probability(best_alignment,
                                                        aMol_c[0],
                                                        aMol_c[1]);

  plot_file << best_alignment.NumMatches() << " "
            << RMSD(best_alignment, aMol_c[0], aMol_c[1]) << " "
            << best_alignment.FindFurthestResPair(worst_i,
                                                  worst_j,
                                                  aMol_c[0],
                                                  aMol_c[1]) << " "
            << ((probability == 0.0) ? 0.0 : -log10(probability)) << "\n";
} //SearchOrNW::WritePlotFile()



