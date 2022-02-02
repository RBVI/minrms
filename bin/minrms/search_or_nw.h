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

#ifndef _SEARCH_OR_NW_H
#define _SEARCH_OR_NW_H
#include <sys/time.h> //required for "struct timeval" and "gettimeofday()"
#include <biopolymer.h>
#include <pair_alignment.h>
#include <dynNW.h>
#include "or_gen.h"
#include "search_or.h"

namespace minrms {

class SearchOrNW:public SearchOr
{
public:

  //Search()
  //Search() attempts to find the best alignments possible
  //(using inter-molecular criteria) between two structures,
  //when the two structures are rotated and translated (oriented)
  //in many different ways.  Alignments are generated between the
  //structures at each of these orientations and the best alignment(s)
  //is(are) retained.
  //    The set of orientations that are sampled
  //over is specified by the "og", "start" and "stop" arguments in the
  //following way:  Search() samples all orientations in the range
  //        [og[start], og[stop])
  //(that is, beginning at og[start], and stopping just before og[stop]).
  //   This function will take care of the details of finding
  //sub-optimal alignments if requested (ie. settings.num_solutions>1).
  void Search(OrientationGenerator const& og,
              long start,
              long stop);

  //Refine()
  //Refine() improves the alignment(s) found by Search(), by
  //using these alignments to move the two molecules
  //to a more optimal relative position.  Alignments are generated
  //at this new orientation, and the process is repeated until
  //either the alignment remains the same after two successive iterations,
  //or settings.refine_max_iters iterations have been performed.
  void Refine();

  //WriteMSF()
  //WriteMSF() writes the alignment(s) calculated so far to a file
  //in MSF format.
  void WriteMSF(bool  show_number_lines = true,
                bool  compress_using_lower_case = false) const;


  //WriteChimeraInfo()
  //WriteChimeraInfo() writes additional information specific to
  //Chimera, to a few files.  Unlike WriteMSF(), the files created
  //by this function are not intended for use outside Chimera.
  void WriteChimeraInfo() const;

  void WriteMidasGraphics(bool show_residue_markers = false,
                          bool connect_the_dots = true) const;

  void WriteTransform() const;





  //The following class contains settings required exclusivly by the the SearchOrNW class.
  struct Settings: public virtual SearchOr::Settings
  {
    //
    //      "nw_which_criteria" selects the method used to decide how
    //                          many matches will be in the alignment.
    //                          It must be one of the following:

    enum NwWhichCriteria{SINGLE_ANSWER_FIXED_RMSD, //These values select which
                         SINGLE_ANSWER_FIXED_N,    //criteria (or method) will use
                         SINGLE_ANSWER_MAX_DIST,   //for selecting the gap penalty
                         SINGLE_ANSWER_MANUAL_GAP};//that determins the alignment.

    NwWhichCriteria nw_which_criteria;

    //        SINGLE_ANSWER_FIXED_RMSD <--> stop adding matches when you exceed
    //                                      a given RMSD.
    //                                      (specified by "nw_fixed_rmsd")
    //        SINGLE_ANSWER_FIXED_N,   <--> stop adding matches when you exceed
    //                                      a given number of matches.
    //                                      (specified by "nw_fixed_n")
    //        SINGLE_ANSWER_MAX_DIST   <--> stop adding matches when it is
    //                                      impossible without exceeding a
    //                                      maxumum cutoff distance
    //                                      (specified by "nw_fixed_max_dist")
    //        SINGLE_ANSWER_NW         <--> interpret "nw_fixed_gap" directly
    //                                      as the 'raw' gap penalty to use
    //                                      when generating alignments according
    //                                      to Needleman & Wunsch's method.
    //                                      (sa_parameter should have dimensions
    //                                       of distance squared.)

    Real  nw_fixed_rmsd;
    int   nw_fixed_n;
    Real  nw_max_dist;
    Real  nw_gap;

    //      "nw_max_iters" is only necessary when "which_criteria" is set to
    //                     SINGLE_ANSWER_FIXED_RMSD, or SINGLE_ANSWER_FIXED_N.
    //                          In order to get an alignment with the desired
    //                     rmsd, or number-of-matches, this function will
    //                     have to try several different "gap_penalties",
    //                     iterating to find a root.  This parameter specifies
    //                     the maximum number of times Needleman & Wunsch
    //                     (see DynNW::Solve()) will be executed, before we
    //                     give up and accept the closest answer so far.
    int   nw_max_iters;

    static const long nw_DEFAULT_MAX_ITERS = 20;


#ifdef COMPILER_BUG2
    Settings();

    //copy constructor:
    //(The alpha compiler produces buggy
    // binaries unless I define one explicitely)
    Settings(Settings const &s);
#endif //#ifdef COMPILER_BUG2

  }; //struct SearchOrNW::Settings

  SearchOrNW(Settings const& s);


private:

  Settings settings;

  SearchOrNW(const SearchOrNW &); // dissable, no implicit copying allowed
  SearchOrNW& operator =(const SearchOrNW &);//no explicit copying either

  DynNW dynNW;  //This performs the minimization of RMSD
                //at a fixed rotation and translation,

  int numRes1, numRes2;   //The number of residues in each filtered sequence.
                          //(This is saved in a variable to make the code more
                          // readable.  It's easier to write "numRes1"
                          // than "settings.aMol_f[0].size()").
  int max_num_matches;    //The maximum number of matches possible between
                          //these two proteins = MIN(numRes1, numRes2).
                          //(This is saved in a variable for ease of notation.)
  int range_of_n;         // Number of different values of n being considered.
                          // (This is just n_max-n_min+1)  Like the others,
                          // this variable just exists for notational conveniance.


  PairwiseAlignment best_alignment;
  Matrix3x4         best_orientation;

  Vect3 *aRotateUsingTheseCoords1; //(This is temporary space for coordinate data.
  Vect3 *aRotateUsingTheseCoords2; // This is so I don't have to keep re-allocating it.)


  static const char *BASE_FILE_NAME; //The beginning of the name of
                                     //any file generated by this
                                     //algorighm ("alignNW" as of 9/8/99)



  //The next function generates an alignment at the current
  //orientation between the structures aMol_c[0] and aMol_c[1].
  //It automatically picks the correct gap parameter to generate
  //an alignment with the desired properties,
  //according to the prefences specified by the user
  //in the settings.nw_which_criteria variable.
  //
  //For a return-value, it returns the RMSD of the alignment generated.
  //
  //On the outside chance it is impossible to align the two proteins
  //at this orientation (while satisfying the constraints used),
  //then FindBestAlignmentForSingleOrientation() returns
  //an impossible value: ORIENTATION_TRIVIALLY_REJECTED.
  Real
  FindBestAlignmentForSingleOrientation(PairwiseAlignment& alignment,
                                        ConstMatrix3x4 orientation_orig,
                                //the following four args are only necessary
                                //if settings.refine_while_searching == true
                                        Matrix3x4 orientation_final,
                                        Superimpose *pSuperimpose = NULL,
                                        Vect3 *aRotateUsingTheseCoords1= NULL,
                                        Vect3 *aRotateUsingTheseCoords2= NULL);

  static const Real ORIENTATION_TRIVIALLY_REJECTED;

  //The following function is invoked by WriteChimeraInfo()
  void WritePlotFile() const;


}; // class SearchOrNW

} //namespace minrms

#endif // #ifndef _SEARCH_OR_NW_H

