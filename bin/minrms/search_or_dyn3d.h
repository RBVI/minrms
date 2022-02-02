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

#ifndef _SEARCH_OR_DYN3D_H
#define _SEARCH_OR_DYN3D_H

#include <fstream>
#include <set>
using namespace std;

#include <sys/time.h> //required for "struct timeval" and "gettimeofday()"
#include <biopolymer.h>
#include <pair_alignment.h>
#include <superimpose.h>
#include <dyn3d.h>
#include "search_or.h"
#include "or_gen.h"

using namespace std;


namespace minrms {


class SearchOrDyn3d:public SearchOr
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
  //     This function will take care of the details of finding
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


  // *********** OUTPUT MEMBER-FUNCTIONS ***************


  // WriteMSF() creates a list of msf-files, each one corresponding two the
  //best alignments so-far, each one with a different number of matches
  //Each msf-file is named "align<n>.msf"  where <n> is the number of matches.
  //It prints only the alignments whose number-of-matches lies in the range
  //[min_n,max_n].
  //Note: It ignores alignments which are either not initialized
  //      (for any reason) or are so poor their RMSD exceeds the cutoff:
  //      cutoff = sqrt( Dyn3d::UPPER_BOUND_OF_COST )
  void WriteMSF(bool  show_number_lines = true,
                bool  compress_using_lower_case = false) const;


  //WriteChimeraInfo()
  //WriteChimeraInfo() writes additional information specific to
  //Chimera, to a few files.  Unlike WriteMSF(), the files created
  //by this function are not intended for use outside Chimera.
  //The file that is created is called <labelfilename> which stores
  //the names of a lot of other files that this program generates.
  //This file is used by one of Conrad's chimera delegate for displaying
  //the data.
  //  All of the char* arguments are strings storing the names
  //  of the files (or if the filename depends on the rank or number of matches
  //  made in the alignments, then the base filename is specified)
  //  that will be written out.
  //  The integer "rank" argument indicates whether we want chimera to display
  //  the 1st-best scores (then pass 0), or the 2nd-best scores (then pass 1),
  //  or the 3rd-best scores (pass 2), etc...
  void WriteChimeraInfo() const;


  //Basically, this next function dumps the state of the aTransforms[] array.
  //Note: This function is now presently unnecessary, since we've started storing
  //      this information in the MSF-files themselves, when we call "WriteMSF()".
  //      So, like "DisplayUsingMidas", this function is almost never called.
  //It creates a list of files, each one storing the transformation from
  //the second sequence's original position to it's optimal position
  //as it appears in the alignment with n matches.
  //Note: It prints only the transform-matrices for the alignments whose number
  //of matches lies in the range [n_start,n_end].  It ignores alignments which
  //not initialized (for any reason) or are so poor that their RMSD exceeds the
  // cutoff = sqrt( Dyn3d::UPPER_BOUND_OF_COST )
  void WriteTransform() const;


  //WriteMidasGraphics(), creates a list of Midas' ".gfx" formated files the user
  //can use to view the alignments visually.
  //This function is almost never used anymore, but I left it in anyway.
  //WriteMidasGraphics() creates a set of graphics-object files for midas
  // (called BASE_FILE_NAME+<n>+".gfx", where <n> is  an integer specifying how
  //  many matches were made).
  //Each graphics file just stores the green line-segments connecting pairs
  //of matched residues.
  //You can specify/limit which alignments you are interested in, by
  //specifying an interval for n, the number of matches made in the alignment.

  void WriteMidasGraphics(bool show_residue_markers = false,
                          bool connect_the_dots = true) const;



  // *********************************************************************
  // ****     The Settings class contains input data and parameters   ****
  // ****  unique to using the Dyn3d dynamic programming algorithm    ****
  // ****  (see "dyn3d.h"), which uses a score-pyramid, rather than   ****
  // ****  a score matrix.  This is all data which should be          ****
  // ****  provided directly by the user in the beginning,            ****
  // ****  and will not be modified for the life of the class.        ****
  // ****                                                             ****
  // ****  This data is needed by the the SearchOrDyn3d class         ****
  // ****  during instantiation.                                      ****
  // ****                                                             ****
  // *********************************************************************
  struct Settings: public virtual SearchOr::Settings
  {
    int   n_max;//The maximum number of matches made in an alignment
    static const int NO_MAX_NUM_MATCHES = -1; //this dissables placing an
                                              //upper bound on the number of
                                              //matches made in the alignment.

    Real  dyn3d_max_rmsd; //Places an upper bound on the
                          //maximum RMSD of an alignment
                          //To dissable this behavior set
                          //dyn3d_max_rmsd = dyn3d_NO_MAX_RMSD
    static const Real dyn3d_NO_MAX_RMSD;

    Real  dyn3d_max_dist;//The maximum distance allowed between any pair
                    //of residues matched in an alignment.
                    //(Note, that: if some cases, alignments will be chosen
                    // which are not optimal in RMSD in order to
                    // comply with the cutoff_max_dist cutoff.)
                    //When it is no longer possible to match additional
                    //pairs of residues without violating the
                    //the maximum distance constraint, we stop calculating.
                    //To dissable this behavior, set dyn2d_max_dist
                    //to "DistanceMetric::NO_MAX_PHYS_DISTANCE"


    // Important Optimization:
    //After an orientation has been calculated, it has to pass a test
    //to see if it is suitable for consideration.
    //Later on, I will call OrientationFilter::NeedlemanWunschFilter() on each
    //orientation, to see if it is suitable.  The next two parameters will
    //determine the behavior of OrientationFilter::NeedlemanWunschFilter().
    //
    //To dissable this filtering behavior, set orientation_filter_nw = false.
    bool     orientation_filter_nw; //enable orientation filtering?
    Real     orientation_filter_nw_max_dist;
    int      orientation_filter_nw_min_num_matches;

    //Search() prints out a prediction on the amount of time remaining
    //at regular intervals to the user.
    //The following parameters determine the time between the intervals.
    //(Note: This happens in addition to the statistics that are printed
    // when an update occurs.)
    //To dissable printing statistics at regular time intervals, set
    // display_stats_interval to DISSABLE_DISPLAY_STATS_INTERVAL.
    long display_stats_interval;
    static const long DISSABLE_DISPLAY_STATS_INTERVAL = -1;

#ifdef COMPILER_BUG2
    Settings();

    //copy constructor:
    //(The alpha compiler produces buggy
    // binaries unless I define one explicitely)
    Settings(Settings const &s);
#endif //#ifdef COMPILER_BUG2

  }; //struct SearchOrDyn3d::Settings





  SearchOrDyn3d(Settings const& s); //constructor requires the settings above.
  ~SearchOrDyn3d();





private:


  // ***********************************************************
  // ***                                                     ***
  // ***         Implementation Details:                     *** 
  // ***                                                     ***
  // ***********************************************************

  Settings settings;

  SearchOrDyn3d(const SearchOrDyn3d &); // dissable, no implicit copying allowed
  SearchOrDyn3d& operator =(const SearchOrDyn3d &);//no explicit copying either



  Dyn3d dyn3d;           //This class is responsible for generating
                          //all alignments.  It generates an alignment of
                          //minimal RMSD, once the relative position
                          //between the two structures is known and fixed.
                          //See <dyn3d.h> for details.

  int numRes1, numRes2;   //The number of residues in each filtered sequence.
                          //(This is saved in a variable to make the code more
                          // readable.  It's easier to write "numRes1"
                          // than "settings.aMol_f[0].size()").
  int max_num_matches;    //The maximum number of matches possible between
                          //these two proteins = MIN(numRes1, numRes2).
                          //(This is saved in a variable for ease of notation.)
  int n_min, n_max;       //Shorthand for "settings.n_min" and "settings.n_max"
                          //These symbols are used so often, the code is easier
                          //to read if I give them names in shorthand.
  int range_of_n;         // Number of different values of n being considered.
                          // (This is just n_max-n_min+1)  I calculate it
                          // once and refer to it many times later.

  PairwiseAlignment alignment_temp;//(temporary work-space.
                                   // I want to avoid calling new and delete
                                   // more than I have to.
                                   // Instead of allocating and de-allocating
                                   // PairwiseAlignment variables frequently,
                                   // I keep one around that is pre-allocated.)

  Vect3 *aRotateUsingTheseCoords1; //(Again, this is temporary space for coordinate data.
  Vect3 *aRotateUsingTheseCoords2; // This is so I don't have to keep re-allocating it.)

  long  *aNumLayersHistogram; //Enables minrms to print out the distribution of
                              //how many matches it was possible to make (N)
                              //at each orientation, given the constraints
                              //specified by the user:
                              //    dyn3d_max_rmsd, and dyn3d_max_dist
                              //This gives a measure of the quality of
                              //the orientations considered and the efficiency.
                              //(Note: As of 8/31/1999, I've dissabled this
                              //       feature because it was not useful
                              //       to most users. -Andrew)



  // ***************************************************************
  // ****    The following 4 members are used exclusively for   ****
  // **** handling the search for multiple alternate solutions  ****
  // **** (ie. 2nd-best, 3rd-best alignments)                   ****
  // ***************************************************************

  Real         **aaSumSqdDists; // Stores the 'num_solutions' lowest-sum-squared
                           //distances of alignments
                           //separately, for each number of matches made,
                           //from 1, to the length of the shorter sequence.
                           //  This multidimensional array is indexed according
                           // to: aaSumSqdDists[r][n],
                           // where n+1 is the number of matched-pairs-of-
                           //           residues in the alignment, and
                           //   and r+1 is the rank (1st-best, 2nd-best, etc.)
                           //           of the alignment.  (lower r is better)
  Matrix3x4     **aaTransforms; //Stores the correspoding rotation and
                                //translation that get's applied to
                                //molecule aMol[MOL2] (the 2nd protein)
                                //for this lowest-cost allignment.
                              //This array is indexed the same way as aaSumSqdDists[][]
  PairwiseAlignment **aaAlignments; //Stores the actual alignments associated
                                    //with the cost stored in aaSumSqdDists[][]
                              //This array is indexed the same way as aaSumSqdDists[][]

  bool alt_method_requires_alt_solution_table; //Create temporary files to speed
                                               //up search for alternate solutions?


  //This stores the prefix of the name of all the files generated by
  //DisplayUsingMidas() ".gfx", MakeMSF() ".msf",
  // and WriteTransformFiles() ".trans"
  // The filenames will be of the form 
  //   BASE_FILE_NAME + <n> + extension(either ".gfx", ".msf", or ".trans")
  // where <n> indicates how many pairs of residues got matched in that
  // alignment.
  static const char * BASE_FILE_NAME;  //Right now (9/8/1999) it is set to "align"


  // *****************************************
  // ****    Private Member Functions:    ****
  // *****************************************

  // FirstPass() finds the alignments with the lowest rmsd for each
  // number of matches possible, over the set of orientations indicated
  // by the "og", "start", and "stop" arguments.  While alternate alignments
  // are not calculated, FirstPass() makes sure all,
  // promissing-looking suboptimal solutions are recorded in a backup table
  // so they can be considered as candidates for alternate-best alignments.
  // 
  // The orientations that are sampled over are orientations in the range:
  //        [og[start], og[stop])
  // (that is, starting at og[start], and ending just before og[stop])
  void FirstPass(OrientationGenerator const& og,
                 long start,
                 long stop);


  //Copies the sum_sqd_dist, orientation, and alignment
  //into the aaSumSqdDists[][], aaTransforms[][], and aaAlignments[][] arrays.
  void RecordSolution(int r,
                      int n,
                      Real sum_sqd_dist,
                      ConstMatrix3x4 orientation,
                      PairwiseAlignment const& alignment);

  //Copies the alignment from the dyn table into the alignment member.
  //If (settings.refine_while_searching == true), then it also
  //also calculates the "reduced" RMSD of this alignment
  //(after optimal rotation and translation).
  void ExtractDynTableSolution(int               n,
                               PairwiseAlignment &alignment,
                               Real&             sum_sqd_dist_orig,
                               Real&             sum_sqd_dist_reduced,
                               Matrix3x4         orientation_reduced,
                               Superimpose       *pSuperimpose = NULL);

  //When calculating an alignment with n matches at a particular relative
  //orientation, this function checks to see if any of the distances
  //between residues in the alignment generated were larger
  //than the value set by DistanceMetric::SetMaxDistance(),
  //then this function will return true.
  bool MaxDistanceExceeded(int n);

  //As long as settings.dyn3d_max_rmsd != settings.dyn3d_NO_MAX_RMSD
  //(in which case it returns false),
  //MaxRMSDExceeded() returns true whenever an rmsd constraint was violated.
  //(This very simple function basically just returns
  // sqrt(sum_sqd_dist/n) >= settings.dyn3d_max_rmsd.
  bool MaxRMSDExceeded(Real sum_sqd_dist, int n);


  //WritePlotFile() prints out a file with two collumns containing
  //number-of-matches-made vs. RMSD of best alignment matching that number.
  //It prints only the costs of alignments whose number-of-maches lie in
  //the range [min_n,max_n], as well as  the Levitt&Gerstein P_str value.
  //WritePlotFile() is not called by the user, but invoked by
  //WriteChimeraInfo().
  //Note: It ignores alignments which are either not initialized
  //      (for any reason) or are so poor their RMSD exceeds the cutoff:
  //      cutoff = sqrt( Dyn3d::UPPER_BOUND_OF_COST )
  void WritePlotFile(string plot_filename_suffix) const;







  // **************************************************************
  // ***                                                        ***
  // ***     The following code is necessary to implement the   ***
  // *** calculation of "alternate" or sub-optimal alignments:  ***
  // ***                                                        ***
  // ***                                                        ***
  // ***   This is where the code gets especially ugly.         ***
  // ***                                                        ***
  // ***   If you are working on minrms, you may want to just   ***
  // ***   delete all of these functions.  They are only        ***
  // ***   used for generating alternate alignments, and        ***
  // ***   very few people seem to use that feature anyway.     ***
  // ***      -Andrew 9/8/1999                                  ***
  // ***                                                        ***
  // **************************************************************




  //   Note:
  //    Looking for 2nd-best and 3rd-best solutions is something
  //we can not do on the fly (ie. the first time the alignments are calculated
  //because we also require that all the alternate alignments,
  //(ie, 1st-best and 2nd-best alignments etc.) be "sufficiently different"
  //from eachother (as measured by the RMSD of their orientations)
  //   The problem is, we need to know the 1st-best solution in advance of
  //choosing the 2nd-best solution (an explanation of why is provided below)
  //and if we calculate solutions on the fly, the pseudo-random order we
  //search orientation space may not present the solutions to us in that order.
  //  Suppose, for example, you were to try and keep track of the 1st-best and
  //2nd-best solutions encountered so-far on the fly.  Half-way through the
  //search, if you discard a candidate alignment because it's too similar
  //to the 1st-best alignment so-far.  Then suppose later on (3/4 of the way),
  //you find an alignment that's even better than the 1st-best
  //solution so-far, you may have inadvertently thrown way your 2nd-best
  //solution.
  //  So what we do is keep track of a very large large pool of previously
  //calculated alignments, saving them in a temporary file.
  //In this file, I save each "promising-looking" (see "RMSD_tolerance")
  //alignment so they can be examined later.
  //This is done separately for each different number-of-matches (n),
  //hence there are multiple temporary files.

  // . . . So, For a given number of matches made, n,
  // history of solutions calculated by FirstPass()
  // is saved to a file, sorted by RMSD, and later recovered into a data structure.
  // This list of solutions is scanned, and the 1st, 2nd, 3rd-lowest
  // sufficiently different RMSD alignments in preference of the order
  // they appear in the list.  For this implementation, I used suggestions
  // from Nick, Greg, and Eric. Andrew, 11/14/1998
  //
  //Implementation Low-Level Details:
  //  To save space, I don't store the entire alignment from each solution
  //  in the file, just the orientation that created it, and it's RMSD
  //  -see the "OrientationUsageInstance" data structure.
  //  If an alignment is found to be 2nd-best or 3rd-best, then it can
  //  be reconstructed from this orientation again by calling
  //  Dyn3d::Solve(), or
  //  Dyn3d::Solve_IncrN()
  //    During the process of reconstructing the alignments we want to keep,
  //  we use some careful book-keeping.
  //  Keep in mind that a single orientation generates many alignments
  //  with differing numbers-of-matched-pairs-of-residues (n), and many of
  //  these alignments could potentially be ones we are interested in
  //  saving, either as best, 2nd-best, or 3rd-best (etc.) solutions.
  //  To figure out where to save all of these reconstructed alignments
  //  we use a data structure, that takes an orientation and returns a list
  //  of which alignments which are wanted, and where to save them.
  //  (see the "OrientationOwners" structure, and the "CalcAltSolutions()"
  //   member function).



  //The following struct is used to keep track of whatever critical information
  //about alternate alignments that can be saved in a file.
// DEC compiler doesn't like this as private...
public:
  struct OrientationUsageInstance
  {
    long  orig_orientation_id;//stores an integer identifying the
             //orientation that generated the alignment.
             //(Remember, a sequence of orientations is created by the
             //"OrientationGenerator" class, either stored in advance,
             // or calculated on the fly.  This number identifies where
             //this orientation lies in the sequence.)

    Real rmsd;//the rmsd of the alignment at this orientation with n matches

    Matrix3x4 orientation_reduced;//This is either uninitialized,
                           //or it stores the orientation after
                           //optimal superposition has been applied.
                           //(It is amost definitely not the same orientation
                           // as the one indicated by orig_orientation_id.)

    bool operator < (SearchOrDyn3d::OrientationUsageInstance const& B) const
    { return (this->rmsd < B.rmsd); }
  };
private:




  //SortAltList() uses the SolutionInfo data structure:
public:
  struct SolutionInfo
  {
    PairwiseAlignment alignment;
    Real sum_sqd_dist;
    Matrix3x4 transform;
    bool operator < (SearchOrDyn3d::SolutionInfo const& B) const
    { return (this->sum_sqd_dist < B.sum_sqd_dist); }
  };
private:


  //The following struct is used to store information about when a
  //particular orientation produced an alignment that was one of the
  //best, or alternate-best alignments.
  //   If the same orientation was needed for multiple alignments with
  //differing values of n, or r, we don't want to end up recalculating
  //a complete the alignments from scratch for each time.
  //   Instead, these structs will be stored associatively indexed by
  //orientation-number.  This container will store only the orientations that
  //are optimal.  At each orientation, it stores the number of matches in the
  //alignment that needs it, and it's rank (1st best, 2nd best, 3rd best, etc...)
  //   This is so that all the minimum-rmsd alignments at this orientation
  //will only have to be calculated once, because once they are,
  //we will be able to know where to send the results
  //(ie. to all of that orientation's "owners").
  //This data structure is used by ReadAltOrientationFile() and CalcAltSolutions().
public:
  struct OrientationOwner
  {
    int n; //the number of matches made in the alignment needing this orientation
    int r; //the rank number (ie, best, second-best, 3rd-best etc)
           //(numbering begins at 0, not 1, so first-best has r=0)
    bool operator < (SearchOrDyn3d::OrientationOwner const& B) const
    { return ((this->n < B.n) || ((this->n == B.n) && (this->r < B.r))); }
  };
private:


  fstream **apAltOrientationFiles; //An array of temporary files used for storing
                        //potential orientations that produce 2nd-best,
                        //3rd-best, etc. alignments.
  long *aNumAltOrientations; //An array storing the number of orientations
                             //stored in each of the apAltOrientationFiles[].
  //  After the search over orientation space is finished,
  //each file should contain a list of all promissing-looking alignments.
  //Then, the alignments in a file are loaded into a large array,
  //and sorted according to RMSD. (see "ReadAltOrientationFiles()" and "ReadAltOrientationFile()")


  //The next three functions are used to find alternate solutions.
  //Depending on the value of settings.alt_method,
  //one of these three functions get's invoked to calculate alternate
  //solutions immediately after FirstPass() is invoked.

  void CalcAltSolutionsUsing_ALT_METHOD_PAIRS(OrientationGenerator const& og,
                                              long start,
                                              long stop);
  void CalcAltSolutionsUsing_ALT_METHOD_3D(OrientationGenerator const& og,
                                           long start,
                                           long stop);
  void CalcAltSolutionsUsing_ALT_METHOD_3D_ORIG(OrientationGenerator const& og,
                                                long start,
                                                long stop);

  //     Unlike "CalcAltSolutions_ALT_METHOD_PAIR()", this function is
  //a lot more efficient, because it does not have to search the entire
  //space of orientations.  It uses a database which tells it which
  //orientations produce the alternate alignments, to speed up this search.
  //     The other command line options have exactly the same meaning as they
  //do in "FindBestAlignmentsOverMultipleOrientations()"
  // CalcAltSolutionsUsingOrientation() is invoked by both:
  //      CalcAltSolutionsUsing_ALT_METHOD_3D() and
  //      CalcAltSolutionsUsing_ALT_METHOD_3D_ORIG()
  void CalcAltSolutionsUsingOrientation(OrientationGenerator const& og);



  // MultiPass() is a general-purpose function almost identical
  // to "FirstPass()", which runs the search over the set of orientations
  // "use_these_orientations" multiple times, this time,
  // looking for alternate-best alignments,
  // (instead of the "best" alignments, which is what "FirstPass()" does.)
  //   This is not the most efficient way of generating alternate alignments,
  // In fact it is _very_slow_, however it is also very flexible.
  // (The "CalcAltSolutions_ALT_METHOD_PAIR()" function is implemented
  //  using this function.)

  void MultiPass(OrientationGenerator const& og,
                 set<long> const& use_these_orientations,
                 bool (*pDissimilar)(PairwiseAlignment const &,
                                     PairwiseAlignment const &));

  //Find the "r"th best solutions over the set of
  void MultiPass(int r,
                 OrientationGenerator const& og,
                 set<long> const& use_these_orientations,
                 bool (*pDissimilar)(PairwiseAlignment const &,
                                     PairwiseAlignment const &));



  void InitAltSolutionTable();
  void CleanUpAltSolutionTable();

  //KeepTrackOfAltSolutions()
  // ...sends information about the quality of an alignemnt with n matches
  //which was generated at orientation indicated by "orig_orientation_id"
  //to a table (..er file, actually) storing data on each solution calculated.
  //Later on this table will be examined for potential alternate solutions.
  //-If (settings.refine_while_searching == true), it also saves the orientation
  // that minimizes the RMSD of this alignment in the file.  This last piece
  // of information is passed via the "orientation_reduced" argument.
  //-The last argument, "alignment" is ignored for now (9/8/1999)
  // Later on you may want to save the actual alignment as well in this file
  // as well.  (I don't do this now because it would make the file too big!)
  void KeepTrackOfAltSolutions(int n,
                               long orig_orientation_id,
                               Real sum_sqd_dist,
                               Matrix3x4 orientation_reduced,
                               PairwiseAlignment const& alignment);


  void ReadAltSolutionFiles(OrientationGenerator const& og,
                            multimap<long, OrientationOwner>& orientation_owners);

  void ReadAltSolutionFile(int n,
                           fstream& alt_file,
                           int num_alt_orientations,
                           OrientationGenerator const& og,
                           multimap<long, OrientationOwner>&
                           orientation_owners);

  //the next particular function is called by ReadAltSolutionFile().
  bool AltSolutionsSufficientlyDissimilar(int n,
                                          int r,
                                          int r_prev);

  void
  GetOrientationsOfLowRmsdFromAltTable(OrientationGenerator const& og,
                                       set<long>& desirable_orientations);


  void CalcAltSolutionsFromAltSolutionTable(OrientationGenerator const& og,
                                            multimap<long, OrientationOwner>&
                                            orientation_owners);


  //The following function is only ever called by ReadAltSolutionFile().
  //See ReadAltSolutionFile() (in search_or_dyn3d.cc) for usage context.
  bool SufficientlyDifferentOrientationFromRthAlternate(
                           OrientationUsageInstance const& oui,
                           long      *aAltOrigOrientationIDs,
                           Matrix3x4 *aAltReducedOrientations,
                           int  r,
                           OrientationGenerator const& og);


  //  After iterative-refinement has been applied to the alternate alignments,
  //  things may have changed.  
  //  They may be out of order, ie the 1st-best alignment may have an
  //  RMSD that exceeds that of the 2nd-best alignment.
  //  They will have to be re-sorted.
  //  Additionally, several of the alignments may have become very similar
  //  to eachother.
  //  I use the following function to sort the list of solutions by rmsd,
  //  detect redundancy, and remove redundant alignments from
  //  the list of solutions.
  void CheckAlternatesForConsistency();

  //CheckAlternatesForConsistance() calls SortAltList()
  void SortAltList(int n);




  // *** Output functions for alternate alignments:

  void WriteMSF(int rank, //limit to only printing the "rank"-th best solutions
                bool  show_number_lines = true,
                bool  compress_using_lower_case = false) const;

  void WriteTransform(int rank//limit to only printing
                              //the "rank"-th best solutions
                      ) const;

  void WriteMidasGraphics(int rank,//limit to only printing the
                                   //"rank"-th best solutions
                          bool show_residue_markers = false,
                          bool connect_the_dots = true) const;

  void WriteChimeraInfo(int rank) const; //limit to only printing the
                                         //"rank"-th best solutions

  void WritePlotFile(int rank,  //limit to only printing the
                                //"rank"-th best solutions
                     string plot_filename_suffix) const;





  // ***************************************************************
  // ***************************************************************
  // ****    End of the section of this class used for          ****
  // **** handling the search for multiple alternate solutions  ****
  // **** (ie. 2nd-best, 3rd-best alignments) Whew!             ****
  // ***************************************************************
  // ***************************************************************




#ifdef CREATE_PAIR_HISTOGRAM
  //The next function prints out a big two-dimensional array of
  //whitespace-delimited long-integers.  Each one being the number of times
  //in the set of alignments (that were accepted in the final results)
  //that a particular residue (i) from the first sequence (A)
  //was matched with a particular residue (j) from the second sequence (B)
  //  The first two lines of this file are:
  //"num rows equals length of fixed sequence1:"
  //(which is the number of rows in the array)
  //(sequence 1 is fixed and does not reorient itself with each alignment)
  //"num collumns equals length of moveable sequence 2:"
  //(which is the number of collumns).
  //(sequence 2 can move and reorients itself to minimize the RMS_D)
  //
  //The files are named filename_base<r>, where <r> is the rank
  //of the solution set the histogram represents.
  //"r" ranges from [1 .. num_solutions].  If num_solutions > 1, multiple
  //pair-histogram files are generated one for each alternate solution named:
  //"<filename_base>1.txt", "<filename_base>2.txt", "<filename_base>3.txt" ...
  //containing only pairs of residues matched in the 1st, 2nd, and 3rd best
  //alignments, respectively.  Additionally, a single file named
  //"<filename_base>.txt" is created which counts the number of times any pair
  //of residues were matched in any of the solutions (1st, 2nd, 3rd best etc.).
  //
  void PrintPairHistogram(string filename_base) const;

private:
  // Invoked by PrintPairHistogram(string filename_base)
  void PrintPairHistogram(int start_rank,
                          int end_rank,
                          string hist_filename) const;


#endif //#ifdef CREATE_PAIR_HISTOGRAM


#ifdef CREATE_PROFILE_HISTOGRAM
public:
  // **************************************************
  // ******   "CREATE_PROFILE_HISTOGRAM" is      ******
  // ******  only for assessing the efficiency   ******
  // ******  of the algorithm.                   ******
  // **************************************************

  //Prints the contents of the "aNumLayersHistogram[]" array to a file.
  //See "aNumLayersHistogram[]" above for a description of this variable.
  void PrintNumLayersHistogram(const char *filename) const;
#endif //#ifdef CREATE_PROFILE_HISTOGRAM



}; // class SearchOrDyn3d

} //namespace minrms

#endif  // #ifndef _SEARCH_OR_DYN3D_H



