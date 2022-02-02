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
#include <cassert>
#include <cstdio> //needed for "sprintf()"
#include <algorithm>  //needed for call to "random_shuffle()"
using namespace std;

#include <global_utils.h>
#include <short_string.h>
#include <simple_numeric_utils.h> //required for definition of "AccumDist()"
#include <vect_3d.h> //required for referrence to Vect3
#include <fast_rot_metric.h>
#include <apply_trans.h>
#include <superimpose.h>
#include <biopolymer.h>
#include <residue_names.h>
#include <pair_alignment.h>
#include "intervals.h"
#include "pair_align_settings.h"
#include "or_gen.h"

using namespace minrms;

const Real
OrientationGenerator::Settings::fm_NO_CUTOFF_RMSD = -1.0f;



OrientationGenerator::Settings::Settings()
  :PairAlignSettings(),
   fm_fragment_size(fm_DISSABLE_FRAGMENT_MATCHING),
   fm_cutoff_rmsd(fm_NO_CUTOFF_RMSD),
   fm_interval_pair_list(0),
   input_msf_files(0),
   input_matrix_files(0),
   vMsfInputLabels(2)
{
  //cerr << "SearchOrDyn3d::Settings, this = " << this << endl;
  //cerr << "default constructor called." << endl;
}


#ifdef COMPILER_BUG2
 //copy constructor:
 //(The alpha compiler produces buggy
 // binaries unless I define one explicitly)
OrientationGenerator::Settings::Settings(Settings const& s)
  :PairAlignSettings(s),
   fm_fragment_size(s.fm_fragment_size),
   fm_cutoff_rmsd(s.fm_cutoff_rmsd),
   fm_interval_pair_list(s.fm_interval_pair_list),
   input_msf_files(s.input_msf_files),
   input_matrix_files(s.input_matrix_files),
   vMsfInputLabels(s.vMsfInputLabels)
{
  //cerr << "OrGen::Settings, this = " << this << endl;
  //cerr << "s.aMol = " << s.aMol << ";  aMol = " << aMol << endl;
  //cerr << "s.aMol_f = " << s.aMol_f << ";  aMol_f = " << aMol_f << endl;
}
#endif //#ifdef COMPILER_BUG2



Real OrientationGenerator::
SumSqdDistBetweenOrientations(ConstMatrix3x4 o1, ConstMatrix3x4 o2) const
{
  Matrix3x4 o1_in_pa, o2_in_pa;
  Mult_Mat3x4_by_Mat3x4_equals_Mat3x4(o1,
                                      aPrincipleAxisAndCentroid[1],
                                      o1_in_pa);

  Mult_Mat3x4_by_Mat3x4_equals_Mat3x4(o2,
                                      aPrincipleAxisAndCentroid[1],
                                      o2_in_pa);
  return FastRotSumSqdDist(o1_in_pa,
                           o2_in_pa,
                           aPrincipleAxisEigenvalues[1],
                           (settings.aMol_f[1].size()
                            *
                            Biopolymer::Residue::NumBackboneAtoms()));
}





OrientationGenerator::OrientationGenerator(Settings const& s)
  :settings(s)
{
  // *** Setup:
  // *** Fill aaRotateUsingTheseCoords[][] and aPrincipleAxisAndCentroid[]
  // *** (These are needed later on.)
  //
  //Allocate space for the arrays that will store the
  //positions of all the quickatoms from each molecule.
  aaRotateUsingTheseCoords[0] = new Vect3[ settings.aMol_f[0].size() *
                                     Biopolymer::Residue::NumBackboneAtoms()];
  aaRotateUsingTheseCoords[1] = new Vect3[ settings.aMol_f[1].size() *
                                     Biopolymer::Residue::NumBackboneAtoms()];
  CHECK_ALLOC(aaRotateUsingTheseCoords[0] && aaRotateUsingTheseCoords[1]);

  //Now, copy the coordinates into the aaRotateUsingTheseCoords[] arrays.
  for(int m = 0; m < 2; ++m) //loop over both molecules
  {
    int match = 0;
    for(Biopolymer::const_iterator ps = settings.aMol_f[m].begin();
        ps < settings.aMol_f[m].end();
        ++ps)
    {
      for(int q=0; q < Biopolymer::Residue::NumBackboneAtoms(); ++q)
      {
        Biopolymer::Residue::const_iterator pa;
        pa = ps->GetBackboneAtom(q);
        assert(pa != ps->end());
        aaRotateUsingTheseCoords[m][match][0] = (*pa).second.xyz[0];
        aaRotateUsingTheseCoords[m][match][1] = (*pa).second.xyz[1];
        aaRotateUsingTheseCoords[m][match][2] = (*pa).second.xyz[2];
        ++match;

      }
    }
  } // for(int m = 0; m < 2; ++m)   loop over both molecules

  //Later on, we may want to invoke the SumSqdDistBetweenOrientations()
  //which measures the differences between the various orientations.
  //In order to compute this quickly, we need to precompute the principle
  //axis of the two structures.

  for(int m = 0; m < 2; ++m) //loop over both molecules
  {
    FindPrincipleAxis(aaRotateUsingTheseCoords[m],
                      settings.aMol_f[m].size() *
                      Biopolymer::Residue::NumBackboneAtoms(),
                      aPrincipleAxisEigenvalues[m],
                      aPrincipleAxisAndCentroid[m]);
  }


  //Now, we are ready to compute the orientations.
  //Again, we compute the entire list of orientations beforehand,
  //and store them in a big array.
  PrecomputeAllOrientations();

  //Now, (this part is optional) we randomize the order in the sequence
  //of orientations. (The order of the orientations is not important.)
  //We do this because if the orientations are randomly selected, it makes it
  //easier to make predictions about how much longer the search
  //over orientations based on how much time has elapsed and
  //how many orientations remain.
  random_shuffle(vO.begin(), vO.end());
} //OrientationGenerator::OrientationGenerator()




OrientationGenerator::~OrientationGenerator()
{
  if (aaRotateUsingTheseCoords[0]) delete [] aaRotateUsingTheseCoords[0];
  if (aaRotateUsingTheseCoords[1]) delete [] aaRotateUsingTheseCoords[1];
} //OrientationGenerator::~OrientationGenerator()




void
OrientationGenerator::PrecomputeAllOrientations()
{
  cout << "Generating superpositions:" << endl;
  if (settings.fm_fragment_size != Settings::fm_DISSABLE_FRAGMENT_MATCHING)
  {
    cout << " (fragment size = " << settings.fm_fragment_size << ")" << endl;
    GenerateOrientationsUsingFragmentMatching();
  }
  if (settings.input_msf_files.size() > 0)
  {
    cout << " Reading msf-files." << endl;
    GenerateOrientationsReadingMSFfiles();
  }
  if (settings.input_matrix_files.size() > 0)
  {
    cout << " Reading matrix-files." << endl;
    GenerateOrientationsReadingMatrixFiles();
  }
  if (vO.size() == 0)
  {
    //If no orientation-generating command-line-options were requested, 
    //use a single default orientation: the original orientation
    //of the second sequence as read from the pdb-file.
    //(The transformation matrix that corresponds to this orientation
    // is the identity.)
    cout << "Default superposition used.  Fragment-matching dissabled.\n"
            " (No rotation or translation will be applied to either structure)."
         << endl;

    OrGenEntry identity;
    ResetToIdentity3x4(identity.transform);
    vO.push_back(identity);
  }
  cout << "Number of orientations that will be considered "
          "in final caclulation: " << vO.size() << endl;
} //OrientationGenerator::PrecomputeAllOrientations()






void OrientationGenerator::GenerateOrientationsUsingFragmentMatching()
{
  Biopolymer const& m1 = settings.aMol_f[0];  //Shorthand for making 
  Biopolymer const& m2 = settings.aMol_f[1];  //the code more readable.

  // Now, we implement part of the "windowing optimization."
  //  If the user places a lower bound in the number of correspondences made
  //  in the alignments, "n_min", then we should never need to match generate
  //  orientations by matching certain blocks together, since the residues
  //  contained those blocks could never be matched together in any alignment
  //  containing at least n_min matches.  (See paper)
  // The variables
  // "windowsize1&2" represents the distance in each direction one
  // can stray from the main diagonal (i=j), when matching blocks together,
  //  WITHOUT preventing the possibility of generating an alignment
  //  that with at least "n_min" corresponding pairs of residues.
  int windowsize1 = m1.size() - settings.n_min + 1;
  int windowsize2 = m2.size() - settings.n_min + 1;
  assert((windowsize1 >= 1) && (windowsize2 >= 1));

  /*
  Each pair-of-intervals in "settings.fm_interval_pair_list"
  defines limits the set of fragments that can be matched together.
  Fragments from inside interval of residues can only be matched with fragments
  from the that interval's partner in the other molecule.
  The set of all possible pairs of fragments from either molecule can be
  represented by the set of pairs of integers (I,J), which lie in the range:
  0 <= I < (m1.size() + 1 - fragment_size), and
  0 <= J < (m2.size() + 1 - fragment_size).
  where I denotes the position of the first residue the fragment from molecule1
    and J denotes the position of the first residue the fragment from molecule2
  Each pair of intervals in the list defines a rectantular window within this
  set of (I,J) pairs.  The union of all of these subsets of (I,J) pairs
  is looped over.  For each I,J pair, a pair of fragments is matched
  and new orientation is generated.

                              J /|\
                                 |  
                     m2.size() ->|____________________
                          endB3->|         .---.      |
                                 |         |   |      |
                        startB3->|         `---'      |
                          endB2->|               .--. |
                                 |               |  | |
                        startB2->|               `--' |
                          endB1->| .-----.            |
                                 | |     |            |
                                 | |     |            |
                        startB1->| |_____|            |
                                 |                    |
                                 |____________________|___\
                                   |     | |   | |  | |   / I
                             startA1 endA1 |   | |  | |
                                     startA2endA2|  | m1.size()
                                           startA3 endA3

  Figure: The three rectangles denote pairs of intervals from
          "settings.fm_interval_pair_list".  Pairs of fragments
          which do not lie entirely within these rectangles
          are not matched together to generate orientations.

  --- Problem:---
  The intervals in "settings.fm_interval_pair_list" do not have to be disjoint.
  Suppose, one of these pairs of intervals overlaps with another pair
  of intervals.  If you loop through each pair of fragments in
  each pair of intervals, you may be matching the same fragments twice:
  In the figure below, the fragments beginning at I,J (marked with a *)
  gets matched twice: once inside (startA1,endA1),(startB1,endB1),
  and once inside (startA2,endA2),(startB2,endB2).
          
                J /|\        
                   |  
        m2.size()->|____________________
                   |                    |
                   |                    |
            endB2->| .---------.        |
                   | |         |        |
            endB1->| |    .----+----.   |
                   | |    |I,J |    |   |
                   | |    | *  |    |   |
          startB2->| `----+----'    |   |
                   |      |         |   |
          startB1->|      |_________|   |
                   |                    |
                   |____________________|___\
                     |    |    |    |   |   / I
                startA1   |   endA1 |   |
                      startA2     endA2 |
                                       m1.size()

  --- Solution: (Use an OR mask) ---

  We pre-computing a 2-D array of bools called aaMakeMatch[][].
  Every pair of fragments (denoted by I,J) that lies within
  at least one rectangular window has aaMakeMatch[I][J] = true.  Otherwise
  aaMakeMatch[I][J] = false.
  (I left out one detail:
   We filter out pairs of fragments that violate the windowing optimization.
   That is, if n_min > 1, and the pair of fragments denoted by I,J cannot
   be matched in an alignment with at least n_min matches,
   then aaMakeMatch[I][J] = false also.)

  After we have filled in aaMakeMatch[][], we loop through
  the entire aaMakeMatch[][] array, and only match pairs fragments
  with aaMakeMatch[I][J] == true.  This prevents the possibility
  of matching the same pair of fragments twice.
  */

  bool **aaMakeMatch = new bool* [m1.size()];
  CHECK_ALLOC(aaMakeMatch);
  for (int I = 0; I < m1.size(); ++I)
  {
    aaMakeMatch[I] = new bool [m2.size()];
    CHECK_ALLOC(aaMakeMatch[I]);
    for (int J = 0; J < m2.size(); ++J)
    {
      aaMakeMatch[I][J] = false; //initialze to false.
    }
  }

  //Now, find out which pairs of fragments will be matched and
  //set aaMakeMatch[I][J]=true for those fragment pairs
  for(long p = 0; p < settings.fm_interval_pair_list.size(); ++p)
  {
    #if DBG_INTERVAL_LISTS
    DEBUG_MSG(DBG_INTERVAL_LISTS, 
              "Generating orientations from matching blocks in interval-pair #"
              << p);
    long block_match_counter = 0;
    #endif
    
    //Retrieve the current pair of intervals:
    assert(settings.fm_interval_pair_list[p].size() == 2);

    Biopolymer::const_iterator start1
      = m1.find(settings.fm_interval_pair_list[p][0].first);

    Biopolymer::const_iterator end1
      = m1.find(settings.fm_interval_pair_list[p][0].last);

    Biopolymer::const_iterator start2
      = m2.find(settings.fm_interval_pair_list[p][1].first);

    Biopolymer::const_iterator end2
      = m2.find(settings.fm_interval_pair_list[p][1].last);

    assert((start1 != m1.end()) && (end1 != m1.end()) &&
           (start2 != m2.end()) && (end2 != m2.end()));

    //Make sure the intervals make sense.
    if ((end1 < start1) || (end2 < start2))
      ERR("Subset-Match Interval cannot have start > end.");
    else if ((start1 < m1.begin())
             ||
             (start2 < m2.begin()))
      ERR("Subset-Match Interval cannot begin before start of sequence.");
    else if ((m1.end() <= end1)
             ||
             (m2.end() <= end2))
      ERR("ERROR: fragment-size cannot exceed length of sequence!");

    //If interval is skipped because there is not enough space, print
    //a warning to the user, but do not exit.
    //  (Note: Later on I added a function which re-sizes the intervals
    //         to make sure that the intervals are large enough, so
    //         This condition should never occur.  I check for it anyway.
    //         -Andrew 8/26/1999)
    if (start1 > end1 - settings.fm_fragment_size + 1)
      cerr << "-----------------------------------------------------------\n"
           << "Note: interval ["
        //<< (start1 - m1.begin()) << " - "
           << (*start1).id << " - "
        //<< ( end1  - m1.begin())
           << (*end1).id
           << "] from sequence A is too small\n"
           << "to accomodate a subset of size "
           << settings.fm_fragment_size
           << ".\n  Skipping over this pair of intervals.\n"
           << "-----------------------------------------------------------"
           << endl;
    else if (start2 > end2-settings.fm_fragment_size+1)
      cerr << "-----------------------------------------------------------\n"
           << "Note: interval ["
        //<< (start2 - m2.begin()) << " - "
           << (*start2).id << " - "
        //<< ( end2  - m2.begin())
           << (*end2).id
           << "] from sequence B is too small\n"
           << "to accomodate a subset of size "
           << settings.fm_fragment_size
           << ".\n  Skipping over this pair of intervals.\n"
           << "-----------------------------------------------------------"
           << endl;
    else
    {
      for(Biopolymer::const_iterator i = start1;
          i <= end1-settings.fm_fragment_size+1;
          ++i)
      {

        // Now clip the interval to the diagonal strip which exists
        // whenever the windowing optimization is enabled.
        //
        // Instead of start2 and end2, the iterator through the
        // second molecule will iterate from j_lower_bound,
        // to j_upper_bound (defined below).
        //
        // The purpose of this optimization ignores pairs of "fragments"
        // whose residues could not possibly be matched in any alignment 
        // containing at least "n_min" matches.
        // The "windowing" optimization is enabled whenever the user
        // specifies a value of "--cutoff-min-n" larger than 1.
        

        Biopolymer::const_iterator
          j_lower_bound = MAX(start2,
                              (m2.begin()
                               + (i - m1.begin())
                               - (windowsize1 - 1)
                               )
                              );

        Biopolymer::const_iterator
          j_upper_bound = MIN( (m2.begin()
                                + (i - m1.begin())
                                + (windowsize2 - 1)
                                ),
                              (end2 - settings.fm_fragment_size + 1));

        for(Biopolymer::const_iterator j = j_lower_bound;
            j <= j_upper_bound;
            ++j)
        {
          int I = i - m1.begin();
          int J = j - m2.begin();
          assert((0 <= I) && (I < m1.size()));
          assert((0 <= J) && (J < m2.size()));
          if (aaMakeMatch[I][J] == true)
            cerr << //I put this message here because I was curious if this
                    //ever occurs.  -Andrew 8/26/1999
              "Note: two helices/sheets or other intervals appear to overlapp" 
                 << endl;

          aaMakeMatch[I][J] = true;

          #if DBG_INTERVAL_LISTS
          ++block_match_counter;
          #endif
        } //j loops over m2 (second molecule)
      } //i loops over m1 (first molecule)

    } // final else clause for "if (start1 > end1-settings.fm_fragment_size+1)"
    #if DBG_INTERVAL_LISTS
    DEBUG_MSG(DBG_INTERVAL_LISTS, block_match_counter
              << " pairs of blocks matched in last interval pair.");
    #endif
  } //for(long p = 0; p < settings.fm_interval_pair_list.size(); ++p)


  // *** Now that we have filled in the aaMakeMatch[][] array,
  // *** we know which pairs of fragments to match.
  // *** Time to go ahead and match them.

  //(First, the "superimpose" class needs to be initialized.)
  int max_num_matches = MIN(m1.size(), m2.size());
  assert(max_num_matches > 0);
  Superimpose superimpose(max_num_matches
                                 *
                          Biopolymer::Residue::NumBackboneAtoms());

  for(Biopolymer::const_iterator i = m1.begin(); i != m1.end(); ++i)
  {
    for(Biopolymer::const_iterator j = m2.begin(); j != m2.end(); ++j)
    {
      int I = i - m1.begin();
      int J = j - m2.begin();
      assert((0 <= I) && (I < m1.size()));
      assert((0 <= J) && (J < m2.size()));
      if (aaMakeMatch[I][J] == true)
      {
        OrGenEntry orientation;
        Real fragment_rmsd 

          = superimpose.FindTransformMinRMSD(

               settings.fm_fragment_size 
                         *
               Biopolymer::Residue::NumBackboneAtoms(),

               aaRotateUsingTheseCoords[0] +
                 (I * Biopolymer::Residue::NumBackboneAtoms()),

               aaRotateUsingTheseCoords[1] +
                 (J * Biopolymer::Residue::NumBackboneAtoms()),

               orientation.transform);

        //If fm_cutoff_rmsd has been set, we discard orientations whose
        //rmsd between matching fragments is too high.
        //(Note: fm_cutoff_rmsd is rarely used. -AJ 8/26/1999)
        if ((settings.fm_cutoff_rmsd == Settings::fm_NO_CUTOFF_RMSD) ||
            (fragment_rmsd <= settings.fm_cutoff_rmsd))
        {
          vO.push_back(orientation);
        }
      } //if (aaMakeMatch[I][J] == true)
    } //loop over molecule m2   "for(j = m2.begin(); j != m2.end(); ++j)"
  } //loop over molecule m1   "for(i = m1.begin(); i != m1.end(); ++i)"


  //Check for possible error:
  //     If the list of interval pairs was empty, or
  // if all the intervals were too small to accomodate
  // the desired fragment size, print an error and exit.
  if (vO.size() == 0)
      ERR("Error: You have requested to generate orientations by\n"
          "       fitting together blocks of "
          << settings.fm_fragment_size << " residues each\n"
          "       belonging to certain portions of each molecule\n"
          "       (ie. helices and sheets, perhaps)\n"
          "       However, there are no valid regions from which to\n"
          "       to draw the " << settings.fm_fragment_size << " residue blocks\n"
          "       that get chosen to fit togther.\n"
          "       Possible causes:\n"
          "         1) Most likely you have a PDB file that does not contain any\n"
          "       helix or sheet records, and you are using the default parameters.\n"
          "       (By default minrms uses the location of helices and sheets to\n"
          "        limit the number of orientations that are generated.)\n"
          "       You can circumvent this behavior using the \"-HS\"\n"
          "       command line option.\n"
          "         2) This could also happen if your argument to -fm was\n"
          "       larger than the number of selected residues in either\n"
          "       molecule.\n"
          "         3) This could also happen if you have not dissabled\n"
          "       helix and sheet matching (using the \"-HS\" option)\n"
          "       and the residues you have selected do not enclose any\n"
          "       helices or sheets.  Check to see that at least some of\n"
          "       the helices,sheets {and other intervals} you specified\n"
          "       lie in the range of the chain(s) {and other intervals}\n"
          "       that you have selected from each molecule.\n"
          "\n"
          "           (If none of these explanations are correct, please\n"
          "            contact the developer and please include the\n"
          "            input files that generated this error.  Thanks!)\n");


  //Time to clean up.
  for (int I = 0; I < m1.size(); ++I)
    delete [] aaMakeMatch[I];
  delete [] aaMakeMatch;
} // OrientationGenerator::GenerateOrientationsUsingFragmentMatching()





void OrientationGenerator::
GenerateOrientationsReadingMatrixFiles()
{
  for (vector<string>::const_iterator pstr = settings.input_matrix_files.begin();
       pstr != settings.input_matrix_files.end();
       ++pstr)
  {
    OrGenEntry orientation;
    ReadTransFile(pstr->c_str(),
                  orientation.transform);
    vO.push_back(orientation);
  }
} //void OrientationGenerator::GenerateOrientationsReadingMatrixfiles()





void OrientationGenerator::GenerateOrientationsReadingMSFfiles()
{
  assert(settings.vMsfInputLabels.size() == 2);
  assert((settings.vMsfInputLabels[0].size() != 0) &&
         (settings.vMsfInputLabels[1].size() != 0));

  for (vector<string>::const_iterator pstr = settings.input_msf_files.begin();
       pstr != settings.input_msf_files.end();
       ++pstr)
  {
    PairwiseAlignment alignment;

    alignment.ImportMSF(*pstr,
                        settings.vMsfInputLabels[0],
                        settings.vMsfInputLabels[1],
                        settings.aMol[0],
                        settings.aMol[1],
                        settings.msf_use_lower_case);

                                                            
    //Note: I have assumed here that the MSF files will contain all the
    //residues in the entire structure (not limited to
    //residues in the current chain).
    //That's why we import the alignment using molecules
    //aMol[0], and aMol[1] (instead of aMol_f[0] and aMol_f[1]).
    //The alignment read in now is an alignment between aMol[0], and aMol[1].
    //
    //However the structural data we are using in our calculations is
    //from aMol_f[0] and aMol_f[1] which are a subset of aMol[0] and aMol[1].
    //So I re-cast this alignment in terms of aMol_f[0] and aMol_f[1].
    //by invoking the following function:
    TranslateAlignmentBetweenMolecules(alignment,
                                       alignment,
                                       settings.aMol[0],
                                       settings.aMol[1],
                                       settings.aMol_f[0],
                                       settings.aMol_f[1]);

    //Now calculate the orientation that
    //minimzes the RMSD of this alignment.
    OrGenEntry orientation;
    alignment.CalcMinRMSD(settings.aMol_f[0],
                          settings.aMol_f[1],
                          orientation.transform);

    vO.push_back(orientation);
  }
} //void OrientationGenerator::GenerateOrientationsReadingMSFfiles()


