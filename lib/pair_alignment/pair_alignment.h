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

#ifndef _PAIR_ALIGNMENT_H
#define _PAIR_ALIGNMENT_H

// Overview:
//
//   PairwiseAlignment data structure stores alignments between two
//linear-molecules (which could be either proteins or DNA/RNA, for example)
//and has member functions to handle common things like:
// 1) parsing and writing MSF-files and FSSP-files,
// 2) calculating the minimal RMSD of an alignment, (and calculating the
//    the relative translation and rotation between the two molecules
//    necessary to achive this minimal RMSD.  See "Superimpose".)
// 3) Comparing two different alignments between the same two molecules
//    to see how much each alignment differs.
//    (This is used to compare the contents of two MSF-files, for example.)
// 4) Performing statistical analysis on an alignment.
//    (calculating the chance of an alignment being a false-positive.
//     see LevittGerstein98)
//
// ...and less common things like:
// 5) generating graphics-object files that show the pairwise matching
//    (with the relative-orientation chosen to minimize RMSD)
// 6) generating histogram (from a large group of PairwiseAlignments),
//    showing how frequently each pair of residues was matched, and printing
//    out a file.
// 7) the distance between the farthest pair of residues in an alignment.
//The alignment can be sequential or non-sequential, (but certain
// member-functions are not available for non-sequential alignments.)
//
//  Implementation:
//
//   Each "PairwiseAlignment" is just an array of
//  "PairwiseMatch"es which simply stores a pair of short-integers
//   which indicate the indexing into two different Sequences.
//  (Note: Indexing starts at zero, that is, these shorts were intended
//         to have values from  0 to (length_of_sequence-1).)


#include <cctype>
#include <cassert>
using namespace std;

#include <global_utils.h>
#include <vect_3d.h>
#include <biopolymer.h>
#include <superimpose.h>

namespace minrms {

class PairwiseMatch
{
  short residues_in_match[2];
public:
  inline PairwiseMatch() { residues_in_match[0] = -1;
                           residues_in_match[1] = -1; }

  inline PairwiseMatch(short i, short j) {
    residues_in_match[0] = i;
    residues_in_match[1] = j;
  }

  short & operator [](int which_seq) {
    assert((0<=which_seq) && (which_seq<2));
    return residues_in_match[which_seq];
  }

  short const& operator [](int which_seq) const {
    assert((0<=which_seq) && (which_seq<2));
    return residues_in_match[which_seq];
  }

  bool operator == (PairwiseMatch const& m) const
  { return (((*this).residues_in_match[0] == m.residues_in_match[0]) &&
            ((*this).residues_in_match[1] == m.residues_in_match[1])); }

  bool operator != (PairwiseMatch const& m) const
  { return ! ((*this) == m); } //some compilers don't do this for you 
}; 





class PairwiseAlignment
{
  int  num_matches; //number of matches in this alignment
  int  space_available; //Size of the array allocated to store the matches.
                        //The maximum number of matches that can be stored.
                        // (Maybee later I'll use vectors instead so
                        //  there is no maximum, but for now this will do.)
  PairwiseMatch *aMatches;    //the array that stores the matches.

  //This function is only intended to be used for reading MSF-files.
  inline bool IsGap(char c, bool check_for_lower_case = false) {
    return ((c==GAP_CHAR)
            //|| (c=='.') || (c=='-')
            ||
            (check_for_lower_case && islower(c)));
  }

  void DigestMSFContents(vector<string> const& vContent,
                         vector<string> const& vLabels,
                         vector<string> const& vSequences,
                         bool check_for_lower_case = false);

public:
  PairwiseAlignment() {space_available = 0; num_matches = 0; aMatches = NULL;}
  PairwiseAlignment(const PairwiseAlignment &source);
  PairwiseAlignment(int n) { space_available = 0; aMatches = NULL; resize(n); }
  ~PairwiseAlignment() { if (aMatches) delete [] aMatches; }

  typedef PairwiseMatch *iterator;
  typedef PairwiseMatch const *const_iterator;
  iterator begin(){ return aMatches; }
  iterator end(){ return aMatches+num_matches; }
  const_iterator begin() const { return aMatches; }
  const_iterator end() const { return aMatches+num_matches; }



  //resize() is used to change the number of matches in the alignment.
  //Note, resize() destroys the former contents of the alignment (if any).
  void resize(int n);

  //If you know in advance, a reasonable upper bound on the size of the
  //alignment, call reserve() to reserve that memory in advance.
  void reserve(int n);

  //data access and array-subscripting:
  inline int NumMatches() const {return num_matches;}
  inline int size() const {return num_matches;}

  inline PairwiseMatch & operator [](int i)
  { assert((0 <= i) && (i<num_matches)); return aMatches[i];}

  inline PairwiseMatch const& operator [](int i) const
  { assert((0 <= i) && (i<num_matches)); return aMatches[i];}

  PairwiseAlignment& operator =(const PairwiseAlignment &source)//copy operator
  { 
    if (num_matches != source.NumMatches())
      resize(source.NumMatches());
    for(int n=0; n<num_matches; ++n)
      aMatches[n] = source[n];
    return *this;
  }

  bool operator ==(const PairwiseAlignment &source) const //equality operator
  { 
    if (num_matches != source.NumMatches())
      return false;
    for(int n=0; n<num_matches; ++n)
      if (aMatches[n] != source[n]) return false;
    return true;
  }

  bool operator !=(const PairwiseAlignment &source) const
  { return (! (*this == source)); }

  //IsMonotonic() returns whether or not the alignment is sequential.
  //(a.k.a. a "topological" alignment)
  bool IsMonotonic() const;


  //CalcMinRMSD() returns the minimum RMSD between designated atoms belonging
  //to residues which are matched in this alignment.
  //You must pass it two "Biopolymer" arguments (m1 and m2)
  //(This data-structure is defined in "linear_molecule.h")
  //which store the three dimensional positions of the atoms from
  //either molecule.  Explicitely: the RMSD is minimized between
  //all the "BackboneAtoms" from all matched segments.
  //(See "linear_molecule.h")  This is user-configurable,
  //but for proteins, this is typically only the " CA " atom.
  //Afterwards, the rotation and translation that must be applied to the second
  //molecule (m2) in order to minimize the RMSD is saved in a 3x4 matrix
  //called "optimal_transform".
  //
  // (Performance note:
  //  Ordinarily, this function will allocate and deallocate
  //  temporary arrays to store the coordinates of
  //  all the atoms used in the RMSD minimization.
  //  For those who care, instead the user can pass a pointer to
  //  a pre-allocated "Superimpose" variable, if one exists, as well
  //  as passing pre-allocated storage space ("aReservedSpace1/2").
  //  This is useful if you intend to call this function more than once.)

  Real CalcMinRMSD(Biopolymer const& m1,
                   Biopolymer const& m2,
                   Matrix3x4 optimal_transform,
                   Superimpose *pSuperimpose = NULL,
                   Vect3 *aReservedSpace1 = NULL,
                   Vect3 *aReservedSpace2 = NULL
                   ) const;

  static char GAP_CHAR;     //The character used for indicating a gap: (.)
                            //in the MSF files that are read/written.

  static char GAP_CHAR_OFFSET;   //(~) GCG distinguishes between gaps in an
                                 //    alignment, and offsets in the sequences.
                                 //. is used for gaps (in the middle
                                 //  of the alignment)
                                 //~ is used to indicate offsets between
                                 //  the starting/ending points of the
                                 //  two sequences (gaps that appear at the
                                 //  beginning or end of the alignment).  This
                                 //  character (~) never appears in an MSF file
                                 //  generated by this code, but this code will
                                 //  not object to reading an MSF file
                                 //  with (~) in it.
  
  static char aColorNames[][32]; //colors used for the different sequences
  static char unit_name[32]; //Optional. Stores the name of the units for
                             //describing physical distance. (ie. "Angstroms")


  //ImportMSF() is an input function.
  //It reads a standard MSF-file containing at least two sequences.
  //-The "labelA/B" arguments indicate the two sequences you wish to
  // find the alignment between. (MSF files contain a label for each sequence.)
  //-The optional "sequenceA/B" arguments contain strings of characters
  // that store the original sequence content of the two molecules.
  // The amino-acid symbols read from the msf-file are checked against
  // the contents of these two arguments to check that the msf-file has
  // the correct contents.  If they are set to NULL, they are ignored.
  //-The "check_for_lower_case" argument determines whether or not
  // it considers two lower-case letters directly on top of eachother
  // a gap(true) or a match(false).
  void ImportMSF(string filename, //name of msf file that will be created
                 string labelA, //indicates the two sequences you
                 string labelB, //wish to find the alignment between
                 string sequenceA="",//Sequence content. To ensure the
                 string sequenceB="",//sequence in msf-file is correct.
                                            //This is optional.  If you don't
                                            //want to check, pass "".
                 bool check_for_lower_case=false//interprets the msf-file
                                                //using Conrad's
                                                //lower-case-means-no-match
                                                //convention.
                 );


  //The following function is the same, but it reads from any istream.
  void ImportMSF(istream& msf_file,
                 string labelA,
                 string labelB,
                 string sequenceA="",
                 string sequenceB="",
                 bool check_for_lower_case=false);

  //The following alternate version of ImportMSF, accepts two
  //"Biopolymer"s, m1 and m2, which replace the two "char *"s
  //sequenceA, and sequenceB in the version above.
  //In either case these arguments are used to check that the residue-letters
  //belonging to the MSF-file are consistent with the 3-letter residue names
  //in the two molecules in m1 and m2.
  void ImportMSF(string filename,
                 string labelA,
                 string labelB,
                 Biopolymer const& m1, //m1 and m2 contain sequence content to insure
                 Biopolymer const& m2, //the characters in the msf-file are correct.
                 bool check_for_lower_case=false);

  //The following function is the same, but it reads from any istream.
  void ImportMSF(istream& msf_file,
                 string labelA,
                 string labelB,
                 Biopolymer const& m1,
                 Biopolymer const& m2,
                 bool check_for_lower_case=false);

  //-Reads in FSSP files.
  //-Generates an error if the fssp file is not sequential (non-monotonic).
  // (Does not work with non-monotonic FSSP files!
  //  I just didn't want to be bothered adding this feature since we 
  //  haven't needed it {yet}.)
  //-The "protein_labelA/B" arguments specify the name of the two protiens
  // as they appear in the "## EQUIVALENCES:" section.  They must be specified
  // in the same order that they appear in, from left to right, on each line
  // the "EQUIVALENCES" section of the fssp file containing the match between
  // those two proteins.
  //-The m1, and m2 arguments store, (in addition to other, less useful things)
  // the PDB-codes for the residues indicated in the FSSP-file.
  void ImportFSSP(string filename,
                  Biopolymer const& m1,
                  Biopolymer const& m2,
                  string labelA,
                  string labelB);

  //The following function is the same, but it reads from any istream.
  void ImportFSSP(istream& fssp_file,
                  Biopolymer const& m1,
                  Biopolymer const& m2,
                  string labelA,
                  string abelB);

  //ExportMSF generates an msf-file for these alignment.
  void ExportMSF(string filename, //Name of the file you want to create.
                                       //(Note: no ".msf" suffix is added.)
             string seqContentsA, //The contents of the two sequences
             string seqContentsB, //the msf-file refers to.
                 
             string labelA,//These strings will be used to label
             string labelB,//each molecule in the msf-file output.

             string comments="", //Puts your comments before the ".."
                                       //terminator.  To dissable, pass ""

             bool  show_markers=true,  //Draws number-line on top

             bool  compress_using_lower_case=false //Conrad's trick
                                          //of compressing the files by using
                                          //upper&lower case to distinguish
                                          //whether a match has occured.
                                          //(not standard, to my knowledge)
             ) const;

  //The following version of ExportMSF, is similar, except that
  //instead of accepting the "labelA" and "labelB"
  //arguments separately, it takes a single argument, a vector-of-strings
  //containing the two protein labels.  As such, this vector should
  //store two strings (this is enforced with assert()).
  //  Likewise, the same is true of it's "seqContentsA/B" arguments,
  //In this new version, these two arguments have been consolidated
  //into a single vector-of-strings "vSequences".
  void ExportMSF(string filename,
                 vector<string> const& vSequences,
                 vector<string> const& vLabels,
                 string comments = "",
                 bool  show_markers = true,
                 bool  compress_using_lower_case = false) const;



  //The following alternate version of ExportMSF, accepts two
  //"Biopolymer"s, m1 and m2, which contain the
  //sequence data stored in the sequenceA and sequenceB arguments
  //in the version above.
  void ExportMSF(string filename,
                 Biopolymer const& m1,
                 Biopolymer const& m2,
                 string labelA,
                 string labelB,
                 string comments = "",
                 bool  show_markers = true,
                 bool  compress_using_lower_case = false) const;


  // ---------------
  //(The next three functiosn are mostly for internal use, but I leave them
  // public:)
  //The following version of ExportMSF, is similar, except that
  //instead of accepting the "labelA" and "labelB"
  //arguments separately, it takes a single argument, a vector-of-strings
  //containing the two protein labels.  As such, this vector should
  //store two strings (this is enforced with assert()).
  //  Likewise, the same is true of it's linear molecule arguments,
  //"m1" and "m2".  In this new version, these two arguments have been
  //consolidated into a single vector-of-Biopolymers argument, "vMol".
  void ExportMSF(string filename,
                 vector<Biopolymer> const& vMol,
                 vector<string> const& vLabels,
                 string comments = "",
                 bool  show_markers = true,
                 bool  compress_using_lower_case = false) const;
  //The following version of ImportMSF uses vectors-of-strings instead
  //of separate string arguments.
  void ImportMSF(istream &msf_file,
                 vector<string> const& vLabels,
                 vector<string> const& vSequences,
                 bool check_for_lower_case);
  //Likewise, the following version of ImportFSSP uses vectors as well.
  void ImportFSSP(istream &fssp_file,
                  vector<Biopolymer> const& vMol,
                  vector<string> const& vLabels);
  // ---------------


  //ExportGFX() creates a "graphics object" file showing the two proteins
  //rotated and superimposed against eachother, as well as red line-segments
  //showing connections between corresponding residue-centers from
  //either protein according to the alignment.
  //The orientations of the two molecules must be set in advance.
  //(This function will not choose an orientation to makes it "look nice".)
  //In the graphics file, the orientation of the first protein is fixed.
  //You must manually specify the rotation and translation to apply to
  //the coordinates of the second protein using the "transform" argument.
  //There are several customizations you can make (discussed below)
  //Other arguments are as follows:
  void ExportGFX(Biopolymer const& m1, //Stores the coordinates of the
                 Biopolymer const& m2, //atoms from either molecule.

           string filename,     //Name of the graphics object file
                                      //you want to create.
                                      //(Note: no ".gfx" suffix is added.)

           string caption,       //Displays a caption in the file.
                                      //If none desired, pass NULL or "".

           bool show_residue_markers = false, //draw a marker at every
                                      //residue of both sequences

           bool connect_the_dots = false, //draw lines along the backbone
                                          //connecting the residues in
                                          //each sequence.
                 

           //If "use_marker_if_dist_less_than" is specified, matches that are
           //too close together are indicated by some kind of 3-D marker 
           //(either a tetrahedron or a sphere depending on version)
           //as well as a line.  This is because a very short line-segment
           //would probably be very difficult to see.  If you don't want this
           //feature, pass DONT_USE_MARKERS_FOR_SHORT_LENGTHS (the default).
           Real use_marker_if_dist_less_than =
                                     DONT_USE_MARKERS_FOR_SHORT_LENGTHS,

           //  The next argument alows you to mark the furthest
           //  pair of residues from eachother (as measured between
           //  their alpha carbons) in orange (instead of red).
           bool mark_worst_match = 0,

           //------- You probably want to ignore (pass 0) --------
           //-------  the next three arguments (default). --------
           //Use the next 3 arguments if you want to highlight certain
           //interval from both sequences in a particular color.
           unsigned int show_subset_size = 0,
           int subset_start_A = 0,
           int subset_start_B = 0) const;

  //---- Related data members ---
  //You can pass these values to ExportGFX() as arguments.
  //They have special meanings.
  static const Real DONT_USE_MARKERS_FOR_SHORT_LENGTHS; //see above
  static const Real MIDAS_MARKER_SIZE; //Approximate size (in Angstroms)
                                        //of the markers used in midas.
                                        //I think this is 0.4 (see .cc file)


  //FindFurthestAtomPair() finds the pair of residues matched in the alignment
  //with the furthest distance between representative atoms.
  //(By default, the positions of the " CA" atoms are used.)
  //Residues not containing an atom of this type are ignored.
  //If this alignment contains no matches, or none of the matches
  //involve residues containing the representative atoms desired,
  //then FindFurthestAtomPair() prints a warning message, and returns -1.0.
  //  This function assumes the two structures, m1 and m2 have allready
  //been physically superimposed in the way you want them.
  //(It doesn not automatically superimpose them.)
  Real
  FindFurthestAtomPair(Biopolymer::const_iterator& worst_i,
                       Biopolymer::const_iterator& worst_j,
                       Biopolymer const& m1,
                       Biopolymer const& m2,
                       string representative_atom_name = string(" CA")) const;


  //FindFurthestResPair() is just like FindFurthestAtomPair(),
  //except it finds pair of residues with the highest average (RMS)
  //distance between ALL of the BackboneAtoms in either residue,
  //as opposed to the pair of residues with the highest distance
  //between just one atom (which is predicted by FindFurthestAtomPair().)
  //This average becomes just a simple distance if there is only
  //one "BackboneAtom()" in each residue.
  //In this version, it is assumed that all residues refered to in this
  //alignment contain all of the BackboneAtoms they need.
  //(In DEBUG mode, this is enforced with assert().)
  //  Again, this function assumes the two structures, m1 and m2 have allready
  //been physically superimposed in the way you want them.
  //(It doesn not automatically superimpose them.)
  Real
  FindFurthestResPair(Biopolymer::const_iterator& worst_i,
                      Biopolymer::const_iterator& worst_j,
                      Biopolymer const& m1,
                      Biopolymer const& m2) const;

#ifdef CREATE_PAIR_HISTOGRAM
  //Accumulates a histogram of matches from the alignments with a number
  //of matches made, n, in the range from [n_start, n_end] (including
  //the endpoints.  indexing starts at 1, that is, n=1 indicates one match)
  void AccumPairHistogram(long **aaHistogram) const;
#endif
}; // class PairwiseAlignment 



// ***********************************************************
// ***** Utilities for evaulating the quality/closeness  *****
// *****    of a single pairwise structure alignment     *****
// ***********************************************************

Real SumSqdDist(PairwiseAlignment const& a,
                Biopolymer const& m1,
                Biopolymer const& m2);


Real RMSD(PairwiseAlignment const& a,
          Biopolymer const& m1,
          Biopolymer const& m2);


namespace LevittGerstein98
{
  //    Returns the score of this structural alignment, using the simplest
  // form of Levitt & Gerstein's criteria with S_str = M / (1 - (d_ij/d_o)^2).
  // (See "A unified statistical framework for sequence comparison
  //       and structure comparison", Proc Sci May 1998, Vol 95, pp5913-5920)
  // Where M = 20.0,
  //       d_o = 5.0,
  //       the Opening gap-parameter is M/2, and
  //       the Extension gap-parameter is 0.
  //    The answer depends on the relative orientation between the two 
  // molecules.  It will not try to optimize the rotation to maximise S_str.
  //    No weighting factors are used (either side-chain-weighting or
  // exposure-weighting) to sway the outcome.
  Real Sstr_Score(PairwiseAlignment const& a,
                  Biopolymer const& m1,
                  Biopolymer const& m2);

  //    Simply returns the RMSD between the C-alpha carbons only
  // (whose names must be " CA") of the proteins m1 and m2,
  // at their present orientation.  (It does _not_ try to optimize
  // the relative rotate/translate between the two molecules to minimize
  // this RMSD.
  // This function is basically identical to CalcRMSD(),
  // except it uses only the alpha-carbon atoms to calculate the RMSD.
  // Unlike RMSD(), it ignores the BackboneAtoms in the molecules.
  // (Definition: "BackboneAtoms" are usually atoms specified at run-time
  //               in a "atoms_used.txt" file. See <load_pdb.h>)
  Real RMSD_score(PairwiseAlignment const& a,
                  Biopolymer const& m1,
                  Biopolymer const& m2);

  //    The following function attempts to find the "P_str", the probability
  // that two unrelated molecules could get as high a similarity score "S_str",
  // as do the molecules "m1", and "m2" according to alignment "a".
  // That is, it tries to predict the probability that this alignment is
  // a "false positive."
  // This uses Levitt & Gerstein's "S_str" similarity score.
  //    The answer depends on the relative orientation between the two 
  // molecules.  It will not try to optimize the rotation.
  //    Because the answer can be very small, the optional
  // ZextremeValueParameter, is returned to the caller.
  // P_str and Z are related by P_str = 1 - exp(-exp(-Z)).
  // (See "A unified statistical framework for sequence comparison
  //       and structure comparison", Proc Sci May 1998, Vol 95, pp5913-5920)

  Real Sstr_Probability(PairwiseAlignment const& al,
                        Biopolymer const& m1,
                        Biopolymer const& m2,
                        Real *p_str = NULL, //returns the score
                        Real *p_mu_str = NULL, //returns the expected score for this N
                        Real *p_delta_str = NULL, //returns the expected variance in score for this N
                        Real *pZ = NULL);

} //namespace LevittGerstein98




// ***************************************************************
// *****  The following 9 functions are different metrics    *****
// *****  for measuring the similarity between two different *****
// *****  alignments between the same two molecules.         *****
// ***** -Some of them measure the similarity in the actual  *****
// *****  sets of residues matched in the two alignments.    *****
// ***** -Others measure the similarity in the 3-dimensional *****
// *****  supperposition between the two molecules implied   *****
// *****  by the two different alignments.                   *****
// *****  We used these functions to compare the alignments  *****
// *****  produced by MINRMS, to alignments generated by     *****
// *****  other algorithms.                                  *****
// ***************************************************************

//    The NumMatchesInBoth() function compares two pairwise alignments.
// It simply returns the number of identical pairs-of-residues
// which are found in both alignments.
// Note:  Each alignment should be between the same two molecules.
//        (They do not have to have the same number of matches, however.)
// (Running time: quadratic in the number of matches in the alignments.)
int NumMatchesInBoth(PairwiseAlignment const& a1,
                     PairwiseAlignment const& a2);

int NumMatchesInEither(PairwiseAlignment const& a1,
                       PairwiseAlignment const& a2);

//     NumResiduesInBoth() counts the number of residues from
// the first sequence that were matched in both alignments a1, and a2.
// and stores it in the "num_bothA" argument.
// It also counts the number residues in the second sequence from both
// alignments, and stores this in the "num_bothB" argument.
// (Running Time: linear in the lengths of the sequences)
void NumResiduesInBoth(PairwiseAlignment const& a1,
                       PairwiseAlignment const& a2,
                       int& num_bothA,
                       int& num_bothB);

void NumResiduesInEither(PairwiseAlignment const& a1,
                         PairwiseAlignment const& a2,
                         int& num_eitherA,
                         int& num_eitherB);

//    The version below called "NumMatchesInBothAssumingMonotonicity()"
// is identical to "NumMatchesInBoth()", except it also requires
// that the two alignments must be monotonic.
// (That is IsMonotonic(a1) and IsMonotnic(a2) must both be true.)
// Note: that alignments generated by all common dynamic-programming
//       algorithms have this property.
// (Running time: linear in number of matches in the alignments.
//                the original version was quadratic.)
int NumMatchesInBothAssumingMonotonicity(PairwiseAlignment const& a1,
                                         PairwiseAlignment const& a2);

int NumMatchesInEitherAssumingMonotonicity(PairwiseAlignment const& a1,
                                           PairwiseAlignment const& a2);


// The next three functions are difficult to explain.  Here goes:

// Compare3dAll() returns the average distance (in angstroms) between
// two physical superpositions of the same two structures
// implied by alignments: a1 and a2.
// More precisely, with structure mA held fixed,
// it returns the RMS-displacement between positions of the atoms
// in structure mB, after mB was superposed twice with structure mA
// in order to minimize RMSD between the atoms matched 
// in alignments a1, and a2 respectively.
// (Minor detail: Actually, only the positions of the "backbone-atoms" in mB,
//                are considered to calculate the RMS-displacement
//                between the two positions of mB.  Not all the atoms are used.
//                See <linear_molecule.h> for a description of BackboneAtoms.)
Real Compare3dAll(PairwiseAlignment const& a1,
                  PairwiseAlignment const& a2,
                  Biopolymer const& mA,
                  Biopolymer const& mB);

// Compare3dBoth() is just like Compare3dAll() except
// instead of RMS-displacement between all the atoms in structure mB
// after optimal superposition with mA in alilgnments a1 and a2,
// it calculates the RMS-displacement only for the atoms that
// were matched in BOTH alignments, a1 and a2.
Real Compare3dBoth(PairwiseAlignment const& a1,
                   PairwiseAlignment const& a2,
                   Biopolymer const& mA,
                   Biopolymer const& mB);

// Compare3dEither() is just like Compare3dBoth() except
// instead of RMS-displacement between all the atoms in structure mB
// after optimal superposition with mA in alilgnments a1 and a2,
// it calculates the RMS-displacement only for the atoms that
// were matched in EITHER alignment, a1 or a2.
Real Compare3dEither(PairwiseAlignment const& a1,
                     PairwiseAlignment const& a2,
                     Biopolymer const& mA,
                     Biopolymer const& mB);


// **************************************************
// ***** A silly function to extract an amino-  *****
// *****    acid sequences from an MSF-file.    *****
// **************************************************

//ReadMSFsequence() takes an msf_file, a label indicating which sequence
//should be used from the MSF-file, and a destination string (sequence).
//It then reads in all of the amino acid 1-character codes from the
//MSF file, and stores them in the "sequence" argument.
//(It returns whether or not a sequence with the desired label was
// was found in the MSF file.)
bool ReadMSFsequence(istream& msf_file,
                     string  label,
                     string& sequence);


// ********************************************************************
// ***** This next function is not likely to be useful to others. *****
// ********************************************************************

//Translates an alignment between molecules mOrig1 and mOrig2 (which is alOrig)
//to an alignment between similar molecules mNew1  and mNew2
//and stores the result in alignment alNew.
// Details:
//  For every matched pair of residues in alOrig, it looks up the same
//residues in the mNew1 and mNew molecules.  (Right now,
//residues are identified by their PDB chain,seqNum,and insertCode)
//It figures out where the residues with the same ID codes are located
//in the other pair of molecules (mNew1 and mNew2) and generates a
//new alignment for the new pair of molecules.
//An assertion-fail occurs if any of the residue pairs from the alignment
//between mOrig1 and mOrig2 are not found in the new molecules mNew1 and mNew2.
// Note:
//   It's okay if alOrig and alNew are the same.
void TranslateAlignmentBetweenMolecules(PairwiseAlignment const& alOrig,
                                        PairwiseAlignment & alNew,
                                        Biopolymer const& mOrig1,
                                        Biopolymer const& mOrig2,
                                        Biopolymer const& mNew1,
                                        Biopolymer const& mNew2);


} //namespace minrms


#endif  // #ifndef _PAIR_ALIGNMENT_H







