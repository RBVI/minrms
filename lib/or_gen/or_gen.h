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

#ifndef _OR_GEN_H
#define _OR_GEN_H

/*
  "or_gen.h"
  This file declares a class called OrientationGenerator which handles
  all the details of how to rotate and translate the two molecules
  when performing the structure alignment.
  This class generates a large sequence of relative orientations
  between the two molecules.  Each "orientation" is a 3x4 matrix
  describing a rotation and translation to be applied to the coordinates
  of the atoms belonging to the second molecule.  (The first molecule
  is not moved.)

         A note on generating orientations "on-the-fly":

     In this implementation, I precompute all the
  orientations and store them in a vector.  We could,
  calculate orientations on the fly, instead
  (by overloading the derefence operator * for iterators, for example).
  This would save a LOT of space because the array can be huge.
  I don't do this because I wanted to be able to provide
  random access into the sequence of orientations. (The [] operator)
  (I wanted be able to find the i'th orientation in a hurry.)
  Because we use such complicated rules for choosing our orientations
  (helix/sheet matching, windowing, read from MSF-files),
  it's too hard to figure out how to generate the i'th orientation.
  In order to calculate solutions on the fly, we'd have to either:
  1)sacrifice random-access, or 2)use a simpler orientation generation
  method.
 */

#include <vect_3d.h>
#include <biopolymer.h>
#include <intervals.h>
#include <superimpose.h>
#include <pair_align_settings.h>


namespace minrms {


struct OrGenEntry
{
  //To have a container class of fixed-size C arrays,
  //I need to wrap each one in a struct/class.
  Matrix3x4 transform;
}; //struct OrGenEntry


class OrientationGenerator
{
  //(The next three functions are called by the constructor)
  void PrecomputeAllOrientations(); 
  void GenerateOrientationsUsingFragmentMatching();
  void GenerateOrientationsReadingMatrixFiles();
  void GenerateOrientationsReadingMSFfiles();

public:

  typedef vector<OrGenEntry>::const_iterator const_iterator;
  typedef vector<OrGenEntry>::size_type size_type;

  OrGenEntry const& operator [] (long i) const
  { return vO[i]; }

  const_iterator begin() const
  { return vO.begin(); }

  const_iterator end() const
  { return vO.end(); }

  OrGenEntry const& front() const
  { return vO.front(); }

  OrGenEntry const& back() const
  { return vO.back(); }

  size_type size() const
  { return vO.size(); }

  //SumSqdDistBetweenOrientations() returns the sum-squared distance
  //between the 2nd molecule rotated and translated according to the
  //two orientations indicated by o1 and o2.
  Real SumSqdDistBetweenOrientations(ConstMatrix3x4 o1, ConstMatrix3x4 o2) const;

  //The following version of SumSqdDistBetweenOrientations()
  //takes indices instead of of Matrix3x4s.
  //(The indices refer to orientations generated by this OrientationGenerator.
  // Note: The other version of this function will work on any orientations.)

  Real SumSqdDistBetweenOrientations(long o1_id, long o2_id) const
  {
    return SumSqdDistBetweenOrientations((*this)[o1_id].transform,
                                         (*this)[o2_id].transform);
  }

  Real RMSDbetweenOrientations(ConstMatrix3x4 o1, ConstMatrix3x4 o2) const
  {
    Real sum_sqd = SumSqdDistBetweenOrientations(o1,
                                                 o2);
    return sqrt(sum_sqd /
                (settings.aMol_f[1].size()
                 *
                 Biopolymer::Residue::NumBackboneAtoms()));
  }

  Real RMSDbetweenOrientations(long o1_id, long o2_id) const
  {
    return RMSDbetweenOrientations(vO[o1_id].transform,
                                   vO[o2_id].transform);
  }



  //Settings contains all input parameters and data
  struct Settings:public virtual PairAlignSettings
  {
    //Use "fragment matching" to orient the molecules?
    //If so, how many residues should be in each fragment?
    int      fm_fragment_size;

    //To dissable fragment matching set
    //fm_fragment_size = fm_DISSABLE_FRAGMENT_MATCHING
    static const int fm_DISSABLE_FRAGMENT_MATCHING = -1;

    //If fragment matching is enabled,
    //to throw away orientations whose fragment's RMSD is too large, set
    Real     fm_cutoff_rmsd;
    //to use all orientations, set fm_cutoff_rmsd = fm_NO_CUTOFF_RMSD
    static const Real fm_NO_CUTOFF_RMSD;

    // The fm_interval_pair_list member is a list of pairs-of-intervals from
    // either sequence that will be used to limit the set of fragments
    // that get matched together during the process of "fragment matching."
    // As such, each entry in the vector-of-vectors-of-intervals,
    // fm_interval_pair_list[i], should only contain only 2 intervals.
    // (This may be enforced with assert().)
    // Suppose:
    //  fm_interval_pair_list[i] = {[startA, endA], [startB, endB]}
    // These two intervals indicate that fragments which lie between
    // residues: startA and endA in the first protein are superimposed
    // with fragments which lie between residues: startB to endB from
    // the second protein.
    //
    // The orientations generated by superimposing all the fragments from
    // [startA, endA] with all the fragments from [startB, endB].
    // This is done for each pair of intervals in fm_interval_pair_list[].
    //
    // (This is how "helix and sheet matching" is implemented.)
    vector<vector<ClosedInterval> >   fm_interval_pair_list;

    //The following argument is commented out because it
    //redundantly specifies information stored in "fm_interval_pair_list".
    //bool     fm_helix_sheet; //use the helix/sheet matching optimization?

    //Consider additional orientations loaded directly from various files.
    vector<string> input_matrix_files; //names of files to be read
    
    //Consider additional orientations from superimposing
    //residues matched in MSF-files?
    vector<string> input_msf_files;  //The names of the files being read.
    vector<string> vMsfInputLabels;  //The labels that will be used to identify
                                     //which structure corresponds to which
                                     //line of each MSF-file read
                                     //by the program.
                                     //(This may not be useful.)

    Settings();

    #ifdef COMPILER_BUG2
    //copy constructor:
    //(The alpha compiler produces buggy
    // binaries unless I define one explicitly)
    Settings(Settings const& s);
    #endif //#ifdef COMPILER_BUG2

  }; //class OrientationGenerator::Settings


  //The constructor does all the precomputation necessary before invoking
  //begin(), end(), the [] operator, or any other public functions.
  OrientationGenerator(Settings const& settings);

  ~OrientationGenerator();


private:

  vector<OrGenEntry> vO;
  Settings settings;

  OrientationGenerator(const OrientationGenerator &);        //dissable copying
  OrientationGenerator& operator =(const OrientationGenerator &); // "  "
  Vect3 *aaRotateUsingTheseCoords[2];//Stores the atom positions being rotated.
  Matrix3x4 aPrincipleAxisAndCentroid[2];//Used for computing SumSqdDistBetOr()
  Vect3     aPrincipleAxisEigenvalues[2];// "    "      "           "

}; //class OrientationGenerator

} //namespace minrms





#endif //#ifndef _OR_GEN_H

