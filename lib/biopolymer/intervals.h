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


#ifndef _INTERVALS_H
#define _INTERVALS_H

#include <vector>
using namespace std;


#include <pdb++.h>
#include <global_utils.h>
#include <stl_utils.h>
#include "pdb_id.h"
#include "biopolymer.h"



#define DBG_INTERVAL_LISTS 0 //(This just turns off the printing of some
#define DBG_HELIX_SHEET 0    // messages used for debugging. Please ignore.)








namespace minrms 
{

// *******************************************************************
//   The ClosedInterval class was intended to be used to denote
//   intervals of residues from within a molecule.
// "start" and "end" represent the elements that identify the endpoints
// of the interval in question: [start, end].
//    (This subset includes all the members in between
//     start and end, including the endpoints start and end.)

struct ClosedInterval 
{
  PDBresID first;
  PDBresID last;
};


// *******************************************************************
// Filters out the residues in a linear molecule which do not
// belong to any of the intervals specified in the "filters"
// argument (which is of type "ListOfIntervals").
bool
ClipMoleculeToListOfIntervals(Biopolymer& m,
                              vector<ClosedInterval>& filters);

// *******************************************************************
//   BelongsToUnionOfIntervals() tests to see if the residue indicated by key,
//belongs to any of the intervals specified in the list of intervals, vI.
//In other words, this function returns if "key" belongs to the union of
//all the intervals in vI.  This takes time O( i ) where i is the number
//of intervals.
inline bool BelongsToUnionOfIntervals(PDBresID key,
                                      vector<ClosedInterval> &vI,
                                      Biopolymer const& m)
{
  for (long i=0; i < vI.size(); ++i)
  {
    if ( (! comes_after(vI[i].first, key, m))
         &&
         (! comes_after(key, vI[i].last, m)) )
      return true;
  }
  return false;
}



//The following function converts a HELIX, SHEET, or TURN,
//record from a PDB file into a ClosedInterval.
//It returns, the type of record that was read.
//If it is not a HELIX, SHEET, or TURN record, 
//then the argument "interval" remains unmodified.
//(Note: This function includes a hack to circumvent a problem
//       Greg's PDB library has with reading in old-fashioned
//       pdb files.  See below.)
PDB::RecordType
ConvertHelixSheetTurnRecordToInterval(char const *line,
                                      ClosedInterval& interval);


//The current version of PDB.h I am using does not understand the
//old-style helix,sheet,and turn-records.  If such a record is read,
//it will be read as type PDB::UNKOWN.
//The following function converts these kinds of records
//into modern PDB helix,sheet,and turn-records.
void
ConvertOldHelixSheetTurnRecordsToNewFormat(PDB const& source_record,
                                           PDB& dest_record);


// *******************************************************************
// ExtendInterval() takes a closed interval, and pads the interval
// on either side (equally on each side) if it contains fewer than
// "min_size" elements.
// Use:  This is useful for the "helix and sheet matching" feature
//       of minrms where a fixed size block of residues from a helix
//       in one molecule must be compared to a fixed size block from
//       another helix in another molecule.  This function is to address
//       the problem that arises if either helix is too short to
//       accomodate a block of the requested size.
// Details:
//    If it is impossible to pad the interval
// on one side because the interval lies too close to the start or end
// of the molecule, then it will proceed to add extra padding to
// the other side of the interval.  If even this is impossible, 
// because the start and end of the molecule are too close together
// (ie. the molecule contains fewer than "min_size" residues),
// then it will print an error message and exit.
// Return value:
//   ExtendInterval() returns true iff the interval was extended.

bool ExtendInterval(ClosedInterval&        interval,
                    Biopolymer const&  mol,
                    int                    min_size);




// *******************************************************************
//Deletes intervals from a vector of intervals "vI", if any of these
//intervals refer to residues which do not exist in the molecule "mol".
//The function returns whether or not any intervals were deleted.
//The last two arguments are mostly for debugging purposes:
//The optional "pFirstDeletedInterval" allows the caller to find out
//which was the first interval that was deleted, and the optional
//"pHowFarThroughList" parameter allows the caller to find out
//where, in the vector the offending interval resided.
//It returns whether or not any intervals were deleted.
bool DeleteInvalidIntervals(vector<ClosedInterval>& vI,
                            Biopolymer  &mol,
                            ClosedInterval* pFirstDeletedInterval = NULL,
                            int* pHowFarThroughList = NULL);



// *******************************************************************
//   ConvertHelixSheetRecordsToIntervalNtuples()
//is used to impliment minrms' helix-and-sheet-matching optimization.
//   The following function reads a list of n PDB files
//(whose names are stored in the "aaFilenames[]" argument)
//and returns a vector of n-tuples of intervals for ever possible
//combination of helices and sheets from each molecule.
//    (For "minrms", n=2, because as a pairwise protein
//     structural alignment tool, it only needs to consider
//     two proteins at once.)
//For each PDB-file, it finds all the helix intervals for each molecule
//by calling ParseHelixSheetTurnsRecords().
//Then, it generates n-tuples of different helices from each molecule.
//    (Again, in the case of "minrms", which considers only two molecules,
//     these n-tuples are just pairs of helices)
//Each "n-tuple" is a combination of helices from every molecule, 
//a helix from the first molecule, a helix from the second molecule, and so on,
//    (An n-tuple is implemented as a vector of closed intervals (ClosedInterval).)
//There is a separate n-tuple created for every possible combination
//of helix records from every pdb-file.
//    The process is repeated for SHEET records as well.
//All possible n-tuples of sheets are appended to the list
//of all possible n-tuples of helices.
//    The resulting list of n-tuptles of intervals is saved in the "vvI"
//argument (which is implemented as vector of n-tuples).
//
// Other arguments:
//    The "aMolecules[]" argument is required for bullet-proofing.
// Some programs, like minrms, allow the user to remove some of
// the residues from the portion of the molecule they are considering.
// apMolecules[] stores an array-of-pointers-to-molecules corresponding
// to each PDB-file in the aaFilenames[] array.  However these molecules
// only contain residues that the user wants to consider in the calculation.
// We need this information to insure that the only helix and sheet records
// that lie inside the subset residues that the user has selected, 
// will be considered.  The others are discarded quietly.
// Note: Helix or sheet records that lie partially inside and partially outide
//       this subset produce an error message and a call to exit(-1);
//
//    The optional "min_interval_size" argument indicates whether you
// want to automatically extend intervals containing fewer than this
// number of residues, so that they fit the minimum size.

void
ConvertHelixSheetRecordsToIntervalNtuples(int  num_files,
                                          vector<string> const& vFileNames,
                                          vector<vector<ClosedInterval> >& vvI,
                                          Biopolymer const *aMolecules,
                                          int min_interval_size = 0);

} //namespace minrms






#endif //#ifndef _INTERVALS_H


