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


#ifndef _PDB_RES_ID_H
#define _PDB_RES_ID_H

#include <iostream>
using namespace std;



namespace minrms {
// PDBid's are used in records in PDB-files.  They act as residue identifiers.
//They containt a triplet: a chainId    (from which chain the residue belongs
//                                      If there's only one chain,
//                                      it is set to ' ' (space))
//                         a seqNum     (indicates the location of
//                                      the residue is in the chain)
//                         an insertCode (optional, can be used to insert
//                                      extra residues.
//                                      If no there are no extra insertions,
//                                      this field is set to ' ' (space))

struct PDBresID
{
  char    chainId;        // \  These
  int     seqNum;         // -> identify
  char    insertCode;     // /  the residue

  PDBresID(char setChainId = ' ',
           int  setSeqNum = 0,
           char setInsertCode = ' ')
  {
    chainId = setChainId;
    seqNum = setSeqNum;
    insertCode = setInsertCode;
  }

  bool  operator <  (const PDBresID& B) const;//true if residue A comes
                              // before B in the sequence
  bool  operator <=  (const PDBresID& B) const;
  bool  operator >  (const PDBresID& B) const;//true if residue A comes
                              // after B in the sequence
  bool  operator >=  (const PDBresID& B) const;

  bool  operator == (const PDBresID& B) const;// true if residue A
                                                   // and residue B refer
                                                   // to the same location
                                                   // in the sequence
                                                   // (have matching seq_num
                                                   //  and insert_code)
  bool  operator != (const PDBresID& B) const;// oposite of above
}; //class PDBresID


ostream& operator << (ostream&  out_file,
                      PDBresID  res_id);

static const PDBresID NULL_RESIDUE = PDBresID('\0', 0, '\0');

} //namespace minrms

#endif // #ifndef _PDB_ID_H

/*
class PDBidSequence
{
  // --- (In the arrays below, indexing starts at 0, not 1) ---
  vector<char> vPDBids;//One-letter-identifiers for each entry (residue) in the
                 //backbone.  (These are presently used for input-
                 //error-checking when reading-in MSF-files.)

  // FindNoGreaterThan() takes a vector of ResInfos named 'res_list',
  // assumes that res_list is sorted in order of the residue-id#s
  // binary-searches the sub-interval of the list from [lo, high] to find
  // the highest element in that sub-interval that's no greater than 'entry', 
  // and returns the index in the list where that element is stored.
  // (By 'no-greater-than', I mean no further down the amino acid sequence.)
  // It returns (low - 1) if the 'entry' is less than all of the elements
  // in [low, high].
  //Performance:
  //  Worst-case time: log2(high-low), but it only takes unit time, O(1),if the
  //  element we are looking for is at either end of the interval (high or low)
  //  When parsing through PDB-files, this is usually the case, since
  //  PDB-files are usually sequential.)
  long FindNoGreaterThan(// const vector< LoadedResidueInfo > *pvResList,
            //It now uses (*this) instead of (*pvResList)
               ChainSeqInsertTriplet const * pEntry,
               long low,
               long high) const;

public:
  PDBidSequence(long n):vPDBids(n) {}
  long size() const {return vNames.size();}
  char GetPDBid(long index) const; {return vNames[index];}

  //Returns location in the sequence where the
  // chainId, seqNum, insertCode triplet appears.
  //If it was not found, it returns -1.  Otherwise
  //it returns the index in the array where that element appears.
  inline long Find(ChainSeqInsertTriplet const *pEntry) const
  {
    long i = FindNoGreaterThan( pEntry, 0, vResList.size()-1);
    ASSERT( i >= -1 );
    // if it was not found it returns -1;
    if ((i == -1) ||
        (*pEntry != vResList[i].id))
      return -1;
    else
      return i;
  }

} // class PDBidSequence
*/

