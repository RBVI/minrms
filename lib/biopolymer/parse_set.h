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


#ifndef _PARSE_SET_H
#define _PARSE_SET_H

// *** parse_set.h
// *** Minrms allows the user to specify subsets of residues from
// *** a molecule.  This file contains the header for a function
// *** which handles parsing the strings that are input by the
// *** user to specify these subsets.

#include <vector>
using namespace std;

#include "biopolymer.h"
#include "intervals.h"

// ParseChimeraSet()
// converts strings input by the user
// to sets of residues within a molecule,
// using the Chimera/Midas syntax.
//    The output is formatted as
// a collection of disjoint closed intervals.
// The union of these intervals is the set
// read by ParseChimeraSet().
//    The input is an excerpt from a string,
// (specified by the start, and stop string iterators)
// as well as the molecule to which the string refers.
// Once the set of residues has been parsed,
// the start iterator is set to point to the first character
// in the string that follows the set of residues.
// Here are some examples of the string format:
//
// "100"         all residues whose seqNums are 100
// "100-150.A"   all residues in chain A whose seqNums are in [100,150]
// "*-150.A"     all the residues in chain A whose seqNums are up to 150
// "150-*.A"     all the residues in chain A whose seqNums are at least 150
// "*.A"         all the residues in chain A
// ".A"           "   "     "     "    "   "
// "100-150.A-C" residues 100-150 in chains A through C
// "100-150.*"   residues 100-150 in all chains
// "100-150."       "      "   "    "    "   "
// "*.*"         the entire molecule

namespace minrms 
{


void
ParseChimeraSet(string::const_iterator& start,
                string::const_iterator  stop,
                Biopolymer const& m,
                vector<ClosedInterval>& dest);


//The next version reads in all of the text until one
//of the characters in "terminators" is encountered.
//Then, the text up until that point is parsed using the function
//above.
//If not all of the characters up until the terminator are used up during
//the parsing, ParseChimeraSet() returns an error message and exits.
void
ParseChimeraSet(string::const_iterator& start,
                string::const_iterator  stop,
                Biopolymer const& m,
                vector<ClosedInterval>& dest,
                string terminators);

} //namespace minrms 

#endif //#ifndef _PARSE_SET_H




