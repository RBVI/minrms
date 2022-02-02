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


#ifndef _LOAD_PDB_H
#define _LOAD_PDB_H

#include <iostream>
using namespace std;

#include "biopolymer.h"

#define DBG_FILE_PARSER 0  //(This just turns off the printing of some
                           // messages used for debugging. Please ignore.)








namespace minrms {

istream & operator >> (istream&     pdb_file,
                       Biopolymer&  molecule);


//    Deletes residues from the linear molecule that do not contain
// all of the residues specified in "atomsThatMustBePresent"
// which is a NULL or "" terminated array of C-strings.
//    For example, if we are using the positions of the alpha-carbons,
// then remove all residues which do not contain at least one
// alpha-carbon atom.
//  Running time: O(a * n)
//                where a is the total number of atoms in the molecule,
//                  and n is the number of atom-names to check for.
void DeleteResiduesWithMissingAtoms(Biopolymer& molecule,
                                    char const **aAtomsThatMustBePresent);


//Alternatively, rather than use the StoreBackboneAtomName() member function
//manually, one might want to automatically read these names in from a file.
//This file stores a separate atom name on every line.
//(note: whitespace within each line is not skipped!)
void LoadBackboneAtomSymbols(char const *filename);




#ifdef DEBUG
//Print out some qualitative information about a biopolymer.
//(This is useful for debugging only.)
ostream & operator << (ostream&     out_file,
                       Biopolymer&  molecule);
//Print out some qualitative information about an atom.
ostream & operator << (ostream& out_file,
                       Atom&    atom);
#endif //#ifdef DEBUG




} //namespace minrms


#endif // #ifndef _LOAD_PDB_H
