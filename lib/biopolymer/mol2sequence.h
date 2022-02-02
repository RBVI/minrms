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


#ifndef _MOL2SEQUENCE_H
#define _MOL2SEQUENCE_H
#include "biopolymer.h"

namespace minrms {

//   Mol2Sequence() converts a Biopolymer it is given into a C-string
//containing that molecule's 1-letter lookup codes.  Prints a warning if
//that molecule contains a residue whose one-letter lookup name is unkown
//and, the corresponding string will contain the warning character:
//  NameCodeLookup::UNKNOWN_RESIDUE (see "residue_names.cc")
//at that location.
//
// Note:  The storage for the C-string is allocated within the function.
//        It will have to be deallocated later.
string Mol2Sequence(Biopolymer const& m);

//Returns true if the Biopolymer contains a residue
//whose one-letter-lookup-code is unknown.  (Running time O(n))
bool  ContainsUnknownResidues(Biopolymer const& m);

} //namespace minrms

#endif //#ifndef _MOL2SEQUENCE_H
