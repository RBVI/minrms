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


#include "residue_names.h"
#include "mol2sequence.h"

namespace minrms
{

//   Mol2Sequence() converts a Biopolymer it is given into a C-string
//containing that molecule's 1-letter lookup codes.  Prints a warning if
//that molecule contains a residue whose one-letter lookup name is unkown
//and, the corresponding string will contain the warning character:
//  NameCodeLookup::UNKNOWN_RESIDUE (see "residue_names.cc")
//at that location.
//
// Note:  The storage for the C-string is allocated within the function.
//        It will have to be deallocated later.
string
Mol2Sequence(Biopolymer const& m)
{
  string sequence;
  char const *unknown_res_name = NULL;
  for(int s = 0; s < m.size(); ++s) {
        char c = NameCodeLookup::Lookup(m[s].name);
        sequence.push_back(c);
        if ((c == NameCodeLookup::UNKNOWN_RESIDUE) && (! unknown_res_name))
          unknown_res_name = m[s].name;
  }

  if (unknown_res_name)
  {
    cerr << "\n"
      "         ------------------------------------------------------\n"
      "Warning: The molecule contains a residue with an unknown\n"
      "         1-letter lookup code.\n"
      "         The name of the offending residue is:\""
         << unknown_res_name << "\".\n"
      "\n"
      "         To avoid these errors, edit (or create) a file that contains\n"
      "         a line with a 1-letter lookup code (in addition to the other\n"
      "         standard residue types).  (This file is often called \n"
      "         \"rescodes.file\" or \"res_code_dict.txt\") \n"
      "             For example:\n"
      "         "<< unknown_res_name << " Z\n"
      "\n"
      "         (See the \"MinRMS\" documentation for the file format.)\n"
         << endl;
  }

  return sequence;
} //Mol2Sequence()


//Returns true if the Biopolymer contains a residue
//whose one-letter-lookup-code is unknown.  (Running time O(n))
bool
ContainsUnknownResidues(Biopolymer const& m)
{
  for(int s = 0; s < m.size(); ++s) {
        char c = NameCodeLookup::Lookup(m[s].name);
        if (c == NameCodeLookup::UNKNOWN_RESIDUE)
          return true;
  }
  return false;
} //ContainsUnknownResidues()

} //namespace minrms
