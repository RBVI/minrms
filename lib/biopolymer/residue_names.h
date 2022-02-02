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


#ifndef _PEPTIDE_NAMES_H
#define _PEPTIDE_NAMES_H


#include <cstdlib> // defines "sizeof()"  (see below)
#include <vector>
#include <map>
#include <string>
using namespace std;


#include <short_string.h> // defines "VeryShortString"
#include <global_utils.h>
//#include <pdb++.h>

namespace minrms {

  // NameCodeLookup contains a collection functions
  // for looking up translations from 1 letter residue codes,
  // to 4-letter residue codes, and back again.
  //
  //   PDB-files typically indicate residue-type using a 4-character
  //   (actually, usually they only use the first 3-characters) name-codes.
  //   However, MSF files use only 1 character, per residue.
  //
  // Use this table lookup 'A' from "ALA", etc...
  //
  // Do not instantiate NameCodeLookup.  It's just a namespace.)
  // Before using any functions in NameCodeLookup, you MUST invoke
  // NameCodeLookup::Init() first.

  class NameCodeLookup
  {
    //Lookup table
    static map<string, char> nameCodeTableStr2Char;

    //Inverse lookup table
    static map<char, string> nameCodeTableChar2Str;

    NameCodeLookup(); //abstract class.  do not instantiate.

  public:
    static void Init(string res_code_dict = "");
                      // Initialize the name-code lookup table.
                      // Seeks for a file whose name is indicated by the
                      // res_code_dict argument.
                      // This file stores the correspondences
                      // between 4-letter and 1-letter codes.

    //Lookup takes a 3-letter residue name code (like "ALA" or "GLU")
    //and returns its corresponding 1-character code (like 'A' or 'E')
    //If it's unable to find anything for the 3-letter code you gave it,
    //it returns UNKNOWN_RESIDUE.
    static char Lookup(char const *name)
    {
      string temp_str(name);
      map<string, char>::iterator p = nameCodeTableStr2Char.find(temp_str);
      if (p == nameCodeTableStr2Char.end())
        return UNKNOWN_RESIDUE;
      else
        return p->second;
      //return nameCodeTableStr2Char[ temp_str ];
    }

    static char UNKNOWN_RESIDUE;

    static string Lookup(char c) { return nameCodeTableChar2Str[c]; }

  }; // class NameCodeLookup

} //namespace minrms

#endif // #ifndef _PEPTIDE_NAMES_H


