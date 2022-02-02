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


#include <cstring> //defines strcmp(), strlen(), strncpy()
#include <fstream>
using namespace std;

#include "residue_names.h"


using namespace minrms;

char NameCodeLookup::UNKNOWN_RESIDUE = 'x';

map<string, char> NameCodeLookup::nameCodeTableStr2Char;

map<char, string> NameCodeLookup::nameCodeTableChar2Str;

void
NameCodeLookup::Init(string res_code_dict)
{
  // There are two tasks
  // 1) fill the nameCodeTableStr2Char[] table.
  // 2) fill the nameCodeTableChar2Str[] table.

  // ---- First, fill the nameCodeTableStr2Char[] table.                   ---
  // ---- (This is a fast lookup table that goes in the reverse direction. ---
  // ----  It us used by the member function Lookup(char).)                ---

  ifstream input_file(res_code_dict.c_str(), ios::in);

  if ((res_code_dict.size() == 0) || (! input_file))
  {
    //If there was no file supplied by the user, fill it with default values.
    if (res_code_dict.size() != 0)
      cerr
        << "Using default 1-letter residue codes.\n"
                   " (no 3-letter to 1-letter-residue-code dictionary found.)"
        << endl;

    nameCodeTableStr2Char["GLY"] = 'G';
    nameCodeTableStr2Char["ALA"] = 'A';
    nameCodeTableStr2Char["SER"] = 'S';
    nameCodeTableStr2Char["CYS"] = 'C';
    nameCodeTableStr2Char["VAL"] = 'V';

    nameCodeTableStr2Char["THR"] = 'T';
    nameCodeTableStr2Char["ILE"] = 'I';
    nameCodeTableStr2Char["PRO"] = 'P';
    nameCodeTableStr2Char["MET"] = 'M';
    nameCodeTableStr2Char["ASP"] = 'D';

    nameCodeTableStr2Char["ASN"] = 'N';
    nameCodeTableStr2Char["LEU"] = 'L';
    nameCodeTableStr2Char["LYS"] = 'K';
    nameCodeTableStr2Char["GLU"] = 'E';
    nameCodeTableStr2Char["GLN"] = 'Q';

    nameCodeTableStr2Char["ARG"] = 'R';
    nameCodeTableStr2Char["HIS"] = 'H';
    nameCodeTableStr2Char["PHE"] = 'F';
    nameCodeTableStr2Char["TYR"] = 'Y';
    nameCodeTableStr2Char["TRP"] = 'W';

  } //if ((res_code_dict.size() == 0) || (! input_file))
  else
  {
    //Otherwise read in the contents of the file.

    ShortString first_symbol;
    ShortString second_symbol;
    input_file >> first_symbol;
    input_file >> second_symbol;

    if (((! input_file.good()) && (! input_file.eof())) ||
        (strcmp(first_symbol,"unknown") != 0) ||
        (strlen(second_symbol) != 1))

      ERR("\"res_code_dict.txt\" file has invalid format.\n"
          "The first line should be formated as:\n"
          "\"unknown c\", where \"c\" represents the character\n"
          "used to indicate that a residue-name code is not recognized.\n");

    second_symbol[0] = UNKNOWN_RESIDUE;
  
    int row_counter = 2;

    while (input_file)
    {
      string strResName;
      char   chResName;
      input_file >> strResName;
      if (((! input_file.good()) && (! input_file.eof()))
          || 
          (strResName.size() > VERY_SHORT_STRING_LENGTH))
        ERR("\"res_code_dict.txt\" file has invalid format"
            " on row " << row_counter << ".\n"
            "Each line after the first line should be formated as:\n"
            "<residue-name> <1-char-name-code>\n"
            "The <residue-name> cannot exceed "
            << VERY_SHORT_STRING_LENGTH << "characters.\n"
            "(Aborting)");
      input_file >> chResName;
      if ((! input_file.good()) && (! input_file.eof()))
        ERR("\"res_code_dict.txt\" file has invalid format"
            " on row " << row_counter << ".\n"
            "Each line after the first line should be formated as:\n"
            "<residue-name> <1-char-name-code>\n"
            "The <residue-name> cannot exceed "
            << VERY_SHORT_STRING_LENGTH << "characters.\n"
            "(Aborting)");
      nameCodeTableStr2Char[ strResName ] = chResName;
      ++row_counter;

    } //while (input_file)


  } // else clause for "if ((res_code_dict.size() == 0) || (! input_file))"


  // ---- Now, fill the nameCodeTableChar2Str[] array.                               ---
  // ---- (This is a lookup table that goes in the reverse direction.)               ---

  for(map<string, char>::const_iterator p = nameCodeTableStr2Char.begin();
      p != nameCodeTableStr2Char.end();
      p++)
  {
    char const *residue_name  = p->first.c_str();
    char one_letter_code = p->second;
    nameCodeTableChar2Str[one_letter_code] = string(residue_name);
  }


} // NameCodeLookup::NameCodeLookup()





