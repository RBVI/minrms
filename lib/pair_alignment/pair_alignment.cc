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

#include <fstream>
#include <cctype>
#include <cstdio> //needed for call to sprintf()
#include <cstring>//needed for call to strcmp() and strncmp() 
#include <ctime>  //note: this line seems to work, both on SGI and and alpha.
                  //      Not sure if this is more generally portable, however.
using namespace std;


#include <global_utils.h>
#include <stl_utils.h>
#include <parse_utils.h>
#include <biopolymer.h>
#include <mol2sequence.h>
#include <superimpose.h>
#include <fast_rot_metric.h>
#include "pair_alignment.h"

namespace minrms
{


char PairwiseAlignment::GAP_CHAR = '.';
char PairwiseAlignment::GAP_CHAR_OFFSET = '~';

inline bool IsGapChar(char c)
{
  #if 0
  char gap_char = PairwiseAlignment::GAP_CHAR;
  char gap_char_offset = PairwiseAlignment::GAP_CHAR_OFFSET;
  bool return_val = (c == gap_char);
  if (! return_val)
    return_val = (c == gap_char_offset);
  cerr << "c = '" << c
       << "'; gap_char = '" << gap_char
       << "'; gap_char_offset = '" << gap_char_offset
       << "'; return_val = " << return_val << endl;
  #endif //#if 0

  return ((c == PairwiseAlignment::GAP_CHAR) ||
          (c == PairwiseAlignment::GAP_CHAR_OFFSET));
}


char PairwiseAlignment::aColorNames[][32] = {"white", "cyan", "yellow"};
char PairwiseAlignment::unit_name[32] = "Angstrom";
const Real PairwiseAlignment::DONT_USE_MARKERS_FOR_SHORT_LENGTHS = -1.0f;
const Real PairwiseAlignment::MIDAS_MARKER_SIZE = 0.4f;

//The following function is almost identical to the similarly
//named function in "../dyn/distance_metric.h"
//However there are no "maximum_distance", or "UPPER_BOUND"
//cutoffs in this version.
Real SumSqdDistBetweenResidues(Biopolymer::const_iterator i,
                               Biopolymer::const_iterator j)
{
  Real output = 0.0f;

  Biopolymer::Residue::const_iterator pa[2];
  for(int q = 0; 
      q < Biopolymer::Residue::NumBackboneAtoms();
      ++q)
  {
    pa[0] = (*i).GetBackboneAtom(q);
    pa[1] = (*j).GetBackboneAtom(q);
    assert(pa[0] != i->end());
    assert(pa[1] != j->end());
    Vect3 displacement_between_ij;
    SubtractVect3((*pa[0]).second.xyz,
                  (*pa[1]).second.xyz,
                  displacement_between_ij);
    output += DotProduct3(displacement_between_ij, displacement_between_ij);
  }
  return output;
}


PairwiseAlignment::PairwiseAlignment(const PairwiseAlignment &source)
{
  space_available = 0;
  num_matches = 0;
  aMatches = NULL;
  *this = source;
}


void PairwiseAlignment::resize(int n)
{
  assert(n >= 0);
  assert(space_available >= 0);
  assert(aMatches || (space_available == 0));
  if (n > space_available)
  {
    if (aMatches) delete [] aMatches;
    space_available = n;
    aMatches = new PairwiseMatch [space_available];
    if (! aMatches)
      ERR(__FILE__ << ":PairwiseAlignment::resize("
          << n << ") cannot alloc mem.\n"
          "Please report this bug to the developer.");
  }
  num_matches = n;
}



void PairwiseAlignment::reserve(int n)
{
  assert(n >= 0);
  assert(space_available >= 0);
  resize(n);
  resize(0);
}




//ReadNextMSFlabel() reads the first symbol
//(i.e. sequence of characters with no WhiteSpace)
//on this line and stores them in the "dest" argument.
//If there is nothing else on this line, then it skips to the next line.
//It keeps skipping until it either encounters a line with non-
//WhiteSpace text after the first symbol, or an end-of-file is reached.
static void
ReadNextMSFlabel(istream &theFile, string& dest)
{
  if (! theFile)
    return;

  char c = '\n';
  int num_tokens = 0;

  //loop over lines in the file until one contains
  //at least 2 blocks of non-whitespace
  while (num_tokens < 2)
  {
    dest.clear();
    ParseUtils::SkipWhiteSpaceExceptNewlines( theFile );
    //get the first token
    while (theFile.get(c),
           (theFile && (! isspace(c))))
    {
      dest.push_back(c);
      num_tokens = 1;
    }
    //If it was terminated by a '\n', we need it later. put it back in the file.
    //(We only read one line per iteration in this loop.)
    if (c == '\n')
      theFile.putback(c);
    else //else, skip to the next bock of text on this line
      ParseUtils::SkipWhiteSpaceExceptNewlines( theFile );
    if (theFile.eof()) return;

    theFile.get(c);
    if (theFile.eof()) return;
    //If there is another block of text on this line, leave it alone
    //and put back the character from it, we just took out.
    assert((! isspace(c)) || (c == '\n'));
    if (c != '\n')
    {
      theFile.putback(c);
      num_tokens = 2;
    }
    else
    {
      //Otherwise, skip past the '\n' and proceed to the next line
      //(to do this, simply don't putback the character c).
    }
  } //while (num_tokens < 2)
} //ReadNextMSFlabel()



void
PairwiseAlignment::ImportMSF(string filename,
                             string labelA,
                             string labelB,
                             Biopolymer const& m1,
                             Biopolymer const& m2,
                             bool check_for_lower_case)
{
  ifstream msf_file(filename.c_str(), ios::in);
  if (! msf_file)
    ERR("Error opening MSF file: \"" << filename
        << "\" for reading.");
  ImportMSF(msf_file,
            labelA,
            labelB,
            m1,
            m2,
            check_for_lower_case);
}


void
PairwiseAlignment::ImportMSF(istream &msf_file,
                             string labelA,
                             string labelB,
                             Biopolymer const& m1,
                             Biopolymer const& m2,
                             bool check_for_lower_case)
{
  //copy the molecules' sequences of characters into temporary arrays
  vector<string> vSequences(2);
  vSequences[0] = Mol2Sequence(m1);
  vSequences[1] = Mol2Sequence(m2);

  ImportMSF(msf_file,
            labelA,
            labelB,
            vSequences[0],
            vSequences[1],
            check_for_lower_case);
}



void
PairwiseAlignment::ImportMSF(string filename,
                             string labelA,
                             string labelB,
                             string sequenceA,
                             string sequenceB,
                             bool check_for_lower_case)
{
  ifstream msf_file(filename.c_str(), ios::in);
  if (! msf_file)
    ERR("Error opening MSF file: \"" << filename
        << "\" for reading.");
  ImportMSF(msf_file,
            labelA,
            labelB,
            sequenceA,
            sequenceB,
            check_for_lower_case);
}



void
PairwiseAlignment::ImportMSF(istream& msf_file,
                             string   labelA,
                             string   labelB,
                             string   sequenceA,
                             string   sequenceB,
                             bool     check_for_lower_case)
{
  vector<string> vLabels(2);
  vLabels[0].assign(labelA);
  vLabels[1].assign(labelB);
  vector<string> vSequences(2);
  vSequences[0].assign(sequenceA);
  vSequences[1].assign(sequenceB);

  ImportMSF(msf_file, vLabels, vSequences, check_for_lower_case);
}




void
PairwiseAlignment::ImportMSF(istream &msf_file,
                             vector<string> const& vLabels,
                             vector<string> const& vSequences,
                             bool check_for_lower_case)
{
  assert(vLabels.size() == 2);
  assert((vLabels[0] != "") && (vLabels[1] != ""));
  assert(vSequences.size() == 2);

  DEBUG_MSG(DBG_IMPORT_MSF, "Beginning to import MSF file.\n");
  // ********************************************************
  // *****      Skip past all the comments and crap     *****
  // *****  that's at the beginning of every msf-file.  *****
  // ********************************************************
  const char line_separator[4]={'/', '/', '\n', '\0'};//string containing "//"
  char aLastRead[3] = {'\0', '\0', '\0'};
  while(msf_file)
  {
    aLastRead[0] = aLastRead[1];
    aLastRead[1] = aLastRead[2];
    msf_file.get(aLastRead[2]);
    if (strncmp(aLastRead, line_separator, 3) == 0)
      break;
  }

  //check for file read errors and end of file
  if (! msf_file.good())
    ERR("Invalid input MSF file.\n"
        //"\"" << filename << "\":\n"
        "Expecting a separating line containing _only_ the "
        "string: \""
        "/"  //must separate the double forward slash
        "/"  //or compiler thinks it's a comment
        "\\n\".\n"
        "(Where \"\\n\" indicates a new-line (return cairrage), "
        "_not_ \"\\n\")\n"
        "There should be no blank space (tabs or spaces) \n"
        "before or after the \""
        "/"
        "/"
        "\" in your separating line.");

  DEBUG_MSG(DBG_IMPORT_MSF, "Reached the end of the comments section \"/"
            "/\"\n");

  //Okay now we're right before the first line of the MSF-file
  //(give or take some whitespace and/or number-lines)

  vector<string> vContent(2); //All the characters 
                               //(including gap characters '.')
                               //belonging to the relevant
                               //lines of the MSF-file
                               //(ie. those beginning with either label:
                               //     vLabels[0] or vLabels[1]),
                               //excluding the sequence labels themselves
                               //(vLabels[0] or vLabels[1]).

  vector<int> vOrderInMSF(2);
  vOrderInMSF[0] = -1; //Stores the order that the MSF file contains
  vOrderInMSF[1] = -1; //the two sequences being compared.
                       //This is necessary since the relative order
                       //of the two sequences in the MSF-file,
                       //may be different than the order
                       //the to labels were specified
                       //in the vLabels[] array.
                       //       Explicit definition:
                       // If the sequence represented by vLabels[1]
                       // appears in the file before the one
                       // represented by vLabels[0], then vOrderInMSF[0] = 1
                       // and vOrderInMSF[1] = 0.
                       // Otherwise, vOrderInMSF[0] = 0 and vOrderInMSF[1] = 1.

  int iter = 0; //Alternates from 0 to 1 to 0 to 1, each time a
                 //new line is read containing the label for either sequence.

  while(msf_file) //Looping over all lines in the file.
                  //Each pass through this loop reads in a line of
                  //actual residues from one of the two sequences in it.
  {
    // **************************************************************
    // ******  Now, skip to a content line, that is, a line     *****
    // ******  containing the amino acids from either sequence. *****
    // ******  These lines begin with one of the two labels     *****
    // ******  vLabels[0] or vLabels[1].                        *****
    // **************************************************************

    string label; //The label preceding the line we are about to read in.
    int which_label = -1; //Indicates the numeric ID of this label,
                          //if (label == vLabels[0]), then it is 0,
                          //if (label == vLabels[1]), then it is 1.
                          //Initially, it is given an impossible value (-1)
                          //indicating we do not know which label
                          //is comming up yet.
    while(which_label == -1)
    {
      ReadNextMSFlabel(msf_file, label);
      if (msf_file.eof())
      {
        if (vOrderInMSF[0] == -1) //If we have not found either label yet.
          ERR("Error: neither sequence \"" << vLabels[0] << "\"\n"
              "       nor sequence \"" << vLabels[1] << "\"\n"
              "       was found in this msf-file.");
        else if (iter == 1)
          ERR("Error: MSF file ends prematurely or sequence not found.\n"
              "       Expecting a line containing sequence "
              "\"" << vLabels[vOrderInMSF[1]] << "\"");
        else
          break; //else, things are okay.
      }
      DEBUG_MSG(DBG_IMPORT_MSF, " beginning new line in MSF-file. Label = \""
                << label << "\"");
      if (label == vLabels[0])
        which_label = 0;
      else if (label == vLabels[1])
        which_label = 1;
      else //Skip past the next '\n' character
      {
        string dummy_str; //needed to pacify getline
        getline(msf_file, dummy_str);
      }
    }
    if (msf_file.eof())
      break; //break out of this loop.
    assert(msf_file.good() && (which_label != -1));

    // *****************************************************
    // ******   If this is the first "content line",   *****
    // ******  read so far, do some extra setup-work.  *****
    // *****************************************************

    //If this is the first "content-line",
    //make a note of which label appeared first in the file,
    //(vLabels[0] or vLabels[1]).
    //Subsequent lines should alternate between these two labels.
    if (vOrderInMSF[0] == -1)
    {
      vOrderInMSF[0] = which_label;
      vOrderInMSF[1] = (which_label + 1) % 2;
      DEBUG_MSG(DBG_IMPORT_MSF, "Established the order of the labels in the MSF-file:\n"
                "                First sequence is \""
                << vLabels[vOrderInMSF[0]] << "\"\n"
                "               Second sequence is \""
                << vLabels[vOrderInMSF[1]] << "\"");
    }
    //Otherwise, if this is not the first line, we check the order that the
    //labels appear in the MSF file to make sure it is consistent throughout.
    else if (which_label != vOrderInMSF[iter])
      ERR("Error: Invalid MSF file.  Inconsistent sequence ordering\n"
          "       or missing sequence label.\n"
          "       Expected sequence label not found: \""
          << vLabels[vOrderInMSF[iter]] << "\"");

    if (iter == 0)
      DEBUG_MSG(DBG_IMPORT_MSF, "------- Start of new line -------");

    DEBUG_MSG(DBG_IMPORT_MSF, "Reading in residues from sequence \""
              << vLabels[vOrderInMSF[iter]] << "\":");
    // Now, read in the next line
    char c;
    while(ParseUtils::SkipWhiteSpaceExceptNewlines(msf_file),
          msf_file.get(c),
          (msf_file && (c != '\n')))
    {
      vContent[vOrderInMSF[iter]].push_back(c);
      #ifdef DEBUG
      #if DBG_IMPORT_MSF
      cerr << c << flush;
      #endif //#if DBG_IMPORT_MSF
      #endif //#ifdef DEBUG
    }
    #ifdef DEBUG
    #if DBG_IMPORT_MSF
    cerr << endl;
    #endif //#if DBG_IMPORT_MSF
    #endif //#ifdef DEBUG

    //Check to make sure the last two lines are the same length
    if ((iter == 1) && (vContent[0].size() != vContent[1].size()))
      ERR("Invalid MSF file format:\n"
          //"\"" << filename << "\":\n"
          "Two of the lines containing sequences content\n"
          "are not the same length.\n");
    iter = (iter + 1) % 2;;
  } //while(msf_file)


  // *******************************************************
  // ***** At this point, all we have done is skip     *****
  // ***** past all the comments spaces, labels and    *****
  // ***** other crap in this MSF file.                *****
  // ***** The "vContent" array of strings contains   *****
  // ***** what's left afterwards.                     *****
  // *****   The following function does the actual    *****
  // ***** work of figuring out what residues were     *****
  // ***** actually matched together.                  *****
  // *******************************************************
  DigestMSFContents(vContent,
                    vLabels,
                    vSequences,
                    check_for_lower_case);
} //PairwiseAlignment::ImportMSF()









// ********************************************************************
// ********* The following couple member functions are       **********
// ********* called by "ImportMSF()" to help pares MSF-files **********
// ********************************************************************




void
PairwiseAlignment::DigestMSFContents(vector<string> const& vContent,
                                     vector<string> const& vLabels,
                                     vector<string> const& vSequences,
                                     bool check_for_lower_case)
{
  assert(vContent.size() == 2);
  assert(vContent[0].size() == vContent[1].size());
  #ifdef DEBUG
  assert(vLabels.size() == 2);
  assert((vLabels[0] != "") && (vLabels[1] != ""));
  #endif //#ifdef DEBUG

  assert(vSequences.size() == 2);
  //If the sequences are not just empty strings,
  //then check the content of the MSF-files against them.
  bool check_sequences = (vSequences[0].size()!=0) && (vSequences[1].size()!=0);

  vector<PairwiseMatch> vMatchesCopy; //stores a local copy of the matches
                                      //that will end up in the "aMatches[]"
                                      //array.

  vector<int> vCurRes(2);      //Counts the number of characters from either
  vCurRes[0] = vCurRes[1] = 0; //sequence that have been digested so far.
  for (int n=0; //loop over the characters in the vContent[] strings.
       n < vContent[0].size();
       ++n)
  {
    if ((! IsGapChar(vContent[0][n])) && (! IsGapChar(vContent[1][n])))
    {
      vMatchesCopy.push_back(PairwiseMatch(vCurRes[0], vCurRes[1]));
      DEBUG_MSG(DBG_IMPORT_MSF,
                "Read-in match #" << vMatchesCopy.size() << ":\n"
                "The " << vCurRes[0] + 1
                << "th residue of \"" << vLabels[0]
                << "\" (" << vContent[0][n] << "), was matched with\n"
                "the "
                << vCurRes[1] + 1
                << "th residue of \"" << vLabels[1]
                << "\" (" << vContent[1][n] << ")");
    }
    //Check the MSF file against the sequence in the vSequences[] arrays.
    if (check_sequences)
    {
      for(int s=0; s < 2; ++s)
      {
        if ((! IsGapChar(vContent[s][n])) &&
            (vContent[s][n] != vSequences[s][vCurRes[s]]))
        {
          cerr <<
            "Warning: residue #"
               << vCurRes[s]+1 <<
            " from sequence \""
               << vLabels[s] <<
            "\": ("
               << vContent[s][n] <<
            ") differs\n"
            "       from the MSF file differs from the expected residue: ("
               << vSequences[s][vCurRes[s]] <<
            ")." << endl;
        }
      }
    }
    //Now, update the counters:
    if (! IsGapChar(vContent[0][n])) ++vCurRes[0];
    if (! IsGapChar(vContent[1][n])) ++vCurRes[1];
  }//loop over the characters in the vContent[] strings.

  // ***********************************************************
  // *****   Done digesting.  Now check for stupid input.  *****
  // ***********************************************************

  //Check the sequence lengths
  if (check_sequences && (vCurRes[0] != vSequences[0].size()))
    ERR("Number of residues read in MSF file for sequence \""
        << vLabels[0] << "\" (" << vCurRes[0] << ")\n"
        "is not equal to the number of residues in sequence \""
        << vLabels[0] << "\" (" << vSequences[0].size() << ").");
  if (check_sequences && (vCurRes[1] != vSequences[1].size()))
    ERR("Number of residues read in MSF file for sequence \""
        << vLabels[0] << "\" (" << vCurRes[1] << ")\n"
        "is not equal to the number of residues in sequence \""
        << vLabels[0] << "\" (" << vSequences[1].size() << ").");

  //If this is a previously allocated Pairwise alignment, check to make
  //sure that we read in exactly the right amout of matches
  if ((size() != 0) && (vMatchesCopy.size() != size()))
  {
    ERR("Error reading MSF file:\n"
        //"\"" << filename << "\".\n"
        "The number of matches ("
        << vMatchesCopy.size() << ") does not equal the expected number ("
        << size() << ")");
  }
  //Otherwise, allocate space for the alignment
  else
  {
    if (vMatchesCopy.size() == 0)
      cerr <<"Warning: Possible error reading MSF file:\n"
        //"\"" <<filename<< "\".\n"
        "There do not appear to be any matches made in this alignment."
           << endl;
    
    resize(vMatchesCopy.size()); //( <-Note: this also sets num_matches = vMatchesCopy.size() )
  }

  //Now (finally), copy the contents of vCopyOfMatches[], into aMatches[]
  assert(aMatches);
  for(int n = 0; n < vMatchesCopy.size(); ++n)
  {
    aMatches[n] = vMatchesCopy[n];
  }
} // PairwiseAlignment::DigestMSFContents()






bool
ReadMSFsequence(istream& msf_file,
                string label,
                string& sequence)
{
  // ********************************************************
  // *****      Skip past all the comments and crap     *****
  // *****  that's at the beginning of every msf-file.  *****
  // ********************************************************
  bool label_found_in_file = false;
  const char line_separator[4]={'/', '/', '\n', '\0'};//string containing "//"
  char aLastRead[3] = {'\0', '\0', '\0'};
  while(msf_file.good()) //while there may be more characters left to read(),
  {
    aLastRead[0] = aLastRead[1];
    aLastRead[1] = aLastRead[2];
    msf_file.get(aLastRead[2]);
    if (strncmp(aLastRead, line_separator, 3) == 0)
      break;
  }
  //check for file read errors and end of file
  if (! msf_file.good())
    ERR("Invalid input MSF file.\n"
        //"\"" << filename << "\":\n"
        "Expecting a separating line containing _only_ the "
        "string: \""
        "/"  //must separate the double forward slash
        "/"  //or compiler thinks it's a comment
        "\\n\".\n"
        "(Where \"\\n\" indicates a new-line (return cairrage), "
        "_not_ \"\\n\")\n"
        "There should be no blank space (tabs or spaces) \n"
        "before or after the \""
        "/"
        "/"
        "\" in your separating line.");

  //Okay, now we're right before the content section of the MSF-file
  //Start looking for lines beginning with the label symbol.

  assert(label.size() != 0);
  sequence.clear();
  while(msf_file)
  {
    bool next_label_found = false;
    while(! next_label_found)
    {
      string label_in_msf;
      ReadNextMSFlabel(msf_file, label_in_msf);
      if (msf_file.eof())
        break;
      DEBUG_MSG(DBG_IMPORT_MSF, " beginning new line in MSF-file. Label = \""
                << label_in_msf << "\"");
      if (label_in_msf == label)
      {
        next_label_found = true;
        label_found_in_file = true;
      }
      else //else, skip to the beginning of the next line.  
        getline(msf_file, label_in_msf);
        
    } //while(! next_label_found)

    //check for eof()
    if (msf_file.eof())
      break;
    assert(next_label_found);

    //Now copy the remaining characters on this line into sequence
    string current_line;
    getline(msf_file, current_line);
    for(string::const_iterator p = current_line.begin();
        p != current_line.end();
        ++p)
    {
      if ((! isspace(*p)) && (! IsGapChar(*p)))
        sequence.push_back(*p);
    }
  } //while(msf_file)
  return label_found_in_file;
} //ReadMSFsequence()
                   



Real
PairwiseAlignment::CalcMinRMSD(Biopolymer const& m1,
                               Biopolymer const& m2,
                               Matrix3x4 optimal_transform,
                               Superimpose *pSuperimpose,
                               Vect3* aRotateUsingTheseCoords1,
                               Vect3* aRotateUsingTheseCoords2
                               ) const
{
  Vect3* aRotateUsing1 = NULL;
  Vect3* aRotateUsing2 = NULL;
  int max_num_matches = MIN(m1.size(), m2.size());

  Superimpose *pMinRotRMSD;
  if (pSuperimpose)
    pMinRotRMSD = pSuperimpose;
  else
    pMinRotRMSD = new Superimpose(max_num_matches);
  
  if (aRotateUsingTheseCoords1)
    aRotateUsing1 = aRotateUsingTheseCoords1;
  else
    aRotateUsing1 =
      new Vect3[m1.size() * Biopolymer::Residue::NumBackboneAtoms()];

  if (aRotateUsingTheseCoords2)
    aRotateUsing2 = aRotateUsingTheseCoords2;
  else
    aRotateUsing2 =
      new Vect3[m2.size() * Biopolymer::Residue::NumBackboneAtoms()];

  if ((! aRotateUsing1) || (! aRotateUsing2) || (! pMinRotRMSD))
    ERR("Failure to allocate mem at \""
                  << __FILE__ << "\":" << __LINE__
                  << "\n.Please report this error to the developer.");


  int m; //indexes through the matched residue-pairs in the alignment
  int which_pair = 0; //indexes through the actual atoms within each
                      //residue-pair
  for (m=0; m < size(); ++m)
  {
    int i,j;
    i = aMatches[m][0];
    j = aMatches[m][1];
    for(int q = 0; 
        q < Biopolymer::Residue::NumBackboneAtoms();
        ++q)
    {
      Biopolymer::Residue::const_iterator pA1, pA2;
      pA1 = m1[i].GetBackboneAtom(q);
      pA2 = m2[j].GetBackboneAtom(q);

      //Make sure matched residues contain the atoms selected for matching
      if (pA1 == m1[i].end())
        ERR("Error: Cannot calculate the RMSD of this structural alignment.\n"
            "       One of the residues in the alignment lacks an essential atom:\n"
            "       Residue/Residue #" << i+1  << " (" << m1[i].name <<
            ") from the first\n"
            "       structure, is missing a \""
            << Biopolymer::Residue::LookupBackboneAtomSymbol(q) << "\" atom.\n"
            "       The position of this atom is needed to determine the\n"
            "       RMSD between the two structures.\n."
            "       (By default, the \" CA\" atoms from each residue are\n"
            "        required to compute of RMSD of an alignment\n"
            "        This can be changed or customized by supplying\n"
            "        an \"atoms_used.txt\" file. See documentation.)");
      else if (pA2 == m2[j].end())
        ERR("Error: Cannot calculate the RMSD of this structural alignment.\n"
            "       One of the residues in the alignment lacks an essential atom:\n"
            "       Residue/Residue #" << j+1 << " (" << m2[j].name <<
            ") from the second\n"
            "       structure, is missing a \"" <<
            Biopolymer::Residue::LookupBackboneAtomSymbol(q) << "\" atom.\n"
            "       The position of this atom is needed to determine the\n"
            "       RMSD between the two structures.\n."
            "       (By default, the \" CA\" atoms from each residue are\n"
            "        required to compute of RMSD of an alignment\n"
            "        This can be changed or customized by supplying\n"
            "        an \"atoms_used.txt\" file. See documentation.)");

      aRotateUsing1[which_pair][0] = (*pA1).second.xyz[0];
      aRotateUsing1[which_pair][1] = (*pA1).second.xyz[1];
      aRotateUsing1[which_pair][2] = (*pA1).second.xyz[2];
      aRotateUsing2[which_pair][0] = (*pA2).second.xyz[0];
      aRotateUsing2[which_pair][1] = (*pA2).second.xyz[1];
      aRotateUsing2[which_pair][2] = (*pA2).second.xyz[2];

      ++which_pair;
    }
  }

  //Pass these arrays to FindTransformMinRMS()
  Real rmsd = pMinRotRMSD->FindTransformMinRMSD(size()
                                   * Biopolymer::Residue::NumBackboneAtoms(),
                                                aRotateUsing1,
                                                aRotateUsing2,
                                                optimal_transform);
  if (! aRotateUsingTheseCoords1)
    delete [] aRotateUsing1;
  if (! aRotateUsingTheseCoords2)
    delete [] aRotateUsing2;
  if (! pSuperimpose)
    delete pMinRotRMSD;

  return rmsd;
}


Real
SumSqdDist(PairwiseAlignment const& a,
           Biopolymer const& m1,
           Biopolymer const& m2)
{
  Real sum_sqd_dist = 0.0f;

  //Loop over the pairs of matched segments (or "residues", er whatever)
  //in the alignment.
  for(int n=0; n < a.NumMatches(); ++n)
  {
    int i = a[n][0]; //the nth matched segment from molecule #1
    int j = a[n][1]; //the nth matched segment from molecule #2
    assert((0 <= i) && (i < m1.size()) && (0 <= j) && (j < m2.size()));

    sum_sqd_dist += 
      SumSqdDistBetweenResidues(m1.begin()+i, m2.begin()+j);
  }
  return sum_sqd_dist;
}


Real RMSD(PairwiseAlignment const& a,
           Biopolymer const& m1,
           Biopolymer const& m2)
{
  Real ssd = SumSqdDist(a, m1, m2);
  return sqrt(ssd / (Biopolymer::Residue::NumBackboneAtoms()*a.NumMatches()));
}


bool
PairwiseAlignment::IsMonotonic() const
{
  int i;
  int j;
  int prev_i = -1;
  int prev_j = -1;
  for(int m=0; m < size(); ++m)
  {
    i = aMatches[m][0];
    j = aMatches[m][0];
    if ((prev_i >= i) || (prev_j >= j))
      return false;
    prev_i = i;
    prev_j = j;
  }
  return true;
}




static void ReadInChainSeqInsertFromFSSP(PDBresID &id,
                                         istream &fssp_file)
{
  ParseUtils::SkipWhiteSpaceExceptNewlines( fssp_file );
  string sym1, sym2;
  string residueNumberAndInsertCode;
  fssp_file >> sym1;
  ParseUtils::SkipWhiteSpaceExceptNewlines( fssp_file );
  if (! fssp_file)
    ERR("Error FSSP file ends prematurely, in the middle of an EQUIVALENCE line.\n");
  char c;
  fssp_file.get(c);     // If the first character of the next symbol is numeric,
  fssp_file.putback(c); // then this next symbol contains the sequence# and
  if (fssp_file && isdigit(c)) //insert code for the next residue.
  {
    if (sym1.size() != 1) ERR("Error in EQUIVALENCES section:\n"
                               "Expected a 1-letter chain identifier code.");
    id.chainId = sym1[0];
    residueNumberAndInsertCode = sym2;
    fssp_file >> residueNumberAndInsertCode;
    if (! fssp_file)
      ERR("Error FSSP file ends prematurely, in the middle of an EQUIVALENCE line.\n");
  }
  else
  {
    id.chainId = ' ';
    residueNumberAndInsertCode = sym1;
  }

  //Now figure out the sequenceNumber and the ChainId from
  //residueNumberAndInsertCode.
  //This can be messy to do because the two pieces of information are
  //merged into the symbol, which may contain a ')' character which
  //will need to be removed beforehand.
  if (residueNumberAndInsertCode.size() == 0)
  {
    assert(! fssp_file);
    ERR("Error in EQUIVALENCES section:\n"
        "FSSP file ends prematurely.");
  }
  char last_char = *(residueNumberAndInsertCode.end() - 1);
  if (last_char == ')')
  {
    //Then, just chop it off and try again.
    //It is a peculiarity of FSSP files that they enclose their
    //equivalences in parenthesis, so you have to be careful to remove
    //parenthesis from the strings that you read in.
    //Chop off last char:
    residueNumberAndInsertCode.erase(residueNumberAndInsertCode.end()-1);

    if (residueNumberAndInsertCode.size() < 1)
      ERR("Error: in FSSP file.  Missing residue identifier before ')'.");
    last_char = *(residueNumberAndInsertCode.end() - 1);
  }

  //Is there an insertCode?
  if (isdigit(last_char))
  {
    //then, there was no insert code, so pass ' ' by default.
    id.insertCode = ' ';
  }
  else
  {
    //Chop off last char:
    residueNumberAndInsertCode.erase(residueNumberAndInsertCode.end()-1);

    if (residueNumberAndInsertCode.size() < 1)
      ERR("Error: in FSSP file.  Incomplete reside identifier.  No seqRes#");
    id.insertCode = last_char;
  }
  //What's left is the sequence number:
  id.seqNum = atoi(residueNumberAndInsertCode.c_str());

  DEBUG_MSG(DBG_IMPORT_FSSP,
            "Residue " << id << " read from FSSP file.");

} //ReadInChainSeqInsertFromFSSP()



void PairwiseAlignment::ImportFSSP(string filename,
                                   Biopolymer const& m1,
                                   Biopolymer const& m2,
                                   string labelA,
                                   string labelB)
{
  assert(filename != "");

  DEBUG_MSG(DBG_IMPORT_FSSP,"Importing fssp-file: \""
            << filename << "\"");

  ifstream fssp_file(filename.c_str(), ios::in);
  if (! fssp_file)
    ERR("Error opening FSSP file: \"" << filename
                  << "\" for reading.");
}



void PairwiseAlignment::ImportFSSP(istream &fssp_file,
                                   Biopolymer const& m1,
                                   Biopolymer const& m2,
                                   string labelA,
                                   string labelB)
{
  vector<Biopolymer> vMol(2);
  //I know this is kinda wasteful of space, but
  //I don't want to open up the pandora's box of using pointers.
  vMol[0] = m1;
  vMol[1] = m2;

  vector<string> vLabels(2);
  vLabels[0] = labelA;
  vLabels[1] = labelB;
  
  ImportFSSP(fssp_file, vMol, vLabels);
}


void PairwiseAlignment::ImportFSSP(istream &fssp_file,
                                   vector<Biopolymer> const& vMol,
                                   vector<string> const& vLabels)
{
  assert(vMol.size() == 2);
  assert(vLabels.size() == 2);
  assert(vLabels[0] != "");
  assert(vLabels[1] != "");

  //First, allocate temporary space to store the matches, as they are read
  //from the file.  ( If resize() has been previously called, then the
  //contents will eventually be copied into this previously allocated
  //"aMatches[]" array.  Otherwise, Alloc() will be called with the
  //appropriate size = number-of-matches, as read from the msf-file, and
  //the contents copied over. )
  int max_num_matches = MIN( vMol[0].size(), vMol[1].size() );
  assert(max_num_matches > 0);

  PairwiseMatch *aCopyOfMatches = new PairwiseMatch [max_num_matches];
  if (! aCopyOfMatches)
    ERR("Error alocating memory: \"" << __FILE__ <<"\":"<< __LINE__);

  //*** Time to read in the file.  First skip past the "## EQUIVALENCES:" line
  bool found_equivalences_section = false;
  while(fssp_file.good() && (! found_equivalences_section))
  {
    string line;
    getline(fssp_file, line);
    found_equivalences_section =
      (strncmp(line.c_str(),"## EQUIVALENCES:", 16) == 0);
  }

  if (! found_equivalences_section)
    ERR("Format error in FSSP-file:\n"
        //"\"" << filename << "\""
        " lacks a line beginning with the keyword: \n"
        "## EQUIVALENCES:");

  DEBUG_MSG(DBG_IMPORT_FSSP,"Entering \"## EQUIVALENCES:\" SECTION");

  //skip past the next line
  for(char c = '\0'; (c != '\n'); fssp_file.get(c))
    {}
  
  //okay, ready to begin

  int N = 0; //counts a running total of the number of matches made

  vector<int> vOrderInFSSP(2);
  vOrderInFSSP[0] = -1; //Stores the order that the FSSP file contains
  vOrderInFSSP[1] = -1; //the two sequences being compared.
                       //This is necessary since the relative order
                       //of the two sequences in the FSSP-file,
                       //may be different than the order
                       //the to labels were specified
                       //in the vLabels[] array.
                       //       Explicit definition:
                       // If the structure represented by vLabels[1]
                       // appears in the FSSP file before the one
                       // represented by vLabels[0], then vOrderInFSSP[0] = 1
                       // and vOrderInFSSP[1] = 0.
                       // Otherwise, vOrderInFSSP[0] = 0 and vOrderInFSSP[1] = 1.

  //Now read in the matched pairs of residues.
  while(fssp_file.good())
  {
    ParseUtils::SkipWhiteSpaceExceptNewlines( fssp_file );
    
    string str;
    fssp_file >> str; //skip the first token(usually a numeric marker "1:")

                     //If the first symbol on this line is a "##", then
    if (str == "##") //we are leaving the ## EQUIVALENCES section.
      break;         //Consequently, stop reading.  We're done.

    DEBUG_MSG(DBG_IMPORT_FSSP,"--- Reading in a new line begining with: \""
              << str << "\" ---");
    ParseUtils::SkipWhiteSpaceExceptNewlines( fssp_file );
    string testsymbol1, testsymbol2;
    fssp_file >> testsymbol1; //read in the first structure name
    fssp_file >> testsymbol2; //read in the second structure name

    DEBUG_MSG(DBG_IMPORT_FSSP,"two identifiers from FSSP are \""
              << testsymbol1 << "\" and \"" << testsymbol2 << "\"");


     //If the name of the first two symboles on this line, match the
    //two labels, then this line of text contains the alignment
    //data we are seeking.
    if (((testsymbol1 == vLabels[0]) && (testsymbol2 == vLabels[1]))
        ||
        ((testsymbol2 == vLabels[0]) && (testsymbol1 == vLabels[1])))
    {
      if (vOrderInFSSP[0] == -1)   //If this is the first such line, establish
      {                                 //which label comes first in the file:
        if (testsymbol1 == vLabels[0])  //vLabels[0] or vLabels[1]
        {
          vOrderInFSSP[0] = 0;
          vOrderInFSSP[1] = 1;
        }
        else
        {
          vOrderInFSSP[0] = 1;
          vOrderInFSSP[1] = 0;
        }
      }
      else  //If it is not the first such line, then check to see
      {     //that the ordering of labels is consistent.
        if (testsymbol1 != vLabels[vOrderInFSSP[0]])
          ERR("Format error in FSSP-file:\n"
              "  the first structure's label on this line in: \""
              << testsymbol1 <<
              "\n  does not match the first label given\n"
              "in earlier lines: \""
              << vLabels[vOrderInFSSP[0]] << "\"");
      }

      assert(testsymbol1 == vLabels[ vOrderInFSSP[0] ]);
      assert(testsymbol2 == vLabels[ vOrderInFSSP[1] ]);

      DEBUG_MSG(DBG_IMPORT_FSSP,"Block of contiguous matched residues found.\n");
      ParseUtils::SkipNonWhiteSpace( fssp_file ); //skip past the starting number from mol1
      ParseUtils::SkipWhiteSpaceExceptNewlines( fssp_file );
      ParseUtils::SkipNonWhiteSpace( fssp_file ); //skip past the "-" string
      ParseUtils::SkipNonWhiteSpace( fssp_file ); //skip past the ending number from mol1
      ParseUtils::SkipWhiteSpaceExceptNewlines( fssp_file );
      ParseUtils::SkipNonWhiteSpace( fssp_file ); //skip past the "<=>" string
      ParseUtils::SkipNonWhiteSpace( fssp_file ); //skip past the starting number from mol2
      ParseUtils::SkipWhiteSpaceExceptNewlines( fssp_file );
      ParseUtils::SkipNonWhiteSpace( fssp_file ); //skip past the "-" string
      ParseUtils::SkipNonWhiteSpace( fssp_file ); //skip past the ending number from mol2
      ParseUtils::SkipWhiteSpaceExceptNewlines( fssp_file );

      // Now, skip to the '('
      char c;
      do {
        fssp_file.get(c);
        if ((c == '\n') || (!fssp_file))
          ERR("ERROR:  While reading in the \"## EQUIVALENCES\" section\n"
                        "        of an FSSP-file, a '(' character was omited.\n");
      } while (c != '(');

      string str;
      string a_start_name;
      string   a_end_name;
      string b_start_name;
      string   b_end_name;
      PDBresID a_start_id, a_end_id, b_start_id, b_end_id;

      fssp_file >> a_start_name;
      ReadInChainSeqInsertFromFSSP(a_start_id,
                                   fssp_file);
      fssp_file >> str;
      if (str != "-")
        ERR("Error, \"-\" expected in EQUIVALENCE section of FSSP file.\n");
      fssp_file >> a_end_name;
      ReadInChainSeqInsertFromFSSP(a_end_id,
                                   fssp_file);
      fssp_file >> str;
      if (str != "<=>")
        ERR("Error, \"<=>\" expected in EQUIVALENCE section of FSSP file.\n");

      fssp_file >> b_start_name;
      ReadInChainSeqInsertFromFSSP(b_start_id,
                                   fssp_file);
      fssp_file >> str;
      if (str != "-")
        ERR("Error, \"-\" expected in EQUIVALENCE section of FSSP file.\n");
      fssp_file >> b_end_name;
      ReadInChainSeqInsertFromFSSP(b_end_id,
                                   fssp_file);

      ParseUtils::SkipWhiteSpaceExceptNewlines( fssp_file );
      if (! fssp_file)
        ERR("Error: incomplete line in the \"## EQUIVALENCES\" section\n"
            "       of FSSP file:\n"
            //"\"" << filename << "\"\n"
            "       -or- error in format.  Aborting...");

      DEBUG_MSG(DBG_IMPORT_FSSP,"["
                << a_start_id << "," << a_end_id << "] from protein \""
                << vLabels[ vOrderInFSSP[0] ] << 
                "\" matches  ["
                << b_start_id << "," << b_end_id << "] from protein \""
                << vLabels[ vOrderInFSSP[0] ] << "\"");
      Biopolymer::const_iterator a_start = vMol[ vOrderInFSSP[0] ].find(a_start_id);
      Biopolymer::const_iterator a_end   = vMol[ vOrderInFSSP[0] ].find(a_end_id);
      Biopolymer::const_iterator b_start = vMol[ vOrderInFSSP[1] ].find(b_start_id);
      Biopolymer::const_iterator b_end   = vMol[ vOrderInFSSP[1] ].find(b_end_id);

      assert((vMol[ vOrderInFSSP[0] ].begin() <= a_start) && (a_start <= vMol[ vOrderInFSSP[0] ].end()));
      assert((vMol[ vOrderInFSSP[0] ].begin() <= a_end) && (a_end <= vMol[ vOrderInFSSP[0] ].end()));
      assert((vMol[ vOrderInFSSP[1] ].begin() <= b_start) && (b_start <= vMol[ vOrderInFSSP[1] ].end()));
      assert((vMol[ vOrderInFSSP[1] ].begin() <= b_end) && (b_end <= vMol[ vOrderInFSSP[1] ].end()));

      if (a_start == vMol[ vOrderInFSSP[0] ].end())
        ERR("Residue " << a_start_id << " not found in \"" << vLabels[ vOrderInFSSP[0] ] << "\"");
      if (a_end == vMol[ vOrderInFSSP[0] ].end())
        ERR("Residue " << a_end_id << " not found in \"" << vLabels[ vOrderInFSSP[0] ] << "\"");
      if (b_start == vMol[ vOrderInFSSP[1] ].end())
        ERR("Residue " << b_start_id << " not found in \"" << vLabels[ vOrderInFSSP[1] ] << "\"");
      if (b_end == vMol[ vOrderInFSSP[1] ].end())
        ERR("Residue " << b_end_id << " not found in \"" << vLabels[ vOrderInFSSP[1] ] << "\"");

      //*** check for errors ***
      // Verify the names of the residues loaded.
      // Holm and Sander use a special indexing method that discard residues
      // lacking certain backbone-atoms.  We want to make sure that the
      // index-numbers correspond to the residues we think they do.
      // The names of the residues are stored in the second half of the line
      // enclosed in () parenthesis.


      if (a_start_name != (*a_start).name)
        ERR("Error: inconsistency in numbering of residues between\n"
            "       FSSP file and order they appear in the PDB-file.\n"
            
            "       One possible solution to this error is to remove all residues\n"
            "       that are missing any of the following atoms\n"
            "       from molecule \"" << vLabels[ vOrderInFSSP[0] ] << "\":\n"
            "       \" N\", \" C\", \" CA\", \" O\" (backbone atoms)\n"
            "Details : Residue #" << (*a_start).id << 
            "from molecule \"" << vLabels[ vOrderInFSSP[0] ] << "\"'s"
            " name = \"" << (*a_start).name << "\"\n"
            "But the " << a_start-vMol[ vOrderInFSSP[0] ].begin()+1 << "th residue is named \""
            << a_start_name
            << "\" in the FSSP-file.\n"
            "Aborting...");

      if (a_end_name != (*a_end).name)
        ERR("Error: inconsistency in numbering of residues between\n"
            "       FSSP file and order they appear in the PDB-file.\n"
            "       One possible solution to this error is to remove all residues\n"
            "       that are missing any of the following atoms\n"
            "       from molecule \"" << vLabels[ vOrderInFSSP[0] ] << "\":\n"
            "       \" N\", \" C\", \" CA\", \" O\" (backbone atoms)\n"
            "Details : Residue #" << (*a_end).id << 
            "from molecule \"" << vLabels[ vOrderInFSSP[0] ] << "\"'s"
            " name = \"" << (*a_end).name << "\"\n"
            "But the " << a_end-vMol[ vOrderInFSSP[0] ].begin()+1 << "th residue is named \""
            << a_end_name
            << "\" in the FSSP-file.\n"
            "Aborting...");



      if (b_start_name != (*b_start).name)
        ERR("Error: inconsistency in numbering of residues between\n"
            "       FSSP file and order they appear in the PDB-file.\n"
            "       One possible solution to this error is to remove all residues\n"
            "       that are missing any of the following atoms\n"
            "       from molecule \"" << vLabels[ vOrderInFSSP[1] ] << "\":\n"
            "       \" N\", \" C\", \" CA\", \" O\" (backbone atoms)\n"
            "Details : Residue #" << (*b_start).id << 
            "from molecule \"" << vLabels[ vOrderInFSSP[1] ] << "\"'s"
            " name = \"" << (*b_start).name << "\"\n"
            "But the " << b_start-vMol[ vOrderInFSSP[1] ].begin()+1 << "th residue is named \""
            << b_start_name
            << "\" in the FSSP-file.\n"
            "Aborting...");

      if (b_end_name != (*b_end).name)
        ERR("Error: inconsistency in numbering of residues between\n"
            "       FSSP file and order they appear in the PDB-file.\n"
            "       One possible solution to this error is to remove all residues\n"
            "       that are missing any of the following atoms\n"
            "       from molecule \"" << vLabels[ vOrderInFSSP[1] ] << "\":\n"
            "       \" N\", \" C\", \" CA\", \" O\" (backbone atoms)\n"
            "Details : Residue #" << (*b_end).id << 
            "from molecule \"" << vLabels[ vOrderInFSSP[1] ] << "\"'s"
            " name = \"" << (*b_end).name << "\"\n"
            "But the " << b_end-vMol[ vOrderInFSSP[1] ].begin()+1 << "th residue is named \""
            << b_end_name
            << "\" in the FSSP-file.\n"
            "Aborting...");

      #if 0 //commenting out. No longer complain if alignment is not sequential
      if ((a_start > a_end) || (b_start > b_end))
        ERR("Error reading FSSP file:\n"
            "Invalid Intervals: "
            << (*a_start).id << " - " << (*a_end).id << "  and  "
            << (*b_start).id << " - " << (*b_end).id << "\n"
            "   Cannot read in FSSP files containing residues in reverse order.\n"
            "   Alignments must be sequential (sometimes called \"monotonic\").");
      #endif //#if 0


      if ((a_end - a_start) != (b_end - b_start))
        ERR("Error reading FSSP file:\n"
            "Intervals: "
            << (*a_start).id << " - " << (*a_end).id << "  and  "
            << (*b_start).id << " - " << (*b_end).id << "\n"
            "  do not have the same size"
            "("
            << (a_end-a_start) + 1
            << ", and "
            << (b_end-b_start) + 1
            << ", respectively)");


      for(int counter = 0; counter <= (a_end-a_start); ++counter)
      {
        ++N;
        if (N > max_num_matches)
          ERR("Error: FSSP file contains too many matches (" 
              << N << ")\n"
              "       for the two proteins in question, whose sizes are:\n"
              "       " << vMol[ vOrderInFSSP[0] ].size() << " and "
              << vMol[ vOrderInFSSP[1] ].size() );

        //now, store the nth match in this alignment
        aCopyOfMatches[N-1][vOrderInFSSP[0]] = (a_start-vMol[ vOrderInFSSP[0] ].begin()) + counter;
        aCopyOfMatches[N-1][vOrderInFSSP[1]] = (b_start-vMol[ vOrderInFSSP[1] ].begin()) + counter;
        DEBUG_MSG(DBG_IMPORT_FSSP,"match#" << N << ": residue "
                  << aCopyOfMatches[N-1][vOrderInFSSP[0]]+1
                  << " from \"" << vLabels[ vOrderInFSSP[0] ] << "\" <--matches-with--> "
                  << aCopyOfMatches[N-1][vOrderInFSSP[1]]+1
                  << " from \"" << vLabels[ vOrderInFSSP[1] ] << "\".");
      }
    } // if (strcmp(testsymbol2, vLabels[ vOrderInFSSP[1] ]) == 0)

    // skip past the remaining text in this line

    #if 0
    char c;
    do {
      fssp_file.get(c);
    } while (c != '\n');
    #endif 

    getline(fssp_file, testsymbol1);
    //also skip any whitespace following it
    ParseUtils::SkipWhiteSpace( fssp_file );
  } // while(fssp_file.good())  Reading in the alignment


  //If this is a previously allocated Pairwise alignment, check to make
  //sure that we read in exactly the correct number-of-matches
  if (size() != 0)
  {
    if (N < size())
      ERR("Error reading FSSP file:\n"
          //"\"" << filename << "\".\n"
          << "The number of matches ("
          << N << ") is less than the expected number ("
          << size() << ")");
    if (N > size())
      ERR("Error reading FSSP file: \n"
          //"\"" << filename << "\".\n"
          << "The number of matches ("
          << N << ") exceeds the expected number ("
          << size() << ")");
  }
  //Otherwise, allocate space for the alignment
  else
  {
    resize(N); //( <-Note: this also sets num_matches = n )
  }
  if (N == 0)
    ERR("Error.  Alignment in FSSP file contains 0 pairs of matched residues.\n");

  //Now (finally), copy the contents of aCopyOfMatches[], into aMatches[]
  assert(aMatches);
  assert(aCopyOfMatches);
  for(--N; (N>=0); --N)
  {
    aMatches[N][0] = aCopyOfMatches[N][0];
    aMatches[N][1] = aCopyOfMatches[N][1];
  }

  //Check, one last time to make sure we have monotonicity:
  if (! IsMonotonic())
    ERR("Error reading FSSP file:\n"
        "   Cannot read in FSSP files containing residues in reverse order.\n"
        "   Alignments must be sequential (sometimes called \"monotonic\").");
    
  delete [] aCopyOfMatches;
} //PairwiseAlignment::ImportFSSP()




void
PairwiseAlignment::ExportMSF(string filename,
                             vector<Biopolymer> const& vMol,
                             vector<string> const& vLabels,
                             string comments,
                             bool  show_markers,
                             bool  compress_using_lower_case) const
{
  assert(vMol.size() == 2);
  assert(vLabels.size() == 2);
  PairwiseAlignment::ExportMSF(filename,
                               vMol[0],
                               vMol[1],
                               vLabels[0],
                               vLabels[1],
                               comments,
                               show_markers,
                               compress_using_lower_case);
}

void
PairwiseAlignment::ExportMSF(string filename, //required for contents of the MSF-file
                             Biopolymer const& m1,
                             Biopolymer const& m2,
                             string labelA,
                             string labelB,
                             string comments,
                             bool  show_markers,
                             bool  compress_using_lower_case) const
{
  if (compress_using_lower_case
      && 
      (ContainsUnknownResidues(m1) || ContainsUnknownResidues(m2)))
  {
    compress_using_lower_case = false;
    cerr <<
      "Warning: At least one of the molecules contains residues with\n"
      "         unknown 1-letter lookup codes.\n"
      "         The \"compress_using_lower_case\" option will be dissabled\n"
      "         This means that the resulting MSF-file generated will not\n"
      "         mix upper and lower case to merge regions of mutual gaps\n"
      "         All letters will be left in their original cases, and all\n"
      "         gaps will be displayed speparately."
         << endl;
  }

  //copy the molecules' sequences of characters into temporary arrays
  vector<string> vSequences(2);
  vSequences[0] = Mol2Sequence(m1);
  vSequences[1] = Mol2Sequence(m2);

  ExportMSF(filename,
            vSequences[0],
            vSequences[1],
            labelA,
            labelB,
            comments,
            show_markers,
            compress_using_lower_case);
}





void
PairwiseAlignment::ExportMSF(string filename,
                             string sequenceA,
                             string sequenceB,
                             string labelA,
                             string labelB,
                             string comments,
                             bool  show_markers,
                             bool  compress_using_lower_case) const
{
  vector<string> vLabels(2);
  vLabels[0].assign(labelA);
  vLabels[1].assign(labelB);

  vector<string> vSequences(2);
  vSequences[0].assign(sequenceA);
  vSequences[1].assign(sequenceB);

  ExportMSF(filename,
            vSequences,
            vLabels,
            comments,
            show_markers,
            compress_using_lower_case);
}



void
PairwiseAlignment::ExportMSF(string filename,
                             vector<string> const& vSequences,
                             vector<string> const& vLabels,
                             string comments,
                             bool  show_markers,
                             bool  compress_using_lower_case) const
{
  assert(filename != "");
  assert(vSequences.size() == 2);
  assert(vLabels.size() == 2);

  if (! IsMonotonic())
    ERR("Error writing MSF file:\n"
        "   Cannot output an MSF-file for this alignment because\n"
        "   the alignment pairs of residues matched in reverse order.\n"
        "\n"
        "   MSF-files are only meaningful in the context of alignments\n"
        "   that are sequential (sometimes called \"monotonic\").");

  ofstream msf_file(filename.c_str(), ios::out);
  if (! msf_file) {
    ERR("Can't open file for writing: \"" << filename << "\"");
  }

  int seq;

  const int cairrage_return_interval = 50;
  //const int marker_interval          = 25;
  const int space_interval           = 10;

  int labelSize = MAX(vLabels[0].size(), vLabels[1].size()) + 2;
  int padLengthA = labelSize - vLabels[0].size();
  int padLengthB = labelSize - vLabels[1].size();

  if (compress_using_lower_case)
  {
    char non_alphabetic_residue = '\0';
    for (int m = 0; m < 2; ++m)
    {
      for (int i = 0; i != vSequences.size(); ++i)
      {
        char c = vSequences[m][i];
        if ((c < 'A') || ('Z' < c))
          non_alphabetic_residue = c;
      } //for (int i = 0; i != vSequences.size(); ++i)
    } //for (int m = 0; m < 2; ++m)

    if (non_alphabetic_residue != '\0')
    {  
      cerr << "Warning: There are non-upper-case or non-alphabetic\n"
           << "         residue names (" << non_alphabetic_residue << ") present, so\n"
           << "         the lower-case mode is not available.  Using\n"
           << "         regular MSF file format instead."
           << endl;
      compress_using_lower_case = false;
    }
  } 

  int max_output_length = vSequences[0].size() + vSequences[1].size();
  int *aResCountedSoFar = new int [2];
  char **aaOutput = new char * [2];
  for (seq=0; seq < 2; ++seq)
  {
    aaOutput[seq] = new char [max_output_length];
  }


  for (seq=0; seq < 2; ++seq)
  {
    aResCountedSoFar[seq] = 0;
  }

  //now, caption/comment the msf-file
  msf_file << comments << "\n"; //add one blank line
  msf_file << filename << "   MSF: "
           << 
    (vSequences[0].size() - size()) +
    (vSequences[1].size() - size()) +
    size();

  msf_file << "  Type: ?  ";

  time_t absolute_time;
  time(&absolute_time);
  tm *pt = localtime(&absolute_time);
    
  char months[12][10] = { "January",
                          "February",
                          "March",
                          "April",
                          "May",
                          "June",
                          "July",
                          "August",
                          "September",
                          "October",
                          "November",
                          "December" };
  msf_file << months[pt->tm_mon] << " "
           << pt->tm_mday << ", "
           << pt->tm_year+1900 << "  "
           << pt->tm_hour << ":" << pt->tm_min
           << "  Check: 0 "
           << "..\n\n";

  //now, specify the name, length, and 'weights' of the two proteins
  msf_file << "Name: " << vLabels[0];
  for(int i=0; i < padLengthA; ++i)
    msf_file << ' ';//pad to uniform length
  msf_file << "Len: " << vSequences[0].size() <<"   Check: 0   Weight: 1.00\n"
           << "Name: " << vLabels[1];
  for(int i=0; i < padLengthB; ++i)
    msf_file << ' ';//pad to uniform length
  msf_file << "Len: " << vSequences[1].size() <<"   Check: 0   Weight: 1.00\n"
           << "\n"
           << "//\n"
           << "\n";

  //----------Now, output the alignment's contents to the msf-file--------
  int counter = 0;
  //    int cur_n = n;
  int cur_i, cur_j;
  int m; //indexes through the allignment
  //    int cur_i = (*pLS)[0].NumRes();
  //    int cur_j = (*pLS)[1].NumRes()

  //trace forwards through the alignment
  cur_i = 0;
  cur_j = 0;
  int previous_match_i = -1;
  int previous_match_j = -1;
  for(m=0; m <= size(); ++m)
  {
    if (m < size()) {
      cur_i = aMatches[m][0];
      cur_j = aMatches[m][1];
    }
    else {
      //Special case: m == size()
      //              Of course, we want to print all the residues.  This
      //              means we have to contine to print the residues in the
      //              two sequences that come after the last match.
      //              If m==size() then this lets us know that we want
      //              to keep printing until the end of each sequence, but
      //              we should not try to draw a match at position
      //              cur_i and cur_j
      //              (since they correspond to the
      //                (*pLS)[0].NumRes()+1 th residue from sequence A, and
      //                (*pLS)[1].NumRes()+1 th residue from sequence B
      //               which don't exist)
      cur_i = vSequences[0].size();
      cur_j = vSequences[1].size();
    }
    
    //okay, now start filling in the output-array:
    int i = previous_match_i + 1;//start at the i following the last match
    int j = previous_match_j + 1;//start at the j following the last match
    //First, draw the gaps in each sequence
    if (compress_using_lower_case)
    { //if using this mode, draw as many as possible directly overlapping
      while ((i < cur_i) && (j < cur_j))                 
      {
        aaOutput[0][ counter ] = tolower(vSequences[0][aResCountedSoFar[0]]);
        aaOutput[1][ counter ] = tolower(vSequences[1][aResCountedSoFar[1]]);
        ++(aResCountedSoFar[0]);
        ++(aResCountedSoFar[1]);
        ++counter;
        ++i;
        ++j;
      }
    } // if (compress_using_lower_case)
    
    //now if there are any remaining residues that were not paired off
    //above, ..or if we are not in "compress_using_lower_case" mode,
    //then output these remaining residues.
    while (i < cur_i)
    {
      aaOutput[0][ counter ] = vSequences[0][aResCountedSoFar[0]];
      aaOutput[1][ counter ] = PairwiseAlignment::GAP_CHAR;
      ++(aResCountedSoFar[0]);
      ++counter;
      ++i;
    }
    while (j < cur_j)
    {
      aaOutput[0][ counter ] = PairwiseAlignment::GAP_CHAR;
      aaOutput[1][ counter ] = vSequences[1][aResCountedSoFar[1]];
      ++(aResCountedSoFar[1]);
      ++counter;
      ++j;
    }

    //Now, (finally) draw the match
    if (m < size()) //(only do it if it's not the m==size() case (see above))
    {
      aaOutput[0][counter] = vSequences[0][aResCountedSoFar[0]];
      aaOutput[1][counter] = vSequences[1][aResCountedSoFar[1]];
      ++(aResCountedSoFar[0]);
      ++(aResCountedSoFar[1]);
      ++counter;
    }
        
    previous_match_i = cur_i;
    previous_match_j = cur_j;
  } // for(m=0; m <= n; ++m)

  assert( aResCountedSoFar[0] == vSequences[0].size() );
  assert( aResCountedSoFar[1] == vSequences[1].size() );

  // *****************************************************************
  //Now, write out the contents of the aaOutput array to the file which
  //we just filled
  //     We have to write a separate line for every sequence,
  //as well as a separate line for the markers displayed (for each sequence)
  //and we have to interrupt the whole thing by return-cairrages every
  //'cairrage_return_interval' characters.
  //    Indeed, this is as Conrad says: "a pain in the _______."

  for (seq=0; seq < 2; ++seq)
    aResCountedSoFar[seq] = 1;
  int output_length = vSequences[0].size() + vSequences[1].size() - size();
  assert(output_length == counter);
  counter = 0; //initialize counter to point to the first
  //elements in the aaOutput arrays
  while (counter < output_length)
  {
    int count_copy;
    if (show_markers)
    { // First, draw the markers

      for (int i = 0; i < labelSize; ++i) msf_file << ' '; //add space for
                                                           //protein label

#if 0
      //The last time I checked, the code below worked,
      //although I've made a few changes to the surrounding code.  It used
      //to generate evenly spaced markers at every "marker_interval"
      //characters.  I commented it out, though because the
      //resulting file does not conform to the MSF-file standard.
      //I don't want to delete this code yet.

      int num_spaces_crossed;
      for (count_copy = counter;
           ((count_copy < counter + cairrage_return_interval)
            && (count_copy < output_length));
           )
      {
        if (((count_copy+1) % marker_interval == 0) ||
            (count_copy == counter) ||
            (count_copy == counter+cairrage_return_interval-1))
        {
          char label[64];
          int labelSize = sprintf(label, "%d", count_copy+1);
          num_spaces_crossed = 0;
          for(int d=0; d<labelSize; ++d) {
            msf_file << label[d];  //copy the label to the file
            ++count_copy;
            if (count_copy % space_interval == 0)
              ++num_spaces_crossed;
          }
          count_copy -= num_spaces_crossed;
        }
        else
        {
          msf_file << ' ';
          ++count_copy;
        }

        if ((count_copy % space_interval == 0) && (num_spaces_crossed == 0))
          msf_file << ' ';
      } // for (...count_copy > counter - cairrage_return_interval...)
      msf_file << endl;
#endif //#if 0

      char label_begin[64];
      char   label_end[64];
      int labelSize_begin = sprintf(label_begin, "%d", counter+1);
      int num_chars_in_line = MIN(cairrage_return_interval,
                                  (output_length - counter));
      int   labelSize_end = sprintf(label_end,   "%d", 
                                    counter+num_chars_in_line);

      //calculate how many spaces to place between the two labels
      int num_spaces = ((counter+cairrage_return_interval+1)/space_interval)
                                             -
                                 ((counter + 1)/space_interval);
      int space_between_labels = num_chars_in_line + num_spaces
                                              -
                                 (labelSize_begin + labelSize_end);
      //Now, make sure there is at least one space
      space_between_labels = MAX(space_between_labels, 1);

      //draw the label for the beginning of this line
      msf_file << label_begin;

      //now draw the spaces
      for(int i = 0; i < space_between_labels; ++i)
        msf_file << ' ';

      //draw the end-label for this line
      msf_file << label_end << "\n";
    } // if (show_markers)
    
    //*********Now draw Sequence A: 
#if 0
    // --- I've commented out the code that tries to draw individual
    // --- counters for each sequence, instead of a single global counter.
    // --- It kindof worked, but it was buggy.

    // first, the sequence markers (if show_sequence_markers is selected)
    if (show_sequence_markers)
    { // Now, draw the markers for this sequence
        
      for (int i = 0; i < labelSize; ++i) msf_file << ' '; //add space for
                                                             //protein label
      for (count_copy = counter;
           ((count_copy < counter + cairrage_return_interval)
            && (count_copy < output_length));
           )
      {
        if ((aResCountedSoFar[0] % marker_interval == 0)
            &&
            (! IsGapChar(aaOutput[0][count_copy])))
        {
          char label[64];
          int labelSize = sprintf(label, "%d", aResCountedSoFar[0]);
          num_spaces_crossed = 0;
          for(int d=0; d<labelSize; ++d)
          {
            msf_file << label[d];
            if ((count_copy+1) % space_interval == 0)
              ++num_spaces_crossed;
            else if (! IsGapChar(aaOutput[0][count_copy]))
              ++(aResCountedSoFar[0]);
            ++count_copy;
          }
          count_copy -= num_spaces_crossed;
        }
        else
        {
          msf_file << '-'; //just print out a space
          //only increment the number of residues used-up if it's not a gap
          if (! IsGapChar(aaOutput[0][count_copy])
              ++(aResCountedSoFar[0]);
              ++count_copy;
        }
          
        if ((count_copy % space_interval == 0) && (num_spaces_crossed == 0))
          msf_file << ' ';
      } // for (...count_copy > counter - cairrage_return_interval...)
      msf_file << endl;
    }// if (show(sequence_markers))
#endif //#if 0

    // Now, (finally) write out the data for protein A
    msf_file << vLabels[0]; //first, put the protein label (identifier)
                                 //at the beginning of this line.
    for (int i = 0; i < padLengthA; ++i) msf_file << ' '; //pad-out label to
                                                          //uniform length
    // print out protein A's residue-letters
    for (count_copy = counter;
         ((count_copy < counter + cairrage_return_interval)
          && (count_copy < output_length));
         ++count_copy)
    {
      msf_file << aaOutput[0][count_copy];
      if ((count_copy+1) % space_interval == 0)
        msf_file << ' ';
    } // for (...count_copy > counter - cairrage_return_interval...)
    msf_file << endl;


    //*********Now draw Sequence B: 
    msf_file << vLabels[1]; //first, put the protein label (identifier)
                            //at the beginning of this line.
    for (int i = 0; i < padLengthB; ++i) msf_file << ' '; //pad-out label to
                                                            //uniform length
    // now, print out sequence B's residue-letters
    for (count_copy = counter;
         ((count_copy < counter + cairrage_return_interval)
          && (count_copy < output_length));
         ++count_copy)
    {
      msf_file << aaOutput[1][count_copy];
      if ((count_copy+1) % space_interval == 0)
        msf_file << ' ';
    } // for (...count_copy > counter - cairrage_return_interval...)
    msf_file << endl;

#if 0
    // --- I've commented out the code that tries to draw individual
    // --- counters for each sequence, instead of a single global counter.
    // --- It kindof worked, but it was buggy.

    // Then draw the sequence markers (if show_sequence_markers is selected)
    if (show_sequence_markers)
    { // Now, draw the markers for this sequence
      for (int i = 0; i < labelSize; ++i) msf_file << ' '; //add space for
                                                           //protein label
      for (count_copy = counter;
           ((count_copy < counter + cairrage_return_interval)
            && (count_copy < output_length));
           )
      {
        if ((aResCountedSoFar[1] % marker_interval == 0)
            &&
            (! IsGapChar(aaOutput[1][count_copy])))
        {
          char label[64];
          int labelSize = sprintf(label, "%d", aResCountedSoFar[1]);
          num_spaces_crossed = 0;
          for(int d=0; d<labelSize; ++d)
          {
            msf_file << label[d];
            if ((count_copy+1) % space_interval == 0)
              ++num_spaces_crossed;
            else if (! IsGapChar(aaOutput[1][count_copy]))
              ++(aResCountedSoFar[1]);
            ++count_copy;
          }
          count_copy -= num_spaces_crossed;
        }
        else
        {
          msf_file << '-'; //just print out a space
          //only increment the number of residues used-up if it's not a gap
          if (! IsGapChar(aaOutput[1][count_copy]))
            ++(aResCountedSoFar[1]);
          ++count_copy;
        }

        if ((count_copy % space_interval == 0) && (num_spaces_crossed == 0))
          msf_file << ' ';
      } // for (...count_copy > counter - cairrage_return_interval...)
      msf_file << endl;
    }// if (show(sequence_markers))
#endif //#if 0

    msf_file << "\n"; //add a small vertical gap before next round
    counter += cairrage_return_interval;
  } //while (counter < output_length)
  msf_file << flush;

  delete [] aResCountedSoFar;
  for (seq=0; seq < 2; ++seq)
  {
    delete [] aaOutput[seq];
  }

} //PairwiseAlignment::ExportMSF()







void
PairwiseAlignment::ExportGFX(Biopolymer const& m1,
                             Biopolymer const& m2,
                             string filename,
                             string caption,
                             bool show_residue_markers,
                             bool connect_the_dots,
                             Real use_marker_if_dist_less_than,
                             bool mark_worst_match,
                             unsigned int show_subset_size,
                             int subset_start_A,
                             int subset_start_B
                             ) const
{
  Biopolymer const *(apMol[2]);
  apMol[0] = &m1;
  apMol[1] = &m2;

  assert(filename != "");
  ofstream midas_file(filename.c_str(), ios::out);
  if (! midas_file) {
    ERR("Can't open file for writing: \""
        << filename << "\""
        );
  }

  // ***    First, Caption the image:
  // *** Find a good place to store the text label for this object.
  // *** To do this, I find the min/max xyz values of all the coordinates.
  // *** This defines a bounding box.  I then place the caption at the center
  // *** of one of the outer faces.
  // ***  (Midas & Chimera both require a 3-D (not 2-D) position where you
  // ***   would like the text to be displayed.  Unfortuneatly, at certain
  // ***   viewing angles this will overlap with the image being displayed.)
  Real min_coord[3];
  Real max_coord[3];

  bool first_atom = true;

  for (int which_protein = 0;
       which_protein < 2;
       ++which_protein)
  {
    for (Biopolymer::const_iterator pS = apMol[which_protein]->begin();
         pS != apMol[which_protein]->end();
         ++pS)
    {
      for (Biopolymer::Residue::const_iterator pA = pS->begin();
           pA != pS->end();
           ++pA)
      {
        if (first_atom)
        {
          min_coord[0] = max_coord[0] = (*pA).second.xyz[0];
          min_coord[1] = max_coord[1] = (*pA).second.xyz[1];
          min_coord[2] = max_coord[2] = (*pA).second.xyz[2];
          first_atom = false;
        }
        else
        {
          min_coord[0] = MIN(min_coord[0], (*pA).second.xyz[0]);
          min_coord[1] = MIN(min_coord[1], (*pA).second.xyz[1]);
          min_coord[2] = MIN(min_coord[2], (*pA).second.xyz[2]);

          max_coord[0] = MAX(max_coord[0], (*pA).second.xyz[0]);
          max_coord[1] = MAX(max_coord[1], (*pA).second.xyz[1]);
          max_coord[2] = MAX(max_coord[2], (*pA).second.xyz[2]);
        }
      } //loop over atoms
    } //loop over segments
  } //loop over molecules

  if (first_atom) { //make sure not min&max_coord[] are not uninitialized
    min_coord[0] = max_coord[0] = 0.0f;
    min_coord[1] = max_coord[1] = 0.0f;
    min_coord[2] = max_coord[2] = 0.0f;
  }

  Real caption_location[3];
  caption_location[0] = (max_coord[0] + min_coord[0])/2;
  caption_location[1] = max_coord[1];
  caption_location[2] = (max_coord[2] + min_coord[2])/2;

  //then caption the image
  midas_file << ".color cyan\n"
             << ".cmov "
             << caption_location[0] << " "
             << caption_location[1] << " "
             << caption_location[2] << "\n";
  midas_file << caption << "\n";



  //Draw "markers" at the location of every "BackboneAtom" in both molecules
  if (show_residue_markers)
  {
    for (int which_protein = 0;
         which_protein < 2;
         ++which_protein)
    {
      midas_file << ".color " << aColorNames[which_protein] << "\n";
      for (Biopolymer::const_iterator pS = apMol[which_protein]->begin();
           pS != apMol[which_protein]->end();
           ++pS)
      {
        for (int q = 0;
             q < Biopolymer::Residue::NumBackboneAtoms();
             ++q)
        {
          Biopolymer::Residue::const_iterator pA = pS->GetBackboneAtom(q);
          if (pA == pS->end())
            ERR("Error: The atoms specified in the \"atoms_used.txt\"\n"
                "       file need to be present in every segment of \n"
                "       the linear molecule (every residue in the\n"
                "       protein, if this is a protein).\n"
                "          One of the segments has one or more missing\n"
                "       atoms.  Aborting...\n");
          midas_file << ".marker "
                     << (*pA).second.xyz[0]
                     << " "
                     << (*pA).second.xyz[1]
                     << " "
                     << (*pA).second.xyz[2]
                     << "\n";
        } //loop over all "backbone" atoms in the segment
      } //loop over all segments of the second protein
    } //loop over both proteins.
  } // if (show_residue_markers)

  if (connect_the_dots && (Biopolymer::Residue::NumBackboneAtoms() > 0))
  {
    int subset_start[2];
    subset_start[0] = subset_start_A;
    subset_start[1] = subset_start_B;

    for (int which_protein = 0;
         which_protein < 2;
         ++which_protein)
    {
      midas_file << ".color " << PairwiseAlignment::aColorNames[which_protein] << "\n";

      for (Biopolymer::const_iterator pS = apMol[which_protein]->begin();
           pS != apMol[which_protein]->end();
           ++pS)
      {
        for (int q = 0;
             q < Biopolymer::Residue::NumBackboneAtoms();
             ++q)
        {
          //Now draw line segments connecting together the "BackboneAtoms"
          //from every segment
          Biopolymer::Residue::const_iterator pA = pS->GetBackboneAtom(q);
          if (pA == pS->end())
            ERR("Error: The atoms specified in the \"atoms_used.txt\"\n"
                "       file need to be present in every segment of \n"
                "       the linear molecule (every residue in the\n"
                "       protein, if this is a protein).\n"
                "          One of the segments has one or more missing\n"
                "       atoms.  Aborting...\n");


          if (pS == apMol[which_protein]->begin())
            midas_file << ".m ";
          else 
            midas_file << ".d ";
          midas_file << (*pA).second.xyz[0]
                     << " "
                     << (*pA).second.xyz[1]
                     << " "
                     << (*pA).second.xyz[2]
                     << "\n";

          if (show_subset_size > 0)
          {
            if ((pS - apMol[which_protein]->begin())
                == subset_start[which_protein])
              midas_file << ".color magenta\n";
            if ((pS - apMol[which_protein]->begin())
                == (subset_start[which_protein] + show_subset_size - 1))
              midas_file << ".color "
                         << PairwiseAlignment::aColorNames[0]
                         << "\n";
          }
        } //loop over all "BackboneAtoms"
      } //loop over all segments
    } //loop over both proteins
  } // if (connect_the_dots)


  if ((show_subset_size>0) && (! connect_the_dots))
  {
    for (int which_protein = 0;
         which_protein < 2;
         ++which_protein)
    {
      //draw line segments in a different color
      //connecting together the first "BackboneAtom" (atom #0)
      //from every segment in the interval

      midas_file << ".color magenta\n";
      int subset_start[2];
      subset_start[0] = subset_start_A;
      subset_start[1] = subset_start_B;
      // Loop over all Residues in the interval
      // [subset_start,subset_start+subset_size]
      for (int i=0; i < show_subset_size; ++i)
      {
        for (int q = 0;
             q < Biopolymer::Residue::NumBackboneAtoms();
             ++q)
        {
          if ((i == 0) && (q == 0))
            midas_file << ".m ";
          else 
            midas_file << ".d ";
          
          assert(subset_start[which_protein] + i <
                 apMol[which_protein]->size());

          Biopolymer::Residue::const_iterator pA = 
            (*apMol[which_protein])
            [subset_start[which_protein]+i].GetBackboneAtom(q);

          if (pA == (*apMol[which_protein])
              [subset_start[which_protein]+i].end())
            ERR("Error: The atoms specified in the \"atoms_used.txt\"\n"
                "       file need to be present in every segment of \n"
                "       the linear molecule (every residue in the\n"
                "       protein, if this is a protein).\n"
                "          One of the segments has one or more missing\n"
                "       atoms.  Aborting...\n");
          midas_file << (*pA).second.xyz[0]
                     << " "
                     << (*pA).second.xyz[1]
                     << " "
                     << (*pA).second.xyz[2]
                     << "\n";
        } // loop over all "BackboneAtoms" in segment
      } // loop over all Residues in the interval [subset_start,subset_start+subset_size]
    } // loop over all molecules
  } // if ((show_subset_size > 0) && (! connect_the_dots))


  //Then (finally), draw the line segments (or markers) indicating a match

  //"sq_marker_dist" stores the square of the distance less-than-which
  //a midas-marker is used to indicate the match.
  Real sq_marker_dist;
  if (use_marker_if_dist_less_than >= 0)
    sq_marker_dist =
      use_marker_if_dist_less_than * use_marker_if_dist_less_than;
  else
    sq_marker_dist = -1.0f; //otherwise, dissable using markers

  midas_file << ".color red\n";

  Biopolymer::const_iterator worst[2];
  if (mark_worst_match)
    FindFurthestAtomPair(worst[0], worst[1], m1, m2, " CA");
  else {
    worst[0] = m1.end();  //This will dissable the behavior of drawing
    worst[1] = m2.end();  //the "worst" match in the alignment in a 
                          //different color.
  }
    
  //Loop over the pairs of matched segments (or "residues", er whatever)
  //in the alignment.
  for(int n=0; n < size(); ++n)
  {
    Vect3 aPos[2];
    int i[2];
    i[0] = aMatches[n][0]; //the nth matched segment from molecule #1
    i[1] = aMatches[n][1]; //the nth matched segment from molecule #2
    assert((i[0]>=0) && (i[1]>=0));

    //Draw line segments connecting all of the individual "BackboneAtoms"
    //from the two matched segments (segment i with segment j).
    int which_pair = 0;
    for(int q = 0; 
        q < Biopolymer::Residue::NumBackboneAtoms();
        ++q)
    {
      Biopolymer::Residue::const_iterator pA[2];

      for(int which_molecule = 0;
          which_molecule < 2;
          ++which_molecule)
      {
        //sequence A
        pA[which_molecule] =
          (*apMol[which_molecule])[i[which_molecule]].GetBackboneAtom(q);
        if (pA[which_molecule] ==
            (*apMol[which_molecule])[i[which_molecule]].end())
          ERR("Error: The atoms specified in the \"atoms_used.txt\"\n"
              "       file need to be present in every segment of \n"
              "       the linear molecule (every residue in the\n"
              "       protein, if this is a protein).\n"
              "          One of the segments has one or more missing\n"
              "       atoms.  Aborting...\n");

        //store the position of qth atom from segment i of molecule "which_molecule"
        aPos[which_molecule][0] = (*pA[which_molecule]).second.xyz[0];
        aPos[which_molecule][1] = (*pA[which_molecule]).second.xyz[1];
        aPos[which_molecule][2] = (*pA[which_molecule]).second.xyz[2];
      }

      if (i[0]+m1.begin() == worst[0])
      {
        assert(i[1]+m2.begin() == worst[1]);
        midas_file << ".color orange\n";
      }

      Vect3 match_vect;
      SubtractVect3(aPos[1],
                    aPos[0],
                    match_vect);
      Real sq_match_dist = DotProduct3( match_vect, match_vect );
    
      //If the distance between matched residues is greater
      //than the square root of sq_marker_dist, just draw
      //a line connecting the two matched residues.
      midas_file << ".m "
                 << aPos[0][0]
                 << " "
                 << aPos[0][1]
                 << " "
                 << aPos[0][2]
                 << "\n";
      midas_file << ".d "
                 << aPos[1][0]
                 << " "
                 << aPos[1][1]
                 << " "
                 << aPos[1][2]
                 << "\n";

      if (sq_match_dist <= sq_marker_dist)
      {
        //Otherwise, put a marker halfway in between the two residues
        Vect3 marker_position;
        Vect3 match_vect_over_2;

        ScalarMultVect3(0.5f, match_vect, match_vect_over_2);

        AddVect3( aPos[0],
                  match_vect_over_2,
                  marker_position );

        midas_file << ".marker "
                   << marker_position[0]
                   << " "
                   << marker_position[1]
                   << " "
                   << marker_position[2]
                   << "\n";
        
      } //else clause for "if (sq_match_dist > sq_marker_dist)"
      
      if (i[0]+m1.begin() == worst[0])
      {
        assert(i[1]+m2.begin() == worst[1]);
        midas_file << ".color red\n";
      }

      ++which_pair;
    } // loop over the "BackboneAtoms" in each segment"for(q=0; q<NumBackboneAtoms()"
  } // loop over matched segment-pAirs "for(n=0; n < size(); ++n)..."
  midas_file << flush;
} //PairwiseAlignment::ExportGFX()




Real PairwiseAlignment::
FindFurthestAtomPair(Biopolymer::const_iterator& worst_i,
                     Biopolymer::const_iterator& worst_j,
                     Biopolymer const& m1,
                     Biopolymer const& m2,
                     string representative_atom_name) const
{
  Real furthest_dist_sqd = -1.0;
  for (int cur_n=0; cur_n < size(); ++cur_n)
  {
    int cur_i = aMatches[cur_n][0];
    int cur_j = aMatches[cur_n][1];
    assert((cur_i >= 0) && (cur_j >= 0));

    Biopolymer::Residue::const_iterator pCAatom1 =
      m1[cur_i].find(representative_atom_name);
    Biopolymer::Residue::const_iterator pCAatom2 =
      m2[cur_j].find(representative_atom_name);
    //if both residues contain alpha carbons, then find the distance
    //between them.
    if ((pCAatom1 != m1[cur_i].end()) && (pCAatom2 != m2[cur_j].end()))
    {
      Vect3 v_ij;
      SubtractVect3((*pCAatom1).second.xyz,
                    (*pCAatom2).second.xyz,
                    v_ij);
      Real d_ij_sqd = DotProduct3(v_ij, v_ij);

      if ((d_ij_sqd > furthest_dist_sqd) || (cur_n == 0))
      {
        furthest_dist_sqd = d_ij_sqd;
        worst_i = m1.begin()+cur_i;
        worst_j = m2.begin()+cur_j;
      }
    }
  }
  if (furthest_dist_sqd == -1.0)
    cerr << "Warning: Cannot find the furthest pair of matched\""
         << representative_atom_name <<
      "\"\n"
      "       atoms in this alignment.\n"
      "       Not a single pair of matched residues in the alignment\n"
      "       both contain a \"" << representative_atom_name << "\" atom.\n";

  return sqrt(furthest_dist_sqd);
}



Real
PairwiseAlignment::
FindFurthestResPair(Biopolymer::const_iterator& worst_i,
                    Biopolymer::const_iterator& worst_j,
                    Biopolymer const& m1,
                    Biopolymer const& m2) const
{
  worst_i = m1.end();
  worst_j = m2.end();
  Real furthest_dist_sqd = 0.0;

  for (int cur_n=0; cur_n < size(); ++cur_n)
  {
    int cur_i = aMatches[cur_n][0];
    int cur_j = aMatches[cur_n][1];
    assert((cur_i >= 0) && (cur_j >= 0));

    Real c_ij = SumSqdDistBetweenResidues(m1.begin()+cur_i,
                                          m2.begin()+cur_j);

    if ((c_ij > furthest_dist_sqd) || (cur_n == 0))
    {
      furthest_dist_sqd = c_ij;;
      worst_i = m1.begin()+cur_i;
      worst_j = m2.begin()+cur_j;
    }
  }
  assert(furthest_dist_sqd >= 0.0);
  return sqrt(furthest_dist_sqd);
} // PairAlignment::FindFurthestResPair()




int NumMatchesInBoth(PairwiseAlignment const& a1,
                     PairwiseAlignment const& a2)
{
  int count = 0;
  for(int n2 = 0; n2 < a2.NumMatches(); ++n2)
  {
    for(int n1 = 0; n1 < a1.NumMatches(); ++n1)
    {
      if (a1[n1] == a2[n2])
        ++count;
    }
  }
  return count;
}


int NumMatchesInEither(PairwiseAlignment const& a1,
                       PairwiseAlignment const& a2)
{
  return a1.size() + a2.size() - NumMatchesInBoth(a1, a2);
}



// (There's probably an easier way to write this function.
//  But the way I did it below obviously works, but it's not very elegant.)
void NumResiduesInBoth(PairwiseAlignment const& a1,
                       PairwiseAlignment const& a2,
                       int& num_bothA,
                       int& num_bothB)
{

  //Figure out the set of residues that got matched in each
  //alignment.  Maintain separate lists for each sequence and
  //for each alignment.
  vector<vector<bool> > residues_used1(2);
  vector<vector<bool> > residues_used2(2);

  //How large should these arrays be?
  //   Since I don't know the size of the sequences I am dealing with,
  //they should be at least as large as the highest numbered residue
  //that was matched in either alignment, for a given sequence.
  //   Since I do not assume the alignments store their matched
  //residues in monotonically increasing order, I have to
  //scan all the matches in both alignments to find this.
  int max_res[2];
  max_res[0] = 0;
  max_res[1] = 0;
  for(int i = 0; i < a1.size(); ++i)
  {
    for (int s = 0; s < 2; ++s)
    {
      if (a1[i][s]+1 > max_res[s])
        max_res[s] = a1[i][s]+1;
    }
  }
  for(int i = 0; i < a2.size(); ++i)
  {
    for (int s = 0; s < 2; ++s)
    {
      if (a2[i][s]+1 > max_res[s])
        max_res[s] = a2[i][s]+1;
    }
  }

  //might as well reserve the memory needed
  residues_used1[0].reserve(max_res[0]);
  residues_used1[1].reserve(max_res[1]);
  residues_used2[0].reserve(max_res[0]);
  residues_used2[1].reserve(max_res[1]);

  //Initialize residues_used to false;
  for(int s = 0; s < 2; ++s)
  {
    for(int i = 0; i < max_res[s]; ++i)
    {
      residues_used1[s].push_back(false);
      residues_used2[s].push_back(false);
    }
  }

  //Flag only the used residues to true.
  for (int s = 0; s < 2; ++s)
  {
    for(int i = 0; i < a1.size(); ++i)
      residues_used1[s][a1[i][s]] = true;
    for(int i = 0; i < a2.size(); ++i)
      residues_used2[s][a2[i][s]] = true;
  }

#if 0
  #ifdef DEBUG
  #if DBG_IMPORT_MSF
  DEBUG_MSG(DBG_IMPORT_MSF, "checking which residues were used in both alignments.\n");
  //Initialize residues_used to false;
  for(int s = 0; s < 2; ++s)
  {
    for(int i = 0; i < max_res[s]; ++i)
    {
      cerr << ((residues_used1[s][i]) ? 1 : 0);
    }
    cerr << "   ";
  }
  cerr << "\n";
  for(int s = 0; s < 2; ++s)
  {
    for(int i = 0; i < max_res[s]; ++i)
    {
      cerr << ((residues_used2[s][i]) ? 1 : 0);
    }
    cerr << "   ";
  }
  #endif //#if DBG_IMPORT_MSF
  #endif //#ifdef DEBUG
#endif //#if 0

  //Now, (finally) count the number if residues used in both alignments.
  //(separately for each sequence)
  int aCount[2] = {0, 0};
  for(int s = 0; s < 2; ++s)
  {
    for(int i = 0; i < max_res[s]; ++i)
    if (residues_used1[s][i] && residues_used2[s][i])
      ++aCount[s];
  }
  num_bothA = aCount[0];
  num_bothB = aCount[1];
}


void NumResiduesInEither(PairwiseAlignment const& a1,
                         PairwiseAlignment const& a2,
                         int& num_eitherA,
                         int& num_eitherB)
{
  int num_bothA;
  int num_bothB;
  NumResiduesInBoth(a1, a2, num_bothA, num_bothB);
  num_eitherA = a1.size() + a2.size() - num_bothA;
  num_eitherB = a1.size() + a2.size() - num_bothB;
}




int NumMatchesInBothAssumingMonotonicity(PairwiseAlignment const& a1,
                                         PairwiseAlignment const& a2)
{
  assert(a1.IsMonotonic());
  assert(a2.IsMonotonic());

  int count = 0;
  int n1 = 0; //for indexing through a1
  int n2 = 0; //for indexing through a2
  while (n2 < a2.NumMatches())
  {
    while((n1 < a1.NumMatches()) && (a1[n1][0] < a2[n2][0]))
      ++n1;
    if ((n1 < a1.NumMatches()) && (a1[n1] == a2[n2]))
      ++count;
    ++n2;
  }

  #if 0
  assert((a1.NumMatches() > 500) //We don't want to test this assertion
         ||                      //if the alginments are too large
         (a2.NumMatches() > 500) //because NumIdenticalMatches() is very slow.
         ||
         (count == NumIdenticalMatches(a1, a2)));//<-This is what I wanted to assert
  #endif //#if 0
  return count;
}


int NumMatchesInEitherAssumingMonotonicity(PairwiseAlignment const& a1,
                                           PairwiseAlignment const& a2)
{
  return a1.size() + a2.size() -
    NumMatchesInBothAssumingMonotonicity(a1, a2);
}




static Real Compare3dSubset(char mode,
                            //'A' means all residues in mB are used
                            //'B' means only residues matched in both
                            //    alignments are used.
                            //'E' means only residues matched in either
                            //    alignment are used.
                            PairwiseAlignment const& a1,
                            PairwiseAlignment const& a2,
                            Biopolymer const& mA,
                            Biopolymer const& mB);
//This function is defined below.


Real Compare3dAll(PairwiseAlignment const& a1,
                  PairwiseAlignment const& a2,
                  Biopolymer const& mA,
                  Biopolymer const& mB)
{
  return Compare3dSubset('A', a1, a2, mA, mB);
}

Real Compare3dBoth(PairwiseAlignment const& a1,
                   PairwiseAlignment const& a2,
                   Biopolymer const& mA,
                   Biopolymer const& mB)
{
  return Compare3dSubset('B', a1, a2, mA, mB);
}


Real Compare3dEither(PairwiseAlignment const& a1,
                     PairwiseAlignment const& a2,
                     Biopolymer const& mA,
                     Biopolymer const& mB)
{
  return Compare3dSubset('E', a1, a2, mA, mB);
}



static Real Compare3dSubset(char mode,
                            //'A' means all residues in mB are used
                            //'B' means only residues matched in both
                            //    alignments are used.
                            //'E' means only residues matched in either
                            //    alignment are used.
                            PairwiseAlignment const& a1,
                            PairwiseAlignment const& a2,
                            Biopolymer const& mA,
                            Biopolymer const& mB)
{
  Matrix3x4 o1, o2;

  a1.CalcMinRMSD(mA, mB, o1); //orient mB two different ways (o1 & o2) to
  a2.CalcMinRMSD(mA, mB, o2); //minimize the RMSD between alignments a1 & a2.

  int max_num_atoms_B_possible =
    mB.size() * Biopolymer::Residue::NumBackboneAtoms();
  int  num_atoms_B = 0; //The number of atoms from structure mB, whose
  Vect3 *atoms_from_B   //position will be compared at orientations: o1 & o2.
    = new Vect3[max_num_atoms_B_possible];


  //Figure out the set of residues from structure mB that got matched in each
  //alignment: a1 & a2.
  vector<bool> residues_usedB1(mB.size(), false); //initialize to false
  vector<bool> residues_usedB2(mB.size(), false); //initialize to false
  //Flag only the used residues to true.
  for(int i = 0; i < a1.size(); ++i)
    residues_usedB1[a1[i][1]] = true;
  for(int i = 0; i < a2.size(); ++i)
    residues_usedB2[a2[i][1]] = true;

  //Copy the positions of only the atoms from structure mB which
  //which were matched in either alignment a1 or a2
  //into the array "atoms_from_B".
  PairwiseAlignment::const_iterator p1 = a1.begin();
  PairwiseAlignment::const_iterator p2 = a2.begin();
  cerr << "Biopolymer::Residue::NumBackboneAtoms() = "
       << Biopolymer::Residue::NumBackboneAtoms() << endl;
  for(int r=0; r < mB.size(); ++r)
  {
    if ((mode == 'A') 
        ||
        ((mode == 'B') && (residues_usedB1[r] && residues_usedB2[r]))
        ||
        ((mode == 'E') && (residues_usedB1[r] || residues_usedB2[r])))
    {
      for(int q = 0; 
          q < Biopolymer::Residue::NumBackboneAtoms();
          ++q)
      {
        Biopolymer::Residue::const_iterator pa = mB[r].GetBackboneAtom(q);
        if (pa != mB[r].end()) //If this residue has a backboneatom of this type
        {
          atoms_from_B[num_atoms_B][0] = (*pa).second.xyz[0];
          atoms_from_B[num_atoms_B][1] = (*pa).second.xyz[1];
          atoms_from_B[num_atoms_B][2] = (*pa).second.xyz[2];
          ++num_atoms_B;
        }
      }
    } //if residue i is part of match *p1 or *p2, copy
  } //for(int i=0; i < mB.size(); ++i)
  if (num_atoms_B == 0)
    ERR("Error: When calculating RMSD between two positions of the\n"
        "       same structure, an error occured.\n"
        "       The criteria that chooses which atoms belong to the\n"
        "       subset of atoms the RMSD represnts,\n"
        "       has filtered out all of the atoms.\n"
        "       This probably means you are considering only positions\n"
        "       of residues matched in BOTH alignments, and there simply\n"
        "       weren't any (from the second structure, at least).\n"
        "       Aborting.\n");
  Real sum_sqd_dist = SlowRotSumSqdDist(o1, o2, atoms_from_B, num_atoms_B);
  delete [] atoms_from_B;
  assert(sum_sqd_dist >= 0.0);
  return sqrt(sum_sqd_dist / num_atoms_B);
} //Compare3dEither()








Real LevittGerstein98::Sstr_Score(PairwiseAlignment const& a,
                                  Biopolymer const& m1,
                                  Biopolymer const& m2)
{
  int num_gaps = 0;
  int num_gap_extensions = 0;
  static const Real M = 20.0f;
  static const Real d_o = 5.0f;
  static const Real d_o_sqrd = d_o*d_o;
  static const Real GAP_OPENING_PENALTY = 10.0f;
  static const Real GAP_EXTENSION_PENALTY = 0.0f;
  Biopolymer const *(apMol[2]);
  apMol[0] = &m1;
  apMol[1] = &m2;

  Real total_score = 0.0f;

  int prev_i = -1; //(used for keeping track of gap-sizes
  int prev_j = -1; // to assign gap-penalties)
  
  //Loop over the pairs of matched segments (or "residues", er whatever)
  //in the alignment.
  for(int n=0; n < a.NumMatches(); ++n)
  {
    Vect3 pos1;
    Vect3 pos2;
    int i = a[n][0]; //the nth matched segment from molecule #1
    int j = a[n][1]; //the nth matched segment from molecule #2
    assert((i>=0) && (i < (apMol[0])->size()));
    assert((j>=0) && (j < (apMol[1])->size()));
    Biopolymer::Residue::const_iterator pA[2];
    //sequence A
    pA[0] = (*apMol[0])[i].find(" CA");
    //sequence B
    pA[1] = (*apMol[1])[j].find(" CA");
    
    //Check for missing C-alhpa:
    for (int m = 0; m < 2; ++m)
    {
      if (pA[m] == (*apMol[m])[i].end()) {
        cerr << 
          "Error:  Missing alpha-carbon in residue " << (*apMol[m])[i].id << "\n"
          "        of protein #" << m+1 << ".\n"
          "    In order to evaluate an alignment using Levitt & Gerstein's\n"
          "    scoring criteria, each residue involved in the match must\n"
          "    contain a C-alpha. \n" << flush;
        return 0.0f;
      }
    }

    pos1[0] = (*pA[0]).second.xyz[0];//position of qth atom from
    pos1[1] = (*pA[0]).second.xyz[1];//segment i of molecule#1
    pos1[2] = (*pA[0]).second.xyz[2];
    pos2[0] = (*pA[1]).second.xyz[0];//position of qth atom from
    pos2[1] = (*pA[1]).second.xyz[1];//segment i of molecule#1
    pos2[2] = (*pA[1]).second.xyz[2];

    Vect3 match_vect;
    SubtractVect3( pos1,
                   pos2,
                   match_vect );
    Real d_ij_sqrd = DotProduct3( match_vect, match_vect );
    
    total_score += M / (1 + (d_ij_sqrd/d_o_sqrd));

    //------- gap-penalties -------
    //Gap penalties are assigned according to Gotoh's variant
    //of the classic dynamic programming (Needleman and Wunsch)
    //technique.  This seems to be the definition of 
    //"opening and extension gap-penalties" that Gerstein and Levitt use.
    //  (see Gotoh, O. An improved algorithm for matching biological sequences,
    //   J. Mol. Biol., v.162, p.705, 1982)

    //total_gap_size is the cululative gap size in both sequences
    //(i.e. the sum of the gap length sequence 1 and sequence 2 counting
    // from this pair of residues, to the previous pair of residues.
    // If this pair of residues fall right after the previous pair, this is 0)
    assert(i > prev_i);
    assert(j > prev_j);
    int total_gap_size = (i - prev_i - 1) + (j - prev_j - 1);
    if ((total_gap_size > 0) && (prev_i != -1))
    {
      assert(prev_j != -1);
      total_score -= (GAP_OPENING_PENALTY +
                      (total_gap_size - 1)*GAP_EXTENSION_PENALTY);
      num_gap_extensions += (total_gap_size - 1);
      ++num_gaps;
    }

    prev_i = i;
    prev_j = j;
  } // loop over matched segment-pairs "for(n=0; n < size(); ++n)..."
#if 0
  cerr << "num_new_gaps = " << num_gaps << endl;
  cerr << "num_gap_extensions = " << num_gap_extensions << endl;
#endif //#if 0

  return total_score;
} // Real LevittGerstein98::Sstr_Score()



// "There are Lies,
//  there are Damn Lies,
//  ...and then there's statistics."

Real LevittGerstein98::Sstr_Probability(PairwiseAlignment const& al,
                                        Biopolymer const& m1,
                                        Biopolymer const& m2,
                                        //Optionally, the caller
                                        //can request to have the
                                        //following intermediate
                                        //calculations returned
                                        //to him/her, by passing
                                        //the address of a Real
                                        //variable for each
                                        //of these arguments.
                                        Real *p_str,
                                        Real *p_mu_str,
                                        Real *p_delta_str,
                                        Real *pZ)
{
  static const long double c = 18.411424;
  static const long double d = -4.501719;
  static const long double e = 2.637163;
  static const long double f = 21.351857;
  static const long double g = -37.521269;
  static long double ln120 = log(120.0);
  static long double a = 2.0*ln120*c + d;
  static long double b = SQR(ln120)*c + ln120*d + e - a*ln120;

  //The lines of code above were meant to initialize
  //the variables from Gerstein & Levitt's paper.
  //Unfortunately, the values of these variables from that paper are
  //heavilly rounded off, and are thus not as accurate
  //as we would have liked.  Instead, we copied the 
  //values out of Gerstein & Levitt's code directly
  //which we downloaded from their structural alignment server
  //on the web.  The following is an excerpt from that code:
  // $ln60 = log(120);
  // @para = ( 18.411424, -4.501719, 2.637163, 21.351857, -37.521269);
  // $bpar = 2.0*$ln60*$para[0] + $para[1]; 
  // $apar = ($ln60**2)*$para[0] + $ln60*$para[1] + $para[2] - $bpar*$ln60;


  long double mu_str;
  long double delta_str;
  long double lnN = log((float)(al.NumMatches()));
  long double S_str = Sstr_Score(al, m1, m2);

  if (al.NumMatches() < 120)
  {
    //for debugging, I split up the terms in the polynomial.
    #ifdef COMPILER_BUG1
    mu_str =  c * minrms::SQR(lnN);
    #else
    mu_str =  c * SQR(lnN);
    #endif
    mu_str += d * lnN;
    mu_str += e;

    delta_str =  f * lnN;
    delta_str += g;
  }
  else
  {
    mu_str =  a * lnN;
    mu_str += b;

    delta_str = f * ln120 + g;
  }

  #ifdef COMPILER_BUG1
  using namespace minrms;
  #endif

  if (p_str)
    *p_str = S_str;

  if (p_mu_str)
    *p_mu_str = mu_str;

  if (p_delta_str)
    *p_delta_str = delta_str;

  long double Z = (S_str - mu_str) / delta_str;

  if (pZ)
    *pZ = Z;

  //                         Approximation:
  //What we want to calculate is:   1 - exp(-exp(-Z))
  //But, in the high Z limit, we must use an approximation because
  //roundoff error will replace 1-exp(-exp(-Z)) by 0 for large values of Z.
  //
  // We approximate  1 - exp(-exp(-x))  by  exp(-x) for large values of x.
  //
  //        Error analysis:
  //An upper bound for the magnitude of the error can be found by
  //taylor expanding 1 - exp(-exp(-x))  in terms of: exp(-x)
  //1 - exp(-exp(-x))  = 1 - SumOverN { (-exp(-x))^N / N! }
  //                   = 1 - (1 - exp(-x) + (1/2)exp(-2x) + higher order terms)
  //                   ~ exp(-x) + (1/2)exp(-2x)
  //As x increases, exp(-x) approaches zero, and we can throw the
  //higher-order terms away.
  //The difference between the two expressions is therefore:
  // (1/2)exp(-2x)
  //The FRACTIONAL difference between the two expressions is:
  // error / value =  (1/2)exp(-2x) / exp(-x)   =  (1/2)exp(-x)
  //This is more useful because it tells us how large the error is
  //compared to the actual value being calculated.
  //
  //So, to calculate 1-exp(-exp(-Z)) to an accuracy of 0.01% (= 1/10000),
  //or better (ignoring floating point roundoff error),
  //we only switch to exp(-Z) if Z exceeds ln(10000/2)  (~ 8.5).
  // (How does this compare with the floating point error?
  //  The floating point roundoff error for a 32 bit integer with a
  //  23 bit mantissa is fixed at  +/- 1/2^23 after subtraction.
  //  This error is, in turn, about is 0.0006 (= 0.06%)
  //  of the calculated value (exp(-Z))
  //  when Z = ln(10000/2).
  //  So we get no benifit from using a larger value of Z.
  //  The total fractional error is the larger of the two: 0.06% in this case.)

  long double result;
  if (Z < log(10000.0/2))
    result = 1.0f - exp(-exp(-Z));
  else
    result = exp(-Z);

  return (Real)result;
} //Real LevittGerstein98::Sstr_Probability()








void TranslateAlignmentBetweenMolecules(PairwiseAlignment const& alOrig,
                                        PairwiseAlignment & alNew,
                                        Biopolymer const& mOrig1,
                                        Biopolymer const& mOrig2,
                                        Biopolymer const& mNew1,
                                        Biopolymer const& mNew2)
{
  assert(alOrig.NumMatches() == alNew.NumMatches());
  int m; //indexes through the matched residue-pairs in the alignment
  for (m=0; m < alOrig.NumMatches(); ++m)
  {
    
    int locOrig1 = alOrig[m][0];
    int locOrig2 = alOrig[m][1];
    PDBresID segID1 = mOrig1[locOrig1].id;
    PDBresID segID2 = mOrig2[locOrig2].id;
    Biopolymer::const_iterator locNew1 = mNew1.find(segID1);
    Biopolymer::const_iterator locNew2 = mNew2.find(segID2);
    assert(locNew1 != mNew1.end());
    assert(locNew2 != mNew2.end());
    assert(locNew1->id == segID1);
    assert(locNew2->id == segID2);
    
    alNew[m][0] = locNew1 - mNew1.begin(); //store the new integer index
    alNew[m][1] = locNew2 - mNew2.begin(); //into the mNew1/2 molecules
  }
} // TranslateAlignmentBetweenMolecules()

#ifdef CREATE_PAIR_HISTOGRAM
void
PairwiseAlignment::AccumPairHistogram(long **aaHistogram) const
{
  assert(aaHistogram);
  int cur_n, cur_i, cur_j;
  for(cur_n=0; cur_n < size(); ++cur_n)
  {
    cur_i = aMatches[cur_n][0];
    cur_j = aMatches[cur_n][1];
    assert(aaHistogram[cur_i]);
    ++(aaHistogram[cur_i][cur_j]);
  }
}
#endif



} //namespace minrms



