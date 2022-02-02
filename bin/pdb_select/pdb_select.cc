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


#include <iostream>
#include <sstream>
using namespace std;

#include <biopolymer.h>
#include <load_pdb.h>
#include <intervals.h>
#include <parse_set.h>
#include <parse_utils.h>
#include <pdb_id.h>
#include <pdb++.h>

const char g_version_string[] = "0.991";
const char g_date_string[]    = "<00/2/7>";

using namespace minrms;

using namespace ParseUtils;

// LoadStructure takes 3 arguments:
// "input_stream" specifies where the PDB-file will be read from.
// "m_selection"  is a Biopolymer containing the residues that were
//                selected from the PDB file loaded from the "input_stream".
// "subset_str"   is a Midas/Chimera formatted string indicating the
//                residues the user would like to select from the PDB-file.
// The output of this function is the Biopolymer data stored in
// the "m_selection" variable.)
//
// (Note: most of this code was pillaged from the
//  PairAlignSettingsParser::LoadStructurers() member function located in:
//  "../minrms/minrms_parser.cc")


void LoadStructure(istream &input_stream,
                   Biopolymer& m_selection,
                   string subset_str = ".")  //("." selects all the residues)
{
  //The first step is to create the Biopolymer variable "m".
  Biopolymer m;
  input_stream >> m;
  m.Finalize();

  if (m.size() == 0)
    ERR("Error: Format error or PDB-file contains no residues.");

  // *** Okay, now that the original structure has
  // *** been loaded into m,
  // *** filter out the all but the residues that are desired
  // *** for use in the calculation.
  // *** Store this subset of residues in m_selection

  //Now, parse the string containing the desired set of residues.
  vector<ClosedInterval> desired_residues;
  string::const_iterator p = subset_str.begin();
  while (p != subset_str.end())
  {
    vector<ClosedInterval> simple_set;
    ParseChimeraSet(p,
                    subset_str.end(),
                    m,
                    simple_set,
                    ",");

    desired_residues.insert(desired_residues.end(),
                            simple_set.begin(), simple_set.end());
  }
  if (desired_residues.size() == 0)
  {
    ClosedInterval whole_molecule;
    whole_molecule.first = m.front().id;
    whole_molecule.last  = m.back().id;
    desired_residues.push_back(whole_molecule);
  }
  else
  {
    //For debugging purposes, print the set of residues selected.
    cerr << "Residues selected from from structure:" << endl;
    for(vector<ClosedInterval>::const_iterator p = desired_residues.begin();
        p != desired_residues.end();
        ++p)
    {
      cerr << "[" << (*p).first << " - " << (*p).last << "]";
      if (p+1 != desired_residues.end())
        cerr << ", ";
    }
    cerr << endl;
  }

  //We want m_selection to contain a subset of the residues in m.
  //We start by copying all of the residues from m,
  //and delete the ones we don't want later.
  m_selection = m;

  //Make sure the sets we have selected are valid.
  ClosedInterval invalid_interval;
  if (DeleteInvalidIntervals(desired_residues,
                             m_selection,
                             &invalid_interval))
  {
    ERR("Error: Error in format for set selection\n"
        " One of these residues does not belong to the structure:\n"
        << invalid_interval.first.seqNum << invalid_interval.first.insertCode <<
        " from chain " << invalid_interval.first.chainId <<
        ", or\n"
        << invalid_interval.last.seqNum << invalid_interval.last.insertCode <<
        " from chain " << invalid_interval.last.chainId);
  }

  //Now throw out the residues in m_selection that do not
  //belong to the "desired_residues" set.
  ClipMoleculeToListOfIntervals(m_selection,
                                desired_residues);

  //Bullet proofing:
  //Check to make sure the structures are not empty:
  if (m_selection.size() == 0)
    ERR("Error: No residues left in structure after filtering."
        "\"\n"
        "       Either the PDB file contains no residues, or\n"
        "       you have selected a subset of residues from within\n"
        "       this structure that is empty.");
  cerr <<
    "(before filtering) length = "
       << m.size() << " residues.\n"
    " (after filtering) length = "
         << m_selection.size() << " residues." << endl;

} //LoadStructure()




int main(int argc, char **argv)
{
  // ****** Input ******

  // Check syntax
  if (argc != 2)
  {
    cerr << 
      "Error: Wrong number of arguments passed.\n"
      "\n"
      "\n"
      "Syntax:\n"
      "       pdb_select set1,set2,set3,... < orig_pdb_file > new_pdb_file \n"
      "Overview:\n"
      "       This command creates a new_pdb_file containing a subset of\n"
      "       residues from the orig_pdb_file.  The only argument to\n"
      "       pdb_select is a comma-separated list (no spaces) of sets\n"
      "       of residues in Midas/Chimera/MinRMS format. (The selection syntax\n"
      "       is described below, as well as in the Midas/Chimera documentation).\n"
      "       \"new_pdb_file\" will contain residues belonging to the union\n"
      "       of all the sets (set1, set2, set3, etc.).\n"
      "\n"
      "Examples of selection syntax:\n"
      "\n"
      "   set           residues selected:\n"
      " --------       ---------------------\n"
      " \"100\"         all residues whose seqNum is 100 in all chains\n"
      " \"100-150\"     all residues between 100 and 150 in all chains\n"
      " \"100-150.\"       \"      \"   \"    \"    \"   \"   \"   \"    \"\n"
      " \"100-150.*\"      \"      \"   \"    \"    \"   \"   \"   \"    \"\n"
      " \"100A-150\"    all res. between 100A (insert-code:A) and 150 in all chains\n"
      " \"100-150.A\"   all residues between 100 and 150 in chain A\n"
      "\"*-150.A\"      all the residues that come before (and include) 150 in chain A\n"
      "\"150-*.A\"      all the residues that come after (and include) 150 in chain A\n"
      " \"*.A\"         all the residues in chain A\n"
      " \".A\"           \"   \"     \"     \"    \"   \"\n"
      " \"100-150.A-C\" residues 100-150 in chains A through C\n"
      " \"*.*\"         the entire molecule\n"
      " \".\"            \"    \"      \"\n"
      "\n"
      "Notes: If you use the '*' character in any of your sets,\n"
      "       you will have to enclose the first argument in quotes to\n"
      "       circumvent the shell.\n"
      "\n"
      "Details:\n"
      "      -Only the ATOM, HETATM, ANISOU, HELIX, SHEET, and TURN records\n"
      "       are effected.  All other records in the PDB file are blindly\n"
      "       sent to the standard output.\n"
      "      -Helices, sheets, or turns which lie all, or partially\n"
      "       outside the selected sets of residues will be deleted."
      << endl;
    exit(1);
  }

  string subset_str(argv[1]); //Specifies the residues desired.

  //Load the input file into the "input_stream" variable.
  //This stream is a variable in memory that will be read twice.
  //(It's more efficient to do this than to read from the actual
  // file twice.)

  stringstream input_stream;  //read in the entire standard input
  stringstream input_stream_2nd_pass;
  
  char c;
  while (cin.get(c))
  {
    input_stream.put(c);
    input_stream_2nd_pass.put(c);
  }

  // I've allready written a function to parse intervals
  //of residues, but it requires that the PDB file has
  //allready been loaded into a Biopolymer data structure first.
  Biopolymer m_selection;

  //After the next call, m_selection should store only the residues desired.
  LoadStructure(input_stream, m_selection, subset_str);


  // ****** Output ******

  //Now, I will read in the PDB file again, and only output the
  //residues and helix/sheet records if they lie inside the selected set.
  //(As indicated by the residues present in the "m_selection" variable.
  // This is kinda wasteful, because I load the PDB file twice.)
  PDB pdb_record;

  while (input_stream_2nd_pass >> pdb_record)
  {
    char c = input_stream_2nd_pass.get(); // These two lines were necessary at one point,
    input_stream_2nd_pass.putback(c);     // but I don't remember if they still are.

    //Fix some problems with older file formats.
    //(Some of the old HELIX/SHEET/TURN records are not parsed correctly.
    // see "../../lib/linear_molecule/interval.cc")
    ConvertOldHelixSheetTurnRecordsToNewFormat(pdb_record, pdb_record);

    PDBresID res;
    ClosedInterval interval;

    switch (pdb_record.type())
    {
    case PDB::ATOM:
      res.chainId     = pdb_record.atom.residue.chainId;
      res.seqNum      = pdb_record.atom.residue.seqNum;
      res.insertCode  = pdb_record.atom.residue.insertCode;
      if (m_selection.find(res) != m_selection.end())
        cout << pdb_record << endl;
      break;
    case PDB::HETATM:
      res.chainId     = pdb_record.hetatm.residue.chainId;
      res.seqNum      = pdb_record.hetatm.residue.seqNum;
      res.insertCode  = pdb_record.hetatm.residue.insertCode;
      if (m_selection.find(res) != m_selection.end())
        cout << pdb_record << endl;
      break;
    case PDB::ANISOU:
      res.chainId     = pdb_record.anisou.residue.chainId;
      res.seqNum      = pdb_record.anisou.residue.seqNum;
      res.insertCode  = pdb_record.anisou.residue.insertCode;
      if (m_selection.find(res) != m_selection.end())
        cout << pdb_record << endl;
      break;

    case PDB::HELIX:
      interval.first.chainId    = pdb_record.helix.residues[0].chainId;
      interval.first.seqNum     = pdb_record.helix.residues[0].seqNum;
      interval.first.insertCode = pdb_record.helix.residues[0].insertCode;
      interval.last.chainId     = pdb_record.helix.residues[1].chainId;
      interval.last.seqNum      = pdb_record.helix.residues[1].seqNum;
      interval.last.insertCode  = pdb_record.helix.residues[1].insertCode;

      DEBUG_MSG(DBG_HELIX_SHEET, __FILE__ << ":" << __LINE__
                << "Helix found from "
                << interval.first << " to "
                << interval.last);

      if ((m_selection.find(interval.first) != m_selection.end()) &&
          (m_selection.find(interval.last)  != m_selection.end()))
        cout << pdb_record << endl;
      else if (m_selection.find(interval.first) !=
               m_selection.find(interval.last))
      {
        cerr <<
          "Warning: Helix [" << interval.first << " - " << interval.last << "]"
          "         contains some residues that were selected, and others\n"
          "         that were not selected.  (This helix will be deleted.)"
          << endl;
      }
      //else, do not include the helix, and do not print a warning.
      break;

    case PDB::SHEET:
      interval.first.chainId    = pdb_record.sheet.residues[0].chainId;
      interval.first.seqNum     = pdb_record.sheet.residues[0].seqNum;
      interval.first.insertCode = pdb_record.sheet.residues[0].insertCode;
      interval.last.chainId     = pdb_record.sheet.residues[1].chainId;
      interval.last.seqNum      = pdb_record.sheet.residues[1].seqNum;
      interval.last.insertCode  = pdb_record.sheet.residues[1].insertCode;

      DEBUG_MSG(DBG_HELIX_SHEET, __FILE__ << ":" << __LINE__
                << "Strand-of-sheet found from "
                << interval.first << " to "
                << interval.last);

      if ((m_selection.find(interval.first) != m_selection.end()) &&
          (m_selection.find(interval.last)  != m_selection.end()))
        cout << pdb_record << endl;
      else if (m_selection.find(interval.first) !=
               m_selection.find(interval.last))
      {
        cerr <<
          "Warning: Sheet-strand [" << interval.first << " - " << interval.last << "]"
          "         contains some residues that were selected, and others\n"
          "         that were not selected.  (This sheet-record will be deleted.)"
          << endl;
      }
      //else, do not include the sheet, and do not print a warning.

      break;
    case PDB::TURN:
      interval.first.chainId    = pdb_record.turn.residues[0].chainId;
      interval.first.seqNum     = pdb_record.turn.residues[0].seqNum;
      interval.first.insertCode = pdb_record.turn.residues[0].insertCode;
      interval.last.chainId     = pdb_record.turn.residues[1].chainId;
      interval.last.seqNum      = pdb_record.turn.residues[1].seqNum;
      interval.last.insertCode  = pdb_record.turn.residues[1].insertCode;

      DEBUG_MSG(DBG_HELIX_SHEET, __FILE__ << ":" << __LINE__
                << "Turn found from "
                << interval.first << " to "
                << interval.last);


      if ((m_selection.find(interval.first) != m_selection.end()) &&
          (m_selection.find(interval.last)  != m_selection.end()))
        cout << pdb_record << endl;
      else if (m_selection.find(interval.first) !=
               m_selection.find(interval.last))
      {
        cerr <<
          "Warning: turn [" << interval.first << " - " << interval.last << "]"
          "         contains some residues that were selected, and others\n"
          "         that were not selected.  (This turn-record will be deleted.)"
          << endl;
      }
      //else, do not include the turn record, and do not print a warning.
      break;

    default:
      cout << pdb_record << endl;
      break;

    } //switch (pdb_record.type())
  } // while (input_stream_2nd_pass >> pdb_record)
} //main()


