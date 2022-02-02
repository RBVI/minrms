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


#include <cstring> //needed for call to strncmp()
#include <cassert> //defines "assert()"
#include <vector>
#include <fstream>
#include <algorithm> //needed for call to STL's "binary_search()"
using namespace std;

#include <pdb++.h> //Required for "ConvertHelixSheetRecordsToIntervals()"

#include "intervals.h"





//Why not have an RCS-version-string?
//static char rcs_id[] = "@(#)$Header: /home/jewett/babar/src/lib/biopolymer/RCS/intervals.cc,v 1.2 2005/05/19 21:11:34 jewett Exp $";


namespace minrms
{


bool
ClipMoleculeToListOfIntervals(Biopolymer& m,
                              vector<ClosedInterval> &filters)
{
  vector<bool> delete_these(m.size(), false);
  for (long i = 0; i != m.size(); ++i)
  {
    if (! BelongsToUnionOfIntervals(m[i].id, filters, m))
      delete_these[i] = true;
  }
  
  bool clipped = false;

  {
    long i = 0;
    for (Biopolymer::iterator pR = m.begin();
         pR != m.end();
         ++pR)
    {
      if (delete_these[i])
      {
        DEBUG_MSG( DBG_INTERVAL_LISTS, "Residue "
                   << pR->id
                   << " got discarded.");
        --pR; //decrement to avoid invalidating the iterator upon "erase()"
        m.erase(pR+1); //discard the offending residue
        clipped = true;
      }
      ++i;
    }
  }

  if (clipped)
    m.Finalize(); //We must call this, since we modified the contents of m.

  return clipped;
} //ClipMoleculeToListOfIntervals()









PDB::RecordType
ConvertHelixSheetTurnRecordToInterval(char const *line,
                                      ClosedInterval& interval)
{
  assert(line);
  PDB pdb_record = PDB(line);

  //Fix some problems with older file formats.
  ConvertOldHelixSheetTurnRecordsToNewFormat(pdb_record, pdb_record);

  switch (pdb_record.type())
  {
  case PDB::UNKNOWN:
    break;
  case PDB::HELIX:
    interval.first.chainId    = pdb_record.helix.residues[0].chainId;
    interval.first.seqNum     = pdb_record.helix.residues[0].seqNum;
    interval.first.insertCode = pdb_record.helix.residues[0].insertCode;
    interval.last.chainId      = pdb_record.helix.residues[1].chainId;
    interval.last.seqNum       = pdb_record.helix.residues[1].seqNum;
    interval.last.insertCode   = pdb_record.helix.residues[1].insertCode;

    DEBUG_MSG(DBG_HELIX_SHEET, __FILE__ << ":" << __LINE__
              << "Helix found from "
              << interval.first << " to "
              << interval.last);
    break;
  case PDB::SHEET:
    interval.first.chainId    = pdb_record.sheet.residues[0].chainId;
    interval.first.seqNum     = pdb_record.sheet.residues[0].seqNum;
    interval.first.insertCode = pdb_record.sheet.residues[0].insertCode;
    interval.last.chainId      = pdb_record.sheet.residues[1].chainId;
    interval.last.seqNum       = pdb_record.sheet.residues[1].seqNum;
    interval.last.insertCode   = pdb_record.sheet.residues[1].insertCode;

    DEBUG_MSG(DBG_HELIX_SHEET, __FILE__ << ":" << __LINE__
              << "Strand-of-sheet found from "
              << interval.first << " to "
              << interval.last);
    break;
  case PDB::TURN:
    interval.first.chainId    = pdb_record.turn.residues[0].chainId;
    interval.first.seqNum     = pdb_record.turn.residues[0].seqNum;
    interval.first.insertCode = pdb_record.turn.residues[0].insertCode;
    interval.last.chainId      = pdb_record.turn.residues[1].chainId;
    interval.last.seqNum       = pdb_record.turn.residues[1].seqNum;
    interval.last.insertCode   = pdb_record.turn.residues[1].insertCode;

    DEBUG_MSG(DBG_HELIX_SHEET, __FILE__ << ":" << __LINE__
              << "Turn found from "
              << interval.first << " to "
              << interval.last);
    break;
  default:
    break;
  } //switch (pdb_record.type())

  return pdb_record.type();
} //ConvertHelixSheetTurnRecordToInterval()




void
ConvertOldHelixSheetTurnRecordsToNewFormat(PDB const& source_record,
                                           PDB& dest_record)
{
  //The current version of PDB.h I am using does not understand the
  //old-style helix,sheet,and turn-records.  If such a record is read,
  //it will be read as type PDB::UNKOWN.
  //The "case PDB::UNKNOWN" handler checks for
  //this condition and tries to correct for it.

  //first, copy source_record into dest_record.
  dest_record = source_record;

  if (source_record.type() != PDB::UNKNOWN)
    return; //then there are no changes to be made, we're done.

  //otherwise, examine the 'junk' field in the record
  char *line = dest_record.unknown.junk;
  assert(line);

  if ((strncmp(line, "HELIX", 5) == 0) ||
      (strncmp(line, "SHEET", 5) == 0) ||
      (strncmp(line, "TURN", 4) == 0))
  {
    bool hack_failed = false;

    //We can correct for the problem
    //f extra characters at the end are responsible
    //for causing this record to come up as type PDB::UNKOWN.
    //If there aren't any, extra characters, then we're screwed, so we give up.
    if (strlen(line) < 72)
      hack_failed = true; //no extra characters?
    else
    {
      line[71] = '\0'; //Truncate away the extra characters at the
                       //end.  Hopefully this should fool PDB.h into
                       //believing this is a valid HELIX,SHEET,TURN record.

      //Now using the truncated line, overwrite "source_record".
      //Try to trick it into believing this is either
      //a HELIX, SHEET, or TURN record.
      dest_record = PDB(line); //truncated version of line
      //In case this hack fails, print a warning message and abort.
      hack_failed =
        ((dest_record.type() != PDB::HELIX) &&
         (dest_record.type() != PDB::SHEET) &&
         (dest_record.type() != PDB::TURN));
    }
    if (hack_failed)
    {
      ERR("ERROR: What appears to be a HELIX,SHEET,or TURN record of an\n"
          "       unsupported format was found in one of the PDB files:\n"
          "       \"" << source_record << "\"\n"
          "Aborting...");
    }
  }
  //Else, this really is just an unknown record.  Do not change it.
} //ConvertOldHelixSheetTurnRecordsToNewFormat()



bool ExtendInterval(ClosedInterval& interval,
                    Biopolymer const& mol,
                    int min_interval_size)
{
  bool extended = false;

  Biopolymer::const_iterator st = mol.find(interval.first);
  if (st == mol.end())
    ERR("ERROR: \"" << __FILE__ << "\":" << __LINE__ << "\n"
        "       Residue " << interval.last << "\n"
        "       not found in molecule.");

  Biopolymer::const_iterator en = mol.find(interval.last);
  if (en == mol.end())
    ERR("ERROR: \"" << __FILE__ << "\":" << __LINE__ << "\n"
        "       Residue " << interval.last << "\n"
        "       not found in molecule");
  
  if (en < st)
    ERR("Subset-Match Interval cannot have start > end.");

  // If an interval is too small, then print
  // a warning to the user, but do not exit.
  if (st > en-min_interval_size+1)
  {
    cerr << "-----------------------------------------------------------\n"
      "Note: helix/sheet interval ["
         << st->id << " - " << en->id
         << "] from sequence is too\n"
      "small to accomodate a subset of size " << min_interval_size << ".\n";
    int needed_residues = min_interval_size - (en-st+1);
    int pad_left, pad_right;
    assert( needed_residues > 0 );

    if ((needed_residues % 2) == 0)
    {
      //if needed_residues is EVEN
      pad_left = needed_residues/2;
      pad_right = needed_residues/2;
    }
    else
    {
      pad_left = needed_residues/2+1;
      pad_right = needed_residues/2+1;
    }

    // Extend the interval to fit the min_interval_size
    // (being careful not to extend past either end of the molecule)
    // Taken care of the problem when the new interval does not
    // lie in the range of the sequence (ie. in the interval
    // [1,sequence_length].)  If not, shift it to the right or left so
    // that it does.)
    int shift_right = 0;
    int shift_left = 0;
    if (st - mol.begin() >= pad_left)
      st -= pad_left;
    else
    {//calculate the shift necessary
      shift_right = pad_left - (st-mol.begin());
      assert(shift_right > 0);
      st = mol.begin();
    }
    if (mol.end() -1 - en >= pad_right)
      en += pad_right;
    else
    {//calculate the shift necessary
      shift_left = pad_right - (mol.end()-1-en);
      assert(shift_left > 0);
      en = mol.end() -1;
    }

    //If sequence is just to small for interval (on both sides),
    //then the interval is the sequence.
    if (shift_right && shift_left)
    {                             
      st = mol.begin();
      en   = mol.end() -1;
    }
    else
    { //otherwise, apply the shifts
      if (shift_left)
      {
        assert(! shift_right);
        assert(st - mol.begin() >= shift_left);
        st -= shift_left;
      }
      if (shift_right)
      {
        assert(! shift_left);
        assert(mol.end()-1 - en >= shift_right);
        en  += shift_right;
      }
    }

    cerr << "extending the interval to ["
         << st->id << " - "
         << en->id << "]"
         << endl;
    //Finally, store the new "start" and "end" back in the interval-list
    interval.first = st->id;
    interval.last = en->id;

    extended = true;
  } // if (st > en-min_interval_size+1)
  return extended;
} //ExtendInterval()



void
ConvertHelixSheetRecordsToIntervalNtuples(int  num_files,
                                          vector<string> const& vFileNames,
                                          vector<vector<ClosedInterval> >& vvI,
                                          Biopolymer const *aMolecules,
                                          int min_interval_size)
{
  assert(num_files > 0);
  assert(num_files == vFileNames.size());
  assert(aMolecules);

  //The folowing stores all the helices and sheets from each molecule.
  vector<vector<ClosedInterval> > helixLI(num_files);// helixLI[0] stores all the
                                                  // helix intervals in the
                                                  // first molecule,
                                                  // helixLI[1] stores all the
                                                  // helix intervals in the
                                                  // second molecule and so on.
  vector<vector<ClosedInterval> > sheetLI(num_files);// same but stores sheets

  //********   First, fill the helixLI, and sheetLI arrays:   *********
  for (int m=0; m<num_files; ++m)
  {
    DEBUG_MSG(DBG_INTERVAL_LISTS, __FILE__ << ":" << __LINE__
              << "Parsing through pdb file \""
              << vFileNames[m]<<"\"\n"
              " looking for HELIX and SHEET records.");

    ifstream pdb_file(vFileNames[m].c_str(), ios::in);
    if (! pdb_file)
      ERR_INTERNAL("Can't load file: \""
                   << vFileNames[m] << "\"");

    ShortString  line;
    int          line_counter = 0;

    while (pdb_file.getline(line, SHORT_STRING_LENGTH))
    {
      ++line_counter;
      bool start_belongs, end_belongs; //For checking the endpoints to make
                                       //sure the interval is valid
      ClosedInterval interval;
      PDB::RecordType record_type =
        ConvertHelixSheetTurnRecordToInterval(line, interval);

      if (((record_type) == PDB::HELIX) ||
          ((record_type) == PDB::SHEET))
      {
        //print optional debug message
        if (record_type == PDB::HELIX)
          DEBUG_MSG(DBG_HELIX_SHEET, __FILE__ << ":" << __LINE__
                    << " Helix found from "
                    << interval.first << " to "
                    << interval.last);
        else
          DEBUG_MSG(DBG_HELIX_SHEET, __FILE__ << ":" << __LINE__
                    << " Strand of sheet found from "
                    << interval.first << " to "
                    << interval.last);

        //******   Now that we have helix or sheet interval,
        //******   check to make sure it lies within the
        //******   portion of the molecule the user has selected.
        //******   case 1:If it lies entirely inside, use the interval.
        //******   case 2:If it lies entirely outside, ignore the interval.
        //******   case 3:If it lies partially inside and partially outside,
        //******          print an error message and exit.
        start_belongs = aMolecules[m].find(interval.first) != aMolecules[m].end();
        end_belongs = aMolecules[m].find(interval.last) != aMolecules[m].end();

        //******   case 1:
        if (start_belongs && end_belongs) //I am lazy and only check
        {                                 //the endpoints ("start" and "end").
          //If the interval belongs to aMolecule[m],
          //then make sure the interval is large enough to accomodate
          //at least "min_interval_size" residues.
          ExtendInterval(interval,
                         aMolecules[m],
                         min_interval_size);
          if (record_type == PDB::HELIX)
            helixLI[m].push_back(interval);
          else
            sheetLI[m].push_back(interval);
        }

        //******   case 3:
        else if (start_belongs != end_belongs) 
          ERR("Invalid " 
              << ((record_type == PDB::HELIX) ?
                          "HELIX"
                          :
                          "SHEET")
              << " record found on line " << line_counter
              << " of file \"" << vFileNames[m] << "\"" << endl <<
              " START: "
              //<< interval.first << "\n"
              << "#" << interval.first.seqNum << ","
              << "chain:'" << interval.first.chainId << "',"
              << "insert:'" << interval.first.insertCode << "'\n"
              "  END : "
              //<< interval.last << "\n"
              << "#" << interval.last.seqNum << ","
              << "chain:'" << interval.last.chainId << "',"
              << "insert:'" << interval.last.insertCode << "'\n"
              "Part, this interval contains some residues which are invalid\n"
              "(they don't exist, or have been discarded by the user)\n"
              "and other residues which are valid.\n"
              "This could occur if the interval spans multiple chains\n"
              "(one or more of which was _not_ selected by the user).\n"
              "All HELIX/SHEET records must contain residues which\n"
              "are either all valid or all invalid.\n"
              "Check the chain(s)/interval(s) specified right after\n"
              "the pdb-filename to make sure they are consistent with\n"
              "the HELIX/SHEET records in that file.\n"
              "Aborting...");
        //******   case 2: nothing to do
        else
        {
          DEBUG_MSG(DBG_HELIX_SHEET, __FILE__ << ":" << __LINE__ <<
                    " Helix/Sheet lies outside filter-->discarded:\n"
                    "\"" << line << "\"");
        }
      } //if (record_type) == PDB::HELIX or PDB::SHEET) ||
    } //while (pdb_file.getline(line, ShortStringSize)
  } //for (int m=0; m<num_files; ++m)



  int num_proteins = num_files;


  // ***********************************************************************
  // *************   Now that we have the helices and sheets   *************
  // *************   from each molecule, generate all possible *************
  // *************   n-tuples of helices and sheets from       *************
  // *************   these proteins.                           *************
  // ***********************************************************************

  long total_permutations = 1;
  for(int m=0; m < num_proteins; ++m)
  {
    total_permutations *= helixLI[m].size();
    if (total_permutations > 200000)
      ERR("Error: Matching Helices-to-Helices generated more than\n"
          "       200000 interval-pair permutations.  Either your.\n"
          "       proteins have too many helices, or you have too many\n"
          "       proteins (note:this function is general and is not limited\n"
          "       to comparing only 2 proteins).\n"
          "       Aborting...");
  }

  for(long i=0; (i < total_permutations); ++i)
  {
    #if DBG_HELIX_SHEET
    cerr << "DBG_HELIX_SHEET, Adding HELIX-group#" << i << ":";
    #endif //#if DBG_HELIX_SHEET

    long left_over = i;
    long index;
    vector<ClosedInterval> li(num_proteins);
    for(int m = 0; m < num_proteins; ++m)
    {
      index = left_over % helixLI[m].size();
      left_over /= helixLI[m].size();
      //      li[m].vSeqIDs[m] = s;
      li[m] = helixLI[m][ index ];

      #if DBG_HELIX_SHEET
      cerr << "[" << li[m].first << "..." << li[m].last << "]";
      if (m == (num_proteins-1))
        cerr << endl;
      else
        cerr << ",";
      #endif //#if DBG_HELIX_SHEET
    }
    vvI.push_back( li ); //add a new an interval-list corresponding to
                         //this permutation of intervals.
                         //Example: If num_proteins == 2, then this generates
                         //        just a pair of helix intervals: one from
                         //        the first protein, one from the second.
                         //        and the loop itself generates all possible
                         //        pairs of helix intervals.
  }


  //---------------------------------------------------------------------

  //Now, do the same for SHEETS.
  //Generate all possible sets with num_proteins elements, where each
  //element is an interval corresponding to a strand of a SHEET from a
  //different protein.  For example, if there are only two proteins
  //(two pdb-files) then this generates a list of all pairs of
  //sheet-intervals.

  total_permutations = 1;
  for(int m=0; m < num_proteins; ++m)
  {
    total_permutations *= sheetLI[m].size();
    if (total_permutations > 200000)
      ERR("Error: Matching Sheets-to-Sheets generated more than\n"
          "       200000 interval-set permutations.  Either your.\n"
          "       proteins have too many sheet-strands, or you\n"
          "       have too many proteins (this part of the code\n"
          "       was not limited to work with only 2 proteins).\n"
          "       Aborting...");
  }

  for(long i=0; (i < total_permutations); ++i)
  {
    #if DBG_HELIX_SHEET
    cerr << "DBG_HELIX_SHEET, Adding SHEET-group#" << i << ":";
    #endif //#if DBG_HELIX_SHEET

    long left_over = i;
    long index;
    vector<ClosedInterval> li(num_proteins);
    for(int m = 0; m < num_proteins; ++m)
    {
      index = left_over % sheetLI[m].size();
      left_over /= sheetLI[m].size();
      //      li[m].vSeqIDs[m] = s;
      li[m] = sheetLI[m][ index ];

      #if DBG_HELIX_SHEET
      cerr << "[" << li[m].first << "..." << li[m].last << "]";
      if (m == (num_proteins-1))
        cerr << endl;
      else
        cerr << ",";
      #endif //#if DBG_HELIX_SHEET
    }
    vvI.push_back( li ); //add a new an interval-list corresponding to
                         //this permutation of intervals.
  }
} //ConvertHelixSheetRecordsToIntervalNtuples()




//Deletes intervals from a vector of intervals, if any of these
//intervals refer to residues which do not exist in the molecule.
//(When each interval in the list comes from the same molecule,
// use this function.)
bool DeleteInvalidIntervals(vector<ClosedInterval>& vI,
                            Biopolymer  &mol,
                            ClosedInterval* pFirstDeletedInterval,
                            int* pHowFarThroughList)
{
  bool deleted = false;
  vector<ClosedInterval>::iterator pI;
  for (pI = vI.begin();
       pI != vI.end();
       ++pI)
  {
    Biopolymer::const_iterator st = mol.find(pI->first);
    Biopolymer::const_iterator en = mol.find(pI->last);
    //Check the endpoints of the interval to make
    //sure residues with such identifiersw exist in the molecule.
    //We do not check for the residues in between the endpoints
    //because there's no way to do this.
    if ((st == mol.end()) || (en == mol.end()))
    {
      if ((deleted == false) && (pFirstDeletedInterval))
      {
        *pFirstDeletedInterval = *pI; //then return the interval to the caller
                                      //(if they requested this information)
        if (pHowFarThroughList != NULL)
          *pHowFarThroughList = pI - vI.begin();
      }
      deleted = true;
      --pI; //don't want to erase the iterator currently in use
      vI.erase(pI+1); //discard the offending list
                      //of intervals
    }
  }
  return deleted;
} //DeleteInvalidIntervals()




} //namespace minrms

