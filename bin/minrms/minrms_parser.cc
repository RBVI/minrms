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
#include <sstream>
#include <string>
#include <vector>
#include <cassert>
using namespace std;

#include <global_utils.h>
#include <simple_numeric_utils.h>
#include <parse_set.h>
#include <parse_utils.h>
#include <biopolymer.h>
#include <residue_names.h>
#include <load_pdb.h>
#include <intervals.h>
#include "minrms_settings.h"
#include "minrms_parser.h"

using namespace minrms;
using namespace ParseUtils;


// *** InterpretNumMatches() is converts a string of
// *** (for example "0.33") into an integer, in a way that is
// *** appropriate for indicating fractional lengths of sequences.
// *** (the string "0.33" represents a fraction of the length of a sequence.)
// ***-More precisely, if the input string has a numerical value
// *** between 0 and 1.0, then the function returns
// *** floor(0.33 * max_number_res_possible + 0.5) instead of that value.
// ***-But if it exceeds 1.0, then InterpretNumMatches() just returns
// *** that value untouched, (The strings "33" and "33.0" would result in 33.)
// ***-Special (weird) cases:
// *** If the string is "1.0", then it returns max_number_res_possible.
// *** If the string is "1", then it returns 1.
// ***
// *** This function is used a lot in the code that follows.
static int InterpretNumMatches(string num_res_str,
                               int    max_num_res_possible);



void PairAlignSettingsParser::LoadStructures(vector<string>& vArgs)
{
  // *******************************************
  // ***
  // *** The first thing we do, is to load up the
  // *** settings that will modify the structures.
  // ***
  // *******************************************

  //The "atoms_used_filename" contains the name of a file 
  //that controls which atoms get used in the computation.
  string atoms_used_filename = "atoms_used.txt";
  vector<string>::iterator pArg;
  for (pArg = vArgs.begin();
       ((pArg != vArgs.end()) && (*pArg != "-atoms-used"));
       ++pArg)
    {}

  if (pArg != vArgs.end())
  {
    if (pArg + 1 == vArgs.end())
      ERR("Error: missing parameter following the \"-atoms-used\" flag.");
    else
    {
      atoms_used_filename.assign((*(pArg+1)).begin(), (*(pArg+1)).end());
      vArgs.erase(pArg, pArg+2);
    }
  }

  //The "res_code_dict_filename" contains the name of
  //the file that defines the correspondence between
  //three letter lookup codes, and one-letter lookup codes.
  string res_code_dict_filename = "res_code_dict.txt";
  for (pArg = vArgs.begin();
       ((pArg != vArgs.end()) && (*pArg != "-res-code-dict"));
       ++pArg)
    {}

  if (pArg != vArgs.end())
  {
    if (pArg + 1 == vArgs.end())
      ERR("Error: missing parameter following the \"-res-code-dict\" flag.");
    else
    {
      res_code_dict_filename.assign((*(pArg+1)).begin(), (*(pArg+1)).end());
      vArgs.erase(pArg, pArg+2);
    }
  }



  // *****************************************
  // ***
  // ***  Now, load in the two molecules:
  // ***
  // *****************************************

  //The first step is to specify atoms used in calculation.
  LoadBackboneAtomSymbols(atoms_used_filename.c_str());

  // Assign const access to the two molecules.
  aMol = aMol_NC;
  assert(aMol_NC);
  assert(aMol == aMol_NC);
  aMol_f = aMol_f_NC;
  assert(aMol_f_NC);
  assert(aMol_f == aMol_f_NC);

  // The last two arguments should contain the filenames
  // as well as any residue subsets.
  if (vArgs.size() < 3)
    ERR("Error: "
        << vArgs[0] <<
        " expects at least two arguments:  the names\n"
        "       of the pdb-files for the structures you want to compare.");

  // First, figure out the PDB filenames from the argument list:
  vector<string> vPdbFileNamesAndSubsets(2);
  vector<string> vSubsetStrs(2);
  vPdbFileNamesAndSubsets[0] = vArgs[vArgs.size() - 2];
  vPdbFileNamesAndSubsets[1] = vArgs[vArgs.size() - 1];
  vArgs.erase(vArgs.end() - 2, vArgs.end()); //erase the last two arguments.
  for (int m = 0; m < 2; ++m)
  {
    string::const_iterator p = vPdbFileNamesAndSubsets[m].begin();

    //Read the part before the ",".
    //This part is the PDB filename.
    NextToken(p,
              vPdbFileNamesAndSubsets[m].end(),
              vPdbFileNames[m],
              ",");

    //Whatever is left over after the "," specifies
    //the subset of residues to use from this structure.
    //   (The next line circumvents some problem with const access
    //   that the alpha compiler complains about)
    string::const_iterator filename_end = vPdbFileNamesAndSubsets[m].end();
    vSubsetStrs[m].assign(p, filename_end);
    //  And the sgi compiler doesn't like this line
    //vSubsetStrs[m] = string(p, vPdbFileNamesAndSubsets[m].end());

    //Finally: open the PDB-files
    ifstream pdb_file(vPdbFileNames[m].c_str(), ios::in);
    if (! pdb_file)
      ERR("Error: unable to open file \""
          << vPdbFileNames[m] << "\" for reading.");
    pdb_file >> aMol_NC[m];

    #ifdef DEBUG
    #if DBG_LIN_MOL_DATA_STRUCT
    DEBUG_MSG(DBG_LIN_MOL_DATA_STRUCT, "Content of structure \""
              << vPdbFileNames[m] << "\":\n");
    cerr << aMol_NC[m];
    #endif //#if DBG_LIN_MOL_DATA_STRUCT
    #endif //#ifdef DEBUG
  }


  // These last post-processing steps are needed:
  for (int m = 0; m < 2; ++m)
  {
    aMol_NC[m].Finalize();//grant "special fast" access to certain atoms.
  }

  //Specify the translation between
  //1-letter residue codes (like in MSF files), and
  //3-letter residue codes (like in PDB files)
  //converting residue name codes
  //from PDB to MSF format later on.
  NameCodeLookup::Init(res_code_dict_filename);



  // *** Okay, now that the two original structures have
  // *** been loaded into aMol_NC[],
  // *** filter out the all but the residues that are desired
  // *** for use in the calculation.
  // *** Store this subset of residues in aMol_f_NC[]


  for (int m = 0; m < 2; ++m)
  {
    //We want aMol_f_NC[m] to contain a subset of the residues in aMol_NC[m].
    //We start by copying all of the residues from aMol_NC[m],
    //and delete the ones we don't want later.
    aMol_f_NC[m] = aMol_NC[m];

    //Now, parse the string containing the desired set of residues.
    vector<ClosedInterval> desired_residues;
    string::const_iterator p = vSubsetStrs[m].begin();
    while (p != vSubsetStrs[m].end())
    {
      vector<ClosedInterval> simple_set;
      ParseChimeraSet(p,
                      vSubsetStrs[m].end(),
                      aMol_NC[m],
                      simple_set,
                      ",");
      desired_residues.insert(desired_residues.end(),
                              simple_set.begin(), simple_set.end());
    }
    if ((desired_residues.size() == 0) && (aMol_NC[m].size() > 0))
    {
      ClosedInterval whole_molecule;
      whole_molecule.first = aMol_NC[m].front().id;
      whole_molecule.last  = aMol_NC[m].back().id;
      desired_residues.push_back(whole_molecule);
    }
    else
    {
      //For debugging purposes, print the set of residues selected.
      cout << "Residues selected from from structure  \""
           << vPdbFileNames[m] << "\":" << endl;
      for(vector<ClosedInterval>::const_iterator p = desired_residues.begin();
          p != desired_residues.end();
          ++p)
      {
        cout << "[" << (*p).first << " - " << (*p).last << "]";
        if (p+1 != desired_residues.end())
          cout << ", ";
      }
      cout << endl;
    }

    //Make sure the sets we have selected are valid.
    ClosedInterval invalid_interval;
    if (DeleteInvalidIntervals(desired_residues,
                               aMol_f_NC[m],
                               &invalid_interval))
    {
      ERR("Error: Error in format for set selection\n"
          " One of these residues does not belong to \"" << vPdbFileNames[m]
          << "\"\n"
          << invalid_interval.first.seqNum << invalid_interval.first.insertCode <<
          " from chain " << invalid_interval.first.chainId <<
          ", or\n"
          << invalid_interval.last.seqNum << invalid_interval.last.insertCode <<
          " from chain " << invalid_interval.last.chainId);
    }

    //Now clip the structures in aMol_f_NC[m] to contain only
    //the residues specified in "desired_residues".
    ClipMoleculeToListOfIntervals(aMol_f_NC[m],
                                  desired_residues);

    //Bullet proofing:
    //Check to make sure the structures are not empty:
    if (aMol_f_NC[m].size() == 0)
      ERR("Error: No residues found in structure \"" << vPdbFileNames[m] <<
          "\"\n"
          "       Either the PDB file contains no residues, or\n"
          "       you have selected a subset of residues from within\n"
          "       this structure that is empty.");
  } // for (int m = 0; m < 2; ++m)

  cout <<
    "(after filtering) length of molecule "
       << vPdbFileNames[0] << ": "
       << aMol_f_NC[0].size()
       << "\n"
    "                  length of molecule "
       << vPdbFileNames[1] << ": "
       << aMol_f_NC[1].size()
       << endl;

  assert(aMol == aMol_NC);
  assert(aMol_f == aMol_f_NC);
} //PairAlignSettingsParser::LoadStructures()


#include <parse_utils.h>



PairAlignSettingsParser::PairAlignSettingsParser(vector<string>& vArgs)
{
  LoadStructures(vArgs);

  //Initialize with default values:
  vMsfOutputLabels[0] = vPdbFileNames[0];
  vMsfOutputLabels[1] = vPdbFileNames[1];
  n_min = 1;

  //Now find the real values in the argument list:
  for (int i=1; i < vArgs.size(); ++i)
  {
    #ifdef DEBUG
    #if DBG_ARG_PARSER
    DEBUG_MSG(DBG_ARG_PARSER,
              "vArgs[i=" << i << "] = \"" << vArgs[i] << "\"\n"
              "  ---  vArgs[] at beginning of PairAlignSettingsParser loop:  ---");
    DisplayVectorOfStrings(vArgs);
    #endif // #if DBG_ARG_PARSER
    #endif // #ifdef DEBUG

    int num_arguments_deleted = 0;
    if (vArgs[i] == "-minN")
    {
      if ((i + 1 >= vArgs.size()) || (vArgs[i+1][0] == '-'))
        ERR("Error: The \"" << vArgs[i] <<
            "\" argument must be followed by a parameter.");

      int max_num_matches_possible = MIN(aMol_f_NC[0].size(), aMol_f_NC[1].size());
      n_min = InterpretNumMatches(vArgs[i+1], max_num_matches_possible);

      if (n_min <= 0)
        ERR("-minN parameter must be > 0");
      if (n_min > max_num_matches_possible)
        ERR("-minN parameter, (" << n_min << "), must not exceed\n"
            "      the number of residues in the smaller structure, ("
            << max_num_matches_possible << ").");

      num_arguments_deleted = 2;
    } //if (vArgs[i] == "-minN")
    else if (vArgs[i] == "-output-msf-labels")
    {
      if ((i + 1 >= vArgs.size()) || (vArgs[i+1][0] == '-'))
        ERR("Error: The \"" << vArgs[i] <<
            "\" argument must be followed by\n"
            "       two parameters separated by a comma.\n");
      string::const_iterator p = vArgs[i+1].begin();
      NextToken(p,
                vArgs[i+1].end(),
                vMsfOutputLabels[0],
                ",");
      NextToken(p,
                vArgs[i+1].end(),
                vMsfOutputLabels[1],
                ",");
      num_arguments_deleted = 2;
    }

    if (num_arguments_deleted > 0)
    {
      assert(i+num_arguments_deleted <= vArgs.size());
      vArgs.erase(vArgs.begin()+i, vArgs.begin()+i+num_arguments_deleted);
      --i; //rewind by one so that when we invoke ++i
           //at the end of this loop, i's value will not change.
           //This will point us to the next un-read argument.
    }
  } //loop over all arguments:  "for (int i=1; i < vArgs.size(); ++i)"  

} //PairAlignSettingsParser::PairAlignSettingsParser()




OrGenSettingsParser::OrGenSettingsParser(vector<string>& vArgs):
  PairAlignSettingsParser(vArgs)
{
  fm_fragment_size = 4;
  fm_cutoff_rmsd = fm_NO_CUTOFF_RMSD;
  vMsfInputLabels[0] = vMsfOutputLabels[0];
  vMsfInputLabels[1] = vMsfOutputLabels[1];
  bool use_helices_and_sheets = true;

  for (int i=1; i < vArgs.size(); ++i)
  {
    #ifdef DEBUG
    #if DBG_ARG_PARSER
    DEBUG_MSG(DBG_ARG_PARSER,
              "vArgs[i=" << i << "] = \"" << vArgs[i] << "\"\n"
              "  ---  vArgs[] at beginning of OrGenSettingsParser loop:  ---");
    DisplayVectorOfStrings(vArgs);
    #endif // #if DBG_ARG_PARSER
    #endif // #ifdef DEBUG

    int num_arguments_deleted = 0;
    if (vArgs[i] == "-fm")
    {
      if ((i + 1 < vArgs.size()) && (vArgs[i+1][0] != '-'))
      {
        string::const_iterator p = vArgs[i+1].begin();
        NextInt(p,
                vArgs[i+1].end(),
                fm_fragment_size,
                ",");
        NextReal(p,
                 vArgs[i+1].end(),
                 fm_cutoff_rmsd,
                 ",");
        if (fm_fragment_size == 0)
          fm_fragment_size = fm_DISSABLE_FRAGMENT_MATCHING;
        else if (fm_fragment_size < 3)
        {
          ERR("Error: Invalid parameter for \"-fm\" command:\n"
              "       \"" << fm_fragment_size << "\"\n"
              "       Must be at least 3");
        }
        num_arguments_deleted = 2;
      } //if (i + 1 < vArgs.size())
      else
        num_arguments_deleted = 1;
    } //if (vArgs[i] == "-fm")

    else if (vArgs[i] == "-FM")
    {
      fm_fragment_size = fm_DISSABLE_FRAGMENT_MATCHING;
      num_arguments_deleted = 1;
    }

    else if (vArgs[i] == "-fm-interval")
    {
      vector<ClosedInterval> interval_pair(2);
      vector<ClosedInterval> subset_struct_1;
      vector<ClosedInterval> subset_struct_2;
      string::const_iterator p = vArgs[i+1].begin();
      ParseChimeraSet(p,
                      vArgs[i+1].end(),
                      aMol_f[0],
                      subset_struct_1,
                      ",");
      ParseChimeraSet(p,
                      vArgs[i+1].end(),
                      aMol_f[1],
                      subset_struct_2,
                      ",");

      if ((p != vArgs[i+1].end()) ||
          (subset_struct_1.size() != 1) ||
          (subset_struct_2.size() != 1))
      {
        ERR("Error: Error in argument format:\n"
            "       The \"-fm-interval\" command must be followed by an\n"
            "       argument containing a contiguous subset of residues from\n"
            "       each structure separated by a comma.  The subset on each\n"
            "       side of the comma must be contiguous (a single interval).\n"
            "       (Example: expressions like \"50-100.A-C\" refer to\n"
            "        subsets of residues which are not contiguous.)");
      }

      interval_pair[0] = subset_struct_1[0];
      interval_pair[1] = subset_struct_2[0];

      //Check the endpoints of the interval to make
      //sure residues with such identifiers exist in the molecule.
      //We do not check for the residues in between the endpoints
      //because there's no way to do this.
      for (int m = 0; m < 2; ++m)
      {
        Biopolymer::const_iterator a= aMol[m].find(interval_pair[m].first);
        if (a == aMol[m].end())
          ERR("Error: Invalid interval specified to -fm-interval flag.\n"
              "       residue "
              << interval_pair[m].first.seqNum
              << interval_pair[m].first.insertCode <<
              " from chain "
              << interval_pair[m].first.chainId << "\n"
              " either does not belong to \"" << vPdbFileNames[m] << "\", or\n"
              " does not belong to the subset of residues\n"
              " that will be used for the calculation.\n");

        Biopolymer::const_iterator b= aMol[m].find(interval_pair[m].last);
        if (b == aMol[m].end())
          ERR("Error: Invalid interval specified to -fm-interval flag.\n"
              "       residue "
              << interval_pair[m].last.seqNum
              << interval_pair[m].last.insertCode <<
              " from chain "
              << interval_pair[m].last.chainId << "\n"
              " either does not belong to \"" << vPdbFileNames[m] << "\", or\n"
              " does not belong to the subset of residues\n"
              " that will be used for the calculation.\n");
      }

      //Otherwise, the intervals are okay. Insert them into the list.
      fm_interval_pair_list.push_back(interval_pair);

      cout
        << "Fragments from residues "
        <<"[" << interval_pair[0].first << " - " << interval_pair[0].last
        <<"] in structure \"" << vPdbFileNames[0] << "\"\n"
        << "  will be matched with\n"
        << "fragments from residues "
        <<"[" << interval_pair[1].first << " - " << interval_pair[1].last
        <<"] in structure \"" << vPdbFileNames[1] << "\"" << endl;

      num_arguments_deleted = 2;
    } //else if (vArgs[i] == "-fm-interval")


    else if (vArgs[i] == "-hs")
    {
      use_helices_and_sheets = true;
      num_arguments_deleted = 1;
    }

    else if (vArgs[i] == "-HS")
    {
      use_helices_and_sheets = false;
      num_arguments_deleted = 1;
    }

    else if (vArgs[i] == "-read-matrix")
    {
      if ((i+1 >= vArgs.size()) ||
          (vArgs[i+1][0] == '-'))
        ERR("Error: The \"-read-matrix\" flag must be followed by a single\n"
            "       argument containing a comma-separated list of the names\n"
            "       of ascii files containing 3x4 matrices.\n"
            "       (The filenames cannot contain '-' or ','\n"
            "        and no spaces are allowed between the commas).");

      //Now, read in the list of names of matrix files:
      
      for(string::const_iterator p = vArgs[i+1].begin();
          p != vArgs[i+1].end();
          )
      {
        string matrix_filename;
        NextToken(p,
                  vArgs[i+1].end(),
                  matrix_filename,
                  ",");
        input_matrix_files.push_back(matrix_filename);
      }
      num_arguments_deleted = 2;
    } //else if (vArgs[i] == "-read-matrix")

    else if (vArgs[i] == "-read-msf")
    {
      stringstream read_msf_syntax_error;
      read_msf_syntax_error <<
        "Error: The \"-read-msf\" flag must be followed by two arguments\n"
        "       The first argument specifies two comma-separated labels\n"
        "       that should be used to identify the sequences\n"
        "       corresponding to the two structures.\n"
        "       " << vPdbFileNames[0] << ","
                   << vPdbFileNames[1] << "\n"
        "       The second argument should be a comma-separated list\n"
        "       of the names of msf-files to read in.\n"
        "       (The filenames cannot contain '-' or ','\n"
        "        and no spaces are allowed between the commas).";

      if ((i+2 >= vArgs.size()) ||
          (vArgs[i+1][0] == '-') ||
          (vArgs[i+2][0] == '-'))
        ERR(read_msf_syntax_error.str());

      //Now, read in the two msf sequence labels.
      string::const_iterator p = vArgs[i + 1].begin();
      NextToken(p,
                vArgs[i+1].end(),
                vMsfInputLabels[0],
                ",");
      if (p == vArgs[i+1].end())
        ERR(read_msf_syntax_error.str());
      NextToken(p,
                vArgs[i+1].end(),
                vMsfInputLabels[1],
                ",");
      if (p != vArgs[i+1].end())
        ERR(read_msf_syntax_error.str());

      //Now, read in the list of names of msf files:
      for(p = vArgs[i+2].begin();
          p != vArgs[i+2].end();
          )
      {
        string msf_filename;
        NextToken(p,
                  vArgs[i+2].end(),
                  msf_filename,
                  ",");
        input_msf_files.push_back(msf_filename);
      }
      num_arguments_deleted = 3;
    } //end of chain of if-then's (Example: "if (vArgs[i] == "-fm")...else if"

    //Delete all the arguments we have processed in this iteration (if any)
    if (num_arguments_deleted > 0)
    {
      assert(i+num_arguments_deleted <= vArgs.size());
      vArgs.erase(vArgs.begin()+i, vArgs.begin()+i+num_arguments_deleted);
      --i; //rewind by one so that when we invoke ++i
           //at the end of this loop, i's value will not change.
           //This will point us to the next un-read argument.
    }

  } // loop over arguments "for (int i=1; i < vArgs.size(); ++i)"


  //For debugging purposes, print the sets
  //of residues used for fragment matching.
  #ifdef DEBUG
  #if DBG_INTERVAL_LISTS
  DEBUG_MSG(DBG_INTERVAL_LISTS, "paired regions for fragment matching\n"
            "                   before helices and sheets are considered:\n"
            << vPdbFileNames[0] << "," << vPdbFileNames[1]);
  if (fm_interval_pair_list.size() == 0)
    cerr << "list of regions is empty." << endl;
  else
  {
    for(vector<vector<ClosedInterval> >::const_iterator
          p = fm_interval_pair_list.begin();
        p !=  fm_interval_pair_list.end();
        ++p)
    {
      cerr << "[" << (*p)[0].first << " - " << (*p)[0].last << "],"
           << "[" << (*p)[1].first << " - " << (*p)[1].last << "]" << endl;
    }
    cerr << endl;
  }
  #endif //#if DBG_INTERVAL_LISTS
  #endif //#ifdef DEBUG

  // *** now figure out to do if the user has selected the "-hs" flag.

  if (use_helices_and_sheets)
  {
    cout << "Secondary structure considered.\n"
            " Parsing Helix & Sheet records.\n";
    DEBUG_MSG(DBG_INTERVAL_LISTS,
              "\n"
              "********************************************************\n"
              "********************************************************");
    //The way I did it was to have the helix&sheet records loaded 
    //from the pdb-files at a separate time than when I load the
    //coordinates of the residues.
    //That way, I did not have to rewrite my existing pdb-file reader.)

    vector<vector<ClosedInterval> >  helix_and_sheet_interval_pairs;

    ConvertHelixSheetRecordsToIntervalNtuples(2,
                                              vPdbFileNames,
                                              helix_and_sheet_interval_pairs,
                                              aMol_f,
                                              fm_fragment_size);

    if (helix_and_sheet_interval_pairs.size() == 0)
      ERR("Error: At least one of the PDB files contains no\n"
          "       HELIX or SHEET records.\n"
          "       To circumvent this problem, either dissable\n"
          "       the helix-sheet matching optimization using the\n"
          "       \"-HS\" flag,\n"
          "       or, add helix and sheet records to the pdb file\n"
          "       (use a utility like \"ksdssp\")");
    else
      fm_interval_pair_list.insert(fm_interval_pair_list.end(),
                                   helix_and_sheet_interval_pairs.begin(),
                                   helix_and_sheet_interval_pairs.end());

  } //if (use_helices_and_sheets)
  
  //If no intervals were explicitly selected by the user,
  //and the user dissabled helix and sheet matching,
  //The default is to use all the residues in both molecules
  //for fragment matching.
  if (fm_interval_pair_list.size() == 0)
  {
    vector<ClosedInterval> use_all_of_both(2);
    use_all_of_both[0].first = aMol_f[0][0].id;
    use_all_of_both[0].last  = aMol_f[0][aMol_f[0].size() - 1].id;
    use_all_of_both[1].first = aMol_f[1][0].id;
    use_all_of_both[1].last  = aMol_f[1][aMol_f[1].size() - 1].id;
    fm_interval_pair_list.push_back(use_all_of_both);
  }
} //OrGenSettingsParser::OrGenSettingsParser(vector<string>& vArgs)





SearchOrSettingsParser::SearchOrSettingsParser(vector<string>& vArgs):
  PairAlignSettingsParser(vArgs)
{
  refine_method = REFINE_BEST;
  refine_max_iters = ITERATE_UNTIL_CONVERGENCE;
  refine_while_searching = true;
  num_solutions = 1;
  alt_method = ALT_METHOD_3D;
  alt_min_3d_difference = 0.0;
  alt_max_pairs_common = 0;
  alt_rmsd_tolerance = ALTERNATES_NOT_LIMITED_BY_RMSD;

  for (int i=1; i < vArgs.size(); ++i)
  {
    #ifdef DEBUG
    #if DBG_ARG_PARSER
    DEBUG_MSG(DBG_ARG_PARSER,
              "vArgs[i=" << i << "] = \"" << vArgs[i] << "\"\n"
              "  ---  vArgs[] at beginning of SearchOrSettingsParser loop:  ---");
    DisplayVectorOfStrings(vArgs);
    #endif // #if DBG_ARG_PARSER
    #endif // #ifdef DEBUG

    int num_arguments_deleted = 0;

    if (vArgs[i] == "-r")
    {
      refine_method = REFINE_BEST;
      if ((i + 1 < vArgs.size()) && (vArgs[i+1][0] != '-'))
      {
        refine_max_iters = atoi(vArgs[i+1].c_str());
        num_arguments_deleted = 2;
      }
      else
        num_arguments_deleted = 1;
    }        
    else if (vArgs[i] == "-r-all")
    {
      refine_method = REFINE_ALL;
      if ((i + 1 < vArgs.size()) && (vArgs[i+1][0] != '-'))
      {
        refine_max_iters = atol(vArgs[i+1].c_str());
        num_arguments_deleted = 2;
      }
      else
        num_arguments_deleted = 1;
    }
    else if (vArgs[i] == "-R")
    {
      refine_method = NO_REFINE;
      num_arguments_deleted = 1;
    }
    else if (vArgs[i] == "-ir")
    {
      refine_while_searching = true;
      num_arguments_deleted = 1;
    }  else if (vArgs[i] == "-IR")  { refine_while_searching = false; num_arguments_deleted = 1; }


    else if ((vArgs[i] == "-alt-o") ||
             (vArgs[i] == "-alt-o-orig"))
    {
      if (vArgs[i] == "-alt-o")
        alt_method = ALT_METHOD_3D;
      else
        alt_method = ALT_METHOD_3D_ORIG;

      if (! ((i + 1 < vArgs.size()) && (vArgs[i+1][0] != '-')))
      {
        ERR("Error: the \"-alt-o\" flag must be followed by\n"
            "       at least two parameters (separated by a comma, no whitespace):\n"
            "\n"
            "    1) An integer specifying the desired total number of\n"
            "       solutions including alternate solutions\n"
            "       (2nd-best, 3rd best, etc)\n"
            "       you would like reported.  For example, a value of 3 means\n"
            "       the 3 top alignments will be found\n"
            "    2) A distance (in Angstroms) representing the minium\n"
            "       threshold displacment (in RMSD) between the\n"
            "       two positions of the second molecule after\n"
            "       being rotated and translate to accomodate\n"
            "       two different alignments in question, so \n"
            "       that these two alignments will be sufficiently\n"
            "       different to be reported as alternate solutions.\n"
            "    3) A third argument is optional but recommended.\n"
            "       The third argument is a number indicating the\n"
            "       ratio of the maximum RMSD an alignment can have\n"
            "       compared to the RMSD of the optimal alignment\n"
            "       for that alignment to be considered as a\n"
            "       potential alternate solution.  (Note, as such, this\n"
            "       parameter must be at leat 1.0)\n"
            "          As this parameter's value get's closer to 1.0\n"
            "       the space used up by the potentially massive\n"
            "       temporary files is reduced and the\n"
            "       computation should be faster as well.\n"
            "       However, if you make the number too small, you\n"
            "       may not find any alternate solutions at all.\n"
            "       Recommended values are 1.2 to 2.5\n"
            "       ---\n"
            "       For example if the optimal solution has an\n"
            "       RMSD of 2 Anstroms, and the this argument is\n"
            "       set to \"1.2\" (i.e., a 20-percent increase),\n"
            "       then the 2nd best solution will be reported as\n"
            "       long as it's RMSD does not exceed 2.4 Angstroms\n"
            "       (i.e. 1.2 x 2 Angstroms)\n"
            "   You need to specify at least 2 arguments.\n"
            "   Aborting...");
      } // if ((i + 1 < vArgs.size()) && (vArgs[i+1][0] != '-'))

      assert(vArgs[i+1].size() != 0);
      string::const_iterator p = vArgs[i+1].begin();
      NextInt(p,
              vArgs[i+1].end(),
              num_solutions,

              ",");
      if (num_solutions < 2)
        ERR("Error:  The \"-alt-o\" command-line flag's first\n"
            "        parameter must be at least 2.");
      if (p == vArgs[i+1].end())
        ERR("Error:  The \"-alt-o\" command-line flag's must\n"
            "        be followed by at least 2 parameters\n"
            "        (seperated by a comma, no whitespace)");
      NextReal(p,
               vArgs[i+1].end(),
               alt_min_3d_difference,
               ",");
      if (alt_min_3d_difference < 0.0f)
        ERR("Error:  The \"-alt-o\" command-line flag's second\n"
            "        parameter (min_delta_orientation) must be at least 0.");
      NextReal(p,
               vArgs[i+1].end(),
               alt_rmsd_tolerance,
               ",");
      if ((alt_rmsd_tolerance <= 1.0f)
          &&
          (alt_rmsd_tolerance
           !=
           ALTERNATES_NOT_LIMITED_BY_RMSD))
        ERR("Error:  The \"-alt-o\" command-line flag's third\n"
            "        parameter (RMSD_tolerance) must be greater than 1.0.");

      num_arguments_deleted = 2;
    } //else if ((vArgs[i] == "-alt-o") || (vArgs[i] == "-alt-o-orig"))

    else if (vArgs[i] == "-alt-r")
    {

      alt_method = ALT_METHOD_PAIRS;

      if (! ((i + 1 < vArgs.size()) && (vArgs[i+1][0] != '-')))
      {
        ERR("Error: the \"-alt-r\" flag\n"
            "       needs to be immediately followed by at least one parameter.\n"
            "       (A second argument can also be supplied separated from\n"
            "        the first argument by a comma, no whitespace).\n"
            "       These arguments are:\n"
            "       1) An integer specifying the desired total number\n"
            "          of solutions including alternate solutions\n"
            "          (2nd-best, 3rd best, etc)\n"
            "          you would like reported.  For example, a value of\n"
            "          3 means the 3 top alignments will be found.\n"
            "       2) (optional) An integer specifying the maximum number of pairs\n"
            "          of residues two alignments can have in common\n"
            "          without being too similar to be reported as\n"
            "          1st and 2nd best.  A value of 0 means the 2nd-best\n"
            "          alignment cannot match any pair of residues matched\n"
            "          by the best alignment.  The default value is 0.\n"
            "       3) (optional) A real number larger than 1.0.\n"
            "            If this argument is specified, then a database of\n"
            "          temporary files will be created to limit the orientations\n"
            "          that are considered in the search to ones which produce\n"
            "          alignments with low RMSDs.  This reduces computation time.\n"
            "             This optional argument is a number indicating the\n"
            "          ratio of the maximum RMSD an alignment can have\n"
            "          compared to the RMSD of the optimal alignment\n"
            "          for that alignment to be considered as a\n"
            "          potential alternate solution.  (Note, as such, this\n"
            "          parameter must be at leat 1.0)\n"
            "             As this parameter's value get's closer to 1.0\n"
            "          the space used up by the potentially massive\n"
            "          temporary files is reduced and the\n"
            "          computation should be faster as well.\n"
            "          However, if you make the number too small, you\n"
            "          may not find any alternate solutions at all.\n"
            "          Recommended values are 1.2 to 2.5\n"
            "          ---\n"
            "          For example if the optimal solution has an\n"
            "          RMSD of 2 Anstroms, and the this argument is\n"
            "          set to \"1.2\" (i.e., a 20-percent increase),\n"
            "          then the 2nd best solution will be reported as\n"
            "          long as it's RMSD does not exceed 2.4 Angstroms\n"
            "          (i.e. 1.2 x 2 Angstroms).\n"
            "   You need to specify at least 1 argument.\n"
            "   Aborting...");
      } // if ((i + 1 < vArgs.size()) && (vArgs[i+1][0] != '-'))

      assert(vArgs[i+1].size() != 0);
      string::const_iterator p = vArgs[i+1].begin();
      NextInt(p,
              vArgs[i+1].end(),
              num_solutions,

              ",");
      if (num_solutions < 2)
        ERR("Error:  The \"-alt-r\" command-line flag's first\n"
            "        parameter must be at least 2.");
      //if (p == vArgs[i+1].end())
      //  ERR("Error:  The \"-alt-r\" command-line flag's must\n"
      //      "        be followed by at least 2 parameters\n"
      //      "        (seperated by a comma, no whitespace)");
      NextInt(p,
               vArgs[i+1].end(),
               alt_max_pairs_common,
               ",");
      if (alt_max_pairs_common < 0)
        ERR("Error:  The \"-alt-r\" command-line flag's second\n"
            "        parameter (max_pairs_in_common) must be at least 0.");

      NextReal(p,
               vArgs[i+1].end(),
               alt_rmsd_tolerance,
               ",");
      if ((alt_rmsd_tolerance <= 1.0f)
          &&
          (alt_rmsd_tolerance
           !=
           ALTERNATES_NOT_LIMITED_BY_RMSD))
        ERR("Error:  The \"-alt-o\" command-line flag's third\n"
            "        parameter (RMSD_tolerance) must be greater than 1.0.");

      num_arguments_deleted = 2;
    } //else if (vArgs[i] == "-alt-r")
    


    //Delete all the arguments we have processed in this iteration (if any)
    if (num_arguments_deleted > 0)
    {
      assert(i+num_arguments_deleted <= vArgs.size());
      vArgs.erase(vArgs.begin()+i, vArgs.begin()+i+num_arguments_deleted);
      --i; //rewind by one so that when we invoke ++i
           //at the end of this loop, i's value will not change.
           //This will point us to the next un-read argument.
    }

  } // loop over arguments "for (int i=1; i < vArgs.size(); ++i)"


  //Check for forbidden input at this point.
  if (num_solutions > 1)
  {
    switch (alt_method)
    {
    case ALT_METHOD_3D:
      if (refine_while_searching == false)
      {
        cerr <<
          "     -----------------------------------------------------------\n"
          "Warning: Since you have selected the -alt-o flag, there is no\n"
          "        reason not to use the \"intermediate refinement\" feature.\n"
          "        This feature is allways preferable because it can reduce\n"
          "        the RMSD of the alignments slightly.\n"
          "        (The \"intermediate refinement\" feature is on by default\n"
          "         and currently has been dissabled using the -IR flag.)\n"
          "        Use of the -alt-o flag overrides the -IR flag\n"
          "     -----------------------------------------------------------"
             << endl;
        refine_while_searching = true;
      }
      break;
    case ALT_METHOD_3D_ORIG:
      if (refine_while_searching == true)
      {
        cerr <<
          "     -----------------------------------------------------------\n"
          "Warning: The -alt-o-orig flag dissables the\n"
          "         \"intermediate refinement\" feature.\n"
          "         (The \"intermediate refinement\" feature is on by\n"
          "          default and can be dissabled using the -IR flag.)\n"
          "        MINRMS will behave as if you had ran it with\n"
          "        the -IR flag.  If you need \"intermediate refinement\",\n"
          "        use the -alt-o flag instead.\n"
          "     -----------------------------------------------------------"
             << endl;
        refine_while_searching = false;
      }
      break;
    } //switch (alt_method)
  } //checking for forbidden input when (num_solutions > 1)

} //SearchOrSettingsParser::SearchOrSettingsParser()




SearchOrDyn3dParser::SearchOrDyn3dParser(vector<string>& vArgs):
  PairAlignSettingsParser(vArgs), SearchOrSettingsParser(vArgs)
{
  int max_num_matches_possible = MIN(aMol_f[0].size(), aMol_f[1].size());
  n_max = max_num_matches_possible;
  dyn3d_max_rmsd = dyn3d_NO_MAX_RMSD;
  dyn3d_max_dist = DistanceMetric::NO_MAX_PHYS_DISTANCE;
  orientation_filter_nw = false;
  orientation_filter_nw_max_dist = 0.0;
  orientation_filter_nw_min_num_matches = 1;
  display_stats_interval = 300;

  for (int i=1; i < vArgs.size(); ++i)
  {
    #ifdef DEBUG
    #if DBG_ARG_PARSER
    DEBUG_MSG(DBG_ARG_PARSER,
              "vArgs[i=" << i << "] = \"" << vArgs[i] << "\"\n"
              "  ---  vArgs[] at beginning of SearchOrDyn3dParser loop:  ---");
    DisplayVectorOfStrings(vArgs);
    #endif // #if DBG_ARG_PARSER
    #endif // #ifdef DEBUG

    int num_arguments_deleted = 0;

    if (vArgs[i] == "-maxN")
    {
      if ((i + 1 >= vArgs.size()) || (vArgs[i+1][0] == '-'))
        ERR("Error: The \"" << vArgs[i] <<
            "\" argument must be followed by a parameter.");

      n_max = InterpretNumMatches(vArgs[i+1],
                                  max_num_matches_possible);
      //Bullet-proofing:
      if (n_max <= 0)
        ERR("-n_max parameter must be > 0");
      if (n_max > max_num_matches_possible)
        ERR("-n_max parameter, (" << n_max << "), must not exceed\n"
            "      the number of residues in the smaller structure, ("
            << max_num_matches_possible << ").");
      if (n_max < n_min)
        ERR("-n_min parameter, (" << n_min << "), must not exceed n_max ("
            << n_max << ").");
      num_arguments_deleted = 2;
    }
    else if (vArgs[i] == "-max-rmsd")
    {
      if ((i + 1 >= vArgs.size()) || (vArgs[i+1][0] == '-'))
        ERR("Error: The \"" << vArgs[i] <<
            "\" argument must be followed by a parameter.");
      dyn3d_max_rmsd = atof(vArgs[i+1].c_str());
      num_arguments_deleted = 2;
    }
    else if (vArgs[i] == "-max-dist")
    {
      if ((i + 1 >= vArgs.size()) || (vArgs[i+1][0] == '-'))
        ERR("Error: The \"" << vArgs[i] <<
            "\" argument must be followed by a parameter.");
      dyn3d_max_dist = atof(vArgs[i+1].c_str());
      num_arguments_deleted = 2;
    }



    else if (vArgs[i] == "-of")
    {
      orientation_filter_nw = true;

      if ((i + 1 >= vArgs.size()) || (vArgs[i+1][0] == '-'))
        ERR("Error: The \"" << vArgs[i] <<
            "\" argument must be followed by two parameters,\n"
            "       separated by a comma.");

      string::const_iterator p = vArgs[i+1].begin();
      NextReal(p,
               vArgs[i+1].end(),
               orientation_filter_nw_max_dist,
               ",");

      if (p == vArgs[i+1].end())
        orientation_filter_nw_min_num_matches = n_min;
      else
      {
        string min_n_str;
        NextToken(p,
                  vArgs[i+1].end(),
                  min_n_str,
                  ",");
        if (p != vArgs[i+1].end())
          ERR("Error: Extra irrelevant characters were found trailing\n"
              "       the second argument to the \"-of\" flag.");
      
        orientation_filter_nw_min_num_matches =
          InterpretNumMatches(min_n_str, max_num_matches_possible);
      }

      //Bullet-proofing:
      if (orientation_filter_nw_min_num_matches <= 0)
        ERR("-of's second parameter must be > 0");
      if (orientation_filter_nw_min_num_matches > max_num_matches_possible)
        ERR("-of's second parameter, ("
            << orientation_filter_nw_min_num_matches << "), must not exceed\n"
            "      the number of residues in the smaller structure, ("
            << max_num_matches_possible << ").");

      num_arguments_deleted = 2;
    } //else if (vArgs[i] == "-of")
    else if (vArgs[i] == "-OF")
    {
      orientation_filter_nw = false;
      num_arguments_deleted = 1;
    }

    //Delete all the arguments we have processed in this iteration (if any)
    if (num_arguments_deleted > 0)
    {
      assert(i+num_arguments_deleted <= vArgs.size());
      vArgs.erase(vArgs.begin()+i, vArgs.begin()+i+num_arguments_deleted);
      --i; //rewind by one so that when we invoke ++i
           //at the end of this loop, i's value will not change.
           //This will point us to the next un-read argument.
    }
  } // loop over arguments "for (int i=1; i < vArgs.size(); ++i)"


  //Additional bullet-proofing:

  if (n_max < orientation_filter_nw_min_num_matches)
    ERR("-of's second parameter, ("
        << orientation_filter_nw_min_num_matches <<
        "), must not exceed n_max ("
        << n_max << ").");

} //SearchOrDyn3dParser::SearchOrDyn3dParser()



SearchOrNWParser::SearchOrNWParser(vector<string>& vArgs):
  PairAlignSettingsParser(vArgs),
  SearchOrSettingsParser(vArgs)
{
  nw_which_criteria = SINGLE_ANSWER_MANUAL_GAP;
  nw_fixed_rmsd = 3.0;
  nw_fixed_n = (MIN(aMol_f[0].size(), aMol_f[1].size())+1) / 2;
  nw_max_dist = 6.0;
  nw_gap = 18.0;
  nw_max_iters = nw_DEFAULT_MAX_ITERS;
}



MinrmsParser::MinrmsParser(vector<string>& vArgs):
  PairAlignSettingsParser(vArgs),
  OrGenSettingsParser(vArgs),
  SearchOrSettingsParser(vArgs),
  SearchOrDyn3dParser(vArgs),
  SearchOrNWParser(vArgs)
{
  which_algorithm = ALGO_DYN3D;

  for (int i=1; i < vArgs.size(); ++i)
  {
    #ifdef DEBUG
    #if DBG_ARG_PARSER
    DEBUG_MSG(DBG_ARG_PARSER,
              "vArgs[i=" << i << "] = \"" << vArgs[i] << "\"\n"
              "  ---  vArgs[] at beginning of MinrmsParser loop:  ---");
    DisplayVectorOfStrings(vArgs);
    #endif // #if DBG_ARG_PARSER
    #endif // #ifdef DEBUG

    int num_arguments_deleted = 0;

    if (vArgs[i] == "-sa")
    {
      //Bulletproofing minrms for the alpha-release:
      //Actually, -sa is supported but has some fatal bugs.
      //I don't want to fix the -sa flag right now
      //(or the entire SearchOrNW class that that feature depends on),
      //so I print this instead.  Andrew 1/30/00
      ERR("Error: The " << g_version_string << " version of minrms does not\n"
          "       support the \"-sa\" flag.  Please wait for a later\n"
          "       release.\n");

      which_algorithm = ALGO_NEEDLEMAN_WUNSCH;
      cout <<
        "  (single-answer mode selected.\n"
        "   minrms will use an asymptotically faster algorithm,\n"
        "   that generates only one answer.)" << endl;  

      if ((i + 1 >= vArgs.size()) || (vArgs[i+1][0] == '-'))
        ERR("Error: the \"-sa\" command line argument requires\n"
            "      an additional mode specifier followed by a comma, \n"
            "      and then at least one parameter\n"
            "      The mode specifier must be one of \"rmsd\", \"max-dist\", \"n\", or \"nw\"\n"
            "      The parameter should be a number indicating either the desired\n"
            "      RMSD of the alignment, the desired maximum-distance between\n"
            "      residue-pairs, the desired number of residue pairs.\n"
            "      in the alignment, or the desired gap penalty used, repectively.\n"
            "      In addition,\n"
            "         if \"rmsd\" or \"n\" is chosen, in additional parameter\n"
            "      can also be supplied (also separated by a comma, no spaces)\n"
            "      indicating the number of maximum number of iterations of\n"
            "      attempting different gap parameters in order to achieve the\n"
            "      target RMSD or number of matches desired, respectively.\n"
            "      (default " << nw_DEFAULT_MAX_ITERS << ")");

      string::const_iterator p = vArgs[i+1].begin();

      string mode_str;
      NextToken(p,
                vArgs[i+1].end(),
                mode_str,
                ",");

      if (mode_str == "rmsd")
        nw_which_criteria = SINGLE_ANSWER_FIXED_RMSD;
      else if (mode_str == "n")
        nw_which_criteria = SINGLE_ANSWER_FIXED_N;
      else if (mode_str == "max-dist")
        nw_which_criteria = SINGLE_ANSWER_MAX_DIST;
      else if (mode_str == "nw")
        nw_which_criteria = SINGLE_ANSWER_MANUAL_GAP;
      else
        ERR("Error: Unrecongnized mode following \"-sa\" flag: \""
            << mode_str << "\"\n"
            "       Should be one of: \"rmsd\", \"n\", \"max-dist\", \"nw\"");

      if (p == vArgs[i+1].end())
        ERR("Error: the \"-sa\" flag requires a numeric parameter\n"
            "       following the \"-sa " << mode_str << ",\" flag.");

      switch (nw_which_criteria)
      {
      case SINGLE_ANSWER_FIXED_RMSD:
        NextReal(p,
                 vArgs[i+1].end(),
                 nw_fixed_rmsd,
                 ",");
        break;
      case SINGLE_ANSWER_FIXED_N:
        {
          string n_str;
          NextToken(p,
                    vArgs[i+1].end(),
                    n_str,
                    ",");
          int max_num_matches_possible = MIN(aMol_f[0].size(), aMol_f[1].size());
          nw_fixed_n =
            InterpretNumMatches(vArgs[i+1], max_num_matches_possible);
        }
        break;
      case SINGLE_ANSWER_MAX_DIST:
        NextReal(p,
                 vArgs[i+1].end(),
                 nw_max_dist,
                 ",");
        break;
      case SINGLE_ANSWER_MANUAL_GAP:
        NextReal(p,
                 vArgs[i+1].end(),
                 nw_gap,
                 ",");
        break;
      default:
        assert(0);
        break;
      }

      NextInt(p,
              vArgs[i+1].end(),
              nw_max_iters,
              ",");

      num_arguments_deleted = 2;
    } //if (vArgs[i] == "-sa")


    //Delete all the arguments we have processed in this iteration (if any)
    if (num_arguments_deleted > 0)
    {
      assert(i+num_arguments_deleted <= vArgs.size());
      vArgs.erase(vArgs.begin()+i, vArgs.begin()+i+num_arguments_deleted);
      --i; //rewind by one so that when we invoke ++i
           //at the end of this loop, i's value will not change.
           //This will point us to the next un-read argument.
    }
  } // loop over arguments "for (int i=1; i < vArgs.size(); ++i)"

} //MinrmsParser::MinrmsParser()






static int InterpretNumMatches(string num_matches_str,
                               int    max_num_matches_possible)
{
  assert(max_num_matches_possible > 0);
  Real x = atof(num_matches_str.c_str());
  if ((x < 0.0) || (x > max_num_matches_possible))
    ERR("ERROR: Error in input parameter.\n"
                  "       It is nonsensical to match \""
                  << num_matches_str << "\" pairs of residues.\n"
                  "       from a molecule with " << max_num_matches_possible << " residues.\n"
                  "       This parameter should be either a positive integer,\n"
                  "       from 1 to " << max_num_matches_possible << " or a fraction\n"
                  "       between 0.0 and 1.0.\n"
                  "            (In the later case, the parameter is\n"
                  "       interpreted as the integer which most closely\n"
                  "       represents that fraction of the total number\n"
                  "       of residues in the shorter molecule (" << max_num_matches_possible << ").  In the\n"
                  "       event that the fraction your seeking is 1.0,\n"
                  "       (ie. the length of the entire molecule), in order to be parsed\n"
                  "       correctly, there must be a '.' period present,\n"
                  "       in order to distinguish it from the integer 1.)\n"
                  "       Aborting...\n");
  if (x == 1.0)
  {
    string::const_iterator p;
    for(p = num_matches_str.begin();
        ((*p != '.') && (*p != '\0'));
        ++p)
      {}

    if (*p == '.')
      return max_num_matches_possible;
    else
      return 1;
  }
  else if (x < 1.0)
    return (int)floor(x * max_num_matches_possible+0.5);
  else
    return (int)floor(x);
} // InterpretNumMatches()




