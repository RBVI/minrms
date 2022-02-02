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
#include <cstring>
#include <cassert>
using namespace std;

#include <global_utils.h>
#include <simple_numeric_utils.h>
#include <parse_utils.h>
#include <residue_names.h>
#include <load_pdb.h>
#include <pair_alignment.h>

using namespace minrms;
using namespace ParseUtils;


int
main(int argc, char const **argv)
{
  enum MSFCompareMetric {PAIR_BOTH,
                         PAIR_EITHER,
                         PAIR_BOTH_OVER_PAIR_EITHER,
                         PAIR_BOTH_OVER_MAX_N1_N2,
                         RES_BOTH,
                         RES_EITHER,
                         RES_BOTH_OVER_RES_EITHER,
                         RES_BOTH_OVER_2_MAX_N1_N2,
                         ALL_3D,
                         BOTH_3D,
                         EITHER_3D};

  char help_msg[] =
    "Error in syntax.\n"
    "\n"
    "usage:\n"
    "   msf_compare [metric] msf1 msf2 labelA1 labelB1 [labelA2 labelB2]\n"
    "\n"
    "msf_compare returns a measure of similarity between two different\n"
    "alignments between the same two sequences (or structures).\n"
    "It only compares alignments between two sequences.\n"
    "The two MSF-files msf1, and msf2, can contain many other sequences as\n"
    "well, but msf_compare only calculates the similarity in how they align\n"
    "the two sequences indicated by labelA1/2, and labelB1/2.  (See below)\n"
    "\n"
    " Input format:\n"
    "\n"
    "-The two pairwise alignments are stored in MSF files: msf1 and msf2.\n"
    " They should contain alignments between the same two sequences.\n"
    //" Gaps in the alignment must be indicated by '.' characters in\n"
    //" the MSF file."
    " (See MinRMS documentation for more details on the\n"
    "  MSF file format.)\n"
    "-label arguments:\n"
    " msf_compare compares two alignments between only two sequences,\n"
    " Because MSF files can contain alignments between many sequences,\n"
    " the user needs to specify which pair of sequences from\n"
    " each MSF-file is being aligned.  Then, the alignments can be compared.\n"
    "    MSF-files contain labels at the beginning of every line to\n"
    " distinguish the sequences from eachother.\n"
    " The labelA1 and labelB1 arguments identify the two sequences from\n"
    " msf1 that are being aligned.  Likewise,\n"
    " the labelA2 and labelB2 arguments identify the two sequences from\n"
    " msf2 that are being aligned. (These later arguments are optional and\n"
    " if omitted, labelA1 and labelB1 will be used for both MSF files.)\n"
    "\n"
    "   Output Format:\n"
    "\n"
    "   The value(s) that msf_compare returns to the user\n"
    "depend on the metric selected.  This is explained below:\n"
    "\n"
    "   the \"metric\" argument:\n"
    "\n"
    " ----- metrics: ----- \n"
    "\n"
    "-E (default)\n"
    "   If -E is selected, msf_compare returns the number of identical\n"
    "   residue equivalences common to BOTH alignments.\n"
    "\n"
    "-R\n"
    "   If -R is selected, msf_compare returns two integers:\n"
    "    i)The first integer specifies the number of same residues from\n"
    "      the first sequence (A) that were matched in both alignments.\n"
    "   ii)The second integer specifies the number of same residues from\n"
    "      the second sequence (B) that were matched in both alignments.\n"
    "\n"
    "-3d,pdbA,pdbB\n"
    "   The -3d metric compares the similarity of two structural alignments.\n"
    "   Unlike -E and -R which count the number of identical residues\n"
    "   in the two alignments, the -3d metric measures the\n"
    "   physical distance in angstroms between the two 3-D superpositions\n"
    "   implied by the alignments in msf1 and msf2.  The arguments\n"
    "   \"pdbA\" and \"pdbB\" represent the names of two Brookhaven\n"
    "   protein databank files containing the structures corresponding\n"
    "   to sequenceA and sequenceB.  Once the two structures are known,\n"
    "   they are superimposed together twice:\n"
    "    i)With structure pdbA held fixed, structure pdbB\n"
    "      is rotated and translated in order to minimize the\n"
    "      intermolecular RMS-displacement between the CA atoms\n"
    "      matched with CA atoms in pdbA according to the\n"
    "      alignment in MSF-file msf1.\n"
    "   ii)Next, structure pdbB is superimposed with structure pdbA\n"
    "      to minimize the RMSD of the alignment in msf2.\n"
    "  Output:\n"
    "      msf_compare returns the RMS-displacement (in angstroms)\n"
    "      between the positions of the CA atoms in structure pdbB\n"
    "      rotated an translated in these two ways.\n"
    "\n"
    "  (Note:  Atom(s) other than the CA atoms may be used by supplying\n"
    "          \"atoms_used.txt\" file.  See minrms documentation for details.)\n"
    "\n"
    " ----- other metrics: ----- \n"
    "\n"
    "-e \n"
    "   If -e is selected, msf_compare returns the total number of\n"
    "   equivaliences from EITHER alignment.\n"
    "   (This may not be interesting to most users.)\n"
    "\n"
    "-r\n"
    "   If -r is selected, msf_compare returns two integers:\n"
    "    i)The first integer specifies the number of residues from\n"
    "      the first sequence (A) that were matched in EITHER alignment.\n"
    "   ii)The second integer specifies the number of residues from\n"
    "      the second sequence (B) that were matched in EITHER alignment.\n"
    "\n"
    "-3d-both,pdbA,pdbB\n"
    "   (This is a little difficult to explain.)\n"
    "   The -3d-both metric is identical to -3d metric, except that it returns\n"
    "   the RMS-displacement between the two different positions of the SUBSET\n"
    "   of the structure from pdbB that was matched in both alignments.\n"
    "   That is, with -3d-both selected, msf_compare considers only the CA atoms\n"
    "   from structure pdbB that were among those residues matched in both\n"
    "   alignment: msf1 and msf2.\n"
    "   The positions of atoms that were not matched in both alignments\n"
    "   are not considered in the calculation of displacement.\n"
    "\n"
    "-3d-either,pdbA,pdbB\n"
    "   The -3d-either metric is analogous to the -3d-both metric.\n"
    "   It computes the RMS-displacement between the CA atoms of\n"
    "   belonging to residues in the second structure that were matched\n"
    "   in EITHER alignment.\n"
    //    "   Just change the word \"both\" to the word \"either\"\n"
    //    "   in the explanation above.\n"
    "\n"
    " ----- metrics provided for conveniance: ----- \n"
    "\n"
    "-E-over-e\n"
    "   Returns the fraction of the number of equivalences common to both\n"
    " alignments divided by the number of equivalences from either alignment.\n"
    " (This can be easily calculated from the quantities described above.)\n"
    "\n"
    "-E-over-N\n"
    "   Returns the fraction of the number of equivalences common to both\n"
    " alignments divided by the number of equivalances in the larger alignment\n"
    " That is, it returns E / max(N1,N2), where:\n"
    " N1 = # equivalences in the alignment from msf1, and\n"
    " N2 = # equivalences in the alignment from msf2\n"
    "\n"
    "-R-over-r\n"
    "   Returns the fraction of the size of the subset of residues\n"
    " (from either sequence) that were matched in both alignments, \n"
    " divided by the size of the subset of residues (from either sequence)\n"
    " that were matched in either alignment.  (This can be easily calculated\n"
    " from the quantities described above.)\n"
    "\n"
    "-R-over-2N\n"
    "   Returns the fraction of the size of the subset of residues\n"
    " (from either sequence) that were matched in both alignments, \n"
    " divided by the number of residues matched in the larger of the two\n"
    " alignments.  That is, it returns, R / (2 x max(N1, N2)), where:\n"
    " N1 = # equivalences in the alignment from msf1, and\n"
    " N2 = # equivalences in the alignment from msf2\n";


  if (! (((argc == 5) && (argv[1][0] != '-')) ||
         ((argc == 6) && (argv[1][0] == '-')) ||
         ((argc == 7) && (argv[1][0] != '-')) ||
         ((argc == 8) && (argv[1][0] == '-'))))
    ERR(help_msg);

  int argv_offset = 1;         //<-- Where to start parsing arguments.
                               //    (right after the program-name, argv[0])

  //Set the metric.
  //(and adjust the number of remaining arguments accordingly.)
  MSFCompareMetric metric = PAIR_BOTH;
  string pdb_filename_A, pdb_filename_B;
  if (argv[1][0] == '-')
  {
    argv_offset = 2;           //(skip program-name and "metric" argument)

    if (strcmp(argv[1], "-E") == 0)
      metric = PAIR_BOTH;
    else if (strcmp(argv[1], "-e") == 0)
      metric = PAIR_EITHER;
    else if (strcmp(argv[1], "-E-over-e") == 0)
      metric = PAIR_BOTH_OVER_PAIR_EITHER;
    else if (strcmp(argv[1], "-E-over-N") == 0)
      metric = PAIR_BOTH_OVER_MAX_N1_N2;
    else if (strcmp(argv[1], "-R") == 0)
      metric = RES_BOTH;
    else if (strcmp(argv[1], "-r") == 0)
      metric = RES_EITHER;
    else if (strcmp(argv[1], "-R-over-r") == 0)
      metric = RES_BOTH_OVER_RES_EITHER;
    else if (strcmp(argv[1], "-R-over-2N") == 0)
      metric = RES_BOTH_OVER_2_MAX_N1_N2;
    else if ((strncmp(argv[1], "-3d,", 4) == 0) ||
             (strncmp(argv[1], "-3d-both,", 9) == 0) ||
             (strncmp(argv[1], "-3d-either,", 11) == 0))
    {
      if (strncmp(argv[1], "-3d,", 4) == 0)
        metric = ALL_3D;
      else if (strncmp(argv[1], "-3d-both,", 9) == 0)
        metric = BOTH_3D;
      else if (strncmp(argv[1], "-3d-either,", 11) == 0)
        metric = EITHER_3D;
      const string arg(argv[1]);
      string::const_iterator pc = arg.begin();
      NextToken(pc, arg.end(), pdb_filename_A, string(",")); //<--skip the "-3d,"
      NextToken(pc, arg.end(), pdb_filename_A, ","); //parse in the two
      NextToken(pc, arg.end(), pdb_filename_B, ","); //pdb filenames
      if ((pdb_filename_A.size() == 0) || (pdb_filename_B.size() == 0))
        ERR("Missing the two comma separated pdb-file names\n"
            "that follow the \"-3d\" flag.");
    }
    else
      ERR("Error: unrecognized option: \"" << argv[1] << "\".\n");
  } //if (argv[1][0] == '-')

  //Open the MSF-files:
  char const *msf_filename1 = argv[0+argv_offset];
  char const *msf_filename2 = argv[1+argv_offset];
  //Parse in the MSF-file labels.
  char const *labelA1 = argv[2+argv_offset];
  char const *labelB1 = argv[3+argv_offset];
  char const *labelA2 = labelA1;
  char const *labelB2 = labelB1;
  if ((argc - argv_offset) == 6)
  {
    labelA2 = argv[4+argv_offset];
    labelB2 = argv[5+argv_offset];
  }
  else
    assert((argc - argv_offset) == 4);


  // ************************************************************
  // *****                 Bulletproofing.                  *****
  // ************************************************************


  // ** It is a good idea to check that the MSF-files
  // ** contain the same two sequences.

  //First, we read in all the sequences from the two msf-files;
  //(All the debug statements, and error checking make this code ugly. Alas.)
  string sequenceA1, sequenceB1, sequenceA2, sequenceB2;

  {
    DEBUG_MSG(DBG_IMPORT_MSF, "Comparing sequences from each MSF-file for consistency:\n"
              "               Scanning in sequence \"" << labelA1 << "\"\n"
              "               from file \"" << msf_filename1 << "\"");
    ifstream msf1(msf_filename1, ios::in);
    if (! msf1)
      ERR("Unable to open file: \"" << msf_filename1 << "\" for reading.");
    if (! ReadMSFsequence(msf1, labelA1, sequenceA1))
      //Check to see that the MSF-file contains
      //the labels that the user has specified. If not, print error.
      ERR("Error: sequence \"" << labelA1 <<
          "\" not found in msf-file \"" << msf_filename1 << "\"");
  }
  {
    DEBUG_MSG(DBG_IMPORT_MSF, "Comparing sequences from each MSF-file for consistency:\n"
              "               Scanning in sequence \"" << labelB1 << "\"\n"
              "               from file \"" << msf_filename1 << "\"");
    ifstream msf1(msf_filename1, ios::in);
    if (! ReadMSFsequence(msf1, labelB1, sequenceB1))
      ERR("Error: sequence \"" << labelB1 <<
          "\" not found in msf-file \"" << msf_filename1 << "\"");
  }

  {
    DEBUG_MSG(DBG_IMPORT_MSF, "Comparing sequences from each MSF-file for consistency:\n"
              "               Scanning in sequence \"" << labelA2 << "\"\n"
              "               from file \"" << msf_filename2 << "\"");
    ifstream msf2(msf_filename2, ios::in);
    if (! msf2)
      ERR("Unable to open file: \"" << msf_filename2 << "\" for reading.");
    if (! ReadMSFsequence(msf2, labelA2, sequenceA2))
      ERR("Error: sequence \"" << labelA2 <<
          "\" not found in msf-file \"" << msf_filename2 << "\"");
  }
  {
    DEBUG_MSG(DBG_IMPORT_MSF, "Comparing sequences from each MSF-file for consistency:\n"
              "               Scanning in sequence \"" << labelB2 << "\"\n"
              "               from file \"" << msf_filename2 << "\"");
    ifstream msf2(msf_filename2, ios::in);
    if (! ReadMSFsequence(msf2, labelB2, sequenceB2))
      ERR("Error: sequence \"" << labelB2 <<
          "\" not found in msf-file \"" << msf_filename2 << "\"");
  }



  // After reading the sequence content,
  // check to see that the sequence content is the same.
  if (sequenceA1 != sequenceA2)
    ERR("Error:  Sequences are not the same:\n"
        "        The first sequence appearing in file\n"
        "        \""<< msf_filename1 << "\"\n"
        "        is not the same as the first sequence appearing in file\n"
        "        \""<< msf_filename2 << "\"");

  if (sequenceB1 != sequenceB2)
    ERR("Error: Sequences are not the same:\n"
        "       The second sequence appearing in file\n"
        "       \""<< msf_filename1 << "\"\n"
        "       is not the same as the second sequence appearing in file\n"
        "       \""<< msf_filename2 << "\"");

  // ***************************************
  // ****     Read in msf1 and msf2,    ****
  // ****  calculate and print result.  ****
  // ***************************************

  PairwiseAlignment a1, a2;

  //Finally, perform the calculation and send the result to the standard out.
  switch (metric)
  {

  case PAIR_BOTH:
    a1.ImportMSF(msf_filename1, labelA1, labelB1);
    a2.ImportMSF(msf_filename2, labelA2, labelB2);
    cout << NumMatchesInBoth(a1, a2) << endl;
    break;

  case PAIR_EITHER:
    a1.ImportMSF(msf_filename1, labelA1, labelB1);
    a2.ImportMSF(msf_filename2, labelA2, labelB2);
    cout << NumMatchesInEither(a1, a2) << endl;
    break;

  case PAIR_BOTH_OVER_PAIR_EITHER:
    {
      a1.ImportMSF(msf_filename1, labelA1, labelB1);
      a2.ImportMSF(msf_filename2, labelA2, labelB2);
      int num_matches_both   =  NumMatchesInBoth(a1, a2);
      int num_matches_either =  NumMatchesInEither(a1, a2);
      if (num_matches_either == 0)
        ERR("No matches were made in either alignment.  Aborting.");
      cout << ((float)num_matches_both)/((float)num_matches_either) << endl;
    }
    break;

  case PAIR_BOTH_OVER_MAX_N1_N2:
    {
      a1.ImportMSF(msf_filename1, labelA1, labelB1);
      a2.ImportMSF(msf_filename2, labelA2, labelB2);
      int num_matches_both   =  NumMatchesInBoth(a1, a2);
      int max_num_pairs = MAX(a1.size(), a2.size());
      if (max_num_pairs == 0)
        ERR("No matches were made in either alignment.  Aborting.");
      cout << ((float)num_matches_both)/((float)max_num_pairs) << endl;
    }
    break;

  case RES_BOTH:
    {
      a1.ImportMSF(msf_filename1, labelA1, labelB1);
      a2.ImportMSF(msf_filename2, labelA2, labelB2);
      int num_res_in_both_A, num_res_in_both_B;
      NumResiduesInBoth(a1, a2, num_res_in_both_A, num_res_in_both_B);
      cout << num_res_in_both_A << " " << num_res_in_both_B << endl;
    }
    break;

  case RES_EITHER:
    {
      a1.ImportMSF(msf_filename1, labelA1, labelB1);
      a2.ImportMSF(msf_filename2, labelA2, labelB2);
      int num_res_in_either_A, num_res_in_either_B;
      NumResiduesInEither(a1, a2, num_res_in_either_A, num_res_in_either_B);
      cout << num_res_in_either_A << " " << num_res_in_either_B << endl;
    }
    break;

  case RES_BOTH_OVER_RES_EITHER:
    {
      a1.ImportMSF(msf_filename1, labelA1, labelB1);
      a2.ImportMSF(msf_filename2, labelA2, labelB2);
      int num_res_in_both_A, num_res_in_both_B;
      int num_res_in_either_A, num_res_in_either_B;
      NumResiduesInBoth(a1, a2, num_res_in_both_A, num_res_in_both_B);
      NumResiduesInEither(a1, a2, num_res_in_either_A, num_res_in_either_B);
      int num_res_both   = num_res_in_both_A   + num_res_in_both_B;
      int num_res_either = num_res_in_either_A + num_res_in_either_B;
      if (num_res_either == 0)
        ERR("No matches were made in either alignment.  Aborting.");
      cout << ((float)num_res_both)/((float)num_res_either) << endl;
    }
    break;

  case RES_BOTH_OVER_2_MAX_N1_N2:
    {
      a1.ImportMSF(msf_filename1, labelA1, labelB1);
      a2.ImportMSF(msf_filename2, labelA2, labelB2);
      int num_res_in_both_A, num_res_in_both_B;
      NumResiduesInBoth(a1, a2, num_res_in_both_A, num_res_in_both_B);
      int num_res_both   = num_res_in_both_A   + num_res_in_both_B;
      int max_num_pairs = MAX(a1.size(), a2.size());
      if (max_num_pairs == 0)
        ERR("No matches were made in either alignment.  Aborting.");
      cout << ((float)num_res_both)/((float)(2*max_num_pairs)) << endl;
    }
    break;

  case ALL_3D:
  case BOTH_3D:
  case EITHER_3D:
    {
      ifstream pdbA(pdb_filename_A.c_str(), ios::in);
      if (! pdbA)
        ERR("Error: unable to open file \"" << pdb_filename_A << "\" for reading.");
      ifstream pdbB(pdb_filename_B.c_str(), ios::in);
      if (! pdbB)
        ERR("Error: unable to open file \"" << pdb_filename_B << "\" for reading.");

      LoadBackboneAtomSymbols("atoms_used.txt");

      Biopolymer mA, mB;
      pdbA >> mA; //load the PDB files
      pdbB >> mB;

      mA.Finalize();
      mB.Finalize();
      NameCodeLookup::Init();
      a1.ImportMSF(msf_filename1, labelA1, labelB1, mA, mB); //check the sequences against
      a2.ImportMSF(msf_filename2, labelA2, labelB2, mA, mB); //residues in the PDB files.
      switch (metric)
      {
      case ALL_3D:
        cout << Compare3dAll(a1, a2, mA, mB) << endl;
        break;
      case BOTH_3D:
        cout << Compare3dBoth(a1, a2, mA, mB) << endl;
        break;
      case EITHER_3D:
        cout << Compare3dEither(a1, a2, mA, mB) << endl;
        break;
      default:
        assert(0);
      } //switch (metric)
      break; //case ALL_3D, BOTH_3D, EITHER_3D:
    }
  default:
    assert(0);
    break;
  } //switch (metric)

} //main()

