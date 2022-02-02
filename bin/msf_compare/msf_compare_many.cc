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

#ifdef SUPPORTS_STD_NAMESPACE
#include <fstream>
#include <sstream>
using namespace std;
#else
#include <fstream.h>
#include <sstream.h>
#endif

#include <global_utils.h>
#include <parse_utils.h>
#include <residue_names.h>
#include <load_pdb.h>
#include <pair_alignment.h>

using namespace minrms;
using namespace minrms::ParseUtils;


int
main(int argc, char const **argv)
{
  enum MSFCompareMode {PAIR,
                       PAIR_OVER_MAX_N1_N2,
                       RES_A,
                       RES_B,
                       RES_AB,
                       RES_AB_OVER_2_MAX_N1_N2,
                       ALL_3D,
                       BOTH_3D,
                       EITHER_3D};

  char help_msg[] =
    "Error in syntax.\n"
    "\n"
    "usage:\n"
    "   msf_compare_many [mode] msf1 maxN labelA1 labelB1 [labelA2 labelB2]\n"
    "\n"
    "returns a measure of similarity between two different pairwise\n"
    "alignments between the same two amino acid sequences or structures.\n"
    "\n"
    " Input format:\n"
    "-These pairwise alignments are stored in MSF files: msf1 and\n"
    " a file whose name is \"alignN.msf\", where N denotes a number\n"
    " between 2, and the maxN argument you specify (should be an integer)\n"
    " The two MSF-files must both contain alignments between the same two\n"
    " sequences.\n"
    "   (The two MSF-files may also contain alignments\n"
    "    between other sequences as well, but these\n"
    "    other sequences will be ignored.)\n"
    "-The user must manually specify the sequences to align from\n"
    " each MSF file.  MSF-files contain labels at the beginning\n"
    " of every line to distinguish the various amino-acid sequences\n"
    " from eachother.  The user needs to supply the labels that indicate\n"
    " which two sequences are being compared in each MSF-file.\n"
    " \"labelA1\" indicates the symbol representing the first sequence (A)\n"
    " of interest in the first MSF-file msf1.\n"
    " \"labelB1\" indicates the symbol representing the second sequence (B)\n"
    " of interest in  the first MSF-file msf1.\n"
    "   (Likewise labelA2 and labelB2 represent the labels for sequences A and B\n"
    "    from the second MSF-file, msf2.  They are optional and if they\n"
    "    are omitted, labelA1 and labelB1 are used for both msf-files.)\n"
    "-The the alignment in msf1 is compared against\n"
    " the alignment that is most similar to it, according to the metric\n"
    " chosen, from the set of alignments stored in the full range of\n"
    " alignN.msf files.\n"
    "\n"
    "Output Format:\n"
    "\n"
    "   The metric used for comparison depends on the optional \"mode\" argument.\n"
    "\n"
    "-E (default)\n"
    "   If -E is selected, msf_compare_many returns\n"
    "   a fraction specifying the number of\n"
    "   pairs-of-residues that were common to both alignments.\n"
    "\n"
    "-RA \n"
    "   If -RA is selected, msf_compare_many returns\n"
    "   a fraction specifying the number of distinct residues from\n"
    "   the first sequence that were common to both alignments.\n"
    "\n"
    "-RB \n"
    "   If -RB is selected, msf_compare_many returns\n"
    "   a fraction specifying the number of distinct residues from\n"
    "   the second sequence that were common to both alignments.\n"
    "\n"
    "-R\n"
    "   If -R is selected, msf_compare_many returns\n"
    "   a fraction specifying the number of distinct residues from\n"
    "   the two sequences that were common to both alignments.\n"
    "\n"
    "-E-over-N\n"
    "   Returns the fraction of the number of equivalences common to both\n"
    " alignments divided by the number of equivalances in the larger alignment\n"
    " That is, it returns E / max(N1,N2), where:\n"
    " N1 = # equivalences in the alignment from msf1, and\n"
    " N2 = # equivalences in the alignment from msf2\n"
    "\n"
    "-R-over-2N\n"
    "   Returns the fraction of the size of the subset of residues\n"
    " (from either sequence) that were matched in both alignments, \n"
    " divided by the number of residues matched in the larger of the two\n"
    " alignments.  That is, it returns, R / (2 x max(N1, N2)), where:\n"
    " N1 = # equivalences in the alignment from msf1, and\n"
    " N2 = # equivalences in the alignment from msf2\n"
    "\n"
    "-3d,pdbA,pdbB\n"
    "   Explained in man page for msf_compare.\n"
    "\n"
    "-3d-both,pdbA,pdbB\n"
    "   Explained in man page for msf_compare.\n"
    "\n"
    "-3d-either,pdbA,pdbB\n"
    "   Explained in man page for msf_compare.";


  if (! (((argc == 5) && (argv[1][0] != '-')) ||
         ((argc == 6) && (argv[1][0] == '-')) ||
         ((argc == 7) && (argv[1][0] != '-')) ||
         ((argc == 8) && (argv[1][0] == '-'))))
    ERR(help_msg);

  int argv_offset = 1;         //<-- Where to start parsing arguments.
                               //    (right after the program-name, argv[0])

  //Set the mode.
  //(and adjust the number of remaining arguments accordingly.)
  MSFCompareMode mode = PAIR;
  string pdb_filename_A, pdb_filename_B;
  if (argv[1][0] == '-')
  {
    argv_offset = 2;           //(skip program-name and "mode" argument)

    if (strcmp(argv[1], "-E") == 0)
      mode = PAIR;
    else if (strcmp(argv[1], "-E-over-N") == 0)
      mode = PAIR_OVER_MAX_N1_N2;
    else if (strcmp(argv[1], "-RA") == 0)
      mode = RES_A;
    else if (strcmp(argv[1], "-RB") == 0)
      mode = RES_B;
    else if (strcmp(argv[1], "-R") == 0)
      mode = RES_AB;
    else if (strcmp(argv[1], "-R-over-2N") == 0)
      mode = RES_AB_OVER_2_MAX_N1_N2;
    else if ((strncmp(argv[1], "-3d,", 4) == 0) ||
             (strncmp(argv[1], "-3d-both,", 9) == 0) ||
             (strncmp(argv[1], "-3d-either,", 11) == 0))
    {
      if (strncmp(argv[1], "-3d,", 4) == 0)
        mode = ALL_3D;
      else if (strncmp(argv[1], "-3d-both,", 9) == 0)
        mode = BOTH_3D;
      else if (strncmp(argv[1], "-3d-either,", 11) == 0)
        mode = EITHER_3D;
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
  string msf_filename1(argv[0+argv_offset]);
  int maxN = atoi(argv[1+argv_offset]);
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
    ASSERT((argc - argv_offset) == 4);


  //If necessary, load in the PDB-files into LinearMolecule data structures.
  LinearMolecule mA, mB;
  if ((mode == ALL_3D) || (mode == BOTH_3D) || (mode == EITHER_3D))
  {
    ifstream pdbA(pdb_filename_A.c_str(), ios::in);
    if (! pdbA)
      ERR("Error: unable to open file \"" << pdb_filename_A << "\" for reading.");
    ifstream pdbB(pdb_filename_B.c_str(), ios::in);
    if (! pdbB)
      ERR("Error: unable to open file \"" << pdb_filename_B << "\" for reading.");
    LoadQuickAtomSymbols("atoms_used.txt");
    pdbA >> mA; //load the PDB files
    pdbB >> mB;
    mA.Finalize();
    mB.Finalize();
    NameCodeLookup::Init();
  } //if ((mode == ALL_3D) || (mode == BOTH_3D) || (mode == EITHER_3D))

  //Now, loop over all the alignment files named "alignN.msf" for N=2...maxN
  //and compare them with msf_filename1, according to the desired criteria.
  int   bestN = -1;
  Real best_similarity = 0.0;
  Real best_rmsd_3d = 1000000000.0;
  for(int N = 2; N <= maxN; ++N)
  {
    //Now, figure out the filename of the second MSF-file.
    stringstream msf_filename2_stream(ios::out);
    msf_filename2_stream << "align" << N << ".msf";
    string msf_filename2 = msf_filename2_stream.str();
    bool msf2_found = false;
    {
      ifstream msf2(msf_filename2.c_str(), ios::in);
      msf2_found = msf2;
    }
    if (! msf2_found)
      cerr << "Warning: No file named \"" << msf_filename2 << "\" found."
           << endl;
    else
    {

      // ************************************************************
      // *****                 Bulletproofing.                  *****
      // ************************************************************


      // ** It is a good idea to check that the MSF-files
      // ** contain the same two sequences.

      //First, we read in all the sequences from the two msf-files;
      //(All the debug statements, and error checking make this code ugly. Alas.)
      string sequenceA1, sequenceB1, sequenceA2, sequenceB2;
      

      {
        DEBUG_MSG(DBG_IMPORT_MSF,"Comparing sequences from each MSF-file for consistency:\n"
                  "               Scanning in sequence \"" << labelA1 << "\"\n"
                  "               from file \"" << msf_filename1 << "\"");
        ifstream msf1(msf_filename1.c_str(), ios::in);
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
        ifstream msf1(msf_filename1.c_str(), ios::in);
        if (! ReadMSFsequence(msf1, labelB1, sequenceB1))
          ERR("Error: sequence \"" << labelB1 <<
              "\" not found in msf-file \"" << msf_filename1 << "\"");
      }

      {
        DEBUG_MSG(DBG_IMPORT_MSF, "Comparing sequences from each MSF-file for consistency:\n"
                  "               Scanning in sequence \"" << labelA2 << "\"\n"
                  "               from file \"" << msf_filename2 << "\"");
        ifstream msf2(msf_filename2.c_str(), ios::in);
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
        ifstream msf2(msf_filename2.c_str(), ios::in);
        if (! ReadMSFsequence(msf2, labelB2, sequenceB2))
          ERR("Error: sequence \"" << labelB2 <<
              "\" not found in msf-file \"" << msf_filename2 << "\"");
      }



      //Finally, check to see that the sequence content is the same.
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
      Real similarity;
      Real rmsd_3d;
      int numerator, denominator;
      int num_matches_in_both, num_matches_in_either;
      int num_res_in_both_A, num_res_in_both_B;
      int num_res_in_either_A, num_res_in_either_B;
      switch (mode)
      {
      case PAIR:
      case RES_A:
      case RES_B:
      case RES_AB:
      case PAIR_OVER_MAX_N1_N2:
      case RES_AB_OVER_2_MAX_N1_N2:
        a1.ImportMSF(msf_filename1, labelA1, labelB1);
        a2.ImportMSF(msf_filename2, labelA2, labelB2);
        num_matches_in_both = NumMatchesInBoth(a1, a2);
        num_matches_in_either = NumMatchesInEither(a1, a2);
        NumResiduesInBoth(a1, a2, num_res_in_both_A, num_res_in_both_B);
        NumResiduesInEither(a1, a2, num_res_in_either_A, num_res_in_either_B);
        switch (mode)
        {
        case PAIR:
          numerator = num_matches_in_both;
          denominator = num_matches_in_either;
          break;
        case RES_A:
          numerator = num_res_in_both_A;
          denominator = num_res_in_either_A;
          break;
        case RES_B:
          numerator = num_res_in_both_B;
          denominator = num_res_in_either_B;
          break;
        case RES_AB:
          numerator = num_res_in_both_A + num_res_in_both_B;
          denominator = num_res_in_either_A + num_res_in_either_B;
          break;
        case PAIR_OVER_MAX_N1_N2:
          numerator = num_matches_in_both;
          denominator = MAX(a1.size(), a2.size());
          break;
        case RES_AB_OVER_2_MAX_N1_N2:
          numerator = num_res_in_both_A + num_res_in_both_B;
          denominator = 2 * MAX(a1.size(), a2.size());
          break;
        default:
          ASSERT_NOT_REACHED();
          break;
        }
        similarity = ((float)numerator)/((float)denominator);

        cout << "N = " << N << ", similarity = " << similarity << endl;
        cout << "  num_matches_in_both = " << num_matches_in_both
             << ",   num_matches_in_either = " << num_matches_in_either << endl;
        cout << "  num_res_in_both_A = " << num_res_in_both_A
             << ",   num_res_in_both_B = " << num_res_in_both_B << endl;
        cout << "  num_res_in_either_A = " << num_res_in_either_A
             << ",   num_res_in_either_B = " << num_res_in_either_B << endl;

        if (similarity > best_similarity)
        {
          bestN = N;
          best_similarity = similarity;
        }
        break;

      case ALL_3D:
      case BOTH_3D:
      case EITHER_3D:
        {
          a1.ImportMSF(msf_filename1, labelA1, labelB1, mA, mB); //check the sequences against
          a2.ImportMSF(msf_filename2, labelA2, labelB2, mA, mB); //residues in the PDB files.
          switch (mode)
            {
            case ALL_3D:
              rmsd_3d = Compare3dAll(a1, a2, mA, mB);
              break;
            case BOTH_3D:
              rmsd_3d = Compare3dBoth(a1, a2, mA, mB);
              break;
            case EITHER_3D:
              rmsd_3d = Compare3dEither(a1, a2, mA, mB);
              break;
            default:
              ASSERT_NOT_REACHED();
            } //switch (mode)
          cout << "N = " << N << ", 3D_RMSD = " << rmsd_3d << endl;
          if (rmsd_3d < best_rmsd_3d)
          {
            bestN = N;
            best_rmsd_3d = rmsd_3d;
          }
          break; //case ALL_3D, BOTH_3D, 3D_EITHER:
        }
      default:
        ASSERT_NOT_REACHED();
        break;
      } //switch (mode)
    } //if (! msf2_found)
  } //for(int N = 2; N < maxN; ++N)

  cout << "   -----------------------" << endl;

  switch (mode)
  {
  case PAIR:
  case RES_A:
  case RES_B:
  case RES_AB:
  case PAIR_OVER_MAX_N1_N2:
  case RES_AB_OVER_2_MAX_N1_N2:
    cout << "N = " << bestN << ", highest_similarity = " << best_similarity << endl;
    break;
  case ALL_3D:
  case BOTH_3D:
  case EITHER_3D:
    cout << "N = " << bestN << ", lowest_3D_RMSD = " << best_rmsd_3d << endl;
    break;
  default:
    ASSERT_NOT_REACHED();
    break;
  }
} //main()

