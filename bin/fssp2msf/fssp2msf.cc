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
#include <cassert>
using namespace std;

#include <global_utils.h>
#include <simple_numeric_utils.h>
#include <biopolymer.h>
#include <load_pdb.h>
#include <residue_names.h>
#include <pair_alignment.h>


const char g_version_string[] = "0.5";


using namespace minrms;

int
main(int argc, char **argv)
{
  char help_str[] = 
    "usage:\n"
    "   fssp2msf fssp_file pdb_file1 pdb_file2 labelA labelB [msf_labelA msf_labelB]\n"
    "\n"
    "   fssp2msf converts alignments between a pair of proteins in the\n"
    "FSSP format into an MSF file.  The FSSP file is read from the\n"
    "standard in, and the MSF file is saved as a file with the same name.\n"
    "as the fssp_file, plus an \".msf\" extension.\n"
    "-fssp_file is a file describing a structural alignment between several\n"
    " molecules in FSSP format.  (This is the format used by the DALI\n"
    " server.)  For more information, see: http:/"
    "/www2.ebi.ac.uk/dali/fssp\n"
    "-\"pdb_file1\" and \"pdb_file2\" are PDB files describing the two\n"
    " structures being aligned, and \"labelA\",\"labelB\" are the\n"
    " labels used to identify these structures inside the FSSP file.\n"
    "\n"
    "Optional:\n"
    "-The optional parameters \"msf_labelA\" and \"msf_labelB\",\n"
    " allow the user to specify separate labels that will be used to\n"
    " identify the same sequences in the MSF file that gets created\n"
    " If not specified, labelA and labelB are used.\n"
    "-The 1-letter codes used to represent each residue in the sequence\n"
    " can be customized by supplying a \"res_code_dict.txt\" file.\n"
    " (See minrms documentation on how to do this.)\n"
    "\n"
    "Limitations:\n"
    "   fssp2msf can not consider more than two proteins at a time.\n"
    "It is limited to generating MSF files which compare only two molecules.\n";

  if ((argc != 6) && (argc != 8))
    ERR("Syntax Error: wrong number of arguments passed to fssp2msf\n"
        "\n"
        << help_str);

  string fssp_filename(argv[1]);

  vector<string> pdb_filenames(2);
  pdb_filenames[0].append(argv[2]);
  pdb_filenames[1].append(argv[3]);

  vector<string> fssp_labels(2);
  fssp_labels[0].append(argv[4]);
  fssp_labels[1].append(argv[5]);

  vector<string> msf_labels(2);
  if (argc == 8)
  {
    msf_labels[0].append(argv[6]);
    msf_labels[1].append(argv[7]);
  }
  else //Default Behavior: If no MSF-labels were specified,
  {    //then use the same protein labels in both MSF and FSSP files.
    msf_labels[0] = fssp_labels[0];
    msf_labels[1] = fssp_labels[1];
  }

  //------ Now, load the two pdb-files into  -----
  //------ "Biopolymer" data structures. -----
  // Existing libraries I have written for reading and writing alignments of
  // various formats communicate structural and sequence-content information
  // about each molecule using "Biopolymer" data structures.
  // This is useful because, they contain the name, chainID,
  // seqNumber, and insertCode of every residue in each molecule.
  // This is necessary to keep track of because this is how FSSP
  // files indicate which residues are matched together.


  // *****************************************
  // ***
  // ***  Now, load in the two molecules:
  // ***
  // *****************************************

  //The first step is to specify atoms used in calculation
  //by invoking the LoadBackboneAtomSymbols() function.
  //Since there is no calculation, I'm not sure if I need to
  //do this.  I do it to be safe.
  LoadBackboneAtomSymbols("atoms_used.txt");

  Biopolymer aMolOrig[2];
  for (int m = 0; m < 2; ++m)
  {
    DEBUG_MSG(DBG_FILE_PARSER,
              "--Reading in the " << m << "th molecule.");
    ifstream pdb_file(pdb_filenames[m].c_str(), ios::in);
    if (! pdb_file)
      ERR_INTERNAL("Cannot open file \""
                   << pdb_filenames[m]
                   << "\" for reading.  Exiting...");
    //read in the file
    pdb_file >> aMolOrig[m];

    //This next line may be optional, since we are
    //not doing any calculations.  I do it to be safe.
    aMolOrig[m].Finalize();

  } //for (int m = 0; m < 2; ++m)
  NameCodeLookup::Init();//Needed for converting residue name codes
                         //from PDB to MSF format later on.


#if 0
  This comment used to be true, but it appears it is not
  true anymore.  There used to be an incompatibility between the
  way DALI and my software counted and numbered residues.  This doesnt appear
  to be true anymore for some reason.  Either way, I have added
  a check for these kind of incompatibilities later on,
  so we dont need to account for them here.

  This means there is no reason to maintain separate
  aMolOrig[] and aMolFiltered[] molecules anymore.  Oh well.

  //------ Filter The Molecules: ----
  //-After receiving an email from Lisa Holm that FSSP files ignore
  // residues which do not contain all 4 backbone atoms: N,CA,C,and O.
  // I am not aware that this fact appears in any FSSP documentation.
  //-To make sure we are considering the same set of residues they
  // are using, we pre-filter out residues which do not contain these
  // atoms using the "DeleteResiduesWithMissingAtoms()" function.
#endif //#if 0

  Biopolymer aMolFiltered[2];
  aMolFiltered[0] = aMolOrig[0];
  aMolFiltered[1] = aMolOrig[1];

#if 0 
  COMMENTING OUT: This feature not needed. 

  //DeleteResiduesWithMissingAtoms() expects an array of C-strings
  //containing the names of the atoms that must be present.
  //We must now create this array:
  char const *(aAtomsThatMustBePresent[5]);
  aAtomsThatMustBePresent[0] = " N";
  aAtomsThatMustBePresent[1] = " CA";
  aAtomsThatMustBePresent[2] = " C";
  aAtomsThatMustBePresent[3] = " O";
  aAtomsThatMustBePresent[4] = NULL; //The array is NULL terminated.
#endif //#if 0  

  for (int m = 0; m < 2; ++m)
  {

#if 0  
    COMMENTING OUT: This feature not needed. 

    DEBUG_MSG(DBG_FILE_PARSER,
              "--Deleting residues with missing backbone atoms from " << m << "th molecule.");
    DeleteResiduesWithMissingAtoms(aMolFiltered[m], aAtomsThatMustBePresent);
#endif //#if 0  

    if (aMolFiltered[m].size() < 2)
      ERR("ERROR: Molecule \"" << pdb_filenames[m] << "\"\n"
          "       contains fewer than two residues.\n"
          "       Aborting...\n");
  } //for (int m = 0; m < 2; ++m)

  // ----- Now, load in the alignment: -----
  PairwiseAlignment alignment;
  ifstream fssp_file(fssp_filename.c_str(), ios::in);
  if (! fssp_file)
    ERR("Error: Unable to open file \"" << fssp_filename << "\" for reading.");
  alignment.ImportFSSP(fssp_file,
                       aMolFiltered[0],
                       aMolFiltered[1],
                       fssp_labels[0],
                       fssp_labels[1]);

  //check for nonsensical input
  assert(0 < alignment.NumMatches());
  assert(alignment.NumMatches()
         <=
         MIN(aMolFiltered[0].size(), aMolFiltered[1].size()));

  //We need to create an MSF file which includes all of
  //residues, even the ones with mising bacbone atoms.
  //This function takes care of this:
  TranslateAlignmentBetweenMolecules(alignment,
                                     alignment,
                                     aMolOrig[0],
                                     aMolOrig[1],
                                     aMolFiltered[0],
                                     aMolFiltered[1]);

  // ----- Finally, generate the MSF-file: -----
  stringstream out_filename(ios::out);
  out_filename << fssp_filename << ".msf";
  alignment.ExportMSF(out_filename.str(),
                      aMolOrig[0],
                      aMolOrig[1],
                      msf_labels[0],
                      msf_labels[1],
                      "This MSF-file was converted from FSSP format.\n");
} //main(int argc, char **argv)


