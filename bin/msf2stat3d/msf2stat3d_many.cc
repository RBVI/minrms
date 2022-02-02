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
using namespace std;

#include <global_utils.h>
#include <residue_names.h>
#include <load_pdb.h>
#include <apply_trans.h>
#include <pair_alignment.h>

using namespace minrms;


int
main(int argc, char const **argv)
{
  char help_msg[] =
    "Error in syntax.\n"
    "\n"
    "usage:\n"
    "   msf2stats3d_many maxN labelA labelB pdb_fileA pdb_fileB\n"
    "\n"
    "No description of this program is available.\n"
    "This program was not intended for public use.\n";


  if (argc != 6)
    ERR(help_msg);

  
  string maxNstr(argv[1]);
  int maxN = atoi(maxNstr.c_str());

  string labelA(argv[2]);
  string labelB(argv[3]);
  string pdb_filenameA(argv[4]);
  string pdb_filenameB(argv[5]);

  //First, specify atoms used in calculation.
  LoadQuickAtomSymbols("atoms_used.txt"); 

  //Open the PDB-files:
  Biopolymer mA, mB;
  ifstream pdbA(pdb_filenameA.c_str(), ios::in);
  if (! pdbA)
    ERR("Error: unable to open file \"" << pdb_filenameA << "\" for reading.");
  ifstream pdbB(pdb_filenameB.c_str(), ios::in);
  if (! pdbB)
    ERR("Error: unable to open file \"" << pdb_filenameB << "\" for reading.");

  pdbA >> mA;
  pdbB >> mB;
  mA.Finalize();
  mB.Finalize();
  NameCodeLookup::Init();//Needed for converting residue name codes
                         //from PDB to MSF format later on.

  int  best_N = -1;
  Real best_Pstr = 1.0;

  //Loop over all the msf-files containing more than
  //20 matched pairs of residues.
  //(In our opinion, P_str is not accurately calculated for alignments with
  // fewer matched pairs.)

  for(int N = 20; N <= maxN; ++N)
  {
    //Now, figure out the filename of the second MSF-file.
    stringstream msf_filename_stream(ios::out);
    msf_filename_stream << "align" << N << ".msf";
    string msf_filename = msf_filename_stream.str();
    bool msf_found = false;
    {
      ifstream msf(msf_filename.c_str(), ios::in);
      msf_found = msf;
    }
    if (! msf_found)
      cerr << "Warning: No file named \"" << msf_filename << "\" found."
           << endl;
    else
    {
      PairwiseAlignment a;
      a.ImportMSF(msf_filename, labelA, labelB,
                  mA, mB);//<<--check against the residues in the two PDB-files

      Matrix3x4 RT;
      a.CalcMinRMSD(mA, mB, RT);
      Biopolymer mB_rotated = mB;
      ApplyTransform(RT, mB, mB_rotated);
      Real Pstr = LevittGerstein98::Sstr_Probability(a, mA, mB_rotated);
      cout << "N = " << a.size() << ", P_str = " << Pstr << endl;
      if (Pstr < best_Pstr)
      {
        best_Pstr = Pstr;
        best_N = N;
      }
    }
  }
  cout << "N = " << best_N << ", lowest_P_str = " << best_Pstr << endl;
} //main()

