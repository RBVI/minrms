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
    "   msf2stat3d [-oGg] msf_file labelA labelB pdb_fileA pdb_fileB\n"
    "\n"
    "   msf2stat3d superimposes the structures in pdb_fileA and pdb_fileB\n"
    "in order to minimize the root-mean-squared-displacement (RMSD)\n"
    "between matched CA atoms according to the alignment stored\n"
    "in msf_file.  It then calculates several structural properties\n"
    "which are returned to the user.\n"
    "\n"
    "Input format:\n"
    "\n"
    " -msf_file is an MSF file containing an alignment between\n"
    "  two or more sequences."
    //"  Gaps in the alignment must be\n"
    //"  indicated by '.' characters in the MSF file.
    "  (See the MinRMS documentation\n"
    "  for more details on the MSF file format.)\n"
    " -pdb_fileA and pdb_fileB are Brookhaven protein data bank\n"
    "  files which should contain structures corresponding\n"
    "  to the two sequences in the msf_file prefixed by the\n"
    "  labelA and labelB identifiers, respectively.\n"
    "\n"
    "Sample output format:\n"
    "\n"
    "        ------------------------\n"
    "        N 100\n"
    "     RMSD 3.0\n"
    "    P_str 2.0e-08\n"
    " max_dist 8.0\n"
    "\n"
    "Transform Matrix to apply to structure: pdb1lyz.ent\n"
    "0.55 -0.29 -0.78 61.65 \n"
    "0.78 0.51 0.36 -23.12 \n"
    "0.30 -0.81 0.51 22.11 \n"
    "        ------------------------\n"
    "  Where:\n"
    "  -N is the number of matched pairs of residues in the alignment\n"
    "\n"
    "  -RMSD is the root-mean-squared-distance between\n"
    "   matched pairs of residues.\n"
    "\n"
    "  -P_str is the probability-score based on S_str, as published\n"
    "   by Levitt and Gerstein in Proc. Natl. Acad. Sci. USA Vol 95 (1998).\n"
    "\n"
    "  -max-dist is the maximum distance between any pair of alpha-carbon\n"
    "   atoms that were matched in the alignment.\n"
    "   (Note: This is not effected by the \"atoms_used.txt\" file.  See below.)\n"
    "\n"
    "  -Transformation Matrix Section:\n"
    "   indicates the relative position of the two molecules after optimal\n"
    "   superposition.  It specifies of a rotation and a translation to\n"
    "   apply to the positions of atoms in the indicated structure,\n"
    "   in order to minimize the RMSD between matched atoms in the alignment.\n"
    "   (Note: This matrix can be passed to as an argument rotate_pdb.)\n"
    "\n"
    "   Matrix Format:\n"
    "      When superimposing the structures, only the indicated structure\n"
    "   (in the above example: \"pdb1lyz.ent\"), is moved.\n"
    "   Let X,Y,Z denote the position of an atom in that structure\n"
    "   (\"pdb1lyz.ent\") before optimal superposition,\n"
    "   and let X', Y', Z' be the position of that atom afterwards.\n"
    "   The transformation of coordinates from X,Y,Z to X',Y',Z' is:\n"
    "   X' = M11*X + M12*Y + M13*Z  +  M14\n"
    "   Y' = M21*X + M22*Y + M23*Z  +  M24\n"
    "   Z' = M31*X + M32*Y + M33*Z  +  M34\n"
    "   where Mij is the 3x4 matrix returned to the user.\n" 
    "\n"
    "\n"
    "Optional:\n"
    "\n"
    " -o   The properties calculated by msf2stat3d depend on\n"
    "      the relative orientation between the two structures.\n"
    "      By default, msf2stat3d superimposes the structures to\n"
    "      minimize the RMSD between the positions of corresponding atoms.\n"
    "      The user can dissable this behavior by passing the -o flag,\n"
    "      which leaves the structures in their original orientations.\n"
    "      (By using the rotate_pdb program, or visualization software,\n"
    "       the user can move the the two structures to any position\n"
    "       they desire, and then invoke msf2stat3d -o on the rotated\n"
    "       structures.)\n"
    "\n"
    "Options for visualizing the alignment in Midas or Chimera:\n"
    "\n"
    " -g   If the -g command line option is specified, msf2stat3d creates\n"
    "      a graphics object file (with a \".gfx\" extension)\n"
    "      containing a wireframe 3-D image of the superposition between\n"
    "      the two structures, and the residues that were matched,\n"
    "      viewable in Midas or Chimera.\n"
    "      The two structures are displayed in white and cyan.\n"
    "      Corresponding atoms are connected by red-line-segments\n"
    "      (or small red markers if the line segments are too short\n"
    "       to see clearly).\n"
    "\n"
    " -G   The -G option is the same as the -g option, producing a\n"
    "      Midas graphics object file displaying the structural alignment.\n"
    "      However, if -G is used, only line segments connecting\n"
    "      matched atoms are drawn (in red).\n"
    "      The two superimposed structures are not drawn.\n"
    "      \n"
    "      Explanation of Usage:\n"
    "         The -G option allows the user the freedom to display these\n"
    "      structures by some other means (not necessarily in wireframe).\n"
    "      In Midas and Chimera you can display multiple structures\n"
    "      at once: the two structures being aligned, and a separate\n"
    "      graphics object file displaying the matched atoms.\n"
    "      To superimpose the two structures correctly, one would have\n"
    "      to generate an additional PDB-file, one whose structure was\n"
    "      rotated correctly to the other structure.  To do this, invoke\n"
    "      the rotate_pdb program on one of the two PDB files,\n"
    "      using the transformation-matrix output of msf2stat3d.\n"
    "      The new, rotated PDB-file could then be displayed simultaneously\n"
    "      along with the PDB-file containing the other structure, and the\n"
    "      graphics-object-file generated by msf2stat3d -G.\n"
    "      Users would then have the freedom to use the features of\n"
    "      Midas or Chimera allowing them to highlight certain residues\n"
    "      from the two molecules or display them in alternate ways,\n"
    "      like using ribbons or spheres.\n"
    "      This is not possible using one graphics object file alone\n"
    "      to represent the entire scene.\n"
    "\n"
    "Configuration files:\n"
    "\n"
    "  If an \"atoms_used.txt\" file is located in the directory in which\n"
    "  msf2stat3d is invoked, msf2stat3d will use the positions of the\n"
    "  atoms specified in that file (instead of \" CA\" atoms), in order\n"
    "  to calculate RMSD, and to physically superimpose the two\n"
    "  structures.  This format of this file is documented in the\n"
    "  minrms documentation.\n";



  if (! (((argc == 6) && (argv[1][0] != '-')) ||
         ((argc == 7) && (argv[1][0] == '-'))))
    ERR(help_msg);

  
  int argv_offset = 1;         //<-- Where to start parsing arguments.
                               //    (right after the program-name, argv[0])

  bool create_gfx_file = false;
  bool show_matches_only = false;
  bool optimize_orientation = true;

  if (argv[1][0] == '-')
  {
    argv_offset = 2;           //(skip program-name and "mode" argument)
    for (char const *p = argv[1]+1; *p != '\0'; ++p)
    {
      switch (*p)
      {
      case 'g':
        create_gfx_file = true;
        show_matches_only = false;
        break;
      case 'G':
        create_gfx_file = true;
        show_matches_only = true;
        break;
      case 'o':
      case 'O':
        optimize_orientation = false;
        break;
      default:
        ERR("Error: unrecognized option: \"-" << *p << "\".\n");
      }
    }
  } //if (argv[1][0] == '-')


  string msf_filename(argv[0 + argv_offset]);
  string labelA(argv[1 + argv_offset]);
  string labelB(argv[2 + argv_offset]);
  string pdb_filenameA(argv[3 + argv_offset]);
  string pdb_filenameB(argv[4 + argv_offset]);

  //First, specify atoms used in calculation.
  LoadBackboneAtomSymbols("atoms_used.txt"); 

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
  mA.Finalize();   //Tell mA and mB to grant special fast access
  mB.Finalize();   //to these atoms.
  NameCodeLookup::Init("res_code_dict.txt");//Needed for converting residue name codes
                         //from PDB to MSF format later on.

  //Open the MSF-file:
  PairwiseAlignment a;
  a.ImportMSF(msf_filename, labelA, labelB,
              mA, mB); // <<--check against the residues in the two PDB-files

  //Now figure out how to orient the second molecule.
  Matrix3x4 RT;
  if (optimize_orientation)
    a.CalcMinRMSD(mA, mB, RT); //rotate and translate mB to minimize RMSD
  else
    ResetToIdentity3x4(RT); //don't apply any rotation or translation to mB

  Biopolymer mB_rotated = mB;
  ApplyTransform(RT, mB, mB_rotated);
  cout << "        N " << a.size() << endl;
  cout << "     RMSD " << RMSD(a, mA, mB_rotated) << endl;
  cout << "    P_str " << LevittGerstein98::Sstr_Probability(a, mA, mB_rotated)
       << endl;
  Biopolymer::const_iterator pdummy1, pdummy2;
  cout << " max_dist " << a.FindFurthestAtomPair(pdummy1, pdummy2,
                                                 mA, mB_rotated,
                                                 " CA") //" CA" atoms only
       << endl;
  if (optimize_orientation)
    cout <<
      "Transform Matrix to apply to structure: "
         << pdb_filenameB << "\n"
         << RT[0][0] << " " << RT[0][1] << " " << RT[0][2] << " " << RT[0][3] << "\n"
         << RT[1][0] << " " << RT[1][1] << " " << RT[1][2] << " " << RT[1][3] << "\n"
         << RT[2][0] << " " << RT[2][1] << " " << RT[2][2] << " " << RT[2][3] << "\n"
         << flush;
  else
    cout <<
      "-o option used, so calculations were performed without rotating or\n"
      "translating either structure.  Consequently the RMSD may not be optimal."
         << endl;

  stringstream gfx_filename(ios::out);
  gfx_filename << msf_filename << ".gfx";
  if (create_gfx_file)
    a.ExportGFX(mA,
                mB_rotated,
                gfx_filename.str(),
                "",
                false,
                (! show_matches_only),
                PairwiseAlignment::MIDAS_MARKER_SIZE,
                true);
} //main()







