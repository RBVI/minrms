"README.txt" April 22nd, 1999

The following list contains the description of files copied to 
the "shared" directory:


-----------------------------------------------------------
1) atoms_used.txt

		This file is a configuration file used by all binaries that
use the "linear_molecule" library.  This includes minrms.
This file contains a list of PDB atom-name codes that indicate
that certain atoms that are earmarked for special treatment.
For details, see the that program's documentation.
		(In the case of minrms, the atoms specified in
	this file are the only atoms used to calculate the
	root-mean-squared-distance (RMSD) of the structural
	alignment being made.)
		The "atoms_used.txt" file is usually optional.  It must be present
in the same directory where the program is invoked.  If it is not supplied,
then a default behavior is assumed (in the case of minrms, the alpha-carbons
are used).

Example:
If minrms were invoked in the same directory containing the atoms_used.txt
file shown below:
(note, the quote (") marks were added and should not appear in the file)
" N"
" C"
" CA"
" O"
then the nitrogen, alpha-carbon, carbon, and oxygen atoms
(the atoms along a protein's backbone) would be used to calculate RMSD
between residues in the alignment.


-----------------------------------------------------------
2) res_code_dict.txt

		This file is another configuration file used by all binaries
that use the "linear_molecule" library.  Again, this includes minrms.
This file contains a list of 3-letter PDB residue-name codes and
their 1-letter residue-name-code equivalents.
		In the case of the minrms alignment program, this allows the user
to manually change the letters used in MSF-output-files to represent the
residues that appear in the PDB-files that minrms reads as input.
A sample "res_code_dict.txt" file is shown below.
"unknown x"
"GLY  G"
"ALA  A"
"SER  S"
"CYS  C"
 .   .
 .   .
 .   .
(Of course, the quote marks (") and periods (.) do not appear in the actual
file.)  In this example residues named "GLY" in the PDB-file are represented
by the letter "G" in the MSF-file.  The first line, "unknown x" indicates
that residues whose three-letter-code is not recognized will be represented
by the letter "x" in the MSF-file.
		The "res_code_dict.txt" file is usually optional.
It must be present in the same directory where the program is invoked.
If it is not supplied, then a default translation is used.

