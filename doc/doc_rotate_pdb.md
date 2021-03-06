rotate_pdb
==========

## Syntax

```
   rotate_pdb matrix_file < original.pdb > rotated.pdb
```


## Overview:

rotate_pdb rotates and/or translates the structure stored in a
PDB-file by the amount specified in "matrix_file".
rotate_pdb was intended to be a utility to read minrms and
msf2stat3d output, but it can be used in more general ways.


## Background

Structural alignments from **minrms** have two sections:
1) an alignment between the residues of the two structures
2) the relative orientation (superposition) between the two
   structures which brings matched atoms into close proximity.

Use rotate_pdb to interpret section 2).


## Details

rotate_pdb reads a pdb-file from the standard-in, applies an
affine transformation (usually a rotation plus a translation)
to the positions of it's atoms (and het-atoms), and writes the
resulting pdb-file to the standard-out.
ONLY the ATOM and HETATOM records are effected


(All other records are blindly sent to standard-out.)

## Input File

**rotate_pdb** requires one argument:
An ascii file storing a 3x4 matrix in the following format:
```
M11  M12  M13  M14
M21  M22  M23  M24
M31  M32  M33  M34
```
This 3x4 matrix stores the transformation to apply to the structure.
This matrix is stored at the beginning of every MSF-file
generated by minrms, along with the filename of the structure
to apply the transformation to.
You will have to cut out the matrix from the MSF-file generated
by minrms manually and save it in a new file which you pass
to rotate_pdb.  By doing this you can superimpose the two
structures exactly as they were superimposed by minrms.
(Note: This matrix is also generated by the msf2stat3d program.)

Explanation of Matrix File Format:
The transformation of coordinates from X,Y,Z to X',Y',Z' is:
```
  X' = M11*X + M12*Y + M13*Z  +  M14
  Y' = M21*X + M22*Y + M23*Z  +  M24
  Z' = M31*X + M32*Y + M33*Z  +  M34
```
...where X,Y,Z denote the position of an atom
before rotation and translation, and X', Y', Z'
denote the position of that atom afterwards.