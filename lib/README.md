## lib subdirectory

Each subdirectory contains shared code used by the
programs in the ../bin/ directory.

- **dyn**
  Implementation of
  both the original and
  Needleman & Wunsch
  Dynamic programming
  algorithms.

- **eigen**
  Finds eigenvectors and
  eigenvalues (taken from
  "numerical recipees in C".)
  *(To avoid license issues,
  I should probably replace this
  with "jacobi_pd" eventually.)*

- **fast_rot_metric**
  Calculates the 3-D orientation
  metric in O(1) time.

- **vect_3d**
  Matrices, vectors, and
  simple multiplication
  routines.

- **global_utils**
  DEBUG flags, ASSERT()
  printing DEBUG_MSG() messages,
  some binary file utilities.

- **biopolymer**
  Utils for reading and storing
  PDB files.  Also a random
  access data structure
  allowing fast access to
  positions of the CA atoms.

- **or_gen**
  Superimposes two structures
  using fragment matching,
  MSF-file reading, and
  Needleman & Wunsch filtering.

- **pair_alignment**
  Stores an alignment.
  Import/Export MSF files,
  FSSP files, and GFX files
  Calculates RMSD and P_str.

- **parse_utils**
  General utils for parsing
  argv and argc.  Handles
  configuration files.

- **superimpose**
  Uses Diamond's method (1988)
  to superimpose two structures
  minimizing the RMSD between
  them.

- **pdb++**
  The 2004 C++ friendly version 
  of the PDF parsing portion 
  of the OTF library.
