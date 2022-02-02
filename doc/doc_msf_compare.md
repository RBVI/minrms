msf_compare
===========

## Usage

```
   msf_compare [metric] msf1 msf2 labelA1 labelB1 [labelA2 labelB2]
```


## Description

msf_compare returns a measure of similarity between two different
alignments between the same two sequences (or structures).
It only compares alignments between two sequences.
The two MSF-files msf1, and msf2, can contain many other sequences as
well, but msf_compare only calculates the similarity in how they align
the two sequences indicated by labelA1/2, and labelB1/2.  (See below)


## Input format

- These pairwise alignments are stored in MSF files: msf1 and msf2.
  They should contain alignments between the same two sequences.
  The MSF file format is explained
  [here](http://rothlab.ucdavis.edu/genhelp/chapter_2_using_sequences.html#_Using_Multiple_Sequence_Format_(MSF)
- *label arguments:*
  **msf_compare** compares two alignments between the same two sequences,
  Because, MSF files can contain alignments between many sequences,
  the user needs to specify which pair of sequences from each
  MSF-file is being aligned (even if only two are present).
  Then, the alignments can be compared.
  MSF-files contain labels at the beginning of every line to
  distinguish the sequences from eachother.
  The labelA1 and labelB1 arguments identify the two sequences from
  msf1 that are being aligned.  Likewise,
  the labelA2 and labelB2 arguments identify the two sequences from
  msf2 that are being aligned. (These later arguments are optional and
  if omitted, labelA1 and labelB1 will be used for both MSF files.)


## Output Format

The value(s) that msf_compare returns to the user
depend on the metric selected.  This is explained below:


## Metrics

### -E

***(default)***
If -E is selected, msf_compare returns the number of identical
residue equivalences common to BOTH alignments.

### -R
If -R is selected, msf_compare returns two integers:
i) The first integer specifies the number of residues from
  the first sequence (A) that were matched in both alignments.
ii) The second integer specifies the number of residues from
  the second sequence (B) that were matched in both alignments.

### -3d,pdbA,pdbB

The -3d metric compares the similarity of two structural alignments.
Unlike -E and -R which count the number of identical residues
in the two alignments, the -3d metric measures the
physical distance in angstroms between the two 3-D superpositions
implied by the alignments in msf1 and msf2.
The arguments "pdbA" and "pdbB" represent the names of
two pdb-files containing the structures corresponding to
sequenceA and sequenceB.  Once the two structures are known,
they are superimposed together twice:
i) With structure pdbA held fixed, structure pdbB
  is rotated and translated in order to minimize the
  intermolecular RMS-displacement between the CA atoms
  matched with CA atoms in pdbA according to the
  alignment in MSF-file msf1.
ii) Next, structure pdbB is superimposed with structure pdbA
  to minimize the RMSD of the alignment in msf2.

#### Output

msf_compare returns the RMS-displacement (in angstroms)
between the positions of the CA atoms in structure pdbB
rotated an translated in these two ways.

(Note:  Atom(s) other than the CA atoms may be used by supplying
"atoms_used.txt" file.  See minrms documentation for details.)


## Other metrics

### -e 

If -e is selected, msf_compare returns the total number of
equivaliences from EITHER alignment.
(This may not be interesting to most users.)


### -r

If -r is selected, msf_compare returns two integers:
i) The first integer specifies the number of residues from
  the first sequence (A) that were matched in EITHER alignment.
ii) The second integer specifies the number of residues from
  the second sequence (B) that were matched in EITHER alignment.


### -3d-both,pdbA,pdbB

(This is a little difficult to explain.)
The "-3d-both" metric is identical to "-3d metric", except that it returns
the RMS-displacement between the two different positions of the SUBSET
of the structure from pdbB that was matched in both alignments.
That is, with -3d-both selected, msf_compare considers only the CA atoms
from structure pdbB that were among those residues matched in both
alignment: msf1 and msf2.
The positions of atoms that were not matched in both alignments
are not considered in the calculation of displacement.


### -3d-either,pdbA,pdbB

The -3d-either metric is analogous to the -3d-both metric.
It computes the RMS-displacement between the CA atoms of
belonging to residues in the second structure that were matched
in EITHER alignment.


## Metrics provided for conveniance

### -E-over-e

Returns the fraction of the number of equivalences common to both
alignments divided by the number of equivalences from either alignment.
(This can be easily calculated from the quantities described above.)

### -E-over-N

Returns the fraction of the number of equivalences common to both
alignments divided by the number of equivalances in the larger alignment
That is, it returns E / max(N1,N2), where:
```
N1 = # equivalences in the alignment from msf1, and
N2 = # equivalences in the alignment from msf2
```

### -R-over-r

Returns the fraction of the size of the subset of residues
(from either sequence) that were matched in both alignments, 
divided by the size of the subset of residues (from either sequence)
that were matched in either alignment.  (This can be easily calculated
from the quantities described above.)

### -R-over-2N

Returns the fraction of the size of the subset of residues
(from either sequence) that were matched in both alignments, 
divided by the number of residues matched in the larger of the two
alignments.  That is, it returns, R / (2 x max(N1, N2)), where:
```
N1 = # equivalences in the alignment from msf1, and
N2 = # equivalences in the alignment from msf2
```
