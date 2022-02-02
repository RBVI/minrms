fssp2msf
========

## Usage

```
   fssp2msf fssp_file pdb_file1 pdb_file2 labelA labelB [msf_labelA msf_labelB]
```
fssp2msf converts alignments between a pair of proteins in the
FSSP format into an MSF file.  The FSSP file is read from the
standard in, and the MSF file is saved as a file with the same name.
as the fssp_file, plus an ".msf" extension.
- fssp_file is a file describing a structural alignment between several
  molecules in FSSP format.  (This is the format used by the DALI
  server.)  For more information, see: http://www2.ebi.ac.uk/dali/fssp
- "pdb_file1" and "pdb_file2" are PDB files describing the two
  structures being aligned, and "labelA","labelB" are the
  labels used to identify these structures inside the FSSP file.


## Optional
- The optional parameters "msf_labelA" and "msf_labelB",
  allow the user to specify separate labels that will be used to
  identify the same sequences in the MSF file that gets created
  If not specified, labelA and labelB are used.
- The 1-letter codes used to represent each residue in the sequence
  can be customized by supplying a "res_code_dict.txt" file.
  (See minrms documentation on how to do this.)


## Limitations

fssp2msf can not consider more than two proteins at a time.
It is limited to generating MSF files which compare only two molecules.
