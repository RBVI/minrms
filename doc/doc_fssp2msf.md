fssp2msf
========

# THIS PROGRAM IS DEPRECIATED

***As of 2022-2-02 The
[DALI server](https://web.archive.org/web/20090301064750/http://ekhidna.biocenter.helsinki.fi/dali_server/start)
and the FSSP file format it uses have been depreciated.
There is probably no longer a need for this program.
This program may be removed in the future.
-Andrew***


## Usage

```
fssp2msf fssp_file pdb_file1 pdb_file2 labelA labelB [msf_labelA msf_labelB]
```
fssp2msf converts alignments between a pair of proteins in the
FSSP format into an MSF file.
The MSF file is saved as a file with the same name as the fssp_file,
with an ".msf" extension appended at the end of the name of the file.
- fssp_file is a file describing a structural alignment between
  several molecules in FSSP format.  (This is the format used by the
  [DALI server](https://web.archive.org/web/20090301064750/http://ekhidna.biocenter.helsinki.fi/dali_server/start).)
  ***WARNING: This file format has been depreciated
  and the link points to an archive.***
- The MSF file format is explained
  [here.](http://rothlab.ucdavis.edu/genhelp/chapter_2_using_sequences.html#_Specifying_RSF_Files)
 *(MSF files can be converted to and from other
  more popular alignment file formats (such as FASTA) using
  [aligncopy.](http://emboss.sourceforge.net/apps/cvs/emboss/apps/aligncopy.html))*
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
