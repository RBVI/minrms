msf2sequence
============

## Usage

```
msf2sequence sequence_identifier < alignment.msf > sequence.txt
```

## Explanation

**msf2sequence** extracts an individual sequence from an MSF file
(eg. "alignment.msf") and saves it as a text file (eg. "sequence.txt").
"sequence_identifier" indicates the label next to the sequence
in the MSF-file you want to extract.

The MSF file format is explained
[here.](http://rothlab.ucdavis.edu/genhelp/chapter_2_using_sequences.html#_Specifying_RSF_Files)
*(MSF files can be converted to and from other
more popular alignment file formats (such as FASTA) using
[aligncopy](http://emboss.sourceforge.net/apps/cvs/emboss/apps/aligncopy.html))*
