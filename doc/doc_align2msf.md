align2msf
=========


## Usage Syntax

```
align2msf.py align_server_output pdb_file_1 pdb_file_2 > output.msf
```

**align2msf.py** takes an alignment generated
by the ALIGN server, and converts it to an
[MSF file](http://rothlab.ucdavis.edu/genhelp/chapter_2_using_sequences.html#_Specifying_RSF_Files)
which was named "output.msf" in the example above.
*(MSF files can be converted to and from other
 more popular alignment file formats (such as FASTA) using
[aligncopy.](http://emboss.sourceforge.net/apps/cvs/emboss/apps/aligncopy.html))*

The first argument is a file containing the output from the ALIGN server.
The last two arguments contain PDB files for the structures
compared by the server.  Note: These two PDB files must be
listed in the command line in the same order they were sent
to the ALIGN server.

**align2msf.py** is a python script and it requires python to be installed.
(If the python interpreter is given a different name, for example "python3",
you can either provide a link from "python3"->"python", or insert the word
"python3" at the beginning of the command in the example above.
Eg "python3 align2msf.py ...".)


## Background:

The [pairwise ALIGN server](http://molmovdb.mbb.yale.edu/align)
returns structural alignments between two PDB files according to the
algorithm described in:

"Comprehensive Assessment of Automatic Structural Alignment Against a Manual Standard, the SCOP Classification of Proteins", Gerstein M, Levitt M, Prot. Sci., 7(2):445-56
[https://doi.org/10.1002/pro.5560070226](https://doi.org/10.1002/pro.5560070226)

The server output is somewhat confusing because it contains two
different alignments.  The alignment at the top is stored in familliar
MSF format, however according to Mark Gerstein, you don't want to
use this alignment because it contains some very 'tenuous matches'.

According to Mark Gerstein, the alignment to use is the one
at the bottom, under the "EQUIV" section, which, unfortunately
is not in MSF-format.  **align2msf.py** converts the alignment contained
in the "EQUIV" section into an MSF-file, and sends this MSF-file
to the standard out.
