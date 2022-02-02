This directory contains 7 different tests:
-----------------------------
1) when you run 

msf_compare -P 1.msf 2.msf A B
or just
msf_compare 1.msf 2.msf A B

you get: 2

2) when you run 
msf_compare -p 1.msf 2.msf A B
you get: 6

3) when you run 
msf_compare -R 1.msf 2.msf A B
you get: 3 4

4) when you run 
msf_compare -r 1.msf 2.msf A B
you get: 5 4


-----------------------------
The next three tests consists of two pdb files and two different
MSF-files beteen them.

Tests:
5) When you run 

msf_compare -3d,resemblesZ.pdb,resemblesDiagonal.pdb resemblesZ.msf resemblesX.msf resemblesZ.pdb resemblesDiagonal.pdb

you should get: sqrt(8/3) which equals 1.6329932.
(a terse explanation why is provided at the end of this file.)


6) When you run 

msf_compare -3d_both,resemblesZ.pdb,resemblesDiagonal.pdb resemblesZ.msf resemblesX.msf resemblesZ.pdb resemblesDiagonal.pdb

you should get a score of 2.0
(the residue in the center, which is only aligned in one
 of the alignments, not both, is not included in the RMSD.
 Since it did not change position at all, this brings the RMSD up.)

7) When you run 

msf_compare -3d_either,resemblesZ.pdb,resemblesDiagonal.pdb resemblesZ.msf resemblesX.msf resemblesZ.pdb resemblesDiagonal.pdb

you should get a score of sqrt(8/3) = 1.6329932 again.
The middle residue is included in the RMSD this time,
since it is aligned in one of the 2 alignments.




-----------------------------
-----------------------------
-----------------------------

Two structures are listed in this directory:

-----------------------------
"resemblesZ.pdb" looks like:

      7-6-5
         /
        4
  ^    /
Y |   3-2-1
  |
 -|---->
      X
-----------------------------

-----------------------------
"resemblesDiagonal.pdb" looks like:

          3
         /
        2
  ^    /
Z |   1
  |
 -|---->
      Y
-----------------------------

The two alignments differ by superimposing the straigth line structure
against the Z in the following orientations: / and \
This results in combined superimposed figures that resemble
a Z and a Roman X, respectively.
The difference between the two orientations of the "resemblesDiagonal.pdb"
structure is a 90-degree shift counterclockwise, moving the outlying
residues a distance of 2.0 angstroms, and not moving the center
residue at all.
