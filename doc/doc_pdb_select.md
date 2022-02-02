pdb_select
==========

## Syntax
```
   pdb_select set1,set2,set3,... < orig_pdb_file > new_pdb_file
```

## Overview
This command creates a new_pdb_file containing a subset of
residues from the orig_pdb_file.  The only argument to
pdb_select is a comma-separated list (no spaces) of sets
of residues in Midas/Chimera/MinRMS format. (The selection syntax
is described below, as well as in the Midas/Chimera documentation).
"new_pdb_file" will contain residues belonging to the union
of all the sets (set1, set2, set3, etc.).


## Examples of selection syntax:

```
   set           residues selected:
 --------       ---------------------
 "100"         all residues whose seqNums are 100
 "100-150.A"   all residues in chain A whose seqNums are in [100,150]
 "*-150.A"     all the residues in chain A whose seqNums are up to 150
 "150-*.A"     all the residues in chain A whose seqNums are at least 150
 "*.A"         all the residues in chain A
 ".A"           "   "     "     "    "   "
 "100-150.A-C" residues 100-150 in chains A through C
 "100-150.*"   residues 100-150 in all chains
 "100-150."       "      "   "    "    "   "
 "100-150"        "      "   "    "    "   "
 "*.*"         the entire molecule
 "."            "    "      "
```

### Notes

If you use the '*' character in any of your sets,
you will have to enclose the first argument in quotes to
circumvent the shell.


### Details

- Only the ATOM, HETATM, ANISOU, HELIX, SHEET, and TURN records
  are effected.  All other records in the PDB file are blindly
  sent to the standard output.
- Helices, sheets, or turns which lie all, or partially
  outside the selected sets of residues will be deleted.
