pdb2sequence
============

## Usage
```
pdb2sequence < input_pdb_file > output_sequence_file
```

## Description

**pdb2sequence** reads a PDB-file from the standard in,
converts the 3-letter residue codes in the PDB-file into 1-letter
codes and sends the resulting sequence to the standard out.


## Optional

The translation between 3 letter codes and one letter codes
can be customized by supplying a "res_code_dict.txt" file.
(See the MinRMS documentation for a description of how to do this.)
