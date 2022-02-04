## Nucleic Acids

It is *possible* to use **minrms** to align nucleic acid chains
*(although it's not clear if this is a useful capability).*

To do that, create a custom [atoms_used.txt](./share/README.md) file containing
the PDB codes for the atom type(s) you want to use for structure comparison.
(The lines in that file must contain the name of a PDB atom type
which is common to all RNA or DNA nucleic acids such as " C3'" or " C5'".)
Example "atoms_used.txt" file (containing only one line):
```
" C3'"
```
In that case, the RMSD will be calculated between the positions of C3' atoms
in either RNA chain.
*(Note:
A blank space must precede the atom names on each line, eg " C3", " C5'",...
The quotes (") were added above only to clarify the position of
whitespace in the file and should not appear in the "atoms_used.txt" file.)*

Once you have created a new "atoms_used.txt" file, run **minrms** using the
[-atoms-used atoms_used.txt](https://www.rbvi.ucsf.edu/Research/projects/minrms/docs/minrms.html#customize) and "-HS" arguments.
*(You must include the "-HS" argument to inform minrms that the molecules
being aligned do not contain any helix or sheet secondary structure.)*


## Warning

Minrms has not been carefully tested when applied to RNA molecules, so
please report problems to the issue tracker if these instructions do not work.
