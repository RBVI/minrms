//
// Copyright (c) 2002 The Regents of the University of California.
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//   1. Redistributions of source code must retain the above copyright
//      notice, this list of conditions, and the following disclaimer.
//   2. Redistributions in binary form must reproduce the above
//      copyright notice, this list of conditions, and the following
//      disclaimer in the documentation and/or other materials provided
//      with the distribution.
//   3. Redistributions must acknowledge that this software was
//      originally developed by the UCSF Computer Graphics Laboratory
//      under support by the NIH National Center for Research Resources,
//      grant P41-RR01081.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE REGENTS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
// OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
// BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
// OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
// EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


#include <cassert> //defines assert()
#include <cstring> //defines strcmp()
#include <fstream>
using namespace std;

#include <global_utils.h>
#include <pdb++.h>
#include "load_pdb.h"





namespace minrms
{



istream & operator >> (istream&     pdb_file,
                       Biopolymer&  molecule)
{
  map<PDBresID, Biopolymer::size_type>
    residues_so_far;

  int record_counter = 1;
  int num_models = 0;
  bool alternate_atoms_present = false;

  PDB next_record;
  DEBUG_MSG(DBG_FILE_PARSER, __FILE__ << ":" << __LINE__
                            << " LoadedSequenceInfo::LoadSequence()"
                              << "\n  began parsing through the pdb file.");

  while (pdb_file >> next_record)
  {
    char c = pdb_file.get(); // These two lines were necessary at one point,
    pdb_file.putback(c);     // but I don't remember if they still are.

    switch (next_record.type())
    {
    case PDB::ATOM:
    {
      DEBUG_MSG(DBG_FILE_PARSER,
                "reading " << record_counter << "th record: ATOM");
      if (next_record.atom.altLoc != ' ')      
        alternate_atoms_present = true;

      //Find the residue with the same id as this atom's residue-id.
      PDBresID key;
      key.chainId    = next_record.atom.residue.chainId;
      key.seqNum     = next_record.atom.residue.seqNum;
      key.insertCode = next_record.atom.residue.insertCode;

      map<PDBresID, Biopolymer::size_type>
        ::iterator pResFound = residues_so_far.find(key);
      if (pResFound == residues_so_far.end())
      {
        // Then the residue was not found.
        // We have a new residue.  Insert it into the protein.
        molecule.push_back(key,
                           next_record.atom.residue.name);
        residues_so_far[key] = molecule.size() - 1; //index points to the
                                                     //residue we just added.
         assert(residues_so_far.size() == molecule.size());
        assert(molecule.size() > 0);
      }
      else
      {
        DEBUG_MSG(DBG_FILE_PARSER, 
                  "     matches residue at position ["
                  << (pResFound->second)
                  << "].  Updating residue.");
        assert(molecule.size() >= 1);

        if (pResFound->second != (molecule.size() - 1))
        {
          ERR("Error: Duplicate residue id (" << pResFound->first << "), found\n"
              "       in two disjoint places in the same PDB-file.");
        }
      }
      Biopolymer::iterator pRes = molecule.end() - 1; //points to the
                                                          //last residue.

      assert(next_record.atom.name);

      //The "altLoc" field needs special attention since we 
      //must check for the default case, and change the type
      //from char to short.
      short altLoc;
      if (next_record.atom.altLoc == ' ')
        altLoc = Atom::NO_ALTERNATES;
      else
        altLoc = static_cast<short>(next_record.atom.altLoc);

      //Now that we know which residue we are referring to, inform it
      //that a new atom has arrived.
      pRes->AddNewAtom(
                       Atom(
                            //LookupElement(next_record.atom.element),
                            0, //If anyone cares what the atomic number is,
                               //then I'll change this code later

                            next_record.atom.xyz[0],
                            next_record.atom.xyz[1],
                            next_record.atom.xyz[2],
                            next_record.atom.serialNum,
                            next_record.atom.name,
                            altLoc,
                            0,//LookupCharge(next_record.atom.name),
                            0,//LookupAtomMass(next_record.atom.name),
                               //If anyone cares about charge or atomic mass,
                               //I'll change this code
                            next_record.atom.occupancy,
                            next_record.atom.tempFactor,
                            record_counter
                            )
                       );
      break; // case (next_record.type() == PDB::ATOM)      
    } //case PDB::ATOM
    case PDB::HETATM:
      DEBUG_MSG(DBG_FILE_PARSER,
                "reading " << record_counter << "th record: HETATM");
      break;
    case PDB::HELIX:
      DEBUG_MSG(DBG_FILE_PARSER,
                "reading " << record_counter << "th record: HELIX");
      break;
    case PDB::SHEET:
      DEBUG_MSG(DBG_FILE_PARSER,
                "reading " << record_counter << "th record: SHEET");
      break;
    case PDB::TURN:
      DEBUG_MSG(DBG_FILE_PARSER,
                "reading " << record_counter << "th record: TURN");
      break;
    case PDB::MODEL:
      DEBUG_MSG(DBG_FILE_PARSER,
                "reading " << record_counter << "th record: MODEL");
      ++num_models;
      if (num_models > 1)
        ERR("Error: One of your PDB files contains several \"MODEL\" records.\n"
                      "       This program is unable to process PDB-files containing\n"
                      "       multiple \"MODEL\" records, because in such\n"
                      "       cases, it is ambiguous which version of the\n"
                      "       structure the user wants to use.\n"
                      "       Please re-submit a PDB-file containing only\n"
                      "       the model you want to use.\n"
                      "       Aborting...");
      break;
    default:
      DEBUG_MSG(DBG_FILE_PARSER,
                "reading " << record_counter << "th record: other record type");
      break;
    } //switch (next_record.type())
    ++record_counter; // keep track of the current record# (line#) in the file.
  } // while (! pdb_file.eof())
  if (alternate_atoms_present)
    cerr << 
      "********************************************************************\n"
      "*  WARNING: PDB-file contains \'alternate\' atom records.          *\n"
      "*           The PDB file indicates several possible locations      *\n"
      "*           for an atom.                                           *\n"
      "*             In such cases, this program discards all but the     *\n"
      "*           first instance of any atom in the pdb-file.  So if     *\n"
      "*           there are multiple CA atoms in residue 100B of         *\n"
      "*           chain C, it uses the first one, and discards the rest. *\n"
      "*             This program has no other quarrels with your         *\n"
      "*           PDB-files and should work as long as you don't         *\n"
      "*           mind this behavior.                                    *\n"
      "********************************************************************"
             << endl;

  return pdb_file;
} // istream & operator >> ()


void DeleteResiduesWithMissingAtoms(Biopolymer& molecule,
                                    char const **aAtomsThatMustBePresent)
{
  // loop over all residues in the molecule

  for(Biopolymer::iterator pRes = molecule.begin();
        pRes != molecule.end();
      ++pRes)
  {
    int delete_this_residue = false;

    // loop over atom-names in the list of required atoms
    for(char const **pAtomName = aAtomsThatMustBePresent;
        ! ((*pAtomName == NULL) || (strcmp(*pAtomName,"") == 0));
        ++pAtomName)
    {
      Biopolymer::Residue::iterator pAtom
        = pRes->find(*pAtomName);


      if (pAtom == pRes->end()) //If an atom with a matching name was not found
      {
        DEBUG_MSG(DBG_BIOPOLYMER,
                  "---- Residue #"
                  << pRes - molecule.begin() + 1
                  << "(" << pRes->id << ")"
                  << " did NOT contain the atom: \""
                  << *pAtomName << "\" ----");

        cerr    
          <<"--- Warning: residue ("
          <<"#" << pRes->id.seqNum << ","
          <<"chainId:'" << pRes->id.chainId << "',"
          <<"insert:'" << pRes->id.insertCode << "'"
          <<") lacks a \"" << *pAtomName << "\" atom.\n"
          "--- (The presence of this atom is required in each residue.)\n"
          << flush;
        //proteins[i].erase(proteins[i].begin()+s);
        //--s;
        delete_this_residue = true;
      }
      else {
        DEBUG_MSG(DBG_BIOPOLYMER,
                  "---- Residue #"
                  << pRes - molecule.begin() + 1
                  << "(" << pRes->id << ")"
                  << " contains the atom: \""
                  << *pAtomName << "\" ----");
      }
    } // loop over atom-names in the list of required atoms

    if (delete_this_residue)
    {
      cerr << "--------------This residue will be ignored-------------\n"
        << flush;
      //decrement to avoid invalidation upon deletion and pass pRes+1
      --pRes;
      molecule.erase(pRes +1);
    }
  } // loop over all residues in the molecule

  molecule.Finalize(); //We must call this, since we modified the molecule.
} // void DeleteResiduesWithMissingAtoms()



void LoadBackboneAtomSymbols(char const *filename)
{ 
  assert(filename);
  ShortString atomName;
  ifstream input_file(filename, ios::in);

  if (! input_file)
  {
    cerr << "Using CA positions.\n"
         << " (\"" << filename << "\" file not found.)" << endl;

    Biopolymer::Residue::StoreBackboneAtomName(" CA");
    return;
  }

  int row_counter = 1;
  char c;
  while (c = input_file.get(),
         input_file.putback(c),
         input_file)
  {
    input_file.get(atomName, SHORT_STRING_LENGTH, '\n');
    c = input_file.get(); //skip past the '\n'
    if ((! input_file.good()) && (! input_file.eof()))
      ERR("Bad atom name: \""  << atomName << "\"\n"
                    " File has invalid format or ends prematurely\n"
                    "on line " << row_counter << ".\n"
                    "This file should consist of separate lines containing\n"
                    "strings storing the exact PDB-names of each atom from\n"
                    "every residue to be compared (SPACES are not ignored!)."
                    "Also: no line is to exceed: "
                    << Atom::ATOM_NAME_LENGTH - 1 << "characters in length.\n"
                    << "(Aborting)");
    if (strcmp(atomName,"") != 0) //(store non-blank lines only)
      Biopolymer::Residue::StoreBackboneAtomName(atomName);
    ++row_counter;
  }
} // void LoadBackboneAtomSymbols(char *filename)







#ifdef DEBUG
ostream & operator << (ostream& out_file,
                       Atom&    atom)
{
  out_file << "name: \"" << atom.name << "\"\n";

  if (atom.altLoc != Atom::NO_ALTERNATES)
    out_file << "altLoc: '" << (char)atom.altLoc << "'\n";
  
  out_file << "serial: " << atom.id << "\n"
           << "pos ("
           << atom.xyz[0] << ","
           << atom.xyz[1] << ","
           << atom.xyz[2] << ")\n";
           //<< "element: " << atom.element << "\n";
  return out_file;
} // ostream & operator << (ostream& out_file,





ostream & operator << (ostream&        out_file,
                       Biopolymer&  molecule)
{
  for(Biopolymer::iterator pSeg = molecule.begin();
      pSeg != molecule.end();
      ++pSeg)
  {
    out_file << "=====================================\n"
             << (pSeg - molecule.begin()) + 1
             << "th residue, id(" 
             << pSeg->id
             << "):"
             << endl;
    out_file << "-------------------------------------" << endl;
    for(Biopolymer::Residue::iterator pAtom = pSeg->begin();
        pAtom != pSeg->end();
        ++pAtom)
    {
      out_file << (*pAtom).second;
      if ((*pAtom).second.altLoc == Atom::NO_ALTERNATES)
        out_file << "-------------------------------------" << endl;
      else
        out_file << " . . . ALTERNATE . . ." << endl;
    } //loop over all atoms in Molecule
    out_file << " - BackboneAtoms: -\n";
    for(int i=0;
        i < Biopolymer::Residue::NumBackboneAtoms();
        ++i)
    {
      out_file << i+1 << "th backboneAtom: ";
      Biopolymer::Residue::iterator patom = pSeg->GetBackboneAtom(i);
      if (patom == pSeg->end())
        out_file << "NOT_PRESENT" << endl;
      else
      {
        out_file << "#" << (*patom).second.id;
        out_file << ", \"" << (*patom).second.name << "\"";
        if ((*patom).second.altLoc != Atom::NO_ALTERNATES)
          out_file << ", '"
                   << static_cast<char>((*patom).second.altLoc)
                   << "'";
        out_file << endl;
      }
    }
  } //loop over all residues in Molecule

  return out_file;
} // ostream & operator << ()


#endif //#ifdef DEBUG





#if 0
//The following function returns true if the molecule
//has a residue containing multiple copies of an atom with
//the same "name" _and_ "altLoc" fields.
//It will _not_ return true if there are multiple copies
//of the same atom, but each copy has a different "altLoc" field.
//   (This function was written to safeguard against PDB-files
//with multiple _models_ or nonsensical data.  If a pdb-file has
//an atom with multiple _alternates_, but only one model,
//it will not detect it.)
bool HasDuplicateAtoms(Biopolymer  molecule);
#endif // #if 0






} //namespace minrms
