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

#include <cstring> //needed for "strncpy"
using namespace std;

#include <global_utils.h>
#include <stl_utils.h>
#include "residue_names.h"
#include "biopolymer.h"

namespace minrms
{

const int Biopolymer::Residue::RESIDUE_NAME_LENGTH;
const Biopolymer::Residue::ResidueName
                      Biopolymer::Residue::DEFAULT_RESIDUE_NAME="?";
map<string, int> Biopolymer::Residue::backboneAtomLookup;



//assignment operator
Biopolymer::Residue& Biopolymer::Residue::
operator = (Residue const& source)
{
  id = source.id;
  assert(source.name);
  strncpy(name, source.name, RESIDUE_NAME_LENGTH);
  MmAtoms = source.MmAtoms;
  assert((aBackboneAtomPtrs && source.aBackboneAtomPtrs)
         ||
         ((! aBackboneAtomPtrs)
          && (! source.aBackboneAtomPtrs)
          && (NumBackboneAtoms() == 0))
         );
  if (source.consistentBackboneAtoms)
    //    InitBackboneAtoms(source.aBackboneAtomPtrs, source.end());
    InitBackboneAtoms();
  return *this;
}



//copy constructor
Biopolymer::Residue::
Residue(Residue const& source):MmAtoms()
{
  consistentBackboneAtoms = false;
  aBackboneAtomPtrs = new iterator [NumBackboneAtoms()];
  if (! aBackboneAtomPtrs)
    ERR("ERROR: \"" << __FILE__ << "\":" << __LINE__ << "\n"
                  "       Failure to alloc mem.\n"
                  "       Please report this error to the developer.\n");
  *this = source;
}



//default constructor
Biopolymer::Residue::
Residue():MmAtoms()
{
  consistentBackboneAtoms = false;
  aBackboneAtomPtrs = new iterator [NumBackboneAtoms()];
  if (! aBackboneAtomPtrs)
    ERR("ERROR: \"" << __FILE__ << "\":" << __LINE__ << "\n"
                  "       Failure to alloc mem.\n"
                  "       Please report this error to the developer.\n");
  strncpy(name, DEFAULT_RESIDUE_NAME, RESIDUE_NAME_LENGTH);
  id = 0;
}



//destructor
Biopolymer::Residue::
~Residue()
{
  if (aBackboneAtomPtrs)
    delete [] aBackboneAtomPtrs;
}



void Biopolymer::Residue::
InitBackboneAtoms()
               //iterator *aAssignWhichBackboneAtoms,
               //const_iterator dontAssignIfEqualsThis)
{
  assert(aBackboneAtomPtrs || (NumBackboneAtoms() == 0));
  for(int q = 0; q < NumBackboneAtoms(); ++q)
    aBackboneAtomPtrs[q] = end();

  //now loop over atoms in the residue and assign the pointers to them
  for(multimap<string,Atom>::iterator pa = MmAtoms.begin();
      pa != end();
      ++pa)
  {
    //Okay, now translate this atom.name string into an integer
    int backboneAtomType = LookupBackboneAtomID((*pa).second.name);

    //Update the aBackboneAtomPtrs[] array (if necessary)
    if (backboneAtomType != NOT_A_BACKBONE_ATOM) 
    {
      //      //first, check to see if this kind of atom is "filtered out"
      //      if ((! aAssignWhichBackboneAtoms) ||
      //          (aAssignWhichBackboneAtoms[backboneAtomType] != dontAssignIfEqualsThis))
      //      {
        //Okay, now check to see which atom has the lowest "altLoc" field
        iterator existing_atom = aBackboneAtomPtrs[backboneAtomType];
        if (existing_atom == end()) //if there's not one allready
          aBackboneAtomPtrs[backboneAtomType] = pa;
        else {
          //(check some stuff first)
          //          DEBUG_MSG(DBG_BIOPOLYMER,
          //                    "\"" << (*pa).second.name << "\",");
          //          DEBUG_MSG(DBG_BIOPOLYMER,
          //                "\'" << (*existing_atom).second.altLoc << "\'");
          assert(strcmp((*pa).second.name,
                        (*existing_atom).second.name) == 0);
          assert((*pa).second.altLoc != 
                 (*existing_atom).second.altLoc);
          //Now, if there's one already, then point to whichever
          //one has the lowest "priority" number.
          if ((*pa).second.priority < (*existing_atom).second.priority)
            aBackboneAtomPtrs[backboneAtomType] = pa;
        } //If there's not one allready, figure out if it's the smallest altLoc
        //      } //If it's not filtered out by aAssignWhichBackboneAtoms[], assign it.
    } // if (backboneAtomType != NOT_A_BACKBONE_ATOM) 
  } // for(iterator pa = begin();..
  consistentBackboneAtoms = true;
} // void Biopolymer::Residue::InitBackboneAtoms()



const int Biopolymer::Residue::NOT_A_BACKBONE_ATOM = -1;


int Biopolymer::Residue::
StoreBackboneAtomName(char const *backbone_atom_name)
{
  assert(backbone_atom_name);
  string temp_str(backbone_atom_name);
  if (temp_str.size() > Atom::ATOM_NAME_LENGTH)
    ERR("Atom's name: \"" << temp_str << "\"\n"
                  "is too long. Atom names are not to exceed: \n"
                  << Atom::ATOM_NAME_LENGTH - 1 << "characters in length.\n"
                  << "(Aborting)");

  int num_backbone_atoms = backboneAtomLookup.size();
  DEBUG_MSG(DBG_BIOPOLYMER,
            "Storing backbone_atom name: \""
            << backbone_atom_name << "\""
            "id = "
            << num_backbone_atoms);
  backboneAtomLookup[temp_str] = num_backbone_atoms;

#if 0
  //now update the table that returns the name, given the int.
  //(I want to use arrays not vectors here, and realloc() bad for C++)
  backboneAtomNames = Realloc_CppFriendly(backboneAtomNames,
                                       num_backbone_atoms,
                                       num_backbone_atoms+1);
  backboneAtomNames[num_backbone_atoms] = new char[strlen(backbone_atom_name)];
  if (! backboneAtomNames[num_backbone_atoms])
    ERR("ERROR: \"" << __FILE__ << "\":" << __LINE__ << "\n"
                  "       Unable to alloc mem().");
  strcpy(backboneAtomNames[num_backbone_atoms - 1], backbone_atom_name);
#endif //#if 0

  return num_backbone_atoms;
}


int Biopolymer::Residue::
LookupBackboneAtomID(Atom::AtomName const backbone_atom_name)
{
  assert(backbone_atom_name);
  string temp_str(backbone_atom_name);
  map<string,int>::iterator p = backboneAtomLookup.find(temp_str);
  if (p == backboneAtomLookup.end())
    return NOT_A_BACKBONE_ATOM;
  else
    return (*p).second;
}


//  The next function returns the name/symbol of the i'th backbone-atom.
char const *Biopolymer::Residue::LookupBackboneAtomSymbol(int i)
{
  assert((i >= 0) && (i < NumBackboneAtoms()));
  map<string,int>::iterator p;
  for(p = backboneAtomLookup.begin();
      ((p != backboneAtomLookup.end()) && ((*p).second != i));
      ++p)
    {}
  assert(p != backboneAtomLookup.end());
  return (*p).first.c_str(); //return the name
}




void Biopolymer::
insert(iterator  before_which_residue,
       PDBresID  residue_id,
       Residue::ResidueName const residue_name)
{
  Residue s;
  strncpy(s.name, residue_name, Residue::RESIDUE_NAME_LENGTH);
  s.id = residue_id;
  insert(before_which_residue,
         s);
}

void Biopolymer::
push_back(PDBresID  residue_id,
          Residue::ResidueName const residue_name)
{
  Residue s;
  strncpy(s.name, residue_name, Residue::RESIDUE_NAME_LENGTH);
  s.id = residue_id;
  push_back(s);
}


} //namespace minrms

