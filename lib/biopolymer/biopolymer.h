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


#ifndef _BIOPOLYMER_H
#define _BIOPOLYMER_H


// *****************************************************************
// ****     The Biopolymer class is a data structure for storing
// **** homo or biopolymers.  The goal was:
// ****
// **** 1) To be able to read and write PDB files easily.
// ****    Internally, the data is stored very similar to a PDB-file.
// ****    (PDB files are messy, and the data structure is messy too.)
// ****
// **** 2)  To have fast random-access to find the i'th residue's data.
// ****    There is heirarchy:  Atoms are grouped into residues,
// ****    and residues are grouped into molecules (biopolymers).
// ****    Random access to find the data for the i'th residue requires O(1)
// ****    time, and has low overhead.  STL-like random-access iterators
// ****    can be used to access the residues in the molecule, and also the
// ****    atoms in each residue.  (Residues can also be accessed by their
// ****    PDB information: chainID, sequenceNumber, and insertCode, but
// ****    this requires O(log(n)) time, where n = the length of the polymer)
// ****
// **** 3) To be able to access the the "backbone atoms" quickly.
// ****    For biopolymers, the atoms contained in each residue
// ****    varies from residue to residue.  However there is usually
// ****    a subset of atoms that most if not all the residues have in common.
// ****    (For example, all protein residues have alpha-carbon and N atoms.)
// ****    I call these the "backbone" atoms (regardless of whether they
// ****    are part of the backbone), and rather you can access the i'th
// ****    residue's backbone atom in O(1) time with very low overhead.
// ****    (The user can ignore this feature if they don't need it.)
// ****
// ****    I realize this data structure is not very elegant, and I was a
// **** C++ novice when I wrote it, but it was good enough for my purposes.
// *****************************************************************


#include <cassert> //needed to define assert()
#include <vector>
#include <map> //not <multimap>, at least on the dec-alpha's
#include <string>
#include <algorithm>
using namespace std;


#include <stl_utils.h>
#include "pdb_id.h"
#include "atom.h"


#define DBG_BIOPOLYMER 0 //(This just turns off the printing of some
                            // messages used for debugging. Please ignore.)



namespace minrms {

class Biopolymer
{
public:

  // *****************************************************************
  // **** The Residue data structure:
  // *****************************************************************
  // **** A "residue" stores a list of atoms, as well as some other
  // **** information to make it easier to identify it within the polymer.
  // **** In addition to the atom coordinates, types, masses and names,
  // **** residues have their own names, and their own ID (seq) numbers.
  // ****** Details: ****************************************
  // ****  -Residues mimic the interface of a multimap indexed by atomName.
  // ****   For example, you can find the atom with name " CA" in O(log(n))
  // ****   time, where n=the number of atoms in the residue.
  // ****  -For faster access, the "backbone atoms" have special id numbers.
  // ****   You can lookup a backbone atom using it's id number in O(1) time.
  // ****  -Alternates (alternate atom locations) are
  // ****   stored along with the regular atoms.
  // ****   If an atom has no alternates, it's "altLoc" field
  // ****   stores Atom::NO_ALTERNATES.
  // ****   All alternate atom's are located in consecutive locations
  // ****   in the residue, (not interspersed randomly throughout)
  // ****  -Otherwise, the atoms are not stored in any
  // ****   particular order in the array.
  // *****************************************************************
  // **** Accessing a residue's "Backbone Atoms"                  ****
  // ****                                                         ****
  // **** StoreBackboneAtomName() 
  // **** corresponds strings with ints.
  // **** In the future the int (eg. _CA)
  // **** will be used instead of the string (eg. " CA")
  // **** to speedup the lookup process.
  // ****
  // ****                             int _CA = StoreBackboneAtomName(" CA");
  // ****                            
  // **** Later on this integer can
  // **** be passed to the function:  iterator GetBackboneAtom(_CA)
  // ****                                 |
  // **** which returns an                |
  // **** iterator to an atom             |
  // **** in the residue which            |
  // **** can be accessed using *         |
  // **** (retrieves an STL "pair".       |
  // ****  use its "second" field)        |
  // **** This returns,                  \|/
  // **** the atom you are             (an "Atom" struct)
  // **** seeking.
  // ****
  // ****    -If there are alternate locations for this atom,
  // ****      the record pointed to by FindBackboneAtom()
  // ****      will be the one with the lowest "altLoc" field.
  // ****    -If the residue in question does not contain
  // ****     an atom of the type requested, then
  // ****     GetBackboneAtom() will return end().
  // ******************************************************************
  // ****    Setting up the Name-Lookup for Backbone Atoms:
  // **** Here's an example of how to initialize everything correctly:
  // **** After loading in the atoms into the molecule, you
  // **** must call:
  // ****
  // **** int _N  = Biopolymer::Residue::StoreBackboneAtomName(" N");
  // **** int _CA = Biopolymer::Residue::StoreBackboneAtomName(" CA");
  // **** int _C  = Biopolymer::Residue::StoreBackboneAtomName(" C");
  // **** int _O  = Biopolymer::Residue::StoreBackboneAtomName(" O");
  // ****
  // **** my_molecule.finalize();
  // ****
  // **** (Note: "my_molecule" is of type "Biopolymer".)
  // *****************************************************************
  // **** Chains:
  // **** Also: Presently there is no requirement that residues within a
  // **** biopolymer have to have to belong to the same "chain"
  // **** in the PDB file.  ChainID numbers are used, along with sequenceID
  // **** numbers and insert-codes, to figure out the order of the residues
  // **** in the polymer.  Petides with lower chainID numbers come first.
  // ******************************************************************


  class Residue
  {
    multimap<string, Atom> MmAtoms;
    friend class Biopolymer;

  public:

    Residue(Residue const&);

    Residue();

    ~Residue();

    Residue& operator = (Residue const&);

    typedef VeryShortString ResidueName; //see <short_string.h> for definition
    static const int RESIDUE_NAME_LENGTH = VERY_SHORT_STRING_LENGTH;

    ResidueName name; // A name for this residue indicating its type (optional)
                      // Examples of residue names for residues are "ALA","TRP"

    static const ResidueName DEFAULT_RESIDUE_NAME;//(This name is used whenever
                                                  // the caller does not supply
                                                  // a name.)

    // (Modify this next line as necessary:)
    PDBresID   id; // A unique identifier for each residue (not optional)
                   // indicating the location in the sequence of residues.


    // *****************************************************************
    // **** Data Access:

    typedef multimap<string, Atom>::iterator iterator;
    typedef multimap<string, Atom>::const_iterator const_iterator;

    iterator begin()                  { return MmAtoms.begin(); }
    iterator end()                    { return MmAtoms.end(); }
    const_iterator begin() const      { return MmAtoms.begin(); }
    const_iterator end()   const      { return MmAtoms.end(); }
    long size() const                 { return MmAtoms.size(); }

    iterator AddNewAtom(Atom const& atom)
    {
      string name_str(atom.name);
      consistentBackboneAtoms = false;
      return MmAtoms.insert(pair<string, Atom>(name_str, atom));
    }

    // (I haven't decided if I should make AddNewAtom the only way to 
    // insert elements.  Right now, you can also manually use the
    // STL-wrapper "insert()" function defined below, but I may get
    // rid of it later, so I guess use it at your own peril.
    // Also, it's not quite consistent with the other code, because it
    // requires passing STL-strings instead of char*.)
    iterator insert(pair<string, Atom> name_atom_pair)
    {
      consistentBackboneAtoms = false;
      return MmAtoms.insert(name_atom_pair);
    }

    iterator find(string key){
      return MmAtoms.find(key);
    }

    const_iterator find(string key) const {
      return MmAtoms.find(key);
    }

    iterator lower_bound(string key){
      return MmAtoms.lower_bound(key);
    }

    const_iterator lower_bound(string key) const {
      return MmAtoms.lower_bound(key);
    }

    iterator upper_bound(string key){
      return MmAtoms.upper_bound(key);
    }

    const_iterator upper_bound(string key) const {
      return MmAtoms.upper_bound(key);
    }
 
    pair<iterator,iterator> equal_range(string key) {
      return MmAtoms.equal_range(key);
    }

    pair<const_iterator,const_iterator> equal_range(string key)
      const
    {
      return MmAtoms.equal_range(key);
    }

    //----Backbone Atoms----
    //"Backbone atoms" are atoms which are indexed by id numbers,
    //instead of using a string containing their name.
    //They can be accessed quickly in O(1) time.

  private:
    iterator *aBackboneAtomPtrs;       //A list of pointers (actually iterators)
                                    //to the backboneAtoms that belong to
                                    //this residue only.
    bool consistentBackboneAtoms; //Is the information in aBackboneAtomPtrs accurate?
    static map<string, int> backboneAtomLookup; //stores global string,int pairs
    static char ** aBackboneAtomNames; //stores the inverse lookup(int to string).

    //Now, finally, as mentioned before, in order for
    //Biopolymer::Residue::GetBackboneAtom()
    //to work it is necessary to first call
    //Biopolymer::Finalize() (see "linear_molecule.h")
    //What this function does is to call "Residue::InitBackboneAtoms()" on all
    //its residues.  This enables quick access to the atoms in this residue
    //whose names match one of the names in the BackboneAtoms array.
    //InitBackboneAtoms() should be called only after loading in all the atoms
    //in the entire residue, as well as after specifying which atoms
    //are "backbone" atoms using the "StoreBackboneAtomNames()" function.

    void InitBackboneAtoms();

  public:
    //This next function returns an iterator to the atom specified
    //by "backbone_atom_type".  "backbone_atom_type" should be an
    //integer associated with a particular kind of atom name.
    //This association should be made by setting it equal
    //to the return value of "StoreBackboneAtomName()".
    //If there are alternate atoms present, it returns the one with
    //the lowest "altLoc" field.
    //If this type of atom is not present in the residue, it returns end().

    iterator GetBackboneAtom(int backbone_atom_type) {
      return aBackboneAtomPtrs[backbone_atom_type];
    }

    const_iterator GetBackboneAtom(int backbone_atom_type) const {
      return aBackboneAtomPtrs[backbone_atom_type];
    }

    //The integers that can be passed to GetBackboneAtom() lie in the 
    //range [0, b), where b = NumBackboneAtoms().
    //This way, if there is a set of operations one needs to
    //perform on all of the backbone-atoms, one can just iterate
    //through the backboneAtoms.
    static int  NumBackboneAtoms() { return backboneAtomLookup.size(); }

    //The following function indicates which atoms should be considered
    //"backbone-atoms", according to the atom's name-field.
    //(For example strings like " CA ", " N  " are used in PDB files to
    // identify the alpha-carbon and nitrogen atoms within a residue.)
    //(Generates an error if the name is too long.)
    // It returns the integer associated with that name so it can
    // be more quickly referenced by number instead of by name
    // by calling GetBackboneAtom().
    static int StoreBackboneAtomName(char const *backbone_atom_name);
    //  (If your program "forgets" the number(s) returned by
    //   StoreBackboneAtomName(), you can retrieve them using "LookupBackboneAtomID")
    static int LookupBackboneAtomID(char const *backbone_atom_name);

    //  The next function returns the name/symbol of the i'th backbone-atom.
    //(Not fast. It's linear in the number of backboneAtoms.)
    static char const *LookupBackboneAtomSymbol(int i);

    static const int NOT_A_BACKBONE_ATOM; // = -1,  returned if search fails

  }; //Biopolymer::Residue
  // *****************************************************************
  // **** End of the Residue data structure. 
  // *****************************************************************
 
public:

  // *****************************************************************
  // **** Data Access:
  typedef  vector<Residue>::iterator iterator;
  typedef  vector<Residue>::const_iterator const_iterator;
  typedef  vector<Residue>::size_type size_type;

private:
  vector<Residue> vResidues;
  map<PDBresID, Biopolymer::size_type> mPdbResIdLookup;
  
public:
  Biopolymer():vResidues()  {}
  Biopolymer(iterator start,
                 iterator stop):vResidues(start, stop)  {}
  
  iterator begin()                       { return vResidues.begin(); }
  iterator end()                         { return vResidues.end(); }
  const_iterator begin() const           { return vResidues.begin(); }
  const_iterator end()   const           { return vResidues.end(); }
  Residue& operator[] (size_type i)           { return vResidues[i]; }
  Residue const& operator[](size_type i)const { return vResidues[i]; }
  size_type size() const                 { return vResidues.size(); }
  Residue& front()                       { return vResidues.front(); }
  Residue const& front() const           { return vResidues.front(); }
  Residue& back()                        { return vResidues.back(); }
  Residue const& back() const            { return vResidues.back(); }

  //Note: find(), comes_after() and comes_before(),
  //      only works if Finalize() has been called since
  //      after the last call to insert(), erase(), or push_back().
  iterator find(PDBresID id)
  {
    map<PDBresID, size_type>::iterator
      p = mPdbResIdLookup.find(id);
    if (p == mPdbResIdLookup.end())
      return end();
    else
      return begin() + (*p).second;
  }
  //{ return SortedFind(begin(), end(), id); }

  const_iterator find(PDBresID id) const
  {
    map<PDBresID, size_type>::const_iterator
      p = mPdbResIdLookup.find(id);
    if (p == mPdbResIdLookup.end())
      return end();
    else
      return begin() + (*p).second;
  }
  //{ return SortedFind(begin(), end(), id); }


  //Note: insert() and push_back() do not safeguard against adding two
  //      residues with the same id (chain, seqNum, and insertCode).
  //      To avoid entering duplicate residues, the user must keep
  //      track of which residues have been added so far themselves.
  void insert(iterator        before_which_residue,
              Residue const&  input_residue)
  {
    vResidues.insert(before_which_residue, input_residue);
  }

  void push_back(Residue const& input_residue)
  {
    vResidues.push_back(input_residue);
  }

  //For conveiniance, the next 2 functions identical to the ones above, but
  //you dont have to pass the whole residue, just the "name" and "id" fields.
  //(The other fields can be filled in later.)
  void insert(iterator                   before_which_residue,
              PDBresID         residue_id,
              Residue::ResidueName const residue_name
                                           = Residue::DEFAULT_RESIDUE_NAME);

  void push_back(PDBresID         residue_id,
                 Residue::ResidueName const residue_name
                                            = Residue::DEFAULT_RESIDUE_NAME);

  void erase(iterator which_residue)
  { vResidues.erase(which_residue); }


  //*** You must invoke Finalize() after you have made all modifications
  //*** to the molecule, and after you have specified which atoms
  //*** are "backbone" atoms using the "StoreBackboneAtomNames()" function.
  void Finalize()
  {
    mPdbResIdLookup.clear();

    for (size_type i = 0; i < size(); ++i)
    {
      map<PDBresID, Biopolymer::size_type>::iterator p;
      p = mPdbResIdLookup.find(vResidues[i].id);
      if (p != mPdbResIdLookup.end()) //if residue allready is there.
        ERR("Duplicate residue in molecule identified by:\n"
            "  (chain, seqNum, insertCode) = " << vResidues[i].id);
      mPdbResIdLookup[vResidues[i].id] = i;
      vResidues[i].InitBackboneAtoms();
    }
    TrimContainer(vResidues);//<-delete excess space at end of vResidues
  }

}; // class Biopolymer


//Ordering:

inline bool comes_before(PDBresID a, PDBresID b, Biopolymer const& m)
{
  Biopolymer::const_iterator ra = m.find(a);
  Biopolymer::const_iterator rb = m.find(b);
  assert(ra != m.end());
  assert(rb != m.end());
  return ra < rb;
}

inline bool comes_after(PDBresID a, PDBresID b, Biopolymer const& m)
{
  Biopolymer::const_iterator ra = m.find(a);
  Biopolymer::const_iterator rb = m.find(b);
  assert(ra != m.end());
  assert(rb != m.end());
  return ra > rb;
}

//Returns whether or not the structure contains residues
//that are not in the dictionary used by NameCodeLookup
//(see "residue_names.h")
bool UnrecongizedResidueNames(Biopolymer const& m);



} //namespace minrms

#endif // #ifndef _BIOPOLYMER_H

