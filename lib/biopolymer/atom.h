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


#ifndef _ATOM_H
#define _ATOM_H


#include <global_utils.h>
#include <short_string.h>
#include <vect_3d.h>
using namespace vect_3d;



namespace minrms {

typedef int AtomID;

struct Atom
{
    long       priority;

    //Explanation of "priority":
    //
    //If multiple atoms share the same AtomIDs,
    //the "priority" integer was intended to be used to
    //distinguish them from eachother.
    //
    //Usage:
    //
    //The PDB-file loading function >>() in  "load_pdb.h" uses the "priority"
    //field to store the line number of record in the PDB file that refers
    //to this atom.
    //
    //Behaviour:
    //
    //"priority" is used to distinguish 'alternate' atoms from eachother.
    //When calling Biopolymer::Segment::GetQuickAtom()
    //(see "linear_molecule.h"), the atom that is returned
    //is the one with the lowest "priority" number.


    short      element; //Atomic number of corresponding element
    Real       xyz[3];  //position of nucleus
    short      charge;  //(Ionic?) charge
    short      mass;    //atomic weight
    Real       occupancy, tempFactor;

    typedef VeryShortString AtomName;
    static const int ATOM_NAME_LENGTH = VERY_SHORT_STRING_LENGTH;
    static const AtomName DEFAULT_ATOM_NAME; //set to "?" as of 4/8/1999

    AtomName   name;    //A string identifier (Eg: "CA", "CB", "N", "O", "CG1")
                        //If this atom is part of a "Segment" of a
                        //"LinearMollecule", then this name should be
                        //unique within that segment.
    short      altLoc;  //Is this an "alternative" atom location?
                        //If so, this member will store an integer which
                        //distinguishes the alternates from eachother.
                        //If not, it will store "NO_ALTERNATES"
    static const short NO_ALTERNATES = -1;

    typedef int AtomID;
    AtomID      id;     //Additionally, one can store an identifier for each
                        //atom, seperate of its name.  (This is optional.)
    static const AtomID DONT_USE_ATOM_IDS = -1; //An impossible default value.


    Atom(int    set_element = 0,  //Atomic number of corresponding element
         Real   set_x = 0.0f,
         Real   set_y = 0.0f,
         Real   set_z = 0.0f,
         AtomID set_id = DONT_USE_ATOM_IDS,
         char const *set_name=DEFAULT_ATOM_NAME,
         short  set_altLoc = NO_ALTERNATES,
         int    set_charge = 0,
         int    set_mass = 0,     //atomic weight
         Real   set_occupancy = 0.0f,
         Real   set_temp_factor = 0.0f,
         long   set_priority = 0);

}; // class Atom

} //namespace minrms

#endif // #ifndef _ATOM_H
