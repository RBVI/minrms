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


#include <cstring> //needed for strncpy()
using namespace std;

#include "atom.h"

namespace minrms {

const Atom::AtomName Atom::DEFAULT_ATOM_NAME = "?";



Atom::Atom(int    set_element,  //Atomic number of corresponding element
           Real   set_x,
           Real   set_y,
           Real   set_z,
           AtomID set_id,
           char const *set_name,
           short  set_altLoc,
           int    set_charge,
           int    set_mass,
           Real   set_occupancy,
           Real   set_temp_factor,
           long   set_priority)
{
  id = set_id;
  altLoc = set_altLoc;
  if (set_name)
    strncpy(name, set_name, ATOM_NAME_LENGTH);
  else
    strncpy(name, DEFAULT_ATOM_NAME, ATOM_NAME_LENGTH);
  xyz[0]     = set_x;
  xyz[1]     = set_y;
  xyz[2]     = set_z;
  charge     = set_charge;
  mass       = set_mass;
  occupancy  = set_occupancy;
  tempFactor = set_temp_factor;
  priority   = set_priority;
}

} //namespace minrms
