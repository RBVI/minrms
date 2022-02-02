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

#ifndef _MINRMS_PARSER_H
#define _MINRMS_PARSER_H

// This file contains declarations for the "MinrmsParser" class.
// This class contains all the settings to run Minrms, and
// also provides a function (the constructor) for parsing
// these settings from a list of command-line arguments.
//
//      Minrms is divided into several classes:
// an orientation generator, an alignment generator
// (which is a dynamic programming algorithm),
// and an orientation searcher ("SearchOr").
// (-The orientation searcher is the glue between
//   the first two classes.)
//
//      Similarly, the code that handles storing and reading the
// settings necessary for these classes to operate, is organized
// into separate classes as well.  Because some of the data
// necessary for an OrientationGenerator class, for example,
// is also the same data required for a "SearchOr" class,
// I organized this data in a hierarchy.
//      The "MinrmsParser" class multiply-inherits from
// the all of the parser-classes which handle the settings
// for all of the complenents that make up minrms.
// 
//
//           MinrmsParser class Hierarchy:
//
//            PairAlignSettingsParser
//              /          `-._ 
//             /               `-._
//            /                    `-._
//           /                 SearchOrSettingsParser
//          /                 /                \
// OrGenSettingsParser   SearchOrNW      SearchOrDyn3dSettingsParser
//          \           SettingsParser          /
//           \           |                  _.-'
//            \          |              _.-'
//             \         |          _.-'
//              \        |      _.-'
//               \       |  _.-'
//               MinrmsParser
//
//
// Note: In addition to inheriting from the classes you see
//       here, (each ending in the word "Parser"),
//       each "Parser" class also inherits from the analagous "Settings"
//       class in the MinrmsSettings hierarchy.
//       (See "minrms_settings.h".)

#include "minrms_settings.h"

namespace minrms
{


class PairAlignSettingsParser: public virtual PairAlignSettings
{
  PairAlignSettingsParser();
  PairAlignSettingsParser(const PairAlignSettingsParser &);
  PairAlignSettingsParser& operator =(const PairAlignSettingsParser &);

  void LoadStructures(vector<string>& vArgs);

  Biopolymer   aMol_NC[2];   //The original molecules, exactly from PDB.
  Biopolymer   aMol_f_NC[2]; //These only contain the desired residues.
                             //
                             //(Note to self:  Implementation.
                             // The "NC" stands for "not copied", I think:
                             // All instances of aMol[] are pointers
                             // to the actual stored in aMol_NC.
                             // When the PairAlignSettings variables
                             // are created destroyed, the 
                             // aMol[] arrays are not created or destroyed,
                             // However, when PairAlignSettingsParser
                             // variables are created, the 
                             // aMol_NC[] arrays are created, and the
                             // ::PairAlignSettings.aMol[] array
                             // is set to point to the aMol_NC[] array.
                             // -That's the only difference between
                             // aMol[] and aMol_NC[].
                             // Similar comments apply to
                             // aMol_f_NC[] vs. aMol_f[].

public:
  PairAlignSettingsParser(vector<string>& vArgs); //Parse arguments vArgs
};



class OrGenSettingsParser:public virtual OrientationGenerator::Settings,
                          public virtual PairAlignSettingsParser
{
  OrGenSettingsParser();
  OrGenSettingsParser(const OrGenSettingsParser &);
  OrGenSettingsParser& operator =(const OrGenSettingsParser &);
public:
  OrGenSettingsParser(vector<string>& vArgs); //Parse arguments vArgs
};



class SearchOrSettingsParser:public virtual SearchOr::Settings,
                             public virtual PairAlignSettingsParser
{
  SearchOrSettingsParser();
  SearchOrSettingsParser(const SearchOrSettingsParser &);
  SearchOrSettingsParser& operator =(const SearchOrSettingsParser &);
public:
  SearchOrSettingsParser(vector<string>& vArgs); //Parse arguments vArgs
};



class SearchOrDyn3dParser:public virtual SearchOrDyn3d::Settings,
                          public virtual SearchOrSettingsParser
{
  SearchOrDyn3dParser();
  SearchOrDyn3dParser(const SearchOrDyn3dParser &);
  SearchOrDyn3dParser& operator =(const SearchOrDyn3dParser &);
public:
  SearchOrDyn3dParser(vector<string>& vArgs); //Parse arguments vArgs
};



class SearchOrNWParser:public virtual SearchOrNW::Settings,
                       public virtual SearchOrSettingsParser
{
  SearchOrNWParser();
  SearchOrNWParser(const SearchOrNWParser &);
  SearchOrNWParser& operator =(const SearchOrNWParser &);
public:
  SearchOrNWParser(vector<string>& vArgs); //Parse arguments vArgs
};



class MinrmsParser:public MinrmsSettings,
                   public OrGenSettingsParser,
                   public SearchOrDyn3dParser,
                   public SearchOrNWParser
{
  //dissable copy constructor
  MinrmsParser(const MinrmsParser &);
  //dissable explicit copy
  MinrmsParser& operator =(const MinrmsParser &);
  MinrmsParser();

public:
  MinrmsParser(vector<string>& vArgs);  //Parse arguments vArgs
}; //class MinrmsParser




} //namespace minrms

#endif //#indef _MINRMS_PARSER_H



