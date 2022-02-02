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

#include <global_utils.h>
#include <pair_alignment.h>

using namespace minrms;

int
main(int argc, char const **argv)
{
  char help_msg[] =
    "Syntax error.\n"
    "\n"
    "usage:\n"
    "  msf2sequence sequence_identifier\n"
    "\n"
    "msf2sequence extract an individual sequence from an MSF file\n"
    "(read from the standard in), and sends it to the standard out.\n"
    "\"sequence_identifier\" indicates the label next to the sequence\n"
    "in the MSF-file you want to extract.\n";
  //"Gaps in the alignment must\n"
  //"be represented by '.' characters in the MSF file.  (See MinRMS\n"
  //"documentation for more details on this file format.)\n";

  if (argc != 2)
    ERR(help_msg);

  string sequence;
  bool label_found = ReadMSFsequence(cin,
                                     argv[1],
                                     sequence);
  if (! label_found)
    ERR("Error: Missing sequence in MSF-file:\n"
        "\"" << argv[1] << "\"");

  cout << sequence << endl;
}
