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

#include <vector>
#include <cassert>
using namespace std;

#include <sys/time.h> //required for "struct timeval" and "gettimeofday()"
#include <global_utils.h>
#include <parse_utils.h>
#include <or_gen.h>
#include "search_or.h"
#include "search_or_dyn3d.h"
#include "search_or_nw.h"
#include "minrms_parser.h"

//needed for random-number generation
#include <stdlib.h>//needed for call to getenv() and seed48()
#include <time.h>  //needed for call to time()

//Why not include an RCS-version-string?
//static char rcs_id[] = "@(#)$Header: /usr/local/src/minrms/bin/minrms/RCS/main.cc,v 3.5 2005/02/24 00:49:03 conrad Exp $";

const char g_version_string[] = "1.0.1";
const char g_date_string[]    = "<2019/10/13>";


using namespace minrms;

void InitSys();

int
main(int argc, char **argv)
{
  cout << "minrms v" << g_version_string << " " << g_date_string << "\n"
       << flush;

  #ifdef DEBUG
  cerr <<
    "   --                                  --  \n"
    " ---                                    ---\n"
    " --- DEBUG version (slower performance) ---\n"
    " ---                                    ---\n"
    "   --                                  --  \n\n" << flush;
  #endif

  InitSys(); //Any operating-system-ish initialization goes on here.

  //Convert argc and argv into a form that is easier to work with.
  vector<string> vOrigArgs;

  ParseUtils::
    ConvertArgvToVectorOfStrings(argc, argv, vOrigArgs);

  //If necessary, load up any configuration files that
  //were specified in the argument list.
  //  All we are doing at this stage is replacing any configuration-file-names
  //(following a "--cfg") found in the argument list, vArgs,
  //with the whitespace-delimited tokens contained in these files.
  //(The tokens contained in these files will not
  //actually be parsed until later.)

  char back_slash[2] = {'\\', '\0'};
  const char *default_settings_file_name_c = getenv("MINRMS_SETTINGS_FILE");
  string default_settings_file_name;
  if (default_settings_file_name_c)
    default_settings_file_name.assign(default_settings_file_name_c);

  vector<string> vArgs;
  ParseUtils::
    HandleArgFiles(vOrigArgs,
                   vArgs,
                   "-cfg",
                   "#",
                   back_slash,
                   default_settings_file_name);

  //Now, we are ready to parse the argument list.
  MinrmsParser settings(vArgs);//Parses the command line, retrieves all the
                               //run-time settings, and loads in the molecules.

  //Check to make sure all the arguments were processed.
  //(except for the first argument, which is the executable name)
  if (vArgs.size() != 1)
  {
    assert(vArgs.size() > 1);
    ERR("Error: Unrecognized command line option: \"" << vArgs[1] << "\"\n");
  }


  
  // ****************************************************************
  // ********    We must decide which orientations        ***********
  // ********  we will initially consider in the search.  ***********
  // ****************************************************************
  OrientationGenerator og(settings);//Figures out which orientations will
                                    //be considered.

  //Create the class that will be used to search the
  //space of orientations stored in og.  This class
  //will search for alignments with low RMSD.
  SearchOr *pSearchOr;

  switch (settings.which_algorithm)
  {
  case MinrmsSettings::ALGO_DYN3D:
    pSearchOr = new SearchOrDyn3d(settings);
    if (! pSearchOr)
    ERR("ERROR: Cannot allocate memory necessary for the\n"
        "       dynamic-programming table that is used to\n"
        "       generate the equivalences at a fixed orientation.\n"
        "          (This table grows in size as of O(m^2 n), where\n"
        "       m and n are the lengths of the two sequences.)\n"
        "          Try reducing sequence size, or use various\n"
        "       command-line arguments to reduce the size of the\n"
        "       table by placing an upper-bound and/or lower-bound\n"
        "       on the number of allowed matches made.\n"
        "       Aborting...\n");
    break;
  case MinrmsSettings::ALGO_NEEDLEMAN_WUNSCH:
    pSearchOr = new SearchOrNW(settings);
    CHECK_ALLOC(pSearchOr);
    break;
  default:
    assert(0);
    break;
  }

  // *******************************************************************
  // ******* Search the set of orientation in og to find          ******
  // ******* the lowest RMSD alignments (and alternate alignments ******
  // ******* if requested).                                       ******
  // *******************************************************************
  pSearchOr->Search(og, 0, og.size());

  // ****************************************
  // *******   Refine if requested.   *******
  // ****************************************
  if (settings.refine_method != SearchOr::Settings::NO_REFINE)
    pSearchOr->Refine();

  // ************************************************
  // ******* Print out the results to the user. *****
  // ************************************************
  pSearchOr->WriteMSF(); //saves the MSF files of all the alignments.
  pSearchOr->WriteChimeraInfo(); //saves any other information Chimera needs.


} //main()



// DEC doesn't declare this if "ansi" is on
extern "C" void srand48(long);

void InitSys()
{
  #if 0
  Commented out.  I want minrms to be deterministic. Set to a constant instead.
  // Initialize the random number generator
  struct timeval  the_time;
  struct timezone dummy_tz;
  gettimeofday(&the_time, &dummy_tz);
  srand48( the_time.tv_usec );
  #endif //#if 0

  srand48( 42 ); //why not 42?
}
