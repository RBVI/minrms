#include <cassert> //defines assert()
#include <cstdlib> //necessary for strtol() and strtod()
#include <fstream>
#include <algorithm>
using namespace std;


#include <global_utils.h>
//#include "../lib/global_utils/global_utils.h"
#include <stl_utils.h>
//#include "../lib/global_utils/stl_utils.h"
#include "test_parse_utils_libs.h"


#ifndef DBG_ARG_FILE_PARSER
//Necessary for some of the DEBUG_MSG() macros.  Just ignore this.
//(See "../global_utils/debug.h" for details.)
#define DBG_ARG_FILE_PARSER 0
#endif

#define STRTOLD_UNSUPPORTED
#define STRTOLL_UNSUPPORTED


//Why not includes an RCS-version-string?
//static char rcs_id[] = "@(#)$Header: /remote/local/babar/jewett/src/lib/parse_utils/RCS/parse_utils.cc,v 1.3 2001/05/03 14:20:28 jewett Exp jewett $";


using namespace ParseUtils;



void ParseUtils::
NextInt(string::const_iterator& start,
	string::const_iterator  stop,
	long long& dest)
{

  if ((start != stop) && isspace(*start))
  {
    ERR("Error:  Text containing integer was preceded by unnexpected\n"
        "        whitespace.  (This could be a bug in the program.)");
  }

  //First, I have to copy this text into a c_string to get around strtol()'s
  //weird syntax: strings only grant "const char *" access to their contents,
  //but strtol() requires "char *" arguments.
  string temp_str(start, stop);
  char *ac_temp_str = new char [temp_str.size() + 1];
  strcpy(ac_temp_str, temp_str.c_str());
  char *pstart = ac_temp_str;
  char *pstop;
  long long result = strtol(pstart, &pstop, 10);

  //If there was at least some valid text was read, then go ahead
  //and use the value returned by strtol(), by setting dest = result.
  //(otherwise discard the value).

  if ((pstop == pstart+1) && (*start == '-')) //first, we safeguard against "-"
    pstop = pstart; //because strtol() incorrectly skips past the '-' in "-"

  if (pstop != pstart) //If there was some valid text, use the value.
  {
    dest = result;
    //...and skip past the text that was read.
    start = start + (pstop-pstart);
  }
  delete [] ac_temp_str;
} //NextInt()



void ParseUtils::
NextInt(string::const_iterator& start,
	string::const_iterator  stop,
	long long& dest,
	string terminators)
{
  string token;
  NextToken(start, stop, token, terminators);

  if (token != "")
    //dest = atol(token.c_str()); This works, but it doesn't check the format.
  {
    string::const_iterator got_this_far = token.begin();
    NextInt(got_this_far,
	    token.end(),
	    dest);

    if (got_this_far != token.end())
      ERR("Error: \"" << token << "\" is not a valid integer.");
  }
} //NextInt()




//StrToReal() is a C++-string compatible version of the standard C-library's
//strtod().  Instead of pointers, StrToReal() takes two arguments
//which are string::const_iterators.
//These two iterators point to the beginning and ending characters demarking
//an interval within a larger string where a number is stored.
//This number contained in the string will be returned to the user.
//Afterwards "stop" will point to the character after the last
//character used in the conversion.

long double ParseUtils::
StrToReal(string::const_iterator  start,
          string::const_iterator& stop)
{

  //First, I have to copy this text into a c_string to get around strtod()'s
  //weird syntax: strings only grant "const char *" access to their contents,
  //but strtod() requires "char *" arguments.
  string temp_str(start, stop);
  char *ac_temp_str = new char [temp_str.size() + 1];
  strcpy(ac_temp_str, temp_str.c_str());
  char *pstart = ac_temp_str;
  char *pstop;

  #ifdef STRTOLD_UNSUPPORTED
  long double result = strtod(pstart, &pstop);
  #else
  long double result = strtold(pstart, &pstop);//Useful but not standard ANSI C
  #endif
  //now pstop points past the last valid char

  stop = start + (pstop - pstart);
  delete [] ac_temp_str;

  return result;
} //StrToReal()




void ParseUtils::
NextReal(string::const_iterator& start,
         string::const_iterator  stop,
         long double& dest)
{
  if ((start != stop) && isspace(*start))
  {
    ERR("Error:  Text containing floating point number was preceded by\n"
        "        unnexpected whitespace.  (This could be a bug in the program.)");
  }
  string::const_iterator p = stop;
  long double result = StrToReal(start, p);

  //If there was at least some valid text was read, then go ahead
  //and use the value returned by StrToReal(), by setting dest = result.
  //(otherwise discard the value).

  if ((p == start+1) && (*start == '-')) //first, we safeguard against "-"
    p = start; //because StrToReal() incorrectly skips past the '-' in "-"
  else if ((p == start+1) && (*start == '-')) //safeguard against "."
    p = start; //because StrToReal() incorrectly skips past the '.' in "."
  else if ((p == start+2) && (*start == '-') && (*(start+1) == '.'))
    //safeguards against "-."
    p = start; //because StrToReal() incorrectly skips past "-."

  if (p != start) //If there was some valid text, use the value.
  {
    dest = result;
    //...and skip past the text that was read.
    start = p;
  }

} //NextReal()






void ParseUtils::
NextReal(string::const_iterator& first,
         string::const_iterator  last,
         long double& dest,
         string terminators)
{
  string token;
  NextToken(first, last, token, terminators);

  if (token != "")
    //dest = atof(token.c_str()); This works, but it doesn't check the format.
  {
    string::const_iterator got_this_far = token.begin();
    NextReal(got_this_far,
             token.end(),
             dest);

    if (got_this_far != token.end())
      ERR("Error: \"" << token << "\" is not a valid floating point number.");
  }
} //NextReal()






void ParseUtils::
ConvertArgvToVectorOfStrings(int argc,
			     char **argv,
			     vector<string>& dest)
{
  dest.resize(argc);
  assert(argv);
  for (int i=0; i < argc; ++i)
  {
    assert(argv[i]);
    dest[i].assign(argv[i]);
  }
}


void ParseUtils::
GetlineUntil(istream &in,
             string& dest,
             string terminator_chars)
{
  DeleteContents(dest);
  if (! (in.good() || in.eof()))
    ERR("Error reading file: File ends prematurely.");
  char c;
  while(in.get(c))
  {
    if (BelongsTo(c, terminator_chars.begin(), terminator_chars.end()))
    {
      in.putback(c);
      break;
    }
    else
      dest.push_back(c);
  }
}



#if 0

void ParseUtils::
HandleArgFiles(vector<string> const& source_args,
               vector<string> & dest_args,
               string config_file_flag,
               string comment_indicator,
               string merge_lines_indicator,
               string default_args_filename,
               bool whitespace_delimited)
{
  dest_args = source_args;

  // ***  Read in the arguments from the 
  // ***  default configuration file (if there is one).
  vector<string> default_args;
  if (default_args_filename != "")
  {
    ifstream default_args_file(default_args_filename.c_str(), ios::in);
    if (! default_args_file)
      cerr <<
    "Warning: Unable to open the default arguments file:\""
       << default_args_filename << "\"\n"
    "         No default arguments loaded." << endl;
  }

  //...and insert the arguments from the default configuration file right after
  //right after first argument in the list (the name of the executable).
  //(The first argument in the array of arguments is special.  By convention,
  // the first argument should contain the name of the program,
  // as with argv[0] in C.  We want to preserve this convention.)
  dest_args.insert(dest_args.begin() + 1,
                   default_args.begin(),
                   default_args.end());
  //If there was no config_file_flag specified, we are done.
  if (config_file_flag == "")
    return;

  for (vector<string>::iterator p = dest_args.begin();
       p != dest_args.end();
       ++p)
  {
    //If any arguments do match config_file_flag,
    if ((*p) == config_file_flag)
    {
      //insert arguments from the file indicated.
      int insert_location = p - dest_args.begin();
      ++p;
      if (p == dest_args.end())
        ERR("Error: the \"" << config_file_flag <<
            "\" flag must be followed by the name of\n"
            "      a configuration file.");

      ifstream arg_file((*p).c_str());
      if (! arg_file)
        ERR("Error: Cannot open file \"" << (*p) << "\" for reading.");

      //Read in the arguments from config_file into config_file_args
      vector<string> args_in_file;





      //erase the arguments storing the config_file_flag and the file's name.
      dest_args.erase(p - 1, p + 1);

      //and insert them at this point.
      dest_args.insert(dest_args.begin() + insert_location,
                       args_in_file.begin(),
                       args_in_file.end());

      //Then set p to point to the first argument we inserted.
      p = dest_args.begin() + insert_location;
    }
  } //for (p == dest_args.begin(); p != dest_args.end(); ++p)
} //HandleArgFiles()

#endif //#if 0

