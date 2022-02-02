#ifndef _PARSE_UTILS_H
// ***************************************************************
// ****  This file contains a lot of functions for parsing    ****
// **** strings, specifically strings that were passed to     ****
// **** the command line (through argc/argv).                 ****
// **** In addition, there are some utilities for parsing     ****
// **** streams as well.                                      ****
// ****                                                       ****
// ****  The traditional C's argc and argv array is converted ****
// **** into a vector-of-strings.  There are also             ****
// **** functions for reading in arguments from files.        ****
// **** The names of the files themselves can appear in the   ****
// **** argument list.                                        ****
// ****    (Why this was implemented:                         ****
// ****     These arguments are typically arguments from the  ****
// ****    command line (argc/argv).  The ability to read in  ****
// ****    files was added to allow the user to write         ****
// ****    configuration-files as an alternative to passing   ****
// ****    long argument lists.)                              ****
// ***************************************************************
// **** Note:  I apologize in advance.                        ****
// ****   Some of these utilities may unecessarily duplicate  ****
// ****   the behavior of existing standard-library/STL       ****
// ****   functions I didn't know of at the time of           ****
// ****   when writing this.                                  ****
// ***************************************************************


#include <vector>
#include <string>
#include <algorithm>
#ifdef DEBUG
#include <iostream>
#endif
using namespace std;


//#include <simple_numeric_types.h>






namespace ParseUtils
{


// ****************************************
// **** Functions for parsing streams *****
// ****************************************

//Skips over any characters in the input stream that are present
//in the "chars_to_skip" string.
// (A version which works on strings is included below.)
void Skip(istream &input, string chars_to_skip);

//Skips over any characters with white-space from the input stream
// (A version which works on strings is included below.)
void SkipWhiteSpace(istream& in);

//Same as above but does not skip over return-cairrages.
// (A version which works on strings is included below.)
void SkipWhiteSpaceExceptNewlines(istream& input);

//skips everything _but_ white-space
// (A version which works on strings is included below.)
void SkipNonWhiteSpace(istream& input);

//Similar to GetLine() for strings, except it terminates whenever it
//encounters any of the characters in "terminator_chars".
void GetlineUntil(istream &input,
                  string& dest,
                  string terminator_chars);

// ************************************************
// **** Functions for parsing STL containers: *****
// ************************************************


// ** The next four functions are identical to the
// ** functions with the same name which operate
// ** on streams.
template <class Str, class StrConstIterator>
inline void
Skip(StrConstIterator& first,
     StrConstIterator last,
     Str chars_to_skip)
{
  first = find_first_of(first,
                        last,
                        chars_to_skip.begin(),
                        chars_to_skip.end(),
                        not_equal_to);
}


template <class Str, class StrConstIterator>
inline StrConstIterator
Find(StrConstIterator first,
     StrConstIterator last,
     Str chars_to_find)
{
  return find_first_of(first,
		       last,
		       chars_to_find.begin(),
		       chars_to_find.end());
}


// ***********************************************
// **** Functions for parsing strings-only:  *****
// ***********************************************


//StrToReal() is a C++-string compatible version of the standard C-library's
//strtod().  Instead of pointers, StrToReal() takes two arguments
//which are string::const_iterators.
//These two iterators point to the beginning and ending characters demarking
//an interval within a larger string where a number is stored.
//This number contained in the string will be returned to the user.
//Afterwards "stop" will point to the character after the last
//character used in the conversion.

long double
StrToReal(string::const_iterator  start,
          string::const_iterator& stop);





//NextToken() finds the location the first instances
//of one of the "terminators" in the range [first,last),
//or "last", whichever comes earlier.
//-NextToken() stores all of the characters up until
// the terminator (or "last") in the "dest" argument.
//-NextToken() sets the iterator "first" to point to
// the location right after this terminator
// (or "last", whichever comes earlier).
//-It's unlikely I'll have a use for the following argument:
// The fifth argument, "equals" is the function-object
// used to make the equality comparison.  If it is set to
// "not_equal_to", then instead of stopping after the first
// instance of terminators, it stops after the first character
// that is NOT a terminator.

inline void
NextToken(string::const_iterator& first,
          string::const_iterator  last,
          string&                 dest,
          string                  terminators)
{
  string::const_iterator p = find_first_of(first,
                                           last,
                                           terminators.begin(),
                                           terminators.end());
  dest.assign(first, p);
  if (p != last) ++p; //skip the terminator (if there was one)
  first = p; //advance the read head ("first") to point to the new location.
}


//The following functions retrieve numbers from text strings.
//They are very similar to the strtod() and strtol() functions
//from the C standard library.
//  Parsing begins at the character pointed to by start and continues
//until encountering a character that violates the format for the
//type of data being read, or stop is reached.
//If any meaningful text was read up until that point,
//that text will be converted into the appropriate
//data type (Real, int, or long).
//Otherwise, the value of dest is unmodified.
//  Afterwards, the iterator start is either set to point to the offending
//character that terminated the parsing (if there was one), or stop.
//
//  Note: Unlike strtod() and strtol(), whitespace preceeding the value
//        being read, is not skipped, but generates an error.

void
NextReal(string::const_iterator& start,
         string::const_iterator  stop,
         long double& dest);

inline void
NextReal(string::const_iterator& start,
         string::const_iterator  stop,
         double& dest)
{
  long double x = dest;
  NextReal(start, stop, x);
  dest = static_cast<double>(x);
}

inline void
NextReal(string::const_iterator& start,
         string::const_iterator  stop,
         float& dest)
{
  long double x = dest;
  NextReal(start, stop, x);
  dest = static_cast<float>(x);
}


void
NextInt(string::const_iterator& start,
        string::const_iterator stop,
        long long& dest);

inline void
NextInt(string::const_iterator& start,
        string::const_iterator stop,
        long& dest)
{
  long long i = dest;
  NextInt(start, stop, i);
  dest = static_cast<long>(i);
}

inline void
NextInt(string::const_iterator& start,
        string::const_iterator stop,
        int& dest)
{
  long long i = dest;
  NextInt(start, stop, i);
  dest = static_cast<int>(i);
}


//The following versions of NextReal(), NextInt, and NextLong():
//1) read the text in a string (limited to the range [start,stop))
//   up until one of the terminating characters is reached.
//2) The text up until that point (if there was any) will then
//   be converted to the appropriate data type (Real, int, or long).
//   If there was not any text before a terminating character
//   was reached, the value of dest is left alone.
//   An error message is generated if the text up until the terminator
//   text does not have the correct format for the value being read.
//3) Either way, the iterator "start" is set to point to
//   the location right after the character that terminated the parsing
//   (or to stop, whichever comes first.).

void
NextReal(string::const_iterator& first,
         string::const_iterator last,
	 long double& dest,
	 string terminators);

inline void
NextReal(string::const_iterator& first,
         string::const_iterator last,
	 double& dest,
	 string terminators)
{
  long double x;
  NextReal(first, last, x, terminators);
  dest = static_cast<double>(x);
}

inline void
NextReal(string::const_iterator& first,
         string::const_iterator last,
	 float& dest,
	 string terminators)
{
  long double x;
  NextReal(first, last, x, terminators);
  dest = static_cast<float>(x);
}



void
NextInt(string::const_iterator& first,
        string::const_iterator last,
        long long& dest,
        string terminators);

inline void
NextInt(string::const_iterator& first,
        string::const_iterator last,
        long& dest,
        string terminators)
{
  long long i;
  NextInt(first, last, i, terminators);
  dest = static_cast<long>(i);
}

inline void
NextInt(string::const_iterator& first,
        string::const_iterator last,
        int& dest,
        string terminators)
{
  long long i;
  NextInt(first, last, i, terminators);
  dest = static_cast<int>(i);
}




inline void
SkipWhiteSpace(string::const_iterator& first,
               string::const_iterator  last)
{
  string::const_iterator p;
  for(p = first; p != last; ++p)
  {
    if (! isspace(*p)) break;//Alas, isspace() is not an "STL predicate".  (or
                             //I would have written this function differently.)
  }
  first = p;
}



inline void
SkipNonWhiteSpace(string::const_iterator& first,
                  string::const_iterator  last)
{
  string::const_iterator p;
  for(p = first; p != last; ++p)
    if (isspace(*p)) break;
  first = p;
}


inline void
SkipWhiteSpaceExceptNewline(string::const_iterator& first,
                            string::const_iterator  last)
{
  string::const_iterator p;
  for(p = first; p != last; ++p)
    if ((! isspace(*p)) || ((*p) == '\n')) break;
  first = p;
}


// ***************************************************************
// ****    Functions for manipulating lists of text symbols.  ****
// ***************************************************************


//The next function converts argc/argv into a vector of strings.
void ConvertArgvToVectorOfStrings(int argc,
                                  char **argv,
                                  vector<string>& dest);


// *** ArgFileToArgs()
// ***  takes a file and separates it into a vector of strings,
// *** each string corresponding to one line in the file.
// *** However:
// *** -If a comment_indicator is supplied, text in each line
// ***  following the comment_indicator is deleted.
// *** -If a merge_lines_indicator is supplied, lines ending in
// ***  the merge_lines_indicator are merged with the following lines.
// ***  (Note: the merge_lines_indicator is deleted.)
// *** -If whitespace_delimited is true, then lines containing multiple
// ***  tokens separated by whitespace, are split into separate strings.
// ***  However, any text surrounded by single or double quotes is not split
// ***  into separate arguments, even if whitespace is present between
// ***  the quotes.  The quotes, however are deleted.
// ***  (Note: Nested quotes like "Sean says 'oh fudge'" are supported.)
// *** 
// *** This feature is meant to be used to enable users to pass arguments
// *** from a file instead of using the shell.  However, this function
// *** does not support the features provided by shells, for example:
// *** 
// *** -Macros and environment variables are not supported.
// *** -Characters-sequences like \[   \]   \'   \"  etc...
// ***  (which are typically converted to [  ]  '  " by most shells)
// ***  are not given any special treatment.
void ArgFileToArgs(istream& in,
                   vector<string>& dest,
                   bool whitespace_delimited = true,
                   string comment_indicator = string(),
                   string merge_lines_indicator = string());



// *** HandleArgFiles()
// *** The objective of HandleArgFiles() is to consolidate
// *** all of the arguments from the following sources
// *** into one argument list (a vector of strings called "dest_args").
// *** 1) If a default file is present (default_filename),
// ***    arguments in this file are inserted into dest_args first.
// *** 2) Then, the original list of arguments "source_args"
// ***    are appended to this list.
// *** 3) If an arg_file_flag is declared, then 
// ***    any time one of the arguments in dest_args equals arg_file_flag,
// ***    the following argument is interpreted as the name of a file
// ***    containing more arguments.  The file is opened, and
// ***    these two arguments are replaced by the arguments inside the
// ***    file.  If the arguments in these files match the arg_file_flag,
// ***    the process will be repeated recursively.
// ***    (Warning: This can lead to infinite loops.)
// *** 4), 5), & 6) The parameters
// ***    whitespace_delimited, comment_indicator, and merge_lines_indicator
// ***    determine how argument files are parsed into separate arguments
// ***    and behave the same way they do in ArgFileToArgs().
// ***    (See comment above.)
// *** Note: arguments "source_args" and "dest_args" should not be the same.
// ***
// *** Debugging: If the preprocessor macro "DBG_ARG_FILE_PARSER" is defined
// ***            to be "1" at the time of compilation of this library,
// ***            (this can be managed by including "-DDBG_ARG_FILE_PARSER"
// ***             in the argument list for gcc or g++),
// ***            then a string of diagnostic messages will be printed
// ***            as argument files are parsed.

void HandleArgFiles(vector<string> const& source_args,
                    vector<string> & dest_args,
                    string arg_file_flag = string(),
                    string comment_indicator = string(),
                    string merge_lines_indicator = string(),
                    string default_filename = string(),
                    bool whitespace_delimited = true);




#ifdef DEBUG
// *** It's useful when debugging, to print out the contents
// *** of the argument list.
inline void
DisplayVectorOfStrings(vector<string> const& strings)
{
  for (vector<string>::const_iterator ps = strings.begin();
       ps != strings.end();
       ++ps)
  {
    cerr << "Arg#" << ps - strings.begin() << ": \"" << *ps << "\"" << endl;
  }
}
#endif // #ifdef DEBUG



} //namespace ParseUtils


#endif //#ifndef _PARSE_UTILS_H
