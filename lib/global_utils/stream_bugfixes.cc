#ifdef GETLINE_UNSUPPORTED

#include <cassert> //defines assert()
#include <istream>
#include <string>
using namespace std;

bool
SubstituteForGetline(istream &in,
                     string& line)
{
  line.erase();
  char c;

               // -- A subtle behavior of getline() --
  in.get(c);   //These four lines make sure that the stream fails when we
  if (! in)    //are initially at the end of the file.  The only circumstances
    return false;//in which we should allow the stream to fail is when
  in.putback(c); //we are already at the end of the file.  Otherwise it should
               //read in whatever characters are available and leave
               //the stream in a good state.

  while(in.get(c))
  {
    if (c == '\n')
      break;
    line.push_back(c);
  }
  in.clear();  //This makes the stream "good". (It sets ios_base::goodbit)

  assert(in);  //Is the stream good?  It better be.

  return in;
}
#endif //#ifdef GETLINE_UNSUPPORTED
