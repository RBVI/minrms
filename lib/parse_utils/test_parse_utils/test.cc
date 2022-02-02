// I tried compiling with
// pgcc test.cc
// and
// pgcc -D__USE_ISOC99  test.cc

#include <iostream>
#include <cstdlib> //necessary for strtold()
using namespace std;

int main(int argc, char **argv)
{
  long long i;
  i = strtoll(argv[1]);
  cout << 3*i << endl;

  #if defined __USE_ISOC99
  USE_ISOC99_is_defined;
  #endif

  #if defined __GNUC__
  __GNUC___is_defined;
  #endif

  #if defined __USE_MISC
  __USE_MISC_is_defined;
  #endif
}
