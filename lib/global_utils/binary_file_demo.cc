#include "binary_file.h"

struct MyType
{
  long  i;
  float x;
  MyType(long set_i, float set_x):i(set_i),x(set_x) {}
};

ostream& operator << (ostream& file,
                      MyType& m)
{ file << "i = " << m.i << "\n" << "x = " << m.x << "\n"; return file; }


int
main()
{
  cout << "Testing binary file read/write\n";
  MyType m(12345678, -1234.5678);
  cout << "Initially the struct m, has the contents:\n" << m << endl;

  cout << "...saving to file \"test_binary_file\"..." << endl;

  {
    fstream test_binary_file("test_binary_file");
    binary_write(test_binary_file, m);
  }

  cout << "...reading from file \"test_binary_file\"...\n" << endl;
  MyType n(0, 0.0f);
  {
    fstream test_binary_file("test_binary_file");
    binary_read(test_binary_file, n);
  }
  cout << "Contents after saving and reloading:\n" << n << endl;
  return 0;
}
