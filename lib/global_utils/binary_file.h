#ifndef _BINARY_FILE_H
#define _BINARY_FILE_H

#include<fstream>
using namespace std;


template<class T>
//void binary_write(ofstream& file, T const& entry)
fstream& binary_write(fstream& file, T const& entry)
{
  size_t size = sizeof(T);
  //  char *p = static_cast<char*>(&entry);
  char *p = (char*)&entry;
  file.write(p, size);
  return file;
}

template<class T>
//void binary_read(ifstream& file, T& entry)
fstream& binary_read(fstream& file, T& entry)
{
  size_t size = sizeof(T);
  //  char *p = static_cast<char*>(&entry);
  char *p = (char*)&entry;
  file.read(p, size);
  return file;
}


#endif //#ifndef _BINARY_FILE_H
