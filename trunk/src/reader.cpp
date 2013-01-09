#include "reader.hpp"
#include <iostream>
#include <string>

using namespace std;

Reader::Reader(const string& fileName, const char delimiter, const char newLine, const string& naStr): 
  delimiter_(delimiter),
  newLine_(newLine),
  naStr_(naStr) {

  inStream_.open(fileName);

  

}

Reader::~Reader() {

  if ( inStream_.is_open() ) {
    inStream_.close();
  }

}

ostream& operator<<(ostream& outStream, Reader& reader) {

  string field;

  getline(reader.inStream_,field,reader.delimiter_);

  return( outStream << field );

}
