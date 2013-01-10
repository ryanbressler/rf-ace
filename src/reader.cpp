#include "reader.hpp"
#include <iostream>

using namespace std;

Reader::Reader(const string& fileName, const char delimiter, const string& naStr): 
  delimiter_(delimiter),
  naStr_(naStr) {

  this->init(fileName);

}

Reader::~Reader() {

  if ( inStream_.is_open() ) {
    inStream_.close();
  }

}

void Reader::init(const string& fileName) {

  inStream_.open(fileName.c_str());

  this->setLineFeed("");

  nCols_ = 0;
  nRows_ = 0;

  string line;
  getline(inStream_,line);

  ++nRows_;

  stringstream ss(line);

  string field;
  while ( getline(ss,field,delimiter_) ) {
    ++nCols_;
  }

  while ( getline(inStream_,line) ) {
    ++nRows_;
  }

  this->rewind();

}

bool Reader::endOfLine() {
  
  return( lineFeed_.rdbuf()->in_avail() == 0 );

}

Reader& Reader::nextLine() {
  
  string line;

  getline(inStream_,line);

  this->setLineFeed(line);

  return(*this);

}

Reader& Reader::skipField() {

  string field;

  std::getline(lineFeed_,field,delimiter_);

  return(*this);

}

Reader& Reader::rewind() {

  inStream_.clear();
  inStream_.seekg(ios_base::beg);

  this->setLineFeed("");

  return(*this);

}

void Reader::checkLineFeed() {

  if ( this->endOfLine() ) {
    cerr << "READ ERROR: tried to read from an empty linefeed. Did you forget Reader::nextLine()?" << endl;
    exit(1);
  }

}

void Reader::setLineFeed(const string& str) {
  lineFeed_.clear();
  lineFeed_.str(str);
}
