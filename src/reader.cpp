#include "reader.hpp"
#include <iostream>

using namespace std;

Reader::Reader(const string& fileName, const char delimiter, const char newLine, const string& naStr): 
  delimiter_(delimiter),
  newLine_(newLine),
  naStr_(naStr) {

  this->init(fileName);

}

Reader::~Reader() {

  if ( inStream_.is_open() ) {
    inStream_.close();
  }

}

void Reader::init(const string& fileName) {

  inStream_.open(fileName);

  nCols_ = 0;
  nRows_ = 0;

  string line;
  getline(inStream_,line,newLine_);

  ++nRows_;

  stringstream ss(line);

  string field;
  while ( getline(ss,field,delimiter_) ) {
    ++nCols_;
  }

  while ( getline(inStream_,line,newLine_) ) {
    ++nRows_;
  }

  inStream_.clear();
  inStream_.seekg(ios_base::beg);

}

void Reader::skipLine() {

  string line;
  getline(inStream_,line,newLine_);
  //return( utils::chomp(line) );

}

void Reader::skipField() {

  string field;

  std::getline(inStream_,field,delimiter_);

}


