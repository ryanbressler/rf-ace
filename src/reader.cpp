#include "reader.hpp"
#include <iostream>

using namespace std;

Reader::Reader(const string& fileName, const char delimiter): 
  delimiter_(delimiter) {

  this->init(fileName);

}

Reader::~Reader() {

  if ( inStream_.is_open() ) {
    inStream_.close();
  }

}

void Reader::init(const string& fileName) {

  inStream_.open(fileName.c_str());

  if ( !inStream_.good() ) {
    cerr << "ERROR: failed to open file '" << fileName << "' for reading. Make sure the file exists. Quitting..." << endl;
    exit(1);
  }
  
  this->setLineFeed("");

  nLines_ = 0;

  string line;

  for ( nLines_ = 0; getline(inStream_,line); ++nLines_ ) { }

  this->rewind();

}

bool Reader::endOfLine() const {
  return( lineFeed_.rdbuf()->in_avail() == 0 );
}

bool Reader::nextLine() {
  
  string line;
  
  if ( getline(inStream_,line) ) {
    this->setLineFeed(line);
    return(true);
  } else {
    this->setLineFeed(line);
    return(false);
  }
  
}

bool Reader::skipField() {
  
  string field;
  
  if ( getline(lineFeed_,field,delimiter_) ) {
    return(true);
  } else {
    return(false);
  }
  
}

void Reader::rewind() {

  inStream_.clear();
  inStream_.seekg(ios_base::beg);

  this->setLineFeed("");

}

void Reader::checkLineFeed() const {

  if ( this->endOfLine() ) {
    cerr << "READ ERROR: tried to read from an empty linefeed. Did you forget Reader::nextLine()?" << endl;
    exit(1);
  }

}

void Reader::setLineFeed(const string& str) {
  lineFeed_.clear();
  lineFeed_.str(str);
}
