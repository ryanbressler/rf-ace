#ifndef READER_HPP
#define READER_HPP

#include <cstdlib>
#include <fstream>
#include <string>
#include <sstream>
#include <ios>

#include "utils.hpp"
#include "datadefs.hpp"

class Reader {
public:

  Reader(const std::string& fileName, const char delimiter = '\t');
  ~Reader();
  
  template<typename T> inline friend Reader& operator>>(Reader& reader, T& val) {
    reader.checkLineFeed();
    std::string field;
    std::getline(reader.lineFeed_,field,reader.delimiter_);
    std::stringstream ss( utils::chomp(field) );
    ss >> val;
    return(reader);
  }
  
  bool nextLine();

  bool skipField();

  void rewind();

  bool endOfLine() const;

  size_t nLines() const { return( nLines_ ); }

  void setDelimiter(const char delimiter) { delimiter_ = delimiter; }

#ifndef TEST__
private:
#endif

  void init(const std::string& fileName);

  void checkLineFeed() const;

  void setLineFeed(const string& str);

  std::ifstream inStream_;

  char delimiter_;

  size_t nLines_;

  stringstream lineFeed_;

};

template<> inline Reader& operator>>(Reader& reader, datadefs::num_t& val) {
  reader.checkLineFeed();
  std::string field;
  std::getline(reader.lineFeed_,field,reader.delimiter_);
  field = utils::chomp(field);
  if ( datadefs::isNAN_STR(field) ) {
    val = datadefs::NUM_NAN;
  } else {
    std::stringstream ss( utils::chomp(field) );
    ss >> val;
  }
  return(reader);
}

/*
  template<> inline Reader& operator>>(Reader& reader, datadefs::cat_t& val) {
  reader.checkLineFeed();
  std::string field;
  std::getline(reader.lineFeed_,field,reader.delimiter_);
  field = utils::chomp(field);
  if ( datadefs::isNAN_STR(field) ) {
  val = datadefs::CAT_NAN;
  } else {
  std::stringstream ss( utils::chomp(field) );
  ss >> val;
  }
  return(reader);
  }
*/

template<> inline Reader& operator>>(Reader& reader, string& str) {
  reader.checkLineFeed();
  std::getline(reader.lineFeed_,str,reader.delimiter_);
  str = utils::chomp(str);
  return(reader);
}

#endif
