#ifndef READER_HPP
#define READER_HPP

#include <cstdlib>
#include <fstream>
#include <string>
#include <sstream>
#include <ios>

#include "utils.hpp"

class Reader {
public:

  Reader(const std::string& fileName, const char delimiter = '\t', const std::string& naStr = "NA");
  ~Reader();
  
  template<typename T> inline friend Reader& operator>>(Reader& reader, T& val) {
    reader.checkLineFeed();
    std::string field;
    std::getline(reader.lineFeed_,field,reader.delimiter_);
    std::stringstream ss( utils::chomp(field) );
    ss >> val;
    return(reader);
  }
  
  Reader& nextLine();

  Reader& skipField();

  Reader& rewind();

  bool endOfLine() const;
  bool endOfFile() const;

  size_t nLines() const { return( nLines_ ); }

#ifndef TEST__
private:
#endif

  void init(const std::string& fileName);

  void checkLineFeed() const;

  void setLineFeed(const string& str);

  bool endOfStream(const ios& stream) const;

  std::ifstream inStream_;

  char delimiter_;
  std::string naStr_;

  size_t nLines_;

  stringstream lineFeed_;

};

template<> inline Reader& operator>>(Reader& reader, datadefs::num_t& val) {
  reader.checkLineFeed();
  std::string field;
  std::getline(reader.lineFeed_,field,reader.delimiter_);
  field = utils::chomp(field);
  if ( field == reader.naStr_ ) {
    val = datadefs::NUM_NAN;
  } else {
    std::stringstream ss( utils::chomp(field) );
    ss >> val;
  }
  return(reader);
}

template<> inline Reader& operator>>(Reader& reader, string& str) {
  reader.checkLineFeed();
  std::getline(reader.lineFeed_,str,reader.delimiter_);
  str = utils::chomp(str);
  return(reader);
}

#endif
