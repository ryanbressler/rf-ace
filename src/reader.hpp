#ifndef READER_HPP
#define READER_HPP

#include <cstdlib>
#include <fstream>
#include <string>
#include <sstream>

#include "utils.hpp"

class Reader {
public:

  Reader(const std::string& fileName, const char delimiter = '\t', const char newLine = '\n', const std::string& naStr = "NA");
  ~Reader();
  
  friend std::ostream& operator<<(std::ostream& outStream, Reader& reader) {
    std::string field;
    std::getline(reader.inStream_,field,reader.delimiter_);
    field = utils::chomp(field);
    return( outStream << field );
  }
  
  template<typename T> inline friend void operator>>(Reader& reader, T& val) {
    std::string field;
    std::getline(reader.inStream_,field,reader.delimiter_);
    std::stringstream ss( utils::chomp(field) );
    ss >> val;
  }
  
  void skipLine();
  void skipField();

  size_t nRows() { return( nRows_ ); }
  size_t nCols() { return( nCols_ ); }

#ifndef TEST
private:
#endif

  void init(const std::string& fileName);

  std::ifstream inStream_;

  char delimiter_;
  char newLine_;
  std::string naStr_;

  size_t nCols_;
  size_t nRows_;

};

template<> inline void operator>>(Reader& reader, datadefs::num_t& val) {
  std::string field;
  std::getline(reader.inStream_,field,reader.delimiter_);
  field = utils::chomp(field);
  if ( field == reader.naStr_ ) {
    val = datadefs::NUM_NAN;
  } else {
    std::stringstream ss( utils::chomp(field) );
    ss >> val;
  }
}

#endif
