#ifndef READER_HPP
#define READER_HPP

#include <cstdlib>
#include <fstream>

class Reader {
public:

  Reader(const std::string& fileName, const char delimiter = '\t', const char newLine = '\n', const std::string& naStr = "NA");
  ~Reader();

  friend std::ostream& operator<<(std::ostream& outStream, Reader& reader);

private:

  std::ifstream inStream_;

  char delimiter_;
  char newLine_;
  std::string naStr_;

};

#endif
