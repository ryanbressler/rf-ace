#ifndef GRAYCODE_HPP
#define GRAYCODE_HPP

#include<cstdlib>
#include<vector>

using namespace std;

class GrayCode {

public:

  typedef size_t graycode_t;

  GrayCode(const size_t nMaxBits);
  ~GrayCode();

  graycode_t getCode(const size_t codeIdx);

private:
  
  vector<graycode_t> graycode;

};

#endif
