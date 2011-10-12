#ifndef GRAYCODE_HPP
#define GRAYCODE_HPP

#include<cstdlib>
#include<vector>

using namespace std;

class GrayCode {

public:
  GrayCode();
  ~GrayCode();

private:
  
  typedef char graycode_t;

  vector<graycode_t> graycode;

  //vector<vector<bool> > setBit_;
  //vector<vector<int> > code_;

};

#endif
