#include "graycode.hpp"

GrayCode::GrayCode() {
  
  graycode.clear();
  graycode.resize(2^10);
  
  cout << "GrayCode::ctor" << endl;

}

GrayCode::~GrayCode() {

  cout << "GrayCode::dtor" << endl;

}
