#include <iostream>
#include <cassert>
#include "graycode.hpp"

GrayCode::GrayCode(const size_t nMaxBits) {

  assert( nMaxBits > 0 );
  assert( nMaxBits < sizeof(graycode_t) );

  size_t nCodes = (1 << nMaxBits);

  graycode.resize(nCodes,0);
  
  for ( size_t codeIdx = 0; codeIdx < nCodes; ++codeIdx ) {

    // Gray Code as integer
    graycode_t codeAsInt = static_cast<graycode_t>(codeIdx);
    
    // Integer transformed to the actual Gray Code
    // Reference: http://answers.yahoo.com/question/index?qid=20071013065704AAgrGp3
    graycode[codeIdx] |= codeAsInt ^ ( codeAsInt >> 1 );

  }
  
}

GrayCode::~GrayCode() {

  //cout << "GrayCode::dtor" << endl;

}

GrayCode::graycode_t GrayCode::getCode(const size_t codeIdx) {

  return( graycode[codeIdx] );

}
