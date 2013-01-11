#ifndef NEWTEST_HPP
#define NEWTEST_HPP

#include <cstdlib>
#include <iostream>
#include <sstream>

std::stringstream ERRLOG;

#define newassert(condition) { if(!(condition)){ ERRLOG << " => FAIL: " << #condition << " @ " << __FILE__ << " (" << __LINE__ << ")" << std::endl; N_FAIL++; } else { N_SUCCESS++; } }

extern size_t N_SUCCESS;
extern size_t N_FAIL;

void newtest(const std::string& info, void (*testFunc)(void) ) {

  size_t nOldSuccess = N_SUCCESS;
  size_t nOldFail = N_FAIL;
  size_t nOldTests = N_SUCCESS + N_FAIL;

  std::cout << "TEST: " << info << "..." << std::flush; 
  testFunc();
  std::cout << " DONE [ " << N_SUCCESS - nOldSuccess << " / " << N_SUCCESS + N_FAIL - nOldTests << " OK ] " << std::flush;

  if ( N_FAIL > nOldFail ) {
    std::cout << " !! " << N_FAIL - nOldFail << " FAILURES !! " << std::flush;
  }

  std::cout << std::endl;

  std::string errLine;

  while( std::getline(ERRLOG,errLine) ) {
    std::cerr << errLine << std::endl;
  }

}

#endif
