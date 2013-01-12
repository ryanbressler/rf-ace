#ifndef NEWTEST_HPP
#define NEWTEST_HPP

#define TEST__

#include <cstdlib>
#include <iostream>
#include <sstream>

std::stringstream ERRLOG;

size_t N_SUCCESS = 0;
size_t N_FAIL = 0;

#define newassert(condition) { if(!(condition)){ ERRLOG << " => FAIL: " << #condition << " @ " << __FILE__ << " (" << __LINE__ << ")" << std::endl; N_FAIL++; } else { N_SUCCESS++; } }

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

  ERRLOG.clear();

}

void newtestdone() {
  
  std::cout << "ALL DONE! " << N_SUCCESS + N_FAIL << " tests run: " << N_SUCCESS << " successes and " << N_FAIL << " failures" << std::endl;

}

#endif
