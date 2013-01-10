#include <cstdlib>
#include <iostream>

size_t N_SUCCESS = 0;
size_t N_FAIL = 0;

#include "test_reader.hpp"

using namespace std;

int main() {

  test_reader();

  size_t nTests = N_SUCCESS + N_FAIL;

  cout << "TEST: " << nTests << " tests run: " << N_SUCCESS << " successes and " << N_FAIL << " failures" << endl; 

  return( EXIT_SUCCESS );

}
