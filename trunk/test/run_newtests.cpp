#include <cstdlib>
#include <iostream>

#include "newtest.hpp"
#include "reader_newtest.hpp"
#include "treedata_newtest.hpp"

using namespace std;

int main() {

  reader_newtest();
  treedata_newtest();

  newtestdone();

  return( EXIT_SUCCESS );

}
