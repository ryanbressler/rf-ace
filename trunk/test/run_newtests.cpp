#include <cstdlib>
#include <iostream>

#include "newtest.hpp"
#include "reader_newtest.hpp"
#include "treedata_newtest.hpp"
#include "rface_newtest.hpp"

using namespace std;

int main() {

  newtestinit();

  reader_newtest();
  treedata_newtest();
  rface_newtest();

  newtestdone();

  return( EXIT_SUCCESS );

}
