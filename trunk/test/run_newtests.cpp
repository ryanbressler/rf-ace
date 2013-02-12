#include <cstdlib>
#include <iostream>

#include "newtest.hpp"
#include "reader_newtest.hpp"
#include "treedata_newtest.hpp"
#include "rface_newtest.hpp"
#include "distributions_newtest.hpp"
#include "utils_newtest.hpp"
#include "datadefs_newtest.hpp"
#include "node_newtest.hpp"
//#include "rootnode_newtest.hpp"

using namespace std;

int main() {

  newtestinit();

  reader_newtest();
  treedata_newtest();
  rface_newtest();
  distributions_newtest();
  utils_newtest();
  datadefs_newtest();
  node_newtest();
  //rootnode_newtest();

  newtestdone();

  return( EXIT_SUCCESS );

}
