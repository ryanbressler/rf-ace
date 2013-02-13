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

  cout << endl << "Testing Reader class:" << endl; 
  reader_newtest();

  cout << endl << "Testing Treedata class:" << endl; 
  treedata_newtest();

  cout << endl << "Testing RFACE class:" << endl;
  rface_newtest();

  cout << endl << "Testing Distributions namespace:" << endl;
  distributions_newtest();
  
  cout << endl << "Testing Utils namespace:" << endl;
  utils_newtest();
  
  cout << endl << "Testing Datadefs namespace:" << endl;
  datadefs_newtest();

  cout << endl << "Testing Node class:" << endl;
  node_newtest();
  //rootnode_newtest();

  newtestdone();

  return( EXIT_SUCCESS );

}
