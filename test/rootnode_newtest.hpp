#ifndef ROOTNODE_NEWTEST_HPP
#define ROOTNODE_NEWTEST_HPP

#include <cstdlib>
#include <set>
#incluce <vector>

#include "newtest.hpp"
#include "rootnode.hpp"
#include "datadefs.hpp"

using namespace std;

void rootnode_newtest_getChildLeafTrainData();

void rootnode_newtest() {

  newtest( "Testing extraction of train samples from child leaf nodes", &rootnode_newtest_getChildLeafTrainData );

}

void rootnode_newtest_getChildLeafTrainData() {

  Rootnode rootNode;
  Node nodeL,nodeR,nodeM;

}

#endif
