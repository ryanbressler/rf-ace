#ifndef NODE_NEWTEST_HPP
#define NODE_NEWTEST_HPP

#include <cstdlib>

#include "newtest.hpp"
#include "node.hpp"
#include "datadefs.hpp"

using namespace std;
using datadefs::num_t;

void node_newtest_getChildLeaves();

void node_newtest() {

  newtest( "Testing extraction of child leaf nodes", &node_newtest_getChildLeaves );

}

void node_newtest_getChildLeaves() {

  Node node,nodeL,nodeR,nodeM,nodeLL,nodeLR;

  node.setSplitter("foo",5.0,nodeL,nodeR);
  nodeL.setSplitter("bar",6.0,nodeLL,nodeLR);
  node.missingChild_ = &nodeM;

  nodeLL.setTrainData({1,2,3});
  nodeLR.setTrainData({4,5});
  nodeR.setTrainData({6});
  nodeM.setTrainData({7});

  vector<Node*> childLeaves = node.getChildLeaves();

  set<Node*> childLeavesSet(childLeaves.begin(),childLeaves.end());

  newassert( childLeaves.size() == 4 );
  newassert( childLeavesSet.size() == 4 );
  newassert( childLeavesSet.find(&nodeLL) != childLeavesSet.end() );
  newassert( childLeavesSet.find(&nodeLR) != childLeavesSet.end() );
  newassert( childLeavesSet.find(&nodeR) != childLeavesSet.end() );
  newassert( childLeavesSet.find(&nodeM) != childLeavesSet.end() );

  childLeaves = nodeR.getChildLeaves();
  newassert( childLeaves.size() == 1 );
  newassert( childLeaves[0] == &nodeR );

  childLeaves = nodeM.getChildLeaves();
  newassert( childLeaves.size() == 1 );
  newassert( childLeaves[0] == &nodeM );

  childLeaves = nodeL.getChildLeaves();
  newassert( childLeaves.size() == 2 );
  childLeavesSet = set<Node*>(childLeaves.begin(),childLeaves.end());
  newassert( childLeavesSet.find(&nodeLL) != childLeavesSet.end() );
  newassert( childLeavesSet.find(&nodeLR) != childLeavesSet.end() );
  
  vector<num_t> trainData = nodeLL.getTrainData();
  set<num_t> trainDataSet(trainData.begin(),trainData.end());

  newassert( trainDataSet.find(1) != trainDataSet.end() );
  newassert( trainDataSet.find(2) != trainDataSet.end() );
  newassert( trainDataSet.find(3) != trainDataSet.end() );
  
}

#endif
