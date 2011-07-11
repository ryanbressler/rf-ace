#ifndef GBT_HPP
#define GBT_HPP

#include <cstdlib>
#include <algorithm>
#include "node.hpp"
#include "rootnode.hpp"
#include "treedata.hpp"
#include "datadefs.hpp"

#define LOG_OF_MAX_FLT 70.0f


using namespace std;

class GBT
{
public:

  // Constructs Gradient Boosting Tree "forest".
  // nmaxleaves is typically a small value, such as 6.
  // shrinkage is between 0 and 1, typically something like 0.125
  // subSampleSize is typically 0.4 to 0.6
  GBT(Treedata* treeData, size_t targetIdx, size_t nTrees, size_t nMaxLeaves, num_t shrinkage, num_t subSampleSize);
  ~GBT();
  // produce predictions of a data set
  void predictForest(Treedata* treeData, vector<num_t>& prediction);


private:

  // actual growing
  void growForest();
  void growForestNumerical();
  void growForestCategorical();
  //void growTree(size_t treeidx);

  void predictForestNumerical(  Treedata* treeData, vector<num_t>& prediction);
  void predictForestCategorical(Treedata* treeData, vector<num_t>& categoryPrediction);
  num_t predictSampleByTree(size_t sampleidx, size_t treeidx); // REMOVED Treedata* treeData as it's readily stored in treeData_
  void predictDatasetByTree(size_t treeidx, vector<num_t>& prediction);

  //void recursiveNodeSplit(size_t treeidx, Node* node, vector<size_t>& sampleics);
  void SetNodePrediction( size_t treeidx, Node* node, vector<size_t>& sampleics);
  void transformLogistic(vector<num_t>& prediction, vector<num_t>& probability);

  Treedata* treeData_;
  size_t targetIdx_;      // The index of the true target in treeData_
  size_t nTrees_;
  size_t nMaxLeaves_;
  //size_t nMaxNodes_;
  size_t nodeSize_;		// smallest node size
  float shrinkage_;
  float subSampleSize_;
  size_t numClasses_;
  size_t nLeavesCounter_; 	// counter for constructing each tree

  //Stores the rootnodes for the forest. Access the child nodes by calling rootNodes_[treeIdx]->leftChild() etc.
  //rootNodes_[treeIdx]->hasChildren() returns true if the node has children
  vector<RootNode*> rootNodes_;
};

#endif
