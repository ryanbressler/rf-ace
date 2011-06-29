#ifndef RANDOMFOREST_HPP
#define RANDOMFOREST_HPP

#include<cstdlib>
#include "node.hpp"
#include "treedata.hpp"

using namespace std;

class Randomforest
{
public:

  //Initializes Random Forest with ntrees, mtry, and nodesize
  //NOTE: target feature will be set to 0. Change with Randomforest::select_target()
  Randomforest(Treedata* treedata, size_t targetIdx, size_t nTrees, size_t mTry, size_t nodeSize);
  ~Randomforest();

  void initialize();

  //Selects the target feature that is to be predicted
  void setTarget(size_t targetIdx);

  //Gets the selected target feature
  size_t getTarget();

  //Grow the Random Forest with respect to selected target feature
  void learn(const size_t nPerms, vector<num_t>& pValues, vector<num_t>& importanceValues);
  //void calculate_importance(const num_t alpha, vector<num_t>& importance, num_t& contrast_prc);

  //blacklist_and_kill

private:

  //Grows one tree with respect to selected target feature
  void growTree(size_t treeIdx);

  //Recursive tree-generating node splitter algorithm
  //NOTE: there will be at least two alternative node splitter algorithms in the future
  void recursiveNodeSplit(const size_t treeIdx, const size_t nodeIdx, const vector<size_t>& sampleIcs);

  void percolateSampleIcs(Node& rootNode, vector<size_t>& sampleIcs, map<Node*,vector<size_t> >& trainIcs);
  void percolateSampleIcsAtRandom(size_t featureIdx, Node& rootNode, vector<size_t>& sampleIcs, map<Node*,vector<size_t> >& trainIcs);
  
  void percolateSampleIdx(size_t sampleIdx, Node** nodep);
  void percolateSampleIdxAtRandom(size_t featureIdx, size_t sampleIdx, Node** nodep);

  bool isFeatureInTree(size_t featureIdx, size_t treeIdx);

  void treeImpurity(map<Node*,vector<size_t> >& trainIcs, num_t& impurity);

  //void calculate_importance(const num_t alpha, vector<num_t>& importance, num_t& contrast_prc);
  void calculateImportance(vector<num_t>& importance);

  size_t targetIdx_;

  //Pointer to treedata_ object, stores all the feature information
  Treedata* treedata_;
  
  //Size of the forest
  size_t nTrees_;
  size_t mTry_;
  size_t nodeSize_;

  vector<vector<Node> > forest_; //forest_[i][j] is the j'th node of i'th tree. forest_[i][0] is the rootnode.
  vector<size_t> nNodes_; //Number of used nodes in each tree.
  vector<vector<size_t> > oobMatrix_;
  
  //vector<map<Node*,vector<size_t> > > trainics_;

};

#endif
