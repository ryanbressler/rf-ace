#ifndef RANDOMFOREST_HPP
#define RANDOMFOREST_HPP

#include<cstdlib>
#include "rootnode.hpp"
#include "treedata.hpp"

using namespace std;

class Randomforest {
public:

  Randomforest(Treedata* treedata, 
               size_t targetIdx, 
               size_t nTrees, 
               size_t mTry, 
               size_t nodeSize, 
               bool useContrasts,
               bool isOptimizedNodeSplit);
  ~Randomforest();

  //Selects the target feature that is to be predicted
  //void setTarget(size_t targetIdx);

  //Gets the selected target feature
  //size_t getTarget();

  //void learn();

  //Calculates the feature importance scores for real and contrast features
  vector<size_t> featureFrequency();

  //IMPLEMENTATION MISSING
  vector<num_t> featureImportance();
  
  //Counts the number of nodes in the forest
  size_t nNodes();


private:

  //Resets all internal variables in the forest
  //void initialize();

  //Grows the forest with respect to selected target feature
  void growForest();

  //Grows one tree with respect to selected target feature
  //void growTree(size_t treeIdx);

  //Recursive tree-generating node splitter algorithm
  //void recursiveNodeSplit(const size_t treeIdx, Node* node, const vector<size_t>& sampleIcs);

  //Percolates samples along the trees, starting from the rootNode. Spits out a map<> that ties the percolated train samples to the leaf nodes
  //NOTE: currently, since there's no implementation for a RootNode class, there's no good way to store the percolated samples in the tree, but a map<> is generated instead
  void percolateSampleIcs(Node* rootNode, vector<size_t>& sampleIcs, map<Node*,vector<size_t> >& trainIcs);
  void percolateSampleIcsAtRandom(size_t featureIdx, Node* rootNode, vector<size_t>& sampleIcs, map<Node*,vector<size_t> >& trainIcs);
  
  void percolateSampleIdx(size_t sampleIdx, Node** nodep);
  void percolateSampleIdxAtRandom(size_t featureIdx, size_t sampleIdx, Node** nodep);

  //A check if a feature exists in the tree
  //bool isFeatureInTree(size_t featureIdx, size_t treeIdx);

  //Given the map<>, generated with the percolation functions, a tree impurity score is outputted
  void treeImpurity(map<Node*,vector<size_t> >& trainIcs, num_t& impurity);

  //void forestImpurity(map<Node*,vector<size_t> >& trainIcs, num_t& impurity);

  //Calculates the feature importance scores for real and contrast features
  //vector<num_t> featureImportance();

  //Chosen target to regress on
  size_t targetIdx_;

  //Pointer to treedata_ object, stores all the feature data
  Treedata* treedata_;
  
  //Size of the forests
  size_t nTrees_;
  //size_t mTry_;
  //size_t nodeSize_;

  //Root nodes for every tree
  vector<RootNode*> rootNodes_;

  //A helper variable that stores all splitter features. This will make impurity calculation faster
  map<size_t, set<size_t> > featuresInForest_;

  //Stores the forest. forest_[i][j] is the j'th node of i'th tree. forest_[i][0] is the rootnode
  //NOTE: will be reduced to just RootNodes in the future
  //vector<vector<Node*> > forest_;
  

  //Number of used nodes in each tree
  //vector<size_t> nNodes_;

  //Out-of-box samples for each tree
  vector<vector<size_t> > oobMatrix_;

  //bool useContrasts_;
  
};

#endif
