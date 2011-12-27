#ifndef STOCHASTICFOREST_HPP
#define STOCHASTICFOREST_HPP

#include<cstdlib>
#include "rootnode.hpp"
#include "treedata.hpp"
#include "partitionsequence.hpp"

using namespace std;

class StochasticForest {
public:

  // !! Documentation: instruct the constructor to choose between the two forest implementations
  StochasticForest(Treedata* treeData, const size_t targetIdx, const size_t nTrees);  
  ~StochasticForest();

  // !! Documentation: ... then, on the public side, learning would become collapsed into a single overloaded function learn(...)
  void learnRF(const size_t mTry,
	       const size_t nMaxLeaves,
	       const size_t nodeSize,
	       const bool useContrasts);

  void learnGBT(const size_t nMaxLeaves,
		const num_t shrinkage,
		const num_t subSampleSize);
  
  vector<size_t> featureFrequency();
  
  //Calculates the feature importance scores for real and contrast features
  vector<num_t> featureImportance();
  
  void predict(vector<string>& prediction, vector<num_t>& confidence);
  void predict(vector<num_t>& prediction, vector<num_t>& confidence);
  
  //void predict(Treedata* treeData, vector<string>& prediction, vector<num_t>& confidence);
  //void predict(Treedata* treeData, vector<num_t>& prediction, vector<num_t>& confidence);

  //Counts the number of nodes in the forest
  size_t nNodes();

  void printToFile(const string& fileName);

#ifndef TEST__
private:
#endif

  void growNumericalGBT();
  void growCategoricalGBT();

  // TODO: StochasticForest::transformLogistic() should be moved elsewhere
  void transformLogistic(const size_t numClasses, vector<num_t>& prediction, vector<num_t>& probability);

  // TODO: predictSampleByTree() and percolateSampleIdx families in StochasticForest need to be fused together 
  num_t predictSampleByTree(size_t sampleIdx, size_t treeIdx);
  vector<num_t> predictDatasetByTree(size_t treeIdx);

  void predictWithCategoricalRF(vector<string>& categoryPrediction);
  void predictWithNumericalRF(vector<num_t>& prediction);

  void predictWithCategoricalGBT(vector<string>& categoryPrediction, vector<num_t>& confidence);
  void predictWithNumericalGBT(vector<num_t>& prediction, vector<num_t>& confidence);

  //Percolates samples along the trees, starting from the rootNode. Spits out a map<> that ties the percolated train samples to the leaf nodes
  void percolateSampleIcs(Node* rootNode, const vector<size_t>& sampleIcs, map<Node*,vector<size_t> >& trainIcs);
  void percolateSampleIcsAtRandom(const size_t featureIdx, Node* rootNode, const vector<size_t>& sampleIcs, map<Node*,vector<size_t> >& trainIcs);
  
  void percolateSampleIdx(const size_t sampleIdx, Node** nodep);
  void percolateSampleIdxAtRandom(const size_t featureIdx, const size_t sampleIdx, Node** nodep);

  //Given the map<>, generated with the percolation functions, a tree impurity score is outputted
  void treeImpurity(map<Node*,vector<size_t> >& trainIcs, num_t& impurity);

  //Pointer to treeData_ object, stores all the feature data with which the trees are grown (i.e. training data)
  Treedata* treeData_;

  //Chosen target to regress on
  size_t targetIdx_;
  
  size_t nTrees_;

  //Root nodes for every tree
  vector<RootNode*> rootNodes_;

  //A helper variable that stores all splitter features. This will make impurity calculation faster
  map<size_t, set<size_t> > featuresInForest_;

  //Out-of-box samples for each tree
  vector<vector<size_t> > oobMatrix_;

  enum LearnedModel {NO_MODEL, RF_MODEL, GBT_MODEL};
  
  LearnedModel learnedModel_;
  size_t numClasses_;
  num_t shrinkage_;

  // This object generates a sequence that makes splitting data with categorical splitters fast
  PartitionSequence *partitionSequence_;

};

#endif
