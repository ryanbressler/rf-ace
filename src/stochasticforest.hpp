#ifndef STOCHASTICFOREST_HPP
#define STOCHASTICFOREST_HPP

#include <cstdlib>
#include <fstream>
#include "rootnode.hpp"
#include "treedata.hpp"
#include "options.hpp"
#include "distributions.hpp"

using namespace std;

class StochasticForest {
public:
  
  StochasticForest();
  
  ~StochasticForest();

  void learnRF(Treedata* trainData, const size_t targetIdx, const ForestOptions* forestOptions, const vector<num_t>& featureWeights, vector<distributions::Random>& randoms);
  void learnGBT(Treedata* trainData, const size_t targetIdx, const ForestOptions* forestOptions, const vector<num_t>& featureWeights, vector<distributions::Random>& randoms);

  void loadForest(const string& fileName);
  
  //num_t getError() { return(0.0); }
  //num_t getOobError();

  //void getImportanceValues(Treedata* trainData, vector<num_t>& importanceValues, vector<num_t>& contrastImportanceValues);
  void getMeanMinimalDepthValues(Treedata* trainData, vector<num_t>& depthValues, vector<num_t>& contrastDepthValues);

  void predict(Treedata* treeDataTest, vector<string>& predictions, vector<num_t>& confidence, size_t nThreads = 1);
  void predict(Treedata* treeDataTest, vector<num_t>& predictions, vector<num_t>& confidence, size_t nThreads = 1);

  //vector<num_t> getOobPredictions();
  //vector<num_t> getPermutedOobPredictions(const size_t featureIdx);

  //Counts the number of nodes in the forest
  size_t nNodes();
  size_t nNodes(const size_t treeIdx);

  size_t nTrees();

  RootNode* tree(const size_t treeIdx) { return( rootNodes_[treeIdx] ); }

  inline set<size_t> getFeaturesInForest() { return( featuresInForest_ ); }

  inline string getTargetName() { return( targetName_ ); }
  inline bool isTargetNumerical() { return( isTargetNumerical_ ); }

  void saveForest(ofstream& toFile);

#ifndef TEST__
private:
#endif
  
  void growNumericalGBT(Treedata* trainData, const size_t targetIdx, const ForestOptions* forestOptions, const distributions::PMF* pmf, vector<distributions::Random>& randoms);
  void growCategoricalGBT(Treedata* trainData, const size_t targetIdx, const ForestOptions* forestOptions, const distributions::PMF* pmf, vector<distributions::Random>& randoms);

  // TODO: StochasticForest::transformLogistic() should be moved elsewhere
  void transformLogistic(size_t nCategories, vector<num_t>& prediction, vector<num_t>& probability);
  
  num_t error(const vector<num_t>& data1,
	      const vector<num_t>& data2); 

  datadefs::forest_t forestType_;

  vector<string> categories_;

  vector<num_t> GBTConstants_;
  num_t GBTShrinkage_;

  //num_t shrinkage_;

  // Chosen target to regress on
  string targetName_;
  bool isTargetNumerical_;

  // Root nodes for every tree
  vector<RootNode*> rootNodes_;

  // Container for all features in the forest for fast look-up
  set<size_t> featuresInForest_;
  
};

#endif
