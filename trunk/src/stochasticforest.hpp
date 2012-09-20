#ifndef STOCHASTICFOREST_HPP
#define STOCHASTICFOREST_HPP

#include <cstdlib>
//#include <thread>
#include "rootnode.hpp"
#include "treedata.hpp"
#include "options.hpp"
#include "distributions.hpp"

using namespace std;

class StochasticForest {
public:
  
  StochasticForest(Treedata* treeData, options::General_options* parameters);
  
  // Load an existing forest
  StochasticForest(options::General_options* parameters);

  ~StochasticForest();
  
  num_t getError() { return(0.0); }
  num_t getOobError();

  void getImportanceValues(vector<num_t>& importanceValues,
                           vector<num_t>& contrastImportanceValues);
    
  void getMeanMinimalDepthValues(vector<num_t>& depthValues, vector<num_t>& contrastDepthValues);

  void predict(Treedata* treeDataTest, vector<string>& predictions, vector<num_t>& confidence);
  void predict(Treedata* treeDataTest, vector<num_t>& predictions, vector<num_t>& confidence);

  vector<num_t> getOobPredictions();
  vector<num_t> getPermutedOobPredictions(const size_t featureIdx);

  //Counts the number of nodes in the forest
  size_t nNodes();
  size_t nNodes(const size_t treeIdx);

  size_t nTrees();

  RootNode* tree(const size_t treeIdx) { return( rootNodes_[treeIdx] ); }

  inline set<size_t> getFeaturesInForest() { return( featuresInForest_ ); }

  inline string getTargetName() { return( parameters_->targetStr ); }
  inline bool isTargetNumerical() { return( isTargetNumerical_ ); }

  void printToFile(const string& fileName);

#ifndef TEST__
private:
#endif
  
  void learnRF();
  void learnGBT();

  void growNumericalGBT();
  void growCategoricalGBT();

  // TODO: StochasticForest::transformLogistic() should be moved elsewhere
  void transformLogistic(size_t nCategories, vector<num_t>& prediction, vector<num_t>& probability);
  
  num_t error(const vector<num_t>& data1,
	      const vector<num_t>& data2); 

  // Pointer to treeData_ object, stores all the feature data with which the trees are grown (i.e. training data)
  Treedata* trainData_;

  options::General_options* parameters_;

  // Experimental parameters for making GBT working
  vector<num_t> GBTfactors_;
  num_t GBTconstant_;

  // Chosen target to regress on
  //string targetName_;
  bool isTargetNumerical_;

  // Root nodes for every tree
  vector<RootNode*> rootNodes_;

  // Container for all features in the forest for fast look-up
  set<size_t> featuresInForest_;
  
};

#endif
