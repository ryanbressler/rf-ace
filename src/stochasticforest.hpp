#ifndef STOCHASTICFOREST_HPP
#define STOCHASTICFOREST_HPP

#include<cstdlib>
#include "rootnode.hpp"
#include "treedata.hpp"
#include "options.hpp"

using namespace std;

class StochasticForest {
public:

  StochasticForest(Treedata* treeData, options::General_options& parameters);
  
  // Load an existing forest
  StochasticForest(Treedata* treeData, const string& forestFile);

  ~StochasticForest();
  
  void learnRF();
  void learnGBT();

  num_t getError();
  num_t getOobError();

  void getImportanceValues(vector<num_t>& importanceValues,
                           vector<num_t>& contrastImportanceValues);
    
  void predict(vector<string>& categoryPrediction, vector<num_t>& confidence);
  void predict(vector<num_t>& prediction, vector<num_t>& confidence);

  vector<num_t> getPredictions();
  vector<num_t> getOobPredictions();
  vector<num_t> getPermutedOobPredictions(const size_t featureIdx);

  //Counts the number of nodes in the forest
  vector<size_t> nNodes();
  size_t nNodes(const size_t treeIdx);

  size_t nTrees();

  inline set<size_t> getFeaturesInForest() { return( featuresInForest_ ); }

  inline string getTargetName() { return( targetName_ ); }
  inline bool isTargetNumerical() { return( targetSupport_.size() == 0 ? true : false ); }

  void printToFile(const string& fileName);

#ifndef TEST__
private:
#endif

  // Summarizes predictions across samples and trees in the forest, stored in predictionMatrix
  vector<num_t> getPredictions(const vector<vector<num_t> >& predictionMatrix);

  void growNumericalGBT(const Node::GrowInstructions& GI);
  void growCategoricalGBT(const Node::GrowInstructions& GI);

  // TODO: StochasticForest::transformLogistic() should be moved elsewhere
  void transformLogistic(vector<num_t>& prediction, vector<num_t>& probability);
  
  num_t error(const vector<num_t>& data1,
	      const vector<num_t>& data2); 

  // Pointer to treeData_ object, stores all the feature data with which the trees are grown (i.e. training data)
  Treedata* treeData_;

  options::General_options parameters_;

  // Chosen target to regress on
  string targetName_;
  vector<string> targetSupport_;

  // Root nodes for every tree
  vector<RootNode*> rootNodes_;

  // Container for all features in the forest for fast look-up
  set<size_t> featuresInForest_;
  
};

#endif
