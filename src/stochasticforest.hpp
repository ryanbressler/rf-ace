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
  
  num_t getError();
  num_t getOobError();

  void getImportanceValues(vector<num_t>& importanceValues,
                           vector<num_t>& contrastImportanceValues);
    
  //void predict(vector<string>& categoryPrediction, vector<num_t>& confidence);
  //void predict(vector<num_t>& prediction, vector<num_t>& confidence);

  vector<num_t> getPredictions();

  void predict(Treedata* treeDataTest, vector<string>& predictions, vector<num_t>& confidence);
  void predict(Treedata* treeDataTest, vector<num_t>& predictions, vector<num_t>& confidence);

  vector<num_t> getOobPredictions();
  vector<num_t> getPermutedOobPredictions(const size_t featureIdx);

  //Counts the number of nodes in the forest
  vector<size_t> nNodes();
  size_t nNodes(const size_t treeIdx);

  size_t nTrees();

  inline set<size_t> getFeaturesInForest() { return( featuresInForest_ ); }

  inline string getTargetName() { return( parameters_->targetStr ); }
  inline bool isTargetNumerical() { return( rootNodes_[0]->isTargetNumerical() ); }

  void printToFile(const string& fileName);

#ifndef TEST__
private:
#endif
  
  //void growTreesPerThread(vector<RootNode*>& rootNodes);

  void learnRF();
  void learnGBT();

  //Node::GrowInstructions getGrowInstructions();

  // Summarizes predictions across samples and trees in the forest, stored in predictionMatrix
  //vector<num_t> getPredictions(const vector<vector<num_t> >& predictionMatrix, vector<num_t>& predictions, vector<num_t>& confidence);
  
  void growNumericalGBT();
  void growCategoricalGBT();

  // TODO: StochasticForest::transformLogistic() should be moved elsewhere
  void transformLogistic(vector<num_t>& prediction, vector<num_t>& probability);
  
  num_t error(const vector<num_t>& data1,
	      const vector<num_t>& data2); 

  // Pointer to treeData_ object, stores all the feature data with which the trees are grown (i.e. training data)
  Treedata* treeData_;

  options::General_options* parameters_;

  // Chosen target to regress on
  //string targetName_;
  vector<string> targetSupport_;

  // Root nodes for every tree
  vector<RootNode*> rootNodes_;

  // Container for all features in the forest for fast look-up
  set<size_t> featuresInForest_;
  
};

#endif
