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

  void learnRF(TreeData* trainData, const size_t targetIdx, const ForestOptions* forestOptions, const vector<num_t>& featureWeights, vector<distributions::Random>& randoms);
  void learnGBT(TreeData* trainData, const size_t targetIdx, const ForestOptions* forestOptions, const vector<num_t>& featureWeights, vector<distributions::Random>& randoms);

  void loadForest(const string& fileName);


  void trainForestAndPredictQuantiles(TreeData* trainData,
				      const size_t targetIdx,
				      TreeData* testData,
				      distributions::PMF* pmf,
				      ForestOptions* forestOptions,
				      distributions::Random* random,
				      vector<vector<num_t> >& predictions); 

  //num_t getError() { return(0.0); }
  //num_t getOobError();

  //void getImportanceValues(TreeData* trainData, vector<num_t>& importanceValues, vector<num_t>& contrastImportanceValues);
  void getMDI(TreeData* trainData, vector<num_t>& impurityValues, vector<num_t>& contrastImpurityValues);

  void predict(TreeData* testData, vector<string>& predictions, vector<num_t>& confidence, size_t nThreads = 1);
  void predict(TreeData* testData, vector<num_t>& predictions, vector<num_t>& confidence, size_t nThreads = 1);

  //bool useQuantiles() const;

  void getNumDistributions(TreeData* testData, vector<vector<num_t> >& distributions, distributions::Random* random, const size_t nSamplesPerTree);
  void getCatDistributions(TreeData* testData, vector<vector<cat_t> >& distributions, distributions::Random* random, const size_t nSamplesPerTree);

  //vector<num_t> getOobPredictions();
  //vector<num_t> getPermutedOobPredictions(const size_t featureIdx);

  //Counts the number of nodes in the forest
  //size_t nNodes();
  //size_t nNodes(const size_t treeIdx);

  size_t nTrees();

  //RootNode* tree(const size_t treeIdx) { return( rootNodes_[treeIdx] ); }

  //inline set<size_t> getFeaturesInForest() const { return( featuresInForest_ ); }
  inline string getTargetName() const { assert(rootNodes_.size() > 0); return( rootNodes_[0]->getTargetName() ); }
  inline bool isTargetNumerical() const { assert(rootNodes_.size() > 0); return( rootNodes_[0]->isTargetNumerical() ); }

  void writeForest(ofstream& toFile);

#ifndef TEST__
private:
#endif

  void readForestHeader(ifstream& forestStream);
  
  void growNumericalGBT(TreeData* trainData, const size_t targetIdx, const ForestOptions* forestOptions, const distributions::PMF* pmf, vector<distributions::Random>& randoms);
  void growCategoricalGBT(TreeData* trainData, const size_t targetIdx, const ForestOptions* forestOptions, const distributions::PMF* pmf, vector<distributions::Random>& randoms);

  //num_t error(const vector<num_t>& data1,
  //	      const vector<num_t>& data2); 

  datadefs::forest_t forestType_;

  vector<num_t> GBTConstants_;
  num_t GBTShrinkage_;

  // Root nodes for every tree
  vector<RootNode*> rootNodes_;

  // Container for all features in the forest for fast look-up
  //set<size_t> featuresInForest_;
  
};

#endif
