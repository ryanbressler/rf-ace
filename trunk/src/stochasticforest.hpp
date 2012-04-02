#ifndef STOCHASTICFOREST_HPP
#define STOCHASTICFOREST_HPP

#include<cstdlib>
#include "rootnode.hpp"
#include "treedata.hpp"
//#include "partitionsequence.hpp"

using namespace std;

/*
  TODO: inherit RF and GBT from this class, to simplify the structure
*/
class StochasticForest {
public:

  enum model_t { RF, GBT };

  struct Parameters {
    model_t model;
    size_t  nTrees;
    size_t  mTry;
    size_t  nMaxLeaves;
    size_t  nodeSize;
    num_t   inBoxFraction;
    bool    sampleWithReplacement;
    bool    useContrasts;
    bool    isRandomSplit;
    num_t   shrinkage;

  };

  StochasticForest(Treedata* treeData, const string& targetName, const Parameters& parameters);
  
  // Load an existing forest
  StochasticForest(Treedata* treeData, const string& forestFile);

  ~StochasticForest();
  
  void learnRF();
  void learnGBT();
  
  map<size_t,map<size_t,size_t> > featureFrequency();
  
  vector<num_t> importanceValues();
  vector<num_t> contrastImportanceValues();
  num_t oobError();
    
  void predict(vector<string>& categoryPrediction, vector<num_t>& confidence);
  void predict(vector<num_t>& prediction, vector<num_t>& confidence);

  //Counts the number of nodes in the forest
  vector<size_t> nNodes();
  size_t nNodes(const size_t treeIdx);

  size_t nTrees();

  inline string getTargetName() { return( targetName_ ); }
  inline bool isTargetNumerical() { return( targetSupport_.size() == 0 ? true : false ); }
  inline model_t type() { return( parameters_.model ); }

  void printToFile(const string& fileName);

#ifndef TEST__
private:
#endif

  void growNumericalGBT(const Node::GrowInstructions& GI);
  void growCategoricalGBT(const Node::GrowInstructions& GI);

  // TODO: StochasticForest::transformLogistic() should be moved elsewhere
  void transformLogistic(vector<num_t>& prediction, vector<num_t>& probability);

  // TODO: predictSampleByTree() and percolateSampleIdx families in StochasticForest need to be fused together 
  //num_t predictSampleByTree(size_t sampleIdx, size_t treeIdx);
  //vector<num_t> predictDatasetByTree(size_t treeIdx);

  //Percolates samples along the trees, starting from the rootNode. Spits out a map<> that ties the percolated train samples to the leaf nodes
  //map<Node*,vector<size_t> > percolateSampleIcsByTree(const vector<size_t>& sampleIcs, const size_t treeIdx);
  //map<Node*,vector<size_t> > percolateSampleIcsByTreeAtRandom(const size_t featureIdx, const vector<size_t>& sampleIcs, const size_t treeIdx);
  
  //Node* percolateSampleIdxByTree(const size_t sampleIdx, const size_t treeIdx);
  //Node* percolateSampleIdxByTreeAtRandom(const size_t featureIdx, const size_t sampleIdx, const size_t treeIdx);

  // Calculates prediction error across the nodes provided in the input map<> 
  //num_t predictionError(const map<Node*,vector<size_t> >& trainIcs);
  
  void updateImportanceValues();

  //Pointer to treeData_ object, stores all the feature data with which the trees are grown (i.e. training data)
  Treedata* treeData_;

  Parameters parameters_;

  //Chosen target to regress on
  string targetName_;
  vector<string> targetSupport_;

  //size_t nTrees_;

  //Root nodes for every tree
  vector<RootNode*> rootNodes_;

  vector<num_t> importanceValues_;
  vector<num_t> contrastImportanceValues_;
  num_t oobError_;

  // Maps a tree to the set of features existing in the tree
  //map<size_t, set<size_t> > featuresInForest_;

  //Out-of-box samples for each tree
  //vector<vector<size_t> > oobMatrix_;

  //enum LearnedModel { NO_MODEL, RF, GBT } learnedModel_;
  
  //num_t shrinkage_;
  
};

#endif
