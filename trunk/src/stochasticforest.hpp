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

    void validate();
  };

  StochasticForest(Treedata* treeData, const string& targetName, const Parameters& parameters);
  
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
  
  num_t error(const vector<num_t>& data1,
	      const vector<num_t>& data2); 


  //Pointer to treeData_ object, stores all the feature data with which the trees are grown (i.e. training data)
  Treedata* treeData_;

  Parameters parameters_;

  //Chosen target to regress on
  string targetName_;
  vector<string> targetSupport_;

  //Root nodes for every tree
  vector<RootNode*> rootNodes_;
  
};

#endif
