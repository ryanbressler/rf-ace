#ifndef ROOTNODE_HPP
#define ROOTNODE_HPP

#include <cstdlib>
#include <map>
#include <vector>
#include <set>
#include "node.hpp"
#include "treedata.hpp"

class RootNode : public Node {
public:

  RootNode(Treedata* treeData,
	   const size_t targetIdx);

  /*
    RootNode(Treedata* treeData,
    size_t targetIdx,
    bool sampleWithReplacement,
    num_t sampleSizeFraction,
    size_t maxNodesToStop,
    size_t minNodeSizeToStop,
    bool isRandomSplit,
    size_t nFeaturesForSplit,
    bool useContrasts,
    size_t numClasses,
    const PredictionFunctionType predictionFunctionType);
  */

  ~RootNode();

  void growTree(const GrowInstructions& GI);
  
  size_t nNodes();

  num_t getTrainPrediction(const size_t sampleIdx);
  vector<num_t> getTrainPrediction();

  num_t getOobError();
  size_t nOobSamples();

  set<size_t> getFeaturesInTree() { return( featuresInTree_ ); }

  num_t getImportance(const size_t featureIdx);

private:

  map<Node*,vector<size_t> > percolateSampleIcs(const vector<size_t>& sampleIcs);
  Node* percolateSampleIdx(const size_t sampleIdx);

  map<Node*,vector<size_t> > percolateSampleIcsAtRandom(const size_t featureIdx, const vector<size_t>& sampleIcs);
  Node* percolateSampleIdxAtRandom(const size_t featureIdx, const size_t sampleIdx);

  num_t getPredictionError(const map<Node*,vector<size_t> >& percolatedSampleIcs);

  // Required parameters
  Treedata* treeData_;
  size_t targetIdx_;
  
  // Parameters that are generated only when a tree is grown
  size_t nNodes_;
  vector<size_t> bootstrapIcs_;
  vector<size_t> oobIcs_;
  num_t oobError_;
  set<size_t> featuresInTree_;

  //GrowInstructions GI_;

  

};

#endif
