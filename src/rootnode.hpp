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

  ~RootNode();

  void growTree(const GrowInstructions& GI);
  
  size_t nNodes();

  num_t getTrainPrediction(const size_t sampleIdx);
  num_t getPermutedTrainPrediction(const size_t featureIdx,
				   const size_t sampleIdx);

  vector<num_t> getTrainPrediction();

  vector<size_t> getOobIcs();

  //num_t getOobError();
  size_t nOobSamples();

  set<size_t> getFeaturesInTree() { return(featuresInTree_); }

  //num_t getImportance(const size_t featureIdx);

#ifndef TEST__
private:
#endif

  //map<Node*,vector<size_t> > percolateSampleIcs(const vector<size_t>& sampleIcs);
  Node* percolateSampleIdx(const size_t sampleIdx);

  //map<Node*,vector<size_t> > percolateSampleIcsAtRandom(const size_t featureIdx, const vector<size_t>& sampleIcs);
  Node* percolateSampleIdxAtRandom(const size_t featureIdx, const size_t sampleIdx);

  //num_t getPredictionError(const map<Node*,vector<size_t> >& percolatedSampleIcs);

  // Required parameters
  Treedata* treeData_;
  size_t targetIdx_;
  
  // Parameters that are generated only when a tree is grown
  size_t nNodes_;
  vector<size_t> bootstrapIcs_;
  vector<size_t> oobIcs_;

  set<size_t> featuresInTree_;

  vector<num_t> trainPredictionCache_;

};

#endif
